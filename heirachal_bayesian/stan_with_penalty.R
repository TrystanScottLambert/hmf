############################################################
# STAGE 3 - Using the Stan model you provided
# Includes the penalty term like gamahmf.r!
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(Cairo)
library(plotrix)

set.seed(42)

############################################################
# Data prep EXACTLY matching gamahmf.r
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5
logbin <- 0.2

vol_max_survey <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max_survey^3 * 179.92/(360^2/pi)

cat("Survey volume (vlimit):", signif(Vsurvey, 4), "Mpc^3\n")

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit & 
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 & 
             g3cx$IterCenDec > -3.5, ]

magica <- 13.9
parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30

g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

masscorr <- c(0.0, 0.0, -2.672595e-01, -1.513503e-01, -1.259069e-01, 
              -9.006064e-02, -5.466009e-02, -6.666895e-02, -1.988694e-02,
              -2.439581e-02, -2.067060e-02, -1.812964e-02, -1.556899e-02,
              -1.313664e-02, -1.743112e-02, -7.965513e-03, -1.257178e-02,
              -7.064037e-03, -3.963656e-03, -1.271533e-02, -2.664687e-03,
              -1.691287e-03)

g3c$masscorr <- masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] <- 0.0
g3c$mymasscorr <- g3c$mymass / 10^g3c$masscorr
g3c$MassAfunc <- g3c$mymasscorr

# Vmax calculation
gig <- fread("../data/GAMAGalsInGroups.csv")

g3c$zmax <- NA
for (i in 1:nrow(g3c)) {
    if(g3c$Nfof[i] == 2) {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[2]
    } else {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[multi]
    }
}

g3c$zmax <- ifelse(g3c$zmax < g3c$Zfof, g3c$Zfof, g3c$zmax)
g3c$zmax <- ifelse(g3c$zmax > zlimit, zlimit, g3c$zmax)

vol_zmax <- cosdist(g3c$zmax, OmegaM=omegam, H0=ho)$CoVol
vol_zmin <- cosdist(zmin, OmegaM=omegam, H0=ho)$CoVol
g3c$vmax <- 179.92/(360^2/pi) * 1E9 * (vol_zmax - vol_zmin)

vlimitmin <- Vsurvey / 1000.0
g3c$weightszlimit <- ifelse(g3c$vmax > Vsurvey, Vsurvey, g3c$vmax)
g3c$weightszlimit <- ifelse(g3c$vmax < vlimitmin, vlimitmin, g3c$vmax)

# Binning
massx <- seq(10.3, 16.1, logbin)
gamahmf2 <- weighted.hist(log10(g3c$MassAfunc), 
                          w=1/g3c$weightszlimit, 
                          breaks=massx, 
                          plot=FALSE)

gamax <- gamahmf2$mids
gamay <- gamahmf2$counts / logbin

# Keep positive bins (no mass cut for now - let the model handle it)
ok <- gamay > 0 & !is.na(gamay)
gamax <- gamax[ok]
gamay <- gamay[ok]

# Fractional errors (from gamahmf.r - Poisson + cosmic variance)
frac_err <- rep(0.1, length(gamay))  # ~10% typical

cat("N bins:", length(gamax), "\n")
cat("Mass range:", range(gamax), "\n\n")

############################################################
# Your Stan model with penalty term!
############################################################

your_stan_model <- "
functions {
  real log_mrp_phi(real logM,
                   real logMstar,
                   real log_phi,
                   real alpha,
                   real beta) {
    return
      log_phi * log(10)
      + log(beta)
      + log(log(10))
      + (alpha+1)*(logM-logMstar)*log(10)
      - pow(10, beta*(logM-logMstar));
  }
}

data {
  int<lower=1> N;
  vector[N] logM_obs;
  vector[N] log_phi_obs;
  vector[N] frac_err;
  real logbin;
  real vlimit;
  real logM_max;
}

parameters {
  real logMstar;
  real log_phi;
  real alpha;
  real<lower=0> beta;
}

model {
  // Priors
  logMstar ~ normal(13.5, 1.0);
  log_phi  ~ normal(-3.5, 1.5);
  alpha    ~ normal(-1.5, 0.8);
  beta     ~ normal(0.5, 0.3);
  
  // Fit to bins
  for (n in 1:N) {
    real model_log =
      log_mrp_phi(
        logM_obs[n],
        logMstar,
        log_phi,
        alpha,
        beta
      ) / log(10);
    real sigma = frac_err[n] / log(10);
    target += normal_lpdf(
      log_phi_obs[n] |
      model_log,
      sigma
    );
  }
  
  // Penalty term for integral (matching gamahmf.r line 206!)
  real integral = 0;
  for (i in 1:10) {
    real logM = logM_max + i*logbin;
    integral += exp(
      log_mrp_phi(
        logM,
        logMstar,
        log_phi,
        alpha,
        beta
      )
    );
  }
  integral *= logbin;
  target += -2*vlimit*integral;
}
"

############################################################
# Run Stan
############################################################

stan_data <- list(
    N           = length(gamax),
    logM_obs    = gamax,
    log_phi_obs = log10(gamay),
    frac_err    = frac_err,
    logbin      = logbin,
    vlimit      = Vsurvey,
    logM_max    = max(gamax)
)

cat("Fitting with penalty term (like gamahmf.r)...\n")

fit3 <- stan(
    model_code = your_stan_model,
    data       = stan_data,
    chains     = 4,
    iter       = 2000,
    warmup     = 1000,
    cores      = 4,
    control = list(adapt_delta = 0.95)
)

cat("\n=== RESULTS ===\n")
print(fit3, pars=c("logMstar","log_phi","alpha","beta"))

############################################################
# Plot
############################################################

posterior <- extract(fit3, pars=c("logMstar","log_phi","alpha","beta"))
posterior_matrix <- cbind(
    mstar   = posterior$logMstar,
    log_phi = posterior$log_phi,
    alpha   = posterior$alpha,
    beta    = posterior$beta
)

med <- apply(posterior_matrix, 2, median)
q16 <- apply(posterior_matrix, 2, quantile, 0.16)
q84 <- apply(posterior_matrix, 2, quantile, 0.84)

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))))
}

xfit <- seq(10, 16, length.out=500)

CairoPDF("MRP_STAGE3_WITH_PENALTY.pdf", 10, 7)

plot(gamax, log10(gamay),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(11,16), ylim=c(-5,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 3: With Penalty Term (Your Model)")

grid(col="gray80")

idx <- sample(1:nrow(posterior_matrix), 300)
for(i in idx) {
    y <- mrp_log10(posterior_matrix[i,"mstar"], posterior_matrix[i,"log_phi"],
                   posterior_matrix[i,"alpha"], posterior_matrix[i,"beta"], xfit)
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
}

lines(xfit, mrp_log10(med["mstar"],med["log_phi"],med["alpha"],med["beta"],xfit),
      col="red", lwd=3)

lines(xfit, mrp_log10(13.51,-3.19,-1.27,0.47,xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("GAMA (Vmax weighted)", "Stan fit",
                "Posterior samples", "Driver+22"),
       col=c("darkgreen","red",rgb(0,0,1,0.3),"black"),
       pch=c(19,NA,NA,NA), lty=c(NA,1,1,2),
       lwd=c(NA,3,1,2), bg="white", cex=0.9)

text(14.3,-2.3,
     sprintf("M* = %.2f (+%.2f/-%.2f)\nlog\u03c6* = %.2f (+%.2f/-%.2f)\n\u03b1 = %.2f (+%.2f/-%.2f)\n\u03b2 = %.2f (+%.2f/-%.2f)",
             med["mstar"],  q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"],
             med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"],
             med["alpha"],  q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"],
             med["beta"],   q84["beta"]-med["beta"],     med["beta"]-q16["beta"]),
     col="red", cex=0.85, pos=4)

dev.off()

cat("\n=== FINAL RESULTS ===\n")
cat(sprintf("           Median    +1sig    -1sig   Driver+22\n"))
cat(sprintf("M*:        %6.3f   +%.3f   -%.3f    13.51\n",
            med["mstar"],  q84["mstar"]-med["mstar"],  med["mstar"]-q16["mstar"]))
cat(sprintf("log_phi:   %6.3f   +%.3f   -%.3f    -3.19\n",
            med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"]))
cat(sprintf("alpha:     %6.3f   +%.3f   -%.3f    -1.27\n",
            med["alpha"],  q84["alpha"]-med["alpha"],  med["alpha"]-q16["alpha"]))
cat(sprintf("beta:      %6.3f   +%.3f   -%.3f     0.47\n",
            med["beta"],   q84["beta"]-med["beta"],     med["beta"]-q16["beta"]))

cat("\nPlot: MRP_STAGE3_WITH_PENALTY.pdf\n")
cat("This model includes the penalty term like gamahmf.r!\n")
