############################################################
# HIERARCHICAL MRP MODEL IN STAN (using rstan)
# Fallback version if cmdstanr installation fails
############################################################

library(celestial)
library(Rfits)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(Cairo)

############################################################
# Data
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25

vol_max <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max^3 * 179.92*(pi/180)^2 / (4*pi)

cat("Survey volume:", signif(Vsurvey,4), "Mpc^3\n")

g3cx  <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c   <- g3cx[g3cx$Nfof >= 5 & g3cx$Zfof < zlimit &
                  g3cx$Zfof > 0.01 & g3cx$IterCenDec > -3.5, ]
g3c$MassA <- g3c$MassA * 100 / ho
g3c   <- g3c[is.finite(g3c$MassA) & g3c$MassA > 0, ]

x_obs <- log10(g3c$MassA)
x_obs <- sort(x_obs[x_obs > 12.7])
N     <- length(x_obs)

xlo <- 12.7
xhi <- 15.8
Ng  <- 500

cat("N groups:", N, "\n\n")

############################################################
# STAGE 1 STAN MODEL
############################################################

stage1_stan <- "
data {
  int<lower=0> N;
  vector[N] x;
  real V;
  real xlo;
  real xhi;
  int<lower=1> Ng;
}

transformed data {
  vector[Ng] xgrid;
  real dx;
  dx = (xhi - xlo) / (Ng - 1);
  for(k in 1:Ng)
    xgrid[k] = xlo + (k-1) * dx;
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.1, upper=2> beta;
}

model {
  // Priors
  mstar   ~ normal(13.5, 1.0);
  log_phi ~ normal(-3.5, 1.5);
  alpha   ~ normal(-1.5, 0.8);

  // Integral term
  {
    vector[Ng] phi_grid;
    real Lambda;
    for(k in 1:Ng) {
      real u = beta * (xgrid[k] - mstar);
      phi_grid[k] = beta * log(10) * pow(10, log_phi)
                    * pow(10, (alpha+1) * (xgrid[k] - mstar))
                    * exp(-pow(10, u));
    }
    Lambda = V * sum(phi_grid) * dx;
    target += -Lambda;
  }

  // Sum over observations
  for(i in 1:N) {
    real u = beta * (x[i] - mstar);
    real log_phi_i = log(beta) + log(log(10))
                     + log_phi * log(10)
                     + (alpha+1) * (x[i] - mstar) * log(10)
                     - pow(10, u);
    target += log(V) + log_phi_i;
  }
}
"

############################################################
# Prepare data
############################################################

stan_data <- list(
    N    = N,
    x    = x_obs,
    V    = Vsurvey,
    xlo  = xlo,
    xhi  = xhi,
    Ng   = Ng
)

############################################################
# Run Stan
############################################################

cat("Compiling Stan model (may take a minute)...\n")
fit1 <- stan(
    model_code = stage1_stan,
    data       = stan_data,
    chains     = 4,
    iter       = 3000,
    warmup     = 1000,
    thin       = 1,
    cores      = 4,
    init       = lapply(1:4, function(i) list(
        mstar   = rnorm(1, 13.5, 0.2),
        log_phi = rnorm(1, -3.5, 0.2),
        alpha   = rnorm(1, -1.5, 0.1),
        beta    = runif(1, 0.4, 0.6)
    )),
    control    = list(adapt_delta = 0.9)
)

cat("\n=== POSTERIOR SUMMARY ===\n")
print(fit1, pars=c("mstar","log_phi","alpha","beta"))

############################################################
# Check convergence
############################################################

cat("\n=== CONVERGENCE CHECK ===\n")
cat("Rhat values (should all be ~1.00):\n")
print(summary(fit1)$summary[c("mstar","log_phi","alpha","beta"), "Rhat"])

cat("\nEffective sample size:\n")
print(summary(fit1)$summary[c("mstar","log_phi","alpha","beta"), "n_eff"])

############################################################
# Extract posterior
############################################################

posterior1 <- extract(fit1, pars=c("mstar","log_phi","alpha","beta"))
posterior_matrix <- cbind(
    mstar   = posterior1$mstar,
    log_phi = posterior1$log_phi,
    alpha   = posterior1$alpha,
    beta    = posterior1$beta
)

med <- apply(posterior_matrix, 2, median)
q16 <- apply(posterior_matrix, 2, quantile, 0.16)
q84 <- apply(posterior_matrix, 2, quantile, 0.84)

############################################################
# Plot
############################################################

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
              10^((alpha+1)*(x-mstar)) *
              exp(-10^(beta*(x-mstar))))
}

breaks   <- seq(12.5, 16.2, by=0.2)
hist_plt <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_plt  <- hist_plt$counts / (Vsurvey * 0.2)
ok       <- phi_plt > 0
xfit     <- seq(12, 16, length.out=500)

CairoPDF("MRP_STAGE1_STAN.pdf", 10, 7)

plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 1: Unbinned Poisson Point Process (Stan/HMC)")

grid(col="gray80")

idx <- sample(1:nrow(posterior_matrix), 300)
for(i in idx) {
    y <- tryCatch(
        mrp_log10(posterior_matrix[i,"mstar"], posterior_matrix[i,"log_phi"],
                  posterior_matrix[i,"alpha"],  posterior_matrix[i,"beta"], xfit),
        error = function(e) rep(NA, length(xfit))
    )
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
}

lines(xfit, mrp_log10(med["mstar"],med["log_phi"],
                      med["alpha"],med["beta"],xfit),
      col="red", lwd=3)

lines(xfit, mrp_log10(13.51,-3.19,-1.27,0.47,xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("GAMA (binned for display)",
                "Posterior median","Posterior samples","Driver+22"),
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
            med["beta"],   q84["beta"]-med["beta"],    med["beta"]-q16["beta"]))

cat("\nPlot saved: MRP_STAGE1_STAN.pdf\n")

saveRDS(list(
    fit = fit1,
    posterior = posterior_matrix,
    median = med,
    q16 = q16,
    q84 = q84
), "stage1_stan_results.rds")

cat("Results saved: stage1_stan_results.rds\n")

############################################################
# Trace plots
############################################################

CairoPDF("MCMC_traces_stan.pdf", 12, 8)
traceplot(fit1, pars=c("mstar","log_phi","alpha","beta"), inc_warmup=FALSE)
dev.off()

cat("Trace plots: MCMC_traces_stan.pdf\n")
