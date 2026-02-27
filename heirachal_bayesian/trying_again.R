############################################################
# FINAL: Monte Carlo Eddington Bias + Stan
#
# Exactly replicating gamahmf.r approach:
# 1. Vmax-weighted bins
# 2. Monte Carlo Eddington bias correction (1000 iterations)
# 3. Fit corrected bins in Stan for Bayesian posteriors
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
library(plotrix)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(Cairo)

set.seed(42)

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.015
multi  <- 5
mass_limit <- 12.7
logbin <- 0.2

vol_max_survey <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max_survey^3 * 179.92*(pi/180)^2 / (4*pi)

cat("Survey volume:", signif(Vsurvey,4), "Mpc^3\n")

############################################################
# Data prep - EXACT from gamahmf.r
############################################################

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

xx <- seq(3, 22)
yy <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 
        0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506, 
        0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983, 
        0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)

g3c$log10MassErr <- approx(xx, yy, g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)] <- 0.03
g3c$log10MassErr[g3c$log10MassErr < 0.1] <- 0.1

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
g3c$MassAfunc[g3c$GroupID == 100622] <- 1E9

############################################################
# Vmax
############################################################

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

############################################################
# Vmax-weighted bins (uncorrected)
############################################################

massx <- seq(10.3, 16.1, logbin)
gamahmf_uncorrected <- weighted.hist(log10(g3c$MassAfunc), 
                                     w=1/g3c$weightszlimit, 
                                     breaks=massx, 
                                     plot=FALSE)

############################################################
# MONTE CARLO EDDINGTON BIAS (gamahmf.r lines 362-373)
############################################################

cat("\nRunning Monte Carlo Eddington bias correction (1001 iterations)...\n")

mockcounts <- matrix(0, nrow=1001, ncol=length(gamahmf_uncorrected$mids))

for(iter in 1:1001) {
    if(iter %% 200 == 0) cat("  Iteration", iter, "/", 1001, "\n")
    
    # Perturb masses by their errors
    mockmass <- log10(g3c$MassAfunc) + rnorm(nrow(g3c), 0.0, g3c$log10MassErr)
    
    # Create Vmax-weighted histogram with perturbed masses
    mock_hist <- weighted.hist(mockmass, 
                               w=1/g3c$weightszlimit, 
                               breaks=massx, 
                               plot=FALSE)
    
    mockcounts[iter, ] <- mock_hist$counts
}

# Calculate Eddington bias correction
meancounts <- apply(mockcounts, 2, mean)
edb <- meancounts / gamahmf_uncorrected$counts
edb[is.infinite(edb)] <- 1.0
edb[is.na(edb)] <- 1.0

# Apply correction
gamax <- gamahmf_uncorrected$mids
gamay_corrected <- (gamahmf_uncorrected$counts / logbin) / edb  # Corrected phi
gamay_uncorrected <- gamahmf_uncorrected$counts / logbin  # Uncorrected for comparison

cat("Eddington bias factors:", round(edb, 3), "\n\n")

# Apply mass cut
ok <- gamay_corrected > 0 & gamax > mass_limit
gamax <- gamax[ok]
gamay_corrected <- gamay_corrected[ok]
gamay_uncorrected <- gamay_uncorrected[ok]
edb <- edb[ok]

N_bins <- length(gamax)

cat("N bins (M>12.7):", N_bins, "\n")
cat("EdB correction range:", round(range(edb), 3), "\n\n")

############################################################
# Stan: Fit to CORRECTED bins
############################################################

stan_model <- "
data {
  int<lower=0> N;
  vector[N] x;
  vector<lower=0>[N] y;  // Corrected phi
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.1, upper=2> beta;
}

model {
  mstar   ~ normal(13.5, 1.0);
  log_phi ~ normal(-3.5, 1.5);
  alpha   ~ normal(-1.5, 0.8);
  
  for(i in 1:N) {
    real u = beta * (x[i] - mstar);
    real log_phi_pred = log10(beta) + log10(log(10))
                        + log_phi
                        + (alpha+1) * (x[i] - mstar)
                        - 10^u / log(10);
    
    log10(y[i]) ~ normal(log_phi_pred, 0.15);
  }
}
"

stan_data <- list(
    N = N_bins,
    x = gamax,
    y = gamay_corrected
)

cat("Fitting to EdB-corrected bins in Stan...\n\n")

fit <- stan(
    model_code = stan_model,
    data       = stan_data,
    chains     = 4,        # More chains
    iter       = 3000,     # More iterations
    warmup     = 1500,     # Longer warmup
    cores      = 4,
    control = list(adapt_delta = 0.98, max_treedepth = 12)  # More aggressive
)

print(fit, pars=c("mstar","log_phi","alpha","beta"))

############################################################
# Plot
############################################################

posterior <- extract(fit, pars=c("mstar","log_phi","alpha","beta"))

posterior_matrix <- cbind(
    mstar   = posterior$mstar,
    log_phi = posterior$log_phi,
    alpha   = posterior$alpha,
    beta    = posterior$beta
)

med <- c(
    mstar   = median(posterior$mstar),
    log_phi = median(posterior$log_phi),
    alpha   = median(posterior$alpha),
    beta    = median(posterior$beta)
)

q16 <- c(
    mstar   = as.numeric(quantile(posterior$mstar, 0.16)),
    log_phi = as.numeric(quantile(posterior$log_phi, 0.16)),
    alpha   = as.numeric(quantile(posterior$alpha, 0.16)),
    beta    = as.numeric(quantile(posterior$beta, 0.16))
)

q84 <- c(
    mstar   = as.numeric(quantile(posterior$mstar, 0.84)),
    log_phi = as.numeric(quantile(posterior$log_phi, 0.84)),
    alpha   = as.numeric(quantile(posterior$alpha, 0.84)),
    beta    = as.numeric(quantile(posterior$beta, 0.84))
)

cat("\nDEBUG: q16 =", q16, "\n")
cat("DEBUG: q84 =", q84, "\n")
cat("DEBUG: med =", med, "\n")

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))))
}

############################################################
# Plot 1: Main HMF plot
############################################################

xfit <- seq(12, 16, length.out=500)

CairoPDF("MRP_WITH_EDDINGTON_BIAS.pdf", 10, 7)

plot(gamax, log10(gamay_uncorrected),
     pch=1, col="gray60", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Vmax bins + Eddington bias (Monte Carlo)")

points(gamax, log10(gamay_corrected), pch=19, col="darkgreen", cex=1.2)

grid(col="gray80")

# Add posterior samples
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
       legend=c("Uncorrected (Vmax)", "EdB corrected", 
                "Posterior samples", "Stan fit", "Driver+22"),
       col=c("gray60","darkgreen",rgb(0,0,1,0.3),"red","black"),
       pch=c(1,19,NA,NA,NA), lty=c(NA,NA,1,1,2),
       lwd=c(NA,NA,1,3,2), bg="white")

text(14, -2.5, sprintf(
    "M* = %.2f (+%.2f/-%.2f)\nlog φ* = %.2f (+%.2f/-%.2f)\nα = %.2f (+%.2f/-%.2f)\nβ = %.2f (+%.2f/-%.2f)",
    med["mstar"], q84["mstar"]-med["mstar"], med["mstar"]-q16["mstar"],
    med["log_phi"], q84["log_phi"]-med["log_phi"], med["log_phi"]-q16["log_phi"],
    med["alpha"], q84["alpha"]-med["alpha"], med["alpha"]-q16["alpha"],
    med["beta"], q84["beta"]-med["beta"], med["beta"]-q16["beta"]),
    col="red", pos=4, cex=0.85)

dev.off()

############################################################
# Plot 2: Corner plot with covariances
############################################################

CairoPDF("MRP_COVARIANCE.pdf", 10, 10)

par(mfrow=c(4,4), mar=c(3,3,1,1), oma=c(0,0,2,0))

params <- c("mstar", "log_phi", "alpha", "beta")
labels <- c("M*", "log φ*", "α", "β")

# Driver+22 values for reference
driver_vals <- c(13.51, -3.19, -1.27, 0.47)

# Axis limits matching Driver+22 range
xlims <- list(
    mstar   = c(11, 16),
    log_phi = c(-5, -1),
    alpha   = c(-2.5, 0),
    beta    = c(0, 1.5)
)

for(i in 1:4) {
    for(j in 1:4) {
        if(i == j) {
            # Diagonal: histogram
            hist(posterior_matrix[,params[i]], 
                 breaks=30, col="skyblue", border="white",
                 main="", xlab="", ylab="", xlim=xlims[[params[i]]])
            abline(v=med[params[i]], col="red", lwd=2)
            abline(v=driver_vals[i], col="black", lwd=2, lty=2)
            if(j == 1) mtext(labels[i], side=2, line=2, cex=0.8)
        } else if(i > j) {
            # Lower triangle: 2D density
            plot(posterior_matrix[,params[j]], posterior_matrix[,params[i]],
                 pch=16, col=rgb(0,0,0,0.1), cex=0.3,
                 xlab="", ylab="", xlim=xlims[[params[j]]], ylim=xlims[[params[i]]])
            points(med[params[j]], med[params[i]], pch=3, col="red", cex=2, lwd=2)
            points(driver_vals[j], driver_vals[i], pch=4, col="black", cex=2, lwd=2)
            if(j == 1) mtext(labels[i], side=2, line=2, cex=0.8)
            if(i == 4) mtext(labels[j], side=1, line=2, cex=0.8)
        } else {
            # Upper triangle: correlation coefficient
            corr <- cor(posterior_matrix[,params[j]], posterior_matrix[,params[i]])
            plot.new()
            text(0.5, 0.5, sprintf("ρ = %.3f", corr), cex=1.5)
        }
    }
}

mtext("Parameter Covariances - Posterior", outer=TRUE, cex=1.2, font=2)

dev.off()

cat("\n=== FINAL RESULTS ===\n")
cat(sprintf("           Median    +1sig    -1sig   Driver+22\n"))
cat(sprintf("M*:        %6.3f   +%.3f   -%.3f    13.51\n",
            med["mstar"], q84["mstar"]-med["mstar"], med["mstar"]-q16["mstar"]))
cat(sprintf("log_phi:   %6.3f   +%.3f   -%.3f    -3.19\n",
            med["log_phi"], q84["log_phi"]-med["log_phi"], med["log_phi"]-q16["log_phi"]))
cat(sprintf("alpha:     %6.3f   +%.3f   -%.3f    -1.27\n",
            med["alpha"], q84["alpha"]-med["alpha"], med["alpha"]-q16["alpha"]))
cat(sprintf("beta:      %6.3f   +%.3f   -%.3f     0.47\n",
            med["beta"], q84["beta"]-med["beta"], med["beta"]-q16["beta"]))

cat("\nPlots:\n")
cat("  MRP_WITH_EDDINGTON_BIAS.pdf - HMF with EdB correction\n")
cat("  MRP_COVARIANCE.pdf          - Parameter covariances\n")
cat("\n✓ Vmax weighting (bins)\n")
cat("✓ Eddington bias (Monte Carlo, 1001 iterations, like gamahmf.r)\n")
cat("✓ Bayesian posteriors (Stan, not bootstrap)\n")
cat("✓ Full covariance matrix\n")
cat("\nThis is gamahmf.r's method with Stan for uncertainties!\n")
