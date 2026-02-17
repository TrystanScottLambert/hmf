############################################################
# PROPERLY SPECIFIED JAGS MODEL
# Fixed: Correct prior precisions and better initialization
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(Cairo)
library(rjags)
library(coda)

set.seed(42)  # For reproducibility

############################################################
# Load and prepare data
############################################################

ho <- 67.37
omegam <- 0.3147
zlimit <- 0.25

# CORRECT volume calculation
vol_max <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
vol_sphere <- (4/3) * pi * vol_max^3
area_sr <- 179.92 * (pi/180)^2
Vsurvey <- vol_sphere * area_sr / (4*pi)

cat("Volume:", Vsurvey, "Mpc^3 =", Vsurvey/1e9, "Gpc^3\n")

if(Vsurvey < 1e6) {
    stop("Volume is wrong! Check cosdistCoDist function.")
}

# Load data
g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")

g3c <- g3cx[
    g3cx$Nfof >= 5 &
        g3cx$Zfof < zlimit &
        g3cx$Zfof > 0.01 &
        g3cx$IterCenDec > -3.5,
]

g3c$MassA <- g3c$MassA * 100 / ho
x_obs <- log10(g3c$MassA[is.finite(g3c$MassA) & g3c$MassA > 0])
x_obs <- x_obs[x_obs > 12.7]

cat("N groups:", length(x_obs), "\n")
cat("Mass range:", range(x_obs), "\n\n")

# Create bins
breaks <- seq(12.5, 16, by=0.2)
hist_data <- hist(x_obs, breaks=breaks, plot=FALSE)

counts <- hist_data$counts
x_mids <- hist_data$mids

# Remove empty bins
keep <- counts > 0
counts <- counts[keep]
x_mids <- x_mids[keep]

cat("Bins:", length(counts), "\n")
cat("Total groups:", sum(counts), "\n\n")

############################################################
# JAGS MODEL - CORRECTED PRIORS
############################################################

# Note: JAGS uses PRECISION = 1/variance for normal distributions
# So dnorm(mean, 1/sd^2)

model_string <- "
model {
  # Priors - using precision = 1/variance
  mstar   ~ dnorm(13.5, 1.0)     # precision = 1/1^2 = 1
  log_phi ~ dnorm(-3.5, 0.444)   # precision = 1/1.5^2 = 0.444
  alpha   ~ dnorm(-1.5, 1.5625)  # precision = 1/0.8^2 = 1.5625
  beta    ~ dunif(0.2, 1.5)

  # Transform log_phi to phi_star
  phi_star <- pow(10, log_phi)

  # Likelihood
  for(i in 1:n_bins) {

    # MRP function
    # phi(x) = beta * ln(10) * phi_star * 10^((alpha+1)*(x-mstar)) * exp(-10^(beta*(x-mstar)))

    log_mass_diff[i] <- x_bin[i] - mstar

    power_term[i] <- pow(10, (alpha + 1) * log_mass_diff[i])
    exp_term[i] <- exp(-pow(10, beta * log_mass_diff[i]))

    phi[i] <- beta * 2.302585 * phi_star * power_term[i] * exp_term[i]

    # Expected count in bin
    lambda[i] <- phi[i] * V * bin_width

    # Likelihood
    counts[i] ~ dpois(lambda[i])
  }
}
"

writeLines(model_string, "mrp_corrected.jags")

############################################################
# Initial values - start near the truth
############################################################

# Use Driver et al. parameters as starting point
inits <- list(
    list(mstar = 13.5, log_phi = -3.2, alpha = -1.3, beta = 0.5),
    list(mstar = 13.6, log_phi = -3.3, alpha = -1.4, beta = 0.4),
    list(mstar = 13.4, log_phi = -3.4, alpha = -1.2, beta = 0.6)
)

############################################################
# JAGS data
############################################################

jags_data <- list(
    counts = counts,
    x_bin = x_mids,
    n_bins = length(counts),
    V = Vsurvey,
    bin_width = 0.2
)

cat("JAGS data:\n")
cat("  n_bins =", jags_data$n_bins, "\n")
cat("  V =", jags_data$V, "Mpc^3\n")
cat("  bin_width =", jags_data$bin_width, "dex\n\n")

############################################################
# Run JAGS
############################################################

cat("Compiling model with good initial values...\n")
model <- jags.model("mrp_corrected.jags",
                    data=jags_data,
                    inits=inits,
                    n.chains=3,
                    n.adapt=3000)

cat("\nBurn-in (10000 iterations)...\n")
update(model, 10000, progress.bar="text")

cat("\nSampling (30000 iterations)...\n")
samples <- coda.samples(model,
                        c("mstar", "log_phi", "alpha", "beta"),
                        n.iter=30000,
                        thin=15,
                        progress.bar="text")

posterior <- as.matrix(samples)

cat("\n=== POSTERIOR SUMMARY ===\n")
print(summary(samples))

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
gelman_result <- gelman.diag(samples)
print(gelman_result)

if(max(gelman_result$psrf[,1]) > 1.1) {
    cat("\nWARNING: Gelman-Rubin > 1.1, chains haven't converged!\n")
    cat("Consider running longer or checking model specification.\n")
} else {
    cat("\nGOOD: Chains have converged (R-hat < 1.1)\n")
}

cat("\nEffective sample size:\n")
print(effectiveSize(samples))

############################################################
# Plot results
############################################################

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    phi_star <- 10^log_phi
    phi <- beta * log(10) * phi_star *
        10^((alpha + 1) * (x - mstar)) *
        exp(-10^(beta * (x - mstar)))
    log10(phi)
}

# For plotting, use all bins including empty ones
hist_full <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_emp <- hist_full$counts / (Vsurvey * 0.2)
x_emp <- hist_full$mids

xfit <- seq(12, 16, length.out=500)

# Driver et al. parameters
driver <- list(mstar=13.51, log_phi=-3.19, alpha=-1.27, beta=0.47)

CairoPDF("MRP_FIT_FINAL.pdf", 10, 7)

plot(x_emp, log10(phi_emp),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12, 16), ylim=c(-8, -2),
     xlab=expression("Halo Mass   log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")   [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="MRP Fit to GAMA Data")

grid(col="gray80")

# Plot posterior samples
idx <- sample(1:nrow(posterior), 200)
for(i in idx) {
    lines(xfit,
          mrp_log10(posterior[i,"mstar"],
                    posterior[i,"log_phi"],
                    posterior[i,"alpha"],
                    posterior[i,"beta"],
                    xfit),
          col=rgb(0,0,1,0.03))
}

# Median fit
med <- apply(posterior, 2, median)
lines(xfit,
      mrp_log10(med["mstar"],
                med["log_phi"],
                med["alpha"],
                med["beta"],
                xfit),
      col="red", lwd=3)

# Driver et al.
lines(xfit,
      mrp_log10(driver$mstar, driver$log_phi, driver$alpha, driver$beta, xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("GAMA data", "Your fit", "Posterior samples", "Driver+22"),
       col=c("darkgreen", "red", rgb(0,0,1,0.2), "black"),
       pch=c(19, NA, NA, NA),
       lty=c(NA, 1, 1, 2),
       lwd=c(NA, 3, 1, 2),
       bg="white")

text(12.5, -2.5,
     sprintf("M* = %.2f\nlog(phi*) = %.2f\nalpha = %.2f\nbeta = %.2f",
             med["mstar"], med["log_phi"], med["alpha"], med["beta"]),
     pos=4, col="red", cex=0.9)

dev.off()

cat("\n=== PLOT SAVED ===\n")
cat("File: MRP_FIT_FINAL.pdf\n\n")

############################################################
# Plot traces
############################################################

CairoPDF("MCMC_traces.pdf", 12, 9)
plot(samples)
dev.off()

cat("Trace plots: MCMC_traces.pdf\n\n")

############################################################
# Comparison
############################################################

cat("=== PARAMETER COMPARISON ===\n")
cat("Parameter | Your fit | Driver+22\n")
cat("----------|----------|-----------\n")
cat(sprintf("M*        | %8.2f | %8.2f\n", med["mstar"], driver$mstar))
cat(sprintf("log(phi*) | %8.2f | %8.2f\n", med["log_phi"], driver$log_phi))
cat(sprintf("alpha     | %8.2f | %8.2f\n", med["alpha"], driver$alpha))
cat(sprintf("beta      | %8.2f | %8.2f\n", med["beta"], driver$beta))

# Save results
saveRDS(list(
    posterior = posterior,
    median = med,
    samples = samples,
    convergence = gelman_result
), "final_results.rds")

cat("\n=== SUCCESS! ===\n")
cat("Results saved to: final_results.rds\n")
