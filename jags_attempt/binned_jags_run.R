############################################################
# COMPLETELY FIXED VERSION
# Issues found:
# 1. Volume was 600 million times too small!
# 2. Mass range included groups way too low in mass
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(Cairo)
library(rjags)
library(coda)

############################################################
# Cosmology & survey - CAREFUL VOLUME CALCULATION
############################################################

ho      = 67.37
omegam  = 0.3147
omegal  = 1 - omegam

zlimit  = 0.25
zmin    = 0.01
multi   = 5  # Use multiplicity >= 5 like Driver et al. GAMA5

# Calculate comoving volume to zlimit
vol_max = cosdistCoDist(zlimit, OmegaM=omegam, OmegaL=omegal, H0=ho)
cat("Comoving distance to z=", zlimit, ":", vol_max, "Mpc\n")

# Volume of sphere
vol_sphere = (4/3) * pi * vol_max^3
cat("Volume of sphere to z=", zlimit, ":", vol_sphere, "Mpc^3\n")
cat("Volume of sphere to z=", zlimit, ":", vol_sphere/1e9, "Gpc^3\n")

# Survey area in steradians
area_deg = 179.92  # square degrees
area_sr = area_deg * (pi/180)^2
cat("Survey area:", area_deg, "sq deg =", area_sr, "steradians\n")

# Survey volume = sphere volume × (survey area / full sky)
Vsurvey = vol_sphere * area_sr / (4*pi)

cat("\n=== VOLUME CALCULATION ===\n")
cat("Survey volume:", Vsurvey, "Mpc^3\n")
cat("Survey volume:", Vsurvey/1e9, "Gpc^3\n")
cat("Expected: ~1.26e7 Mpc^3 = 0.0126 Gpc^3\n\n")

# Sanity check
if(Vsurvey < 1e6 | Vsurvey > 1e8) {
    cat("ERROR: Volume is way off! Got", Vsurvey, "but expected ~1.26e7\n")
    cat("Stopping here - please check volume calculation\n")
    quit(status=1)
}

############################################################
# Data
############################################################

g3cx = Rfits_read_table("../data/G3CFoFGroupv10.fits")

cat("=== APPLYING CUTS ===\n")
cat("Starting groups:", nrow(g3cx), "\n")

g3c = g3cx[
    g3cx$Nfof >= multi &                    # Multiplicity cut
        g3cx$Zfof < zlimit &                # Upper redshift
        g3cx$Zfof > zmin &                  # Lower redshift
        g3cx$IterCenDec > -3.5,             # Declination cut
]

cat("After cuts:", nrow(g3c), "groups\n")

# Convert masses
g3c$MassA = g3c$MassA * 100 / ho
g3c = g3c[is.finite(g3c$MassA) & g3c$MassA > 0, ]

x_obs = log10(g3c$MassA)

# IMPORTANT: Only use groups with M > 10^12.7 M_sun
# This is the completeness limit mentioned in Driver et al.
mass_limit = 12.7
x_obs = x_obs[x_obs > mass_limit]
N_total = length(x_obs)

xmin = max(mass_limit, min(x_obs))  # Don't go below completeness limit
xmax = max(x_obs)

cat("\n=== FINAL SAMPLE ===\n")
cat("Number of groups:", N_total, "\n")
cat("Mass limit: 10^", mass_limit, "M_sun\n")
cat("Log mass range:", round(xmin, 2), "to", round(xmax, 2), "\n")
cat("Mass range:", signif(10^xmin, 3), "to", signif(10^xmax, 3), "M_sun\n\n")

# Expected based on Driver et al.
cat("Driver et al. GAMA5 (mult >= 5, M > 10^12.7):\n")
cat("  N = 1732 groups\n")
cat("  Your N =", N_total, "groups\n")

if(N_total < 1000 | N_total > 3000) {
    cat("WARNING: Sample size looks unusual!\n")
}

############################################################
# Create binned data
############################################################

logbin = 0.2
breaks = seq(floor(xmin), ceiling(xmax), by=logbin)
hist_data = hist(x_obs, breaks=breaks, plot=FALSE)

counts = hist_data$counts
x_mids = hist_data$mids
n_bins = length(counts)

# Remove bins with zero counts
keep = counts > 0
counts = counts[keep]
x_mids = x_mids[keep]
n_bins = length(counts)

cat("\n=== BINNED DATA ===\n")
cat("Number of bins:", n_bins, "\n")
cat("Mass range:", round(min(x_mids), 2), "to", round(max(x_mids), 2), "\n")
cat("Counts per bin: min =", min(counts), ", max =", max(counts), "\n")

# Sanity check on phi
phi_test = counts / (Vsurvey * logbin)
cat("log10(phi) range:", round(min(log10(phi_test)), 2), "to",
    round(max(log10(phi_test)), 2), "\n")
cat("Expected: approximately -6 to -3\n\n")

if(min(log10(phi_test)) > -2) {
    cat("ERROR: phi values are way too high!\n")
    cat("Volume is probably still wrong.\n")
    quit(status=1)
}

############################################################
# JAGS MODEL
############################################################

model_string <- "
model {

  # Priors based on Driver et al. GAMA5
  mstar   ~ dnorm(13.5, 1/1.0^2)
  log_phi ~ dnorm(-3.5, 1/1.5^2)
  alpha   ~ dnorm(-1.5, 1/0.8^2)
  beta    ~ dunif(0.2, 1.5)

  phi_star <- exp(log_phi * log(10))

  # Likelihood
  for (i in 1:n_bins) {

    phi_center[i] <-
        beta * log(10) * phi_star *
        pow(10, (alpha + 1) * (x_bin[i] - mstar)) *
        exp(- pow(10, beta * (x_bin[i] - mstar)))

    lambda[i] <- phi_center[i] * V * bin_width

    counts[i] ~ dpois(lambda[i])
  }
}
"

writeLines(model_string, "mrp_fixed.jags")

############################################################
# Prepare JAGS data
############################################################

jags_data <- list(
    counts = counts,
    x_bin = x_mids,
    n_bins = n_bins,
    V = Vsurvey,
    bin_width = logbin
)

cat("JAGS data prepared:\n")
cat("  n_bins =", n_bins, "\n")
cat("  V =", Vsurvey, "Mpc^3\n")
cat("  bin_width =", logbin, "dex\n\n")

############################################################
# Run JAGS
############################################################

cat("Compiling JAGS model...\n")
model = jags.model("mrp_fixed.jags",
                   data=jags_data,
                   n.chains=3,
                   n.adapt=5000)

cat("Burn-in (10000 iterations)...\n")
update(model, 10000)

cat("Sampling (50000 iterations)...\n")
samples = coda.samples(model,
                       c("mstar","log_phi","alpha","beta"),
                       n.iter=50000,
                       thin=25)

posterior = as.matrix(samples)

cat("\n=== POSTERIOR SUMMARY ===\n")
print(summary(samples))

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
cat("Gelman-Rubin (should be < 1.1):\n")
print(gelman.diag(samples))

cat("\nEffective sample size:\n")
print(effectiveSize(samples))

############################################################
# Plot
############################################################

# Recompute full binned data for plotting
hist_full = hist(x_obs, breaks=breaks, plot=FALSE)
phi_emp = hist_full$counts / (Vsurvey * logbin)
x_emp = hist_full$mids

mrp_log10 <- function(mstar, log_phi, alpha, beta, x){
    phi_star = exp(log_phi * log(10))
    phi = beta * log(10) * phi_star *
        10^((alpha + 1) * (x - mstar)) *
        exp(- 10^(beta * (x - mstar)))
    log10(phi)
}

xfit = seq(xmin - 0.5, xmax + 0.5, length.out=600)

# Driver et al. parameters
driver_params = list(
    mstar = 13.51,
    log_phi = -3.19,
    alpha = -1.27,
    beta = 0.47
)

CairoPDF("MRP_FINAL_FIT.pdf", 10, 7)

plot(x_emp, log10(phi_emp),
     pch=19,
     ylim=c(-8, -2),
     xlim=c(12, 16),
     xlab=expression("Halo Mass   log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")   [Mpc"^{-3}*" dex"^{-1}*"]"),
     main=paste0("MRP Fit - FINAL (mult >= ", multi, ", M > 10^", mass_limit, ")"),
     col="darkgreen", cex=1.2)

grid(col="gray90", lty=1)

# Posterior samples
idx = sample(1:nrow(posterior), min(200, nrow(posterior)))
for(i in idx){
    lines(xfit,
          mrp_log10(posterior[i,"mstar"],
                    posterior[i,"log_phi"],
                    posterior[i,"alpha"],
                    posterior[i,"beta"],
                    xfit),
          col=rgb(0,0,1,0.05))
}

# Median fit
med = apply(posterior, 2, median)
lines(xfit,
      mrp_log10(med["mstar"],
                med["log_phi"],
                med["alpha"],
                med["beta"],
                xfit),
      col="red", lwd=3)

# Driver et al. fit
lines(xfit,
      mrp_log10(driver_params$mstar,
                driver_params$log_phi,
                driver_params$alpha,
                driver_params$beta,
                xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("Binned data",
                "Your fit",
                "MCMC samples",
                "Driver+22"),
       col=c("darkgreen", "red", rgb(0,0,1,0.3), "black"),
       pch=c(19, NA, NA, NA),
       lty=c(NA, 1, 1, 2),
       lwd=c(NA, 3, 1, 2),
       cex=0.9,
       bg="white")

med_text = sprintf("M* = %.2f, log(phi*) = %.2f, alpha = %.2f, beta = %.2f",
                   med["mstar"], med["log_phi"], med["alpha"], med["beta"])
mtext(med_text, side=3, line=0.5, cex=0.8, col="red")

dev.off()

cat("\n=== PLOT SAVED ===\n")
cat("File: MRP_FINAL_FIT.pdf\n\n")

############################################################
# Save results
############################################################

results = list(
    median = med,
    quantiles = apply(posterior, 2, quantile, probs=c(0.16, 0.5, 0.84)),
    samples = samples,
    volume = Vsurvey,
    n_groups = N_total,
    mass_range = c(10^xmin, 10^xmax)
)

saveRDS(results, "mrp_final_results.rds")

cat("=== FINAL FITTED PARAMETERS ===\n")
for(par in c("mstar", "log_phi", "alpha", "beta")) {
    q = quantile(posterior[,par], probs=c(0.16, 0.5, 0.84))
    cat(sprintf("%s: %.3f (+%.3f / -%.3f)\n",
                par, q[2], q[3]-q[2], q[2]-q[1]))
}

cat("\n=== COMPARISON WITH DRIVER ET AL. ===\n")
cat("Driver+22 GAMA5:  M*=13.51, log(phi*)=-3.19, alpha=-1.27, beta=0.47\n")
cat("Your fit:         M*=", round(med["mstar"],2),
    ", log(phi*)=", round(med["log_phi"],2),
    ", alpha=", round(med["alpha"],2),
    ", beta=", round(med["beta"],2), "\n", sep="")
