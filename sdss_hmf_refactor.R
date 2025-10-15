# ============================================================
#  Halo Mass Function (HMF) – SDSS DR10 Group Catalogue
#  Cleaned and annotated version
# ============================================================

# ---- Load Required Libraries ----
library(celestial)   # cosmological distance & coordinate utilities
library(ProFound)    # image segmentation / photometry (not heavily used here)
library(Rfits)       # FITS table I/O
library(magicaxis)   # "pretty" astronomical axes for plots
library(data.table)  # fast table reading and manipulation
library(plotrix)     # extra plotting helpers
library(foreign)     # legacy I/O compatibility
library(MASS)        # stats and optimization

# ============================================================
# 1.  MCMC / Corner-Plot Utility
# ============================================================

mymagtri <- function(chains,
                     parameters,
                     best_values,
                     reference_values,
                     filename = NULL,
                     samples = NULL,
                     samptype = "thin",
                     thin = 1) {
  # ---- Basic bookkeeping ----
  n_samples <- dim(chains)[1]
  n_params  <- dim(chains)[2]

  # Choose how to subsample the chain
  if (is.null(samples)) samples <- n_samples
  if (samptype == "thin") {
    idx <- seq(1, n_samples, by = thin)
  } else if (samptype == "random") {
    idx <- sample(seq_len(n_samples), size = samples)
  } else if (samptype == "last") {
    idx <- (n_samples - samples + 1):n_samples
  } else {
    stop("samptype must be 'thin', 'random', or 'last'")
  }

  chains_sub <- chains[idx, , drop = FALSE]

  # ---- Marginal statistics ----
  param_means <- colMeans(chains_sub)
  param_sds   <- apply(chains_sub, 2, sd)

  # ---- Plot layout ----
  if (!is.null(filename)) png(filename, width = 2000, height = 2000, res = 250)
  par(mfrow = c(n_params, n_params),
      mar = c(2, 2, 1.5, 1.5), oma = c(2, 2, 0, 0))

  # ---- Diagonal (1-D histograms) ----
  for (i in seq_len(n_params)) {
    x <- chains_sub[, i]
    hist_x <- hist(x, plot = FALSE, breaks = "FD")
    hist_x$density <- hist_x$counts / max(hist_x$counts)

    plot(hist_x, col = "grey90", border = "grey70",
         main = "", xlab = "", ylab = "",
         axes = FALSE)
    magaxis(1:2, unlog = FALSE, tcl = -0.25)

    # Overlay reference lines
    abline(v = param_means[i], col = "red", lwd = 1.5)
    abline(v = param_means[i] + param_sds[i], col = "red", lty = 2)
    abline(v = param_means[i] - param_sds[i], col = "red", lty = 2)
    abline(v = best_values[i], col = "blue", lwd = 1.5)
    abline(v = reference_values[i], col = "black", lwd = 1.5)
  }

  # ---- Lower-triangle (2-D correlations) ----
  for (i in seq_len(n_params)) {
    for (j in seq_len(n_params)) {
      if (i < j) {
        plot(chains_sub[, j], chains_sub[, i],
             pch = ".", col = "grey30",
             xlab = "", ylab = "", axes = FALSE)
        magaxis(1:2, unlog = FALSE, tcl = -0.25)

        # Overplot key points
        points(reference_values[j], reference_values[i],
               pch = 19, cex = 1.4, col = "black")
        points(best_values[j], best_values[i],
               pch = 19, cex = 1.4, col = "blue")
        points(param_means[j], param_means[i],
               pch = 19, cex = 1.4, col = "red")
      } else if (i > j) {
        # Density map with contours
        magcon(chains_sub[, j], chains_sub[, i],
               xlab = "", ylab = "", axes = FALSE, col = "grey80")
        magaxis(1:2, unlog = FALSE, tcl = -0.25)
      }
    }
  }

  if (!is.null(filename)) dev.off()
}

# ============================================================
# 2.  Halo Mass Function Model / Objective Function
# ============================================================

massfn <- function(x) {
  # Parameters:  M* , Phi , alpha , beta
  m_star <- x[1]
  phi    <- x[2]
  alpha  <- x[3]
  beta   <- x[4]

  # Expected model log10(dn/dlogM)
  model_log <- log10(
    beta * log(10) * exp(-10^(beta * (all_x - m_star))) *
      (phi * ((10^all_x / 10^m_star)^(alpha + 1)))
  )

  # Penalty term discourages sharp kinks or over-normalisation
  penalty <- 2 * volume_sdss *
    sum((diff(model_log, 2) / diff(all_x[1:2])^2)^2, na.rm = TRUE)

  chi2 <- sum(((all_y - model_log) / (all_f / log(10)))^2, na.rm = TRUE) + penalty
  return(chi2)
}

# ============================================================
# 3.  Cosmology and Constants
# ============================================================

# ---- Physical constants ----
h0      <- 67.37
omega_m <- 0.3147
omega_l <- 1 - omega_m
grav_const <- 6.67408e-11
m_sun   <- 1.988e30
parsec  <- 3.0857e16

# Critical density ρ_c = 3H0² / (8πG)
rho_crit <- (3 * (h0 * 1e3 / (parsec * 1e6))^2) / (8 * pi * grav_const) / m_sun * (parsec * 1e6)^3

# ---- Empirical cosmic-variance fit (Driver+2011) ----
cosvar <- function(volume, number_of_objects) {
  # V: survey volume [Mpc³], N: number of objects
  log_v <- log10(volume)
  (10^(1.41 - 0.38 * log_v) + 1.76 /
     sqrt(number_of_objects)) / sqrt(number_of_objects)
}

# ============================================================
# 4.  Input Arguments and Setup
# ============================================================

# The script expects three command-line arguments:
#   1. minimum multiplicity (number of members per group)
#   2. minimum log-mass limit
#   3. minimum redshift
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript hmf_sdss.R <multi> <mlimit> <zmin>")

min_mult  <- as.numeric(args[1])
mass_lim  <- as.numeric(args[2])
z_min     <- as.numeric(args[3])

cat(sprintf("Using multiplicity ≥ %d, logM ≥ %.2f, z ≥ %.3f\n",
            min_mult, mass_lim, z_min))

# ---- Survey area and redshift limits ----
area_deg2 <- 7221
solid_angle <- (area_deg2 / 41253) * 4 * pi          # steradians
z_max <- 0.1

# ---- Comoving-volume function ----
comoving_volume <- function(z) {
  celestial::cosdist(z, h0 = h0, OmegaM = omega_m,
                     OmegaL = omega_l, type = "CoVol")
}

# Volume between z_min and z_max, corrected for survey area
volume_sdss <- solid_angle * (comoving_volume(z_max) - comoving_volume(z_min))
cat(sprintf("Survey comoving volume = %.3e Mpc³\n", volume_sdss))

# ============================================================
# 5.  Load SDSS DR10 Group and Galaxy Catalogues
# ============================================================

cat("Loading SDSS group and galaxy catalogues...\n")
sdss_galaxies <- Rfits_read_table("sdssdr10table1.fits", ext = 2)$imDat
sdss_groups   <- Rfits_read_table("sdssdr10table2.fits", ext = 2)$imDat

# ---- Rename columns for clarity ----
names(sdss_groups)[names(sdss_groups) == "idcl"] <- "group_id"
names(sdss_groups)[names(sdss_groups) == "nrich"] <- "n_rich"
names(sdss_groups)[names(sdss_groups) == "mass"] <- "mass"
names(sdss_groups)[names(sdss_groups) == "zcl"]  <- "z_cl"

names(sdss_galaxies)[names(sdss_galaxies) == "idcl"] <- "group_id"
names(sdss_galaxies)[names(sdss_galaxies) == "rank"] <- "rank"
names(sdss_galaxies)[names(sdss_galaxies) == "redshift"] <- "z"
names(sdss_galaxies)[names(sdss_galaxies) == "absmag_r"] <- "abs_mag_r"

# ============================================================
# 6.  Compute Detection Limits and Volumes
# ============================================================

mag_limit <- 17.77  # SDSS main-sample r-band limit

# Apparent magnitude relation to distance modulus
calc_dmax <- function(abs_mag) {
  10^((mag_limit - abs_mag - 25) / 5) / 1e6  # in Mpc
}

sdss_galaxies$dmax <- calc_dmax(sdss_galaxies$abs_mag_r)
sdss_galaxies$zmax <- celestial::cosdist(sdss_galaxies$dmax, type = "dist2z")

# Merge galaxy- and group-level info
sdss <- merge(sdss_galaxies, sdss_groups, by = "group_id")

# Filter by multiplicity, redshift, and mass limit
sdss <- subset(sdss,
               n_rich >= min_mult & z_cl >= z_min & log10(mass) > mass_lim)

cat(sprintf("%d galaxies across %d groups after cuts\n",
            nrow(sdss), length(unique(sdss$group_id))))

# ---- Compute Vmax (maximum observable volume) per galaxy ----
sdss$vmax <- solid_angle * (comoving_volume(pmin(sdss$zmax, z_max)) -
                              comoving_volume(z_min))
sdss$weight_zlimit <- volume_sdss / pmax(sdss$vmax, 1e-9)

# ============================================================
# 7.  Construct the Halo Mass Function (HMF)
# ============================================================

log_bin_width <- 0.1
mass_bins <- seq(10.3, 16.1, by = log_bin_width)
bin_centres <- mass_bins[-length(mass_bins)] + log_bin_width / 2

# Raw (unweighted) histogram
hmf_raw <- hist(log10(sdss$mass), breaks = mass_bins, plot = FALSE)

# Weighted by 1/Vmax to correct for selection effects
hmf_weighted <- hist(log10(sdss$mass),
                     breaks = mass_bins,
                     weights = 1 / sdss$weight_zlimit,
                     plot = FALSE)

# Convert to number density φ = N / ΔlogM / Volume
phi_raw      <- hmf_raw$counts / diff(mass_bins) / volume_sdss
phi_weighted <- hmf_weighted$counts / diff(mass_bins) / volume_sdss

# Poisson error
phi_err <- sqrt(hmf_raw$counts) / diff(mass_bins) / volume_sdss

# Store globals for fitting function
all_x <- bin_centres
all_y <- log10(phi_weighted)
all_f <- phi_err * volume_sdss  # used as error term in massfn()

# ---- Cosmic-variance estimate ----
cosmic_variance <- cosvar(volume_sdss, hmf_raw$counts)
phi_total_err <- sqrt(phi_err^2 + (phi_weighted * cosmic_variance)^2)

# ---- Diagnostic plot ----
plot(all_x, phi_weighted,
     log = "y", pch = 19, col = "black",
     xlab = expression(log[10](M / M[sun])),
     ylab = expression(phi(M) ~ "[Mpc^{-3}~dex^{-1}]"),
     main = "SDSS Halo Mass Function")
arrows(all_x, phi_weighted - phi_total_err,
       all_x, phi_weighted + phi_total_err,
       code = 3, angle = 90, length = 0.02)

# ============================================================
# Part 3: Monte-Carlo propagation, MRP fit, plotting, and output
# ============================================================

# ---- Derived model (MRP) evaluation helper ----
evaluate_mrp <- function(params, x_vals) {
  # params: numeric vector c(m_star, phi, alpha, beta)
  m_star <- params[1]
  phi    <- params[2]
  alpha  <- params[3]
  beta   <- params[4]

  y_model <- beta * log(10) * exp(-10^(beta * (x_vals - m_star))) *
    (phi * ((10^x_vals / 10^m_star)^(alpha + 1)))
  return(y_model)
}

# ---- Initial literature MRP / model guess (kept from original) ----
beta_mrp_init <- 0.7097976
a_const <- 1.727006e-19
m_star_mrp <- 14.42947
alpha_mrp <- -1.864908

# build a reference MRP curve (mrpx, mrpy) used for plotting comparison
mrp_x <- seq(0, 17, by = 0.001) + log10(100 / H0)
mrp_y <- a_const * beta_mrp_init * 10^((alpha_mrp + 1) * (mrp_x - m_star_mrp)) *
  exp(-10^(beta_mrp_init * (mrp_x - m_star_mrp))) * (H0 / 100)^3

# Normalisation factor such that integrated mass fraction matches Omega_m reference
factor_norm <- sum(10^mrp_x * mrp_y) * 0.001 * M_sun / (1e6 * parsec)^3 / (Omega_m * rho_crit)
phi_mrp <- a_const / factor_norm

# ---- Monte-Carlo propagation of mass errors to number counts ----
n_mc_realizations <- 1001
n_bins <- length(bin_centres)

# Preallocate matrix to store mock counts (each row = realization)
mock_counts <- matrix(0, nrow = n_mc_realizations, ncol = n_bins)

# A quick defensive fix -- ensure per-object mass errors exist
if (!("log10MassErr" %in% names(sdss))) {
  sdss$log10MassErr <- 0.1  # fallback if not defined
}

cat("Running Monte-Carlo resampling for mass uncertainties...\n")
for (i in seq_len(n_mc_realizations)) {
  # perturb the log10(mass) per galaxy by its measurement uncertainty (Gaussian)
  perturbed_logmass <- log10(sdss$mass) + rnorm(nrow(sdss), mean = 0, sd = sdss$log10MassErr)

  # compute weighted histogram (1/Vmax correction) using perturbed masses
  # we mimic original weighted.hist behavior by accumulating counts per bin
  counts <- numeric(n_bins)
  for (j in seq_len(n_bins)) {
    in_bin <- (perturbed_logmass >= mass_bins[j]) & (perturbed_logmass < mass_bins[j + 1])
    if (any(in_bin)) {
      counts[j] <- sum(1 / sdss$weight_zlimit[in_bin])
    } else {
      counts[j] <- 0
    }
  }
  mock_counts[i, ] <- counts
}

# Compute mean / median / 16th / 84th percentiles across realizations
mc_mean_counts  <- apply(mock_counts, 2, mean)
mc_median_counts <- apply(mock_counts, 2, median)
mc_up_counts    <- apply(mock_counts, 2, quantile, probs = 0.84)
mc_down_counts  <- apply(mock_counts, 2, quantile, probs = 0.16)

# Correction factor edb (empirical debiasing factor from MC / observed)
edb <- mc_mean_counts / (hmf_weighted$counts + 1e-12)
edb[is.infinite(edb) | is.na(edb)] <- 1.0

# Monte-Carlo fractional error estimate per bin (use 66% quantile ~ 1σ)
mc_frac_err <- numeric(n_bins)
for (b in seq_len(n_bins)) {
  mc_frac_err[b] <- sqrt(quantile((mc_mean_counts[b] - mock_counts[, b])^2, probs = 0.66)) /
    (mc_mean_counts[b] + 1e-12)
}

# Poisson fractional error
poisson_frac_err <- sqrt(hmf_raw$counts) / pmax(hmf_raw$counts, 1e-12)

# Combine fractional errors in quadrature
frac_err_combined <- sqrt(mc_frac_err^2 + poisson_frac_err^2)
frac_err_combined[is.na(frac_err_combined)] <- 0.9999
frac_err_combined[frac_err_combined >= 1.0] <- 0.9999

# Final number-density and error arrays used for fitting
sdss_x     <- hmf_weighted$mids  # mass bin centres
sdss_y     <- hmf_weighted$counts / (log_bin_width * volume_sdss)  # phi M
sdss_error <- sdss_y * frac_err_combined  # absolute error in phi(M)

# Keep only bins above mass limit and with positive counts
valid_mask <- (!is.na(sdss_y)) & (sdss_y > 0) & (sdss_x > mass_lim)
all_x <- sdss_x[valid_mask]
all_y <- log10(sdss_y[valid_mask])
# all_f: errors used inside massfn(); convert to absolute error on log10 scale
all_f <- sdss_error[valid_mask]

# ============================================================
# Fit the MRP model to the corrected SDSS HMF using optim()
# ============================================================
cat("Running initial MRP fit using optim()...\n")
initial_par <- c(m_star_mrp, phi_mrp, alpha_mrp, beta_mrp_init)
optim_control <- list(maxit = 1000, reltol = 1e-8, parscale = c(1, 1, 1, 0.1))

fit_result <- optim(par = initial_par,
                    fn  = massfn,
                    control = optim_control,
                    method = "Nelder-Mead")

if (fit_result$convergence != 0) {
  warning("optim() did not converge cleanly (convergence != 0). Check diagnostics.")
}
best_par <- fit_result$par
names(best_par) <- c("m_star", "phi", "alpha", "beta")
cat("Best-fit parameters (optim):\n")
print(best_par)

# Evaluate the best-fit MRP on a fine grid for plotting
x_fit <- seq(0, 20, by = 0.01)
y_fit <- evaluate_mrp(best_par, x_fit)

# Compute Omega_M in halos from the fitted MRP (integrated mass fraction)
sdss_omega_matter <- sum((y_fit * 10^x_fit) * 0.01) * M_sun / (1e6 * parsec)^3 / rho_crit
cat(sprintf("Inferred Omega_m (from HMF fit) = %.4f\n", sdss_omega_matter))

# ============================================================
# Monte-Carlo re-fitting to estimate parameter uncertainties
# (Re-fit for many realizations of the HMF using combined errors)
# ============================================================
n_fit_realizations <- 10001

# Preallocate vectors for sampled parameters
m_star_samples  <- numeric(n_fit_realizations)
phi_samples     <- numeric(n_fit_realizations)
alpha_samples   <- numeric(n_fit_realizations)
beta_samples    <- numeric(n_fit_realizations)

cat("Running Monte-Carlo fits (this may take a while)...\n")
for (i in seq_len(n_fit_realizations)) {
  # Perturb the observed phi(M) by a combination of measurement and cosmic variance
  cv_term <- sdss_y * rnorm(length(sdss_y), mean = 0, sd = cosmic_variance)
  perturbed_phi <- sdss_y + sdss_y *
    rnorm(length(sdss_y), mean = 0, sd = frac_err_combined) + cv_term

  # Keep only positive bins and above mass limit
  mask <- (perturbed_phi > 0) & (!is.na(perturbed_phi)) & (sdss_x > mass_lim)
  fit_x <- sdss_x[mask]
  fit_y <- log10(perturbed_phi[mask])
  fit_f <- sdss_error[mask]

  # Move global all_x/all_y/all_f pointers for massfn()
  all_x <- fit_x
  all_y <- fit_y
  all_f <- fit_f

  # run optim on this perturbed dataset (use initial guesses near published MRP)
  mc_fit <- optim(par = initial_par, fn = massfn, control = optim_control, method = "Nelder-Mead")

  m_star_samples[i] <- mc_fit$par[1]
  phi_samples[i]    <- mc_fit$par[2]
  alpha_samples[i]  <- mc_fit$par[3]
  beta_samples[i]   <- mc_fit$par[4]

  # For the first 1000 draws, plot the fit thinly to show spread
  if (i <= 1000) {
    y_monty <- evaluate_mrp(mc_fit$par, x_fit)
    lines(x_fit, log10(y_monty), col = rgb(100 / 255, 149 / 255, 237 / 255, 0.01))
  }
}

# Re-assign the global all_x/all_y/all_f to the nominal (unperturbed) dataset
all_x <- sdss_x[valid_mask]
all_y <- log10(sdss_y[valid_mask])
all_f <- sdss_error[valid_mask]

# ============================================================
# Plot results: data + best-fit + monte-carlo envelope + literature
# ============================================================
png(filename = sprintf("sdss_hmf_multi_%d.png", min_mult), width = 1200, height = 800, res = 150)
par(mar = c(5, 5, 3, 2))

# Plot the corrected HMF (points) with errors
plot(sdss_x, log10(sdss_y),
     pch = 16, col = "purple",
     xlab = expression(log[10](M)),
     ylab = expression(log[10](phi(M))),
     xlim = c(12, 16),
     ylim = c(-8, -2))

# Add error bars from combined fractional errors
arrows(sdss_x, log10(pmax(sdss_y - sdss_error, 1e-20)),
       sdss_x, log10(sdss_y + sdss_error),
       angle = 90, code = 3, length = 0.02, col = "purple")

# Add raw counts (unweighted) as small ticks
points(hmf_raw$mids, log10(hmf_raw$counts / (log_bin_width * volume_sdss)),
       pch = 5, col = "limegreen", cex = 0.8)

# Best-fit line
lines(x_fit, log10(y_fit), col = rgb(100 / 255, 149 / 255, 237 / 255), lwd = 2)

# Add MRP literature curve (mrp_y normalized)
lines(mrp_x - 0.08, log10(mrp_y) - log10(factor_norm) + 0.08, lty = 2, lwd = 2)

# Annotate
legend("topright", legend = c("Corrected HMF", "Raw counts", "Best-fit MRP", "Literature MRP"),
       col = c("purple", "limegreen", rgb(100 / 255, 149 / 255, 237 / 255), "black"),
       pch = c(16, 5, NA, NA), lty = c(NA, NA, 1, 2), bg = "white")

dev.off()

# ============================================================
# Save numeric outputs to disk (CSV + human-readable text)
# ============================================================
cat("Writing HMF tables and parameter samples to disk...\n")
out_table <- data.frame(
  mass_bin_center = sdss_x,
  phi = sdss_y,
  phi_err = sdss_error,
  raw_counts = hmf_raw$counts,
  weighted_counts = hmf_weighted$counts
)
write.csv(out_table, file = sprintf("sdss_hmf_table_multi_%d.csv", min_mult), row.names = FALSE)

# Parameter sample summary and CSV
params_df <- data.frame(m_star = m_star_samples,
                        phi = phi_samples,
                        alpha = alpha_samples,
                        beta = beta_samples)
write.csv(params_df, file = sprintf("sdss_hmf_params_mc_multi_%d.csv", min_mult), row.names = FALSE)

# Summary printed to console
cat("Parameter medians and 16/84 percentiles:\n")
for (name in c("m_star", "phi", "alpha", "beta")) {
  vec <- params_df[[name]]
  cat(sprintf("%s : median = %.4g ; 16%% = %.4g ; 84%% = %.4g\n",
              name, median(vec, na.rm = TRUE),
              quantile(vec, 0.16, na.rm = TRUE),
              quantile(vec, 0.84, na.rm = TRUE)))
}

# ============================================================
# Corner plot (parameter covariance) using mymagtri()
# ============================================================
cat("Creating corner plot of fitted parameters...\n")
# Prepare chain-like matrix: each row is a realization
param_chain <- cbind(m_star_samples, log10(phi_samples), alpha_samples, beta_samples)

# Because mymagtri() expects vectors for best_values and reference_values
best_values <- c(best_par["m_star"], log10(best_par["phi"]), best_par["alpha"], best_par["beta"])
ref_values  <- c(m_star_mrp - 0.075, log10(phi_mrp) + 0.075, alpha_mrp, beta_mrp_init)

# Save a PNG of the corner plot
mymagtri(chains = param_chain,
         parameters = NULL,
         best_values = best_values,
         reference_values = ref_values,
         filename = sprintf("sdss_hmf_corner_multi_%d.png", min_mult),
         samples = nrow(param_chain),
         samptype = "last",
         thin = 1)

cat("Done. All outputs written and plots saved.\n")
