# Refactored Halo Mass Function (HMF) analysis script
# - Split into functions
# - Clear variable names
# - Extensive comments and docstrings
# - Preserves Monte Carlo loop
# - Separate plotting from fitting

###########################
# 0) Libraries
###########################
library(data.table) # fast file IO
library(magicaxis) # enhanced axis and plotting helpers (magplot, magerr, etc.)
library(plotrix) # maghist
library(Cairo) # high-quality PNG output
library(sm) # density / smoothing used in triangle plotting
library(Rfits) # FITS handling if needed
library(MASS) # miscellaneous utilities
library(celestial) # cosmology functions
require(foreign)

###########################
# 1) Configuration / Constants
###########################
# All global constants are provided by setup_constants() to keep namespace clean.
setup_constants <- function() {
  # Cosmology & physical constants
  constants <- list()
  constants$h0 <- 67.37 # Hubble constant [km/s/Mpc]
  constants$omega_matter <- 0.3147 # Omega_m
  constants$omega_matter_err <- 0.0074
  constants$omega_lambda <- 1 - constants$omega_matter
  constants$G <- 6.67408e-11 # gravitational constant [SI]
  constants$msun <- 1.988e30 # solar mass [kg]
  constants$parsec <- 3.0857e16 # metre
  constants$rhocrit <- 3 * (1000 * constants$h0 / (1e6 * constants$parsec))^2 /
    (8 * pi * constants$G)

  # Analysis settings
  constants$fit_bin_width <- 0.001 # resolution for integrals / fits (dex)
  constants$logbin <- 0.2 # used for penalty extension
  constants$alpha_fixed <- -1.864908 # optional fixed alpha for "fix" mode
  constants$RAmid <- c(135, 180, 217, 5) # unused in core func but kept for provenance

  # Survey selection parameters for volume calculation
  constants$zlimit_gama <- 0.25
  constants$zlimit_sdss <- 0.08

  # Survey footprint areas (deg^2) used to compute volumes
  constants$area_gama <- 175.0
  constants$area_sdss <- 7221.0
  constants$volume_reflexii <- 13000000.0 # used in penalty (adopted constant from original script)

  return(constants)
}

###########################
# 2) Utility functions
###########################
# cosmic variance approximation from original script
cos_var_frac <- function(volume_galaxy, n) {
  # returns fractional cosmic variance estimate (not used exactly for each dataset but kept)
  ((219.7 - 52.4 * log10(volume_galaxy) + 3.21 * (log10(volume_galaxy))^2) / n^0.5) / 100.0
}

# Helper that returns a small sequence for extending integrals
extend_mass_tail <- function(log_mass_vector, logbin) {
  max_mass <- max(log_mass_vector)
  seq(max_mass + seq(1, 10, 1) * logbin)
}

###########################
# 3) Data loading & standardization
###########################
#' Load and standardize the observational datasets
#'
#' The function expects files at the specified paths. It returns a list with
#' standardized data.tables: gama, sdss, reflex, reflex2, tpigg, elmo
#'
load_survey_data <- function(base_path = "/Users/sdriver/Drpbx/active/hmf/") {
  d <- list()

  # READ: REFLEX (two files used in original script)
  reflex_path <- file.path(base_path, "reflex.csv")
  reflex2_path <- file.path(base_path, "reflex2.csv")
  if (file.exists(reflex_path)) d$reflex <- fread(reflex_path)
  if (file.exists(reflex2_path)) d$reflex2 <- fread(reflex2_path)

  # 2PIGG
  tpigg_path <- file.path(base_path, "tpigg.dat")
  if (file.exists(tpigg_path)) d$tpigg <- fread(tpigg_path, sep = " ", data.table = FALSE)

  # GAMA
  gama_path <- file.path(base_path, "gamahmfGAMA5.csv")
  if (file.exists(gama_path)) d$gama <- fread(gama_path, sep = ",", data.table = FALSE)

  # SDSS
  sdss_path <- file.path(base_path, "sdsshmf5.csv")
  if (file.exists(sdss_path)) d$sdss <- fread(sdss_path, sep = ",", data.table = FALSE)

  # Elmo/Tempel (SDSS DR10 polygon shading)
  elmo_path <- file.path(base_path, "elmo.csv")
  if (file.exists(elmo_path)) d$elmo <- fread(elmo_path)

  # Standardize and correct cosmology where necessary.
  # NOTE: The corrections replicate the original script's scaling to the chosen H0.
  const <- setup_constants()
  h0 <- const$h0

  if (!is.null(d$reflex2)) {
    d$reflex2$log_mass <- log10(1e14 * d$reflex2$mass) + log10(70 / h0)
    d$reflex2$log_density <- log10(d$reflex2$density) + 3.0 * log10(h0 / 70)
    -14.0 + d$reflex2$log_mass + 1
    d$reflex2$log_min <- log10(d$reflex2$min) + 3.0 * log10(h0 / 70) - 14.0 + d$reflex2$log_mass + 1
    d$reflex2$log_max <- log10(d$reflex2$max) + 3.0 * log10(h0 / 70) - 14.0 + d$reflex2$log_mass + 1
  }

  if (!is.null(d$reflex)) {
    # In original script reflex$x is the log-mass center
    d$reflex$log_mass <- d$reflex$x + log10(70 / h0)
    d$reflex$log_density <- d$reflex$Curve1 + 4.0 * log10(h0 / 70) - 14.0 + d$reflex$log_mass + 1
  }

  if (!is.null(d$tpigg)) {
    d$tpigg$log_mass <- d$tpigg$V1 + log10(100 / h0)
    d$tpigg$log_density <- d$tpigg$V2 + log10((h0 / 100)^3)
    d$tpigg$err_up <- d$tpigg$V3
    d$tpigg$err_down <- d$tpigg$V4
  }

  if (!is.null(d$elmo)) {
    # elmo had its own transformation in original script
    d$elmo$V1 <- d$elmo$V1 + 10.0 + log10(100 / h0)
    d$elmo$V2 <- d$elmo$V2 * (h0 / 100)^3
    d$elmo$V3 <- d$elmo$V3 * (h0 / 100)^3
    d$elmo$V4 <- d$elmo$V4 * (h0 / 100)^3
  }

  if (!is.null(d$sdss)) {
    # apply minimal filtering like the original script
    d$sdss <- d$sdss[d$sdss$V1 > 12.9 & !is.infinite(d$sdss$V4), ]
  }

  if (!is.null(d$gama)) {
    d$gama <- d$gama[d$gama$V1 > 12.7 & !is.infinite(d$gama$V4), ]
  }

  return(d)
}

###########################
# 4) Mass function definition and chi^2 functions
###########################
#' Compute the log10 of the model mass function at given log-mass points
#'
#' @param log_mass numeric vector (dex)
#' @param params numeric vector: c(m_star, phi_star, alpha_slope, beta_cutoff)
#' @return numeric vector of log10(phi)
mass_function_log10 <- function(log_mass, params) {
  m_star <- params[1]
  phi_star <- params[2]
  alpha_slope <- params[3]
  beta_cutoff <- params[4]

  # phi(M) = beta * ln(10) * exp(-10^(beta*(M-M*))) * phi* * (10^M / 10^M*)^(alpha+1)
  model_linear <- beta_cutoff * log(10) * (exp(-10^(beta_cutoff * (log_mass - m_star)))) *
    (phi_star * ((10^log_mass / 10^m_star)^(alpha_slope + 1)))
  return(log10(model_linear))
}

# generic chi-square used by optim()
chi2_mass_function <- function(params, observed_log_mass, observed_log_density,
                               observed_frac_error, constants) {
  # params: (mstar, phi, alpha, beta)
  # Build penalty term for extrapolated tail as original script
  logbin <- constants$logbin
  ext_mass <- extend_mass_tail(observed_log_mass, logbin)

  # penalty uses the linear space model (not log) integrated over the tail
  penalty_term <- 0
  for (vol in c(constants$volume_gama, constants$area_sdss, constants$volume_reflexii)) {
    # This loop replicates original script logic but is conservative if values missing
    # If the volumes aren't defined exactly, skip contributions of that term.
    if (!is.null(vol)) {
      model_tail_linear <- log(10) *
        params[4] *
        (exp(-10^(params[4] * (ext_mass - params[1])))) *
        (params[2] * ((10^ext_mass / 10^params[1])^(params[3] + 1)))
      penalty_term <- penalty_term + 2 * vol * sum(model_tail_linear) * logbin
    }
  }

  model_log10 <- mass_function_log10(observed_log_mass, params)
  # observed_frac_error is fraction (like sigma/log10 scale?) -- replicate original scaling
  # original used ((ally - model) / (allf / log(10)))^2
  denom <- (observed_frac_error / log(10))
  denom[denom == 0] <- 1e-9
  chisq <- sum(((observed_log_density - model_log10) / denom)^2)

  return(chisq + penalty_term)
}

# fixed-alpha variant
chi2_mass_function_fixed_alpha <- function(params, observed_log_mass, observed_log_density,
                                           observed_frac_error, constants) {
  # params: (mstar, phi, beta)
  full_params <- c(params[1], params[2], constants$alpha_fixed, params[3])
  return(chi2_mass_function(
    full_params, observed_log_mass, observed_log_density,
    observed_frac_error, constants
  ))
}

# omegam-constrained variant: add omegam penalty
chi2_mass_function_with_omegam <- function(params, observed_log_mass, observed_log_density,
                                           observed_frac_error, constants) {
  # params: (mstar, phi, alpha, beta)
  xgrid <- rev(seq(0.0, 18, constants$fit_bin_width))
  model_linear <- params[4] * log(10) * ((10^xgrid / 10^params[1])^(params[3] + 1)) *
    (exp(-10^(params[4] * (xgrid - params[1])))) * params[2]
  # integrate to get implied Omega_m from mass function
  implied_omegam <- (sum((model_linear * 10^xgrid) * constants$fit_bin_width) *
                       constants$msun / (1e6 * constants$parsec)^3) / constants$rhocrit
  omegam_penalty <- ((implied_omegam - constants$omega_matter)^2) / (constants$omega_matter_err^2)

  return(omegam_penalty + chi2_mass_function(params, observed_log_mass,
                                             observed_log_density, observed_frac_error, constants))
}

###########################
# 5) Monte Carlo fitter
###########################
#' Run the Monte Carlo fitting loop
#'
#' @param data_list list of standardized datasets
#' @param constants list from setup_constants()
#' @param option what combo of datasets to fit: "GSR", "GS", "GR", "SR", "R", "G", "FIX", "Omega"
#' @param n_iter number of Monte Carlo iterations
#' @return a list with stored parameter chains and derived omegam arrays
run_monte_carlo_fit <- function(data_list, constants, option = "Omega", n_iter = 1001) {
  # Pre-allocate vectors to store chain results
  mstar_chain <- rep(NA_real_, n_iter)
  phi_chain <- rep(NA_real_, n_iter)
  alpha_chain <- rep(NA_real_, n_iter)
  beta_chain <- rep(NA_real_, n_iter)
  omegam_chain <- rep(NA_real_, n_iter)
  omegam_chain_restricted <- rep(NA_real_, n_iter)
  mass_chain <- rep(NA_real_, n_iter)

  # Monte Carlo loop
  for (i in seq_len(n_iter)) {
    # 1) Build perturbed datasets (simulate cosmic variance & measurement noise)
    #    The original script used normal perturbations scaled by cosmic variance estimates.
    cv_reflex <- if (!is.null(data_list$reflex)) data_list$reflex$log_density *
      rnorm(nrow(data_list$reflex), 0.0, 0.05) else NULL
    cv_sdss <- if (!is.null(data_list$sdss)) data_list$sdss$V4 *
      rnorm(nrow(data_list$sdss), 0.0, cos_var_frac(constants$area_sdss, 1)) else NULL
    cv_gama <- if (!is.null(data_list$gama)) data_list$gama$V4 *
      rnorm(nrow(data_list$gama), 0.0, cos_var_frac(constants$area_gama / 3, 3)) else NULL

    # Merge datasets according to option
    if (option %in% c("GSR", "FIX", "Omega")) {
      all_log_mass <- c(data_list$gama$V1, data_list$sdss$V1, data_list$reflex$log_mass)
      combined_density <- c(data_list$gama$V4 + cv_gama, data_list$sdss$V4 +
                              cv_sdss, data_list$reflex$log_density + cv_reflex)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- c(data_list$gama$V8, data_list$sdss$V8,
                               rep(1 / sqrt(20), length(data_list$reflex$log_mass)))
    } else if (option == "GS") {
      all_log_mass <- c(data_list$gama$V1, data_list$sdss$V1)
      combined_density <- c(data_list$gama$V4 + cv_gama, data_list$sdss$V4 + cv_sdss)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- c(data_list$gama$V8, data_list$sdss$V8)
    } else if (option == "GR") {
      all_log_mass <- c(data_list$gama$V1, data_list$reflex$log_mass)
      combined_density <- c(data_list$gama$V4 + cv_gama, data_list$reflex$log_density + cv_reflex)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- c(data_list$gama$V8,
                               rep(1 / sqrt(20), length(data_list$reflex$log_mass)))
    } else if (option == "SR") {
      all_log_mass <- c(data_list$sdss$V1, data_list$reflex$log_mass)
      combined_density <- c(data_list$sdss$V4 + cv_sdss, data_list$reflex$log_density + cv_reflex)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- c(data_list$sdss$V8,
                               rep(1 / sqrt(20), length(data_list$reflex$log_mass)))
    } else if (option == "R") {
      all_log_mass <- c(data_list$reflex$log_mass)
      combined_density <- c(data_list$reflex$log_density + cv_reflex)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- rep(1 / sqrt(20), length(combined_density))
    } else if (option == "G") {
      all_log_mass <- c(data_list$gama$V1)
      combined_density <- c(data_list$gama$V4 + cv_gama)
      combined_density_linear <- 10^combined_density
      combined_frac_error <- c(data_list$gama$V8)
    } else {
      stop("Unknown option passed to run_monte_carlo_fit()")
    }

    # Clean NA / invalids like original script
    valid_idx <- which(!is.na(combined_density_linear) & combined_density_linear > 0)
    observed_log_mass <- all_log_mass[valid_idx]
    observed_log_density <- log10(combined_density_linear[valid_idx])
    observed_frac_error <- combined_frac_error[valid_idx]
    observed_frac_error[is.na(observed_frac_error)] <- 0.9999
    observed_frac_error[is.infinite(observed_frac_error)] <- 0.0
    observed_frac_error <- ifelse(observed_frac_error >= 1.0, 0.9999, observed_frac_error)

    # Choose chi2 function based on option
    if (option == "FIX") {
      # initial guess uses reference MRP parameters adapted from original script
      init_par <- c(14.42947, 1e-20, 0.7097976)
      fit <- optim(par = init_par,
                   fn = chi2_mass_function_fixed_alpha,
                   observed_log_mass = observed_log_mass,
                   observed_log_density = observed_log_density,
                   observed_frac_error = observed_frac_error,
                   constants = constants,
                   control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 0.1)))

      fit_full <- c(fit$par[1], fit$par[2], constants$alpha_fixed, fit$par[3])
    } else if (option == "Omega") {
      init_par <- c(14.42947, 1e-20, -1.864908, 0.7097976)
      fit <- optim(par = init_par,
                   fn = chi2_mass_function_with_omegam,
                   observed_log_mass = observed_log_mass,
                   observed_log_density = observed_log_density,
                   observed_frac_error = observed_frac_error,
                   constants = constants,
                   control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
      fit_full <- fit$par
    } else {
      init_par <- c(14.42947, 1e-20, -1.864908, 0.7097976)
      fit <- optim(par = init_par,
                   fn = chi2_mass_function,
                   observed_log_mass = observed_log_mass,
                   observed_log_density = observed_log_density,
                   observed_frac_error = observed_frac_error,
                   constants = constants,
                   control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
      fit_full <- fit$par
    }

    # compute the model on a grid and derived quantities
    x_grid <- seq(0, 18, constants$fit_bin_width)
    y_model_linear <- fit_full[4] * log(10) *
      (exp(-10^(fit_full[4] * (x_grid - fit_full[1])))) *
      (fit_full[2] * ((10^x_grid / 10^fit_full[1])^(fit_full[3] + 1)))

    # store chains
    mstar_chain[i] <- fit_full[1]
    phi_chain[i] <- fit_full[2]
    alpha_chain[i] <- fit_full[3]
    beta_chain[i] <- fit_full[4]

    # compute implied Omega_m by integrating over x_grid
    implied_omega <- sum((y_model_linear * 10^x_grid) * constants$fit_bin_width) *
      constants$msun / (1e6 * constants$parsec)^3 / constants$rhocrit
    omegam_chain[i] <- implied_omega

    # restricted mass integral (original "bestomegamatter2") for x>12.7
    restricted_idx <- which(x_grid > 12.7)

    implied_omega_restricted <- sum(
      (y_model_linear[restricted_idx] * 10^x_grid[restricted_idx]) * constants$fit_bin_width
    ) *
      constants$msun / (1e6 * constants$parsec)^3 / constants$rhocrit

    omegam_chain_restricted[i] <- implied_omega_restricted

    # store mass scalar (same as omegam_chain but preserved name from original script)
    mass_chain[i] <- implied_omega

    # optional: plot faint realization lines for the first 1000 iters to reproduce original visual
    if (i <= 1000) {
      # draw thin translucent lines on a global device if user opened one
      try(
        {
          lines(x_grid,
                log10(y_model_linear),
                col = rgb(100 / 255, 149 / 255, 237 / 255, alpha = 0.01))
        },
        silent = TRUE
      )
    }
  }

  return(list(mstar = mstar_chain,
              phi = phi_chain,
              alpha = alpha_chain,
              beta = beta_chain,
              omegam = omegam_chain,
              omegam_restricted = omegam_chain_restricted,
              mass = mass_chain))
}

###########################
# 6) Plotting functions
###########################
#' Plot combined mass function and data points
#' This mirrors the plotting behavior in the original script but as a function
plt_mass_func_comp <- function(data_list, fit_params, constants, option = "Omega", out_png = NULL) {
  # if an output filename is provided, open PNG device
  if (!is.null(out_png)) {
    png(filename = out_png, width = 20.0, height = 12.0, units = "cm", res = 240)
  }

  par(fig = c(0, 1, 0, 1), mar = c(3, 3.5, 0.25, 0.25), oma = c(0, 0, 0, 0))
  magplot(0, 0, xlim = c(12.75, 16), ylim = c(-8, -2), unlog = "xy",
          xlab = expression("Halo Mass (M"["\u0298"] ~ ")"),
          ylab = expression("log"[10] ~ "(number density) [Mpc"^-3 ~ "dex"^-1 ~ "]"))

  # theoretical reference (MRP) line from original script
  # We do not recompute the MRP here; user can pass it if needed. We simply plot the fit curve.
  x_grid <- seq(0, 18, constants$fit_bin_width)
  y_fit_linear <- fit_params$beta * log(10) *
    (exp(-10^(fit_params$beta * (x_grid - fit_params$mstar)))) *
    (fit_params$phi * ((10^x_grid / 10^fit_params$mstar)^(fit_params$alpha + 1)))
  lines(x_grid, log10(y_fit_linear), col = rgb(100 / 255, 149 / 255, 237 / 255, 1),
        lwd = 2, lty = 3)

  # plot observational datasets (if available)
  if (!is.null(data_list$reflex)) {
    points(data_list$reflex$log_mass, data_list$reflex$log_density,
           pch = 18, col = "forestgreen", cex = 1.0)
    # magerr expects log-space errors; approximate using original reflexf
    reflexf <- rep(1 / sqrt(20), nrow(data_list$reflex))
    magerr(data_list$reflex$log_mass, data_list$reflex$log_density, yhi = log10(1 + reflexf),
           ylo = log10(1 - reflexf), col = "forestgreen")
  }

  if (!is.null(data_list$tpigg)) {
    points(data_list$tpigg$log_mass, data_list$tpigg$log_density,
           pch = 16, col = "grey50", cex = 1.0)
    magerr(data_list$tpigg$log_mass, data_list$tpigg$log_density, yhi = data_list$tpigg$err_up,
           ylo = -data_list$tpigg$err_down, col = "grey50")
  }

  if (!is.null(data_list$sdss)) {
    points(data_list$sdss$V1, data_list$sdss$V4, pch = 10, col = "purple", cex = 1.0)
    magerr(data_list$sdss$V1, data_list$sdss$V4, yhi = log10(1 + data_list$sdss$V8),
           ylo = log10(1 - data_list$sdss$V8), col = "purple")
  }

  if (!is.null(data_list$gama)) {
    points(data_list$gama$V1, data_list$gama$V4, pch = 16, col = "red", cex = 1.0)
    magerr(data_list$gama$V1, data_list$gama$V4, yhi = log10(1 + data_list$gama$V8),
           ylo = log10(1 - data_list$gama$V8), col = "red")
  }

  # add legend and annotations similar to original
  points(12.75, -5.6, pch = 16, col = "red", cex = 1)
  text(12.75, -5.6, paste0(" GAMA z<", constants$zlimit_gama, " and N>4"),
       col = "red", pos = 4, cex = 1.0)
  points(12.75, -6.0, pch = 10, col = "purple", cex = 1)
  text(12.75, -6.0, paste0(" SDSS DR12, z<", constants$zlimit_sdss, " and N>4"),
       col = "purple", pos = 4, cex = 1.0)
  points(12.75, -6.4, pch = 18, col = "forestgreen", cex = 1)
  text(12.75, -6.4, " REFLEX II, x-ray, z~0.1 (Bohringer et al. 2017)",
       col = "forestgreen", pos = 4, cex = 1.0)
  points(12.75, -6.8, pch = 16, col = "grey50", cex = 1)
  text(12.75, -6.8, " 2PIGG z < 0.12 (Eke et al. 2008)", col = "grey50", pos = 4, cex = 1.0)

  # small inset: histogram of omega_m from Monte Carlo
  par(fig = c(0.6, 0.98, 0.6, 0.98), new = TRUE)
  maghist(fit_params$mass_chain, breaks = seq(-10, 50, 0.01), xlab = expression(Omega[M]),
          grid = FALSE, ylab = "Frequency", majorn = c(3, 2), xlim = c(0.0, 1),
          col = rgb(100 / 255, 149 / 255, 237 / 255, 1.0), verbose = FALSE)
  polygon(c(quantile(fit_params$mass_chain, 0.84), quantile(fit_params$mass_chain, 0.84),
            quantile(fit_params$mass_chain, 0.16), quantile(fit_params$mass_chain, 0.16)),
          c(0, 1e5, 1e5, 0), border = NA, col = rgb(1, 0, 0, 0.25))
  abline(v = constants$omega_matter, col = "black", lwd = 3, lty = 3)
  abline(v = quantile(fit_params$mass_chain, 0.50), col = "red", lwd = 1, lty = 1)

  if (!is.null(out_png)) dev.off()
}

###########################
# 7) Triangle / corner plotting (refactored mymagtri)
###########################
#' Simplified triangle plot for parameter chains
#' This is a cleaned and reduced version of the original mymagtri().
plot_triangle <- function(chain_matrix, param_labels = c("m_star", "log10_phi", "alpha", "beta"),
                          save_to = NULL) {
  # chain_matrix: n x p matrix where p >= 2
  if (!is.null(save_to)) png(save_to, width = 16, height = 16, units = "cm", res = 240)

  par(mfrow = c(ncol(chain_matrix), ncol(chain_matrix)), oma = c(4.1, 4.1, 1.1, 1.1))
  # use lower triangle for scatter, diag for density
  for (i in seq_len(ncol(chain_matrix))) {
    for (j in seq_len(ncol(chain_matrix))) {
      if (i == j) {
        # 1D density
        vals <- chain_matrix[, i]
        if (sd(vals, na.rm = TRUE) == 0) vals <- vals + rnorm(length(vals), sd = 1e-6)
        plot(density(vals, na.rm = TRUE), main = "", xlab = param_labels[i], ylab = "")
        abline(v = mean(vals, na.rm = TRUE), col = "red")
      } else if (i > j) {
        # filled density / contour (approx using magcon from original)
        plot(chain_matrix[, j], chain_matrix[, i], pch = ".", cex = 0.5, xlab = param_labels[j],
             ylab = param_labels[i])
        points(mean(chain_matrix[, j], na.rm = TRUE), mean(chain_matrix[, i], na.rm = TRUE),
               col = "red", pch = 4)
      } else {
        # upper tri: simple scatter with low alpha
        plot(chain_matrix[, j], chain_matrix[, i], pch = ".", cex = 0.5)
      }
    }
  }

  if (!is.null(save_to)) dev.off()
}

###########################
# 8) Main driver
###########################
main <- function(base_path = "/Users/sdriver/Drpbx/active/hmf/", option = "Omega", n_iter = 1001,
                 out_prefix = "/Users/sdriver/Drpbx/active/hmf/") {
  # 1) load constants and data
  constants <- setup_constants()
  data_list <- load_survey_data(base_path)

  # Inject derived constants used elsewhere
  constants$volume_gama <- constants$area_gama / (360^2 / pi) * 1e9 *
    cosdist(constants$zlimit_gama, OmegaM = constants$omega_matter,
            OmegaL = constants$omega_lambda, H0 = constants$h0)$CoVol
  constants$volume_sdss <- constants$area_sdss / (360^2 / pi) * 1e9 *
    cosdist(constants$zlimit_sdss, OmegaM = constants$omega_matter,
            OmegaL = constants$omega_lambda, H0 = constants$h0)$CoVol

  # 2) Run Monte Carlo fitting (this will also attempt to plot faint realization lines)
  cat("Running Monte Carlo fit with", n_iter, "iterations and option =", option, "\n")

  # Open the main figure device to accumulate thin realization lines as in the original script
  png_main <- file.path(out_prefix, paste0("allhmf", option, ".png"))
  png(filename = png_main, width = 20.0, height = 12.0, units = "cm", res = 240)
  par(fig = c(0, 1, 0, 1), mar = c(3, 3.5, 0.25, 0.25), oma = c(0, 0, 0, 0))
  magplot(0, 0, xlim = c(12.75, 16), ylim = c(-8, -2), unlog = "xy",
          xlab = expression("Halo Mass (M"["\u0298"] ~ ")"),
          ylab = expression("log"[10] ~ "(number density) [Mpc"^-3 ~ "dex"^-1 ~ "]"))

  fit_chains <- run_monte_carlo_fit(data_list, constants, option = option, n_iter = n_iter)

  # Close the device so that the faint lines are saved
  dev.off()

  # 3) Fit to the full (unperturbed) dataset for best-fit curve
  # Build final combined dataset without perturbations
  if (option %in% c("GSR", "FIX", "Omega")) {
    all_log_mass <- c(data_list$gama$V1, data_list$sdss$V1, data_list$reflex$log_mass)
    combined_density <- c(data_list$gama$V4, data_list$sdss$V4, data_list$reflex$log_density)
    combined_density_linear <- 10^combined_density
    combined_frac_error <- c(data_list$gama$V8, data_list$sdss$V8,
                             rep(1 / sqrt(20), length(data_list$reflex$log_mass)))
  } else {
    stop("Currently only GSR/FIX/Omega final-fit implemented in main()")
  }

  # Remove NA
  valid_idx <- which(!is.na(combined_density_linear) & combined_density_linear > 0)
  observed_log_mass <- all_log_mass[valid_idx]
  observed_log_density <- log10(combined_density_linear[valid_idx])
  observed_frac_error <- combined_frac_error[valid_idx]

  # final fit
  if (option == "FIX") {
    all_fit <- optim(par = c(14.42947, 1e-20, 0.7097976),
                     fn = chi2_mass_function_fixed_alpha,
                     observed_log_mass = observed_log_mass,
                     observed_log_density = observed_log_density,
                     observed_frac_error = observed_frac_error,
                     constants = constants,
                     control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 0.1)))
    all_fit_full <- c(all_fit$par[1], all_fit$par[2], constants$alpha_fixed, all_fit$par[3])
  } else if (option == "Omega") {
    all_fit <- optim(par = c(14.42947, 1e-20, -1.864908, 0.7097976),
                     fn = chi2_mass_function_with_omegam,
                     observed_log_mass = observed_log_mass,
                     observed_log_density = observed_log_density,
                     observed_frac_error = observed_frac_error,
                     constants = constants,
                     control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
    all_fit_full <- all_fit$par
  } else {
    all_fit <- optim(par = c(14.42947, 1e-20, -1.864908, 0.7097976),
                     fn = chi2_mass_function, observed_log_mass = observed_log_mass,
                     observed_log_density = observed_log_density,
                     observed_frac_error = observed_frac_error,
                     constants = constants,
                     control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
    all_fit_full <- all_fit$par
  }

  # Build a compact list of fit results for downstream plotting and triangle plotting
  fit_result_compact <- list(mstar = all_fit_full[1],
                             phi = all_fit_full[2],
                             alpha = all_fit_full[3],
                             beta = all_fit_full[4],
                             mass_chain = fit_chains$mass)

  # 4) Produce the comparison plot (overwrites the earlier PNG to include best fit and points)
  plt_mass_func_comp(data_list,
    list(mstar = fit_result_compact$mstar,
         phi = fit_result_compact$phi,
         alpha = fit_result_compact$alpha,
         beta = fit_result_compact$beta,
         mass_chain = fit_result_compact$mass_chain),
    constants,
    option = option,
    out_png = file.path(out_prefix,
                        paste0("allhmf_", option, "_final.png"))
  )

  # 5) Triangle plot of Monte Carlo parameter chains
  chain_matrix <- cbind(fit_chains$mstar, log10(fit_chains$phi), fit_chains$alpha, fit_chains$beta)
  colnames(chain_matrix) <- c("m_star", "log10_phi", "alpha", "beta")
  plot_triangle(chain_matrix,
                param_labels = colnames(chain_matrix),
                save_to = file.path(out_prefix, paste0("allcov_", option, ".png")))

  # 6) Print summary results
  cat("Final best-fit parameters (", option, ")\n")
  cat("m_star:", fit_result_compact$mstar, "\n")
  cat("log10(phi):", log10(fit_result_compact$phi), "\n")
  cat("alpha:", fit_result_compact$alpha, "\n")
  cat("beta:", fit_result_compact$beta, "\n")
  cat("Omega_m (median from MC):", signif(median(fit_chains$mass, na.rm = TRUE), 3), "\n")

  invisible(list(constants = constants,
                 data = data_list,
                 fit_chains = fit_chains,
                 final_fit = fit_result_compact))
}

# If this file is sourced, run main() with default arguments
main()

# End of script
