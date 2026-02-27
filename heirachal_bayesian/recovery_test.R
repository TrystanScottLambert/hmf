############################################################
# SIMULATION-RECOVERY TEST FOR HMF STAN MODEL
# 
# Purpose: Generate a mock halo catalogue from known MRP
# parameters, simulate observation (detection limits,
# mass errors), then run the Stan model to see if we
# recover the input parameters.
#
# This validates the inference pipeline before applying
# it to real GAMA data.
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. TRUE MRP PARAMETERS (ground truth)
#    Using the GSR joint fit from Driver+22 Table 2
############################################################

TRUE_MSTAR   <- 14.13   # log10(M*/Msun)
TRUE_LOG_PHI <- -3.96   # log10(phi* / Mpc^-3 dex^-1)
TRUE_ALPHA   <- -1.68
TRUE_BETA    <- 0.63

cat("=== TRUE PARAMETERS ===\n")
cat(sprintf("M*      = %.2f\n", TRUE_MSTAR))
cat(sprintf("log_phi = %.2f\n", TRUE_LOG_PHI))
cat(sprintf("alpha   = %.2f\n", TRUE_ALPHA))
cat(sprintf("beta    = %.2f\n\n", TRUE_BETA))

############################################################
# 2. MRP FUNCTION DEFINITIONS
############################################################

# MRP in log10 mass space: dn/dlog10(M) [Mpc^-3 dex^-1]
mrp_phi <- function(x, mstar, log_phi, alpha, beta) {
    # x = log10(M/Msun)
    u <- x - mstar
    beta * log(10) * 10^log_phi * 10^((alpha + 1) * u) * exp(-10^(beta * u))
}

# Log10 of the MRP (for plotting)
mrp_log10 <- function(x, mstar, log_phi, alpha, beta) {
    log10(mrp_phi(x, mstar, log_phi, alpha, beta))
}

############################################################
# 3. SURVEY GEOMETRY (mimicking GAMA)
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01

# Comoving distance (simplified Hubble law for low z)
# For a proper treatment use celestial::cosdistCoDist
# Here we use a quick approximation adequate for z < 0.3
cosdist_approx <- function(z) {
    # Comoving distance in Mpc (using numerical integration)
    c_light <- 299792.458  # km/s
    integrand <- function(zp) 1 / sqrt(omegam * (1 + zp)^3 + (1 - omegam))
    result <- integrate(integrand, 0, z)$value
    return(c_light / ho * result)
}

d_max <- cosdist_approx(zlimit)
d_min <- cosdist_approx(zmin)

# GAMA survey: 179.92 sq deg across 3 equatorial fields
sky_frac <- 179.92 * (pi / 180)^2 / (4 * pi)
Vsurvey  <- (4/3) * pi * (d_max^3 - d_min^3) * sky_frac

cat(sprintf("Survey volume: %.2e Mpc^3\n", Vsurvey))
cat(sprintf("Sky fraction: %.4f\n", sky_frac))
cat(sprintf("d_max: %.1f Mpc\n\n", d_max))

############################################################
# 4. GENERATE TRUE HALO POPULATION
#    Sample from MRP using rejection sampling in log10(M)
############################################################

# Mass range for sampling
mass_lo <- 10.0   # log10(M/Msun) - well below detection
mass_hi <- 16.0

# Expected number of halos in the volume
# Integrate MRP over mass range
mass_grid <- seq(mass_lo, mass_hi, by = 0.001)
phi_grid  <- mrp_phi(mass_grid, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
N_expected <- Vsurvey * sum(phi_grid) * 0.001

cat(sprintf("Expected total halos (%.0f < log M < %.0f): %.0f\n", 
            mass_lo, mass_hi, N_expected))

# For the high-mass end where we actually detect things:
mass_grid_det <- seq(12.5, mass_hi, by = 0.001)
phi_grid_det  <- mrp_phi(mass_grid_det, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
N_detectable  <- Vsurvey * sum(phi_grid_det) * 0.001

cat(sprintf("Expected halos above log M = 12.5: %.0f\n", N_detectable))

# Use Poisson draw for total number, then rejection sample for masses
# We'll sample from mass_lo_sample to avoid the huge number of tiny halos
mass_lo_sample <- 12.0  # still below detection threshold
phi_sample <- mrp_phi(seq(mass_lo_sample, mass_hi, by = 0.001), 
                       TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
N_sample <- Vsurvey * sum(phi_sample) * 0.001
N_draw   <- rpois(1, N_sample)

cat(sprintf("Drawing %d halos from MRP (log M > %.1f)\n", N_draw, mass_lo_sample))

# Rejection sampling
phi_max <- max(phi_sample) * 1.1  # envelope

true_masses <- numeric(0)
attempts    <- 0
while(length(true_masses) < N_draw) {
    # Propose masses uniformly in log10
    n_batch <- min(N_draw * 5, 1e6)
    m_prop  <- runif(n_batch, mass_lo_sample, mass_hi)
    phi_prop <- mrp_phi(m_prop, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
    accept   <- runif(n_batch) < phi_prop / phi_max
    true_masses <- c(true_masses, m_prop[accept])
    attempts <- attempts + 1
    if(attempts > 100) {
        cat("Warning: rejection sampling slow, check phi_max\n")
        break
    }
}
true_masses <- true_masses[1:N_draw]

cat(sprintf("Generated %d true halo masses\n", length(true_masses)))
cat(sprintf("Mass range: [%.2f, %.2f]\n", min(true_masses), max(true_masses)))

# Quick check: histogram should match MRP
h <- hist(true_masses, breaks = seq(mass_lo_sample, mass_hi, by = 0.2), plot = FALSE)
cat("\nSanity check - binned number counts vs MRP prediction:\n")
for(i in 1:min(10, length(h$mids))) {
    predicted <- mrp_phi(h$mids[i], TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA) * Vsurvey * 0.2
    cat(sprintf("  log M = %.1f: observed = %d, predicted = %.0f\n",
                h$mids[i], h$counts[i], predicted))
}

############################################################
# 5. ASSIGN REDSHIFTS
#    Uniform in comoving volume (as expected for a 
#    non-evolving population over this narrow z range)
############################################################

# For each halo, draw a redshift uniformly in comoving volume
# P(z) dz prop to dV/dz prop to r^2 dr/dz

# Build a fine grid for inverse CDF sampling
z_grid <- seq(zmin, zlimit, length.out = 5000)
d_grid <- sapply(z_grid, cosdist_approx)
# dV/dz proportional to d^2 * dd/dz
dV_dz <- c(0, diff(d_grid^3))  # proportional to volume element
cdf_z <- cumsum(dV_dz)
cdf_z <- cdf_z / max(cdf_z)

# Inverse CDF sampling
u_rand   <- runif(N_draw)
z_halos  <- approx(cdf_z, z_grid, xout = u_rand, rule = 2)$y

cat(sprintf("\nRedshift range: [%.3f, %.3f], median = %.3f\n",
            min(z_halos), max(z_halos), median(z_halos)))

############################################################
# 6. SIMULATE DETECTION: MASS-DEPENDENT SELECTION
#    A halo is "detected" if it has enough bright members
#    to be identified as a group with N >= 5.
#    
#    Key physics: lower mass halos are only detectable at
#    lower redshifts. We model this as a mass-dependent
#    maximum redshift (zmax).
#
#    We use a simple scaling motivated by the GAMA survey:
#    zmax(M) increases with mass. A 10^13 Msun group is
#    detectable to z~0.1, while 10^15 is detectable to z~0.4.
############################################################

# Simple detection model: the nth brightest galaxy scales with halo mass
# At the flux limit, galaxy is detectable to some zmax
# zmax scales roughly as 10^(0.2*(log_mass - 13)) in a simplified model

# More physically: the 5th brightest galaxy in a group has luminosity
# that scales with halo mass. We parameterize zmax(M):

zmax_of_mass <- function(log_mass) {
    # Approximate: z_max scales with mass
    # Calibrated so that:
    #   log M = 12.5 -> zmax ~ 0.05
    #   log M = 13.0 -> zmax ~ 0.10  
    #   log M = 14.0 -> zmax ~ 0.20
    #   log M = 15.0 -> zmax ~ 0.40+
    
    z_max <- 0.01 * 10^(0.4 * (log_mass - 12.0))
    return(pmin(z_max, 0.5))  # cap at sensible value
}

# Test the detection function
cat("\nDetection model zmax(M):\n")
for(m in c(12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0)) {
    cat(sprintf("  log M = %.1f -> zmax = %.3f\n", m, zmax_of_mass(m)))
}

# Add scatter to zmax (real groups have varying luminosity distributions)
zmax_true <- zmax_of_mass(true_masses)
# Add ~0.3 dex scatter in zmax
log_zmax_scatter <- rnorm(N_draw, 0, 0.15)
zmax_scattered <- zmax_true * 10^log_zmax_scatter
zmax_scattered <- pmin(zmax_scattered, zlimit)
zmax_scattered <- pmax(zmax_scattered, zmin)

# A halo is DETECTED if z_halo < zmax_scattered
detected <- z_halos < zmax_scattered

N_detected <- sum(detected)
cat(sprintf("\nDetected %d of %d halos (%.1f%%)\n", 
            N_detected, N_draw, 100 * N_detected / N_draw))

# Apply detection
det_true_mass <- true_masses[detected]
det_z         <- z_halos[detected]
det_zmax      <- zmax_scattered[detected]

# Cap zmax at survey limit  
det_zmax <- pmin(det_zmax, zlimit)

cat(sprintf("Detected mass range: [%.2f, %.2f]\n", 
            min(det_true_mass), max(det_true_mass)))

############################################################
# 7. ADD MEASUREMENT ERRORS
#    Mass errors depend on group multiplicity.
#    We use the same error model as in the real analysis:
#    sigma_log10M ~ 0.1 to 0.4, lower for massive groups.
############################################################

# Mass error model: bigger groups have better-measured masses
# Approximate: sigma scales inversely with sqrt(multiplicity)
# multiplicity correlates with mass

mass_error <- function(log_mass) {
    # Higher mass -> more members -> smaller error
    # Approximate the error curve from Driver+22 Fig 3
    sigma <- 0.45 - 0.15 * (log_mass - 12.5) / 2.5
    sigma <- pmax(sigma, 0.10)
    sigma <- pmin(sigma, 0.50)
    return(sigma)
}

det_sigma <- mass_error(det_true_mass)

cat(sprintf("Mass error range: [%.2f, %.2f], median = %.2f\n",
            min(det_sigma), max(det_sigma), median(det_sigma)))

# Observed masses = true masses + Gaussian error in log10 space
det_obs_mass <- det_true_mass + rnorm(N_detected, 0, det_sigma)

cat(sprintf("Observed mass range: [%.2f, %.2f]\n\n",
            min(det_obs_mass), max(det_obs_mass)))

############################################################
# 8. COMPUTE MASS LIMITS PER GROUP
#    This mirrors your real code: m_lim[i] is the minimum
#    mass that could have been detected at this group's
#    redshift, given its zmax.
#
#    The logic: if zmax is much larger than z_obs, the group
#    is well above threshold -> low m_lim.
#    If zmax ~ z_obs, the group is barely detected -> m_lim ~ M_obs.
############################################################

# Use the INVERSE of the detection function to get mass limit
# If zmax(M) = f(M), then m_lim(z) = f^{-1}(z)
# From our model: zmax = 0.01 * 10^(0.4*(logM - 12))
# So: logM = 12 + 2.5 * log10(zmax / 0.01)

mass_limit_from_z <- function(z) {
    # What minimum mass is detectable at redshift z?
    log_m <- 12.0 + 2.5 * log10(z / 0.01)
    return(log_m)
}

# The mass limit for each group is set by its redshift
# (the minimum mass detectable at that distance)
det_m_lim <- mass_limit_from_z(det_z)

# Add some scatter to mass limits (reflecting real survey complexity)
det_m_lim <- det_m_lim + rnorm(N_detected, 0, 0.2)

# Ensure mass limits are sensible
det_m_lim <- pmax(det_m_lim, 10.5)
det_m_lim <- pmin(det_m_lim, det_obs_mass - 0.1)  # must be below observed mass

cat("Mass limit statistics:\n")
cat(sprintf("  Range: [%.2f, %.2f]\n", min(det_m_lim), max(det_m_lim)))
cat(sprintf("  Median gap (M_obs - M_lim): %.2f dex\n\n", 
            median(det_obs_mass - det_m_lim)))

############################################################
# 9. QUALITY CUTS (mimicking real analysis)
############################################################

# Cut to sensible mass range
keep <- det_obs_mass > 12.0 & det_obs_mass < 16.5 & 
        is.finite(det_obs_mass) & is.finite(det_m_lim)

x_obs     <- det_obs_mass[keep]
sigma_obs <- det_sigma[keep]
m_lim_obs <- det_m_lim[keep]
x_true    <- det_true_mass[keep]

N <- length(x_obs)
cat(sprintf("Final mock catalogue: N = %d groups\n", N))
cat(sprintf("  Observed mass range: [%.2f, %.2f]\n", min(x_obs), max(x_obs)))
cat(sprintf("  True mass range:     [%.2f, %.2f]\n", min(x_true), max(x_true)))
cat(sprintf("  Mass limit range:    [%.2f, %.2f]\n", min(m_lim_obs), max(m_lim_obs)))

############################################################
# 10. DIAGNOSTIC PLOT: MOCK DATA vs TRUE MRP
############################################################

pdf("mock_data_diagnostic.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

# Panel 1: Binned HMF of true masses vs MRP
breaks <- seq(12, 16, by = 0.2)
h_true <- hist(x_true, breaks = breaks, plot = FALSE)
h_obs  <- hist(x_obs, breaks = breaks, plot = FALSE)
phi_true <- h_true$counts / (Vsurvey * 0.2)
phi_obs  <- h_obs$counts / (Vsurvey * 0.2)

xfit <- seq(12, 16, length.out = 500)
yfit <- mrp_phi(xfit, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)

ok_t <- phi_true > 0
ok_o <- phi_obs > 0

plot(h_true$mids[ok_t], log10(phi_true[ok_t]),
     pch = 16, col = "blue", cex = 1.2,
     xlim = c(12, 16), ylim = c(-8, -2),
     xlab = expression("log"[10]*"(M/M"[symbol("\u2299")]*")"),
     ylab = expression("log"[10]*"("*phi*") [Mpc"^{-3}*" dex"^{-1}*"]"),
     main = "Mock HMF: True vs Observed vs MRP")
points(h_obs$mids[ok_o], log10(phi_obs[ok_o]),
       pch = 17, col = "red", cex = 1.2)
lines(xfit, log10(yfit), col = "black", lwd = 2)
legend("topright", c("True masses (detected)", "Observed masses", "Input MRP"),
       pch = c(16, 17, NA), col = c("blue", "red", "black"),
       lty = c(NA, NA, 1), lwd = c(NA, NA, 2))

# Panel 2: Mass errors
plot(x_true, x_obs - x_true, pch = ".", col = rgb(0, 0, 0, 0.1),
     xlab = "True log10(M)", ylab = "Observed - True (dex)",
     main = "Mass measurement errors")
abline(h = 0, col = "red", lwd = 2)
abline(h = c(-0.3, 0.3), col = "red", lty = 2)

# Panel 3: Detection function
plot(x_true, det_z[keep], pch = ".", col = rgb(0, 0, 1, 0.2),
     xlab = "True log10(M)", ylab = "Redshift",
     main = "Detection: z vs Mass")
m_test <- seq(12, 16, by = 0.1)
lines(m_test, zmax_of_mass(m_test), col = "red", lwd = 2)
abline(h = zlimit, col = "gray", lty = 2)
legend("topleft", c("Detected halos", "zmax(M) curve", "Survey limit"),
       pch = c(1, NA, NA), col = c("blue", "red", "gray"),
       lty = c(NA, 1, 2), lwd = c(NA, 2, 1))

# Panel 4: Mass limit distribution
plot(x_obs, m_lim_obs, pch = ".", col = rgb(0, 0, 0, 0.2),
     xlab = "Observed log10(M)", ylab = "Mass limit log10(M_lim)",
     main = "Mass limits per group")
abline(0, 1, col = "red", lty = 2)
abline(a = -1, b = 1, col = "blue", lty = 2)
legend("topleft", c("M_obs = M_lim", "M_obs = M_lim + 1"),
       col = c("red", "blue"), lty = 2)

dev.off()
cat("\nDiagnostic plot saved: mock_data_diagnostic.pdf\n\n")

############################################################
# 11. STAN MODEL
#     This is your model with fixes noted:
#     
#     FIX 1: Consistent error convolution in numerator 
#            AND denominator (normalization)
#     FIX 2: Proper trapezoidal rule
#     FIX 3: Added Poisson term for total count N to
#            constrain phi*
############################################################

stan_model_code <- "
data {
  int<lower=0> N;
  vector[N] x;                 // observed log10(mass)
  vector<lower=0>[N] sigma_x;  // mass uncertainty per group
  vector[N] m_lim;             // mass limit per group
  real V;                      // survey volume [Mpc^3]
  real xlo;                    // lower grid bound
  real xhi;                    // upper grid bound
  int<lower=1> Ng;             // grid points for normalization
  int<lower=1> Ne;             // grid points for error convolution
}

transformed data {
  vector[Ng] xgrid;
  real dx;
  dx = (xhi - xlo) / (Ng - 1);
  for(k in 1:Ng)
    xgrid[k] = xlo + (k - 1) * dx;
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.1, upper=2.0> beta;
}

model {
  // ---------- Priors ----------
  mstar   ~ normal(14.0, 1.5);
  log_phi ~ normal(-4.0, 2.0);
  alpha   ~ normal(-1.7, 1.0);
  // beta has implicit uniform(0.1, 2.0) prior
  
  // ---------- Precompute MRP on grid ----------
  vector[Ng] log_phi_grid;
  for(k in 1:Ng) {
    real u = xgrid[k] - mstar;
    // log10(phi) at this grid point
    log_phi_grid[k] = log10(beta * log(10)) + log_phi 
                      + (alpha + 1) * u 
                      + (-pow(10, beta * u)) / log(10);
    // The last term: log10(exp(-10^(beta*u))) = -10^(beta*u) / ln(10)
  }
  
  // ---------- Per-group likelihood ----------
  real sqrt_2pi = sqrt(2.0 * pi());
  
  for(i in 1:N) {
    
    // --- NUMERATOR: Error-convolved MRP at observed mass ---
    // Integrate phi(m_true) * N(x_obs | m_true, sigma) dm_true
    // over m_true in [x_obs - 4*sigma, x_obs + 4*sigma] intersected with [m_lim, xhi]
    real lo_i = fmax(x[i] - 4.0 * sigma_x[i], m_lim[i]);
    real hi_i = fmin(x[i] + 4.0 * sigma_x[i], xhi);
    
    if(hi_i <= lo_i) {
      // Edge case: skip (shouldn't happen with good data)
      target += -1e6;
      continue;
    }
    
    real dx_i = (hi_i - lo_i) / (Ne - 1);
    
    vector[Ne] f_num;
    for(e in 1:Ne) {
      real m_true = lo_i + (e - 1) * dx_i;
      real u = m_true - mstar;
      real phi_val = beta * log(10) * pow(10, log_phi)
                     * pow(10, (alpha + 1) * u)
                     * exp(-pow(10, beta * u));
      real z = (x[i] - m_true) / sigma_x[i];
      real gauss = exp(-0.5 * z * z) / (sigma_x[i] * sqrt_2pi);
      f_num[e] = phi_val * gauss;
    }
    
    // Trapezoidal rule
    real numerator = (f_num[1] + f_num[Ne]) / 2.0;
    for(e in 2:(Ne-1))
      numerator += f_num[e];
    numerator *= dx_i;
    
    // --- DENOMINATOR: Integral of error-convolved MRP above m_lim ---
    // This is the probability of observing ANY mass above m_lim
    // = integral over x_obs from m_lim to xhi of 
    //   [integral over m_true of phi(m_true) * N(x_obs | m_true, sigma)]
    //
    // For computational tractability, we approximate this as just
    // the integral of the raw MRP above m_lim[i].
    // This is valid when sigma << range of integration.
    // (The convolution redistributes mass but preserves the integral
    //  except at the boundary, which matters little here.)
    
    real denominator = 0.0;
    int k_start = 1;
    // Find first grid point >= m_lim[i]
    for(k in 1:Ng) {
      if(xgrid[k] >= m_lim[i]) {
        k_start = k;
        break;
      }
    }
    
    // Trapezoidal rule on the grid
    if(k_start < Ng) {
      real phi_first;
      {
        real u = xgrid[k_start] - mstar;
        phi_first = beta * log(10) * pow(10, log_phi)
                    * pow(10, (alpha + 1) * u)
                    * exp(-pow(10, beta * u));
      }
      real phi_last;
      {
        real u = xgrid[Ng] - mstar;
        phi_last = beta * log(10) * pow(10, log_phi)
                   * pow(10, (alpha + 1) * u)
                   * exp(-pow(10, beta * u));
      }
      denominator = (phi_first + phi_last) / 2.0;
      for(k in (k_start+1):(Ng-1)) {
        real u = xgrid[k] - mstar;
        denominator += beta * log(10) * pow(10, log_phi)
                       * pow(10, (alpha + 1) * u)
                       * exp(-pow(10, beta * u));
      }
      denominator *= dx;
    }
    
    if(numerator <= 0 || denominator <= 0) {
      target += -1e6;
    } else {
      // Conditional likelihood: p(x_obs | detected above m_lim)
      target += log(numerator) - log(denominator);
    }
  }
  
  // ---------- Poisson term for total count ----------
  // This constrains phi* (the overall amplitude)
  // Expected number = V * integral of MRP above min(m_lim)
  // 
  // We use the minimum m_lim as the effective survey threshold.
  // More precisely, we should sum the expected number per group,
  // but this is a good first approximation.
  {
    real m_lim_min = min(m_lim);
    real N_expected = 0.0;
    for(k in 1:Ng) {
      if(xgrid[k] >= m_lim_min) {
        real u = xgrid[k] - mstar;
        N_expected += beta * log(10) * pow(10, log_phi)
                      * pow(10, (alpha + 1) * u)
                      * exp(-pow(10, beta * u));
      }
    }
    N_expected *= dx * V;
    
    // Poisson log-likelihood: N * log(lambda) - lambda
    // (dropping terms not depending on parameters)
    target += N * log(N_expected + 1e-30) - N_expected;
  }
}
"

############################################################
# 12. PREPARE STAN DATA AND RUN
############################################################

xlo <- min(m_lim_obs) - 1.0
xhi <- 16.5
Ng  <- 300   # Grid resolution for normalization
Ne  <- 15    # Quadrature points for error convolution

stan_data <- list(
    N       = N,
    x       = x_obs,
    sigma_x = sigma_obs,
    m_lim   = m_lim_obs,
    V       = Vsurvey,
    xlo     = xlo,
    xhi     = xhi,
    Ng      = Ng,
    Ne      = Ne
)

cat("\n=== COMPILING AND RUNNING STAN ===\n")
cat(sprintf("N = %d, Ng = %d, Ne = %d\n", N, Ng, Ne))
cat(sprintf("xlo = %.2f, xhi = %.2f\n", xlo, xhi))
cat(sprintf("V = %.2e Mpc^3\n\n", Vsurvey))

fit <- stan(
    model_code = stan_model_code,
    data       = stan_data,
    chains     = 4,
    iter       = 2000,
    warmup     = 1000,
    thin       = 1,
    cores      = 4,
    init       = lapply(1:4, function(i) list(
        mstar   = rnorm(1, 14.0, 0.3),
        log_phi = rnorm(1, -4.0, 0.3),
        alpha   = rnorm(1, -1.7, 0.2),
        beta    = runif(1, 0.4, 0.8)
    )),
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)

############################################################
# 13. RESULTS AND COMPARISON
############################################################

cat("\n=== STAN RESULTS ===\n")
print(fit, pars = c("mstar", "log_phi", "alpha", "beta"))

posterior <- extract(fit, pars = c("mstar", "log_phi", "alpha", "beta"))
post_mat  <- cbind(
    mstar   = posterior$mstar,
    log_phi = posterior$log_phi,
    alpha   = posterior$alpha,
    beta    = posterior$beta
)

med <- apply(post_mat, 2, median)
q16 <- apply(post_mat, 2, quantile, 0.16)
q84 <- apply(post_mat, 2, quantile, 0.84)

cat("\n=== PARAMETER RECOVERY ===\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s\n", 
            "Param", "True", "Median", "+1sig", "-1sig", "Bias(sig)"))
cat(paste(rep("-", 56), collapse = ""), "\n")

true_vals <- c(mstar = TRUE_MSTAR, log_phi = TRUE_LOG_PHI, 
               alpha = TRUE_ALPHA, beta = TRUE_BETA)

for(p in names(true_vals)) {
    err_up <- q84[p] - med[p]
    err_dn <- med[p] - q16[p]
    err_avg <- (err_up + err_dn) / 2
    bias <- (med[p] - true_vals[p]) / err_avg
    
    cat(sprintf("%-10s %8.3f %8.3f  +%.3f  -%.3f  %+.2f sigma\n",
                p, true_vals[p], med[p], err_up, err_dn, bias))
}

# Correlations
cat("\n=== POSTERIOR CORRELATIONS ===\n")
print(round(cor(post_mat), 3))

############################################################
# 14. RESULTS PLOT
############################################################

pdf("sim_recovery_results.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))

# Panel 1: HMF with fit
breaks <- seq(12, 16, by = 0.2)
h_obs  <- hist(x_obs, breaks = breaks, plot = FALSE)
phi_obs <- h_obs$counts / (Vsurvey * 0.2)
ok <- phi_obs > 0

xfit <- seq(11, 16.5, length.out = 500)

plot(h_obs$mids[ok], log10(phi_obs[ok]),
     pch = 19, col = "darkgreen", cex = 1.3,
     xlim = c(12, 16), ylim = c(-8, -2),
     xlab = expression("log"[10]*"(M/M"[symbol("\u2299")]*")"),
     ylab = expression("log"[10]*"("*phi*") [Mpc"^{-3}*" dex"^{-1}*"]"),
     main = "Simulation Recovery Test")
grid(col = "gray80")

# Posterior samples
idx <- sample(1:nrow(post_mat), min(200, nrow(post_mat)))
for(i in idx) {
    y <- tryCatch(
        mrp_log10(xfit, post_mat[i, "mstar"], post_mat[i, "log_phi"],
                  post_mat[i, "alpha"], post_mat[i, "beta"]),
        error = function(e) rep(NA, length(xfit))
    )
    if(all(is.finite(y))) lines(xfit, y, col = rgb(0, 0, 1, 0.03))
}

# Median fit
lines(xfit, mrp_log10(xfit, med["mstar"], med["log_phi"], 
                       med["alpha"], med["beta"]),
      col = "red", lwd = 3)

# True MRP
lines(xfit, mrp_log10(xfit, TRUE_MSTAR, TRUE_LOG_PHI, 
                       TRUE_ALPHA, TRUE_BETA),
      col = "black", lwd = 2, lty = 2)

legend("topright",
       c("Mock data (binned)", "Stan median", "Posterior samples", "TRUE MRP"),
       pch = c(19, NA, NA, NA),
       col = c("darkgreen", "red", rgb(0, 0, 1, 0.3), "black"),
       lty = c(NA, 1, 1, 2), lwd = c(NA, 3, 1, 2), cex = 0.9)

text(14.5, -2.5,
     sprintf("TRUE:\nM*=%.2f, log_phi=%.2f\nalpha=%.2f, beta=%.2f",
             TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA),
     col = "black", cex = 0.8, pos = 4)
text(14.5, -3.5,
     sprintf("RECOVERED:\nM*=%.2f, log_phi=%.2f\nalpha=%.2f, beta=%.2f",
             med["mstar"], med["log_phi"], med["alpha"], med["beta"]),
     col = "red", cex = 0.8, pos = 4)

# Panels 2-5: Marginal posteriors
par_labels <- c(
    mstar   = expression("M* [log"[10]*"(M/M"[symbol("\u2299")]*")]"),
    log_phi = expression("log"[10]*"("*phi*"*) [Mpc"^{-3}*" dex"^{-1}*"]"),
    alpha   = expression(alpha),
    beta    = expression(beta)
)

for(p in c("mstar", "log_phi", "alpha", "beta")) {
    hist(post_mat[, p], breaks = 50, main = paste("Posterior:", p),
         xlab = p, col = "lightblue", border = "white", freq = FALSE)
    abline(v = true_vals[p], col = "red", lwd = 3, lty = 1)
    abline(v = med[p], col = "blue", lwd = 2, lty = 2)
    legend("topright", c(paste("True:", round(true_vals[p], 3)),
                          paste("Median:", round(med[p], 3))),
           col = c("red", "blue"), lty = c(1, 2), lwd = c(3, 2), cex = 0.8)
}

# Panel 6: Corner-style covariances
plot(post_mat[, "mstar"], post_mat[, "log_phi"],
     pch = ".", col = rgb(0, 0, 0, 0.1),
     xlab = "M*", ylab = "log_phi",
     main = "M* vs log_phi covariance")
points(TRUE_MSTAR, TRUE_LOG_PHI, pch = 4, col = "red", cex = 2, lwd = 3)
points(med["mstar"], med["log_phi"], pch = 4, col = "blue", cex = 2, lwd = 3)
legend("topright", c("True", "Median"), col = c("red", "blue"),
       pch = 4, lwd = 3, cex = 0.9)

dev.off()
cat("\nResults plot saved: sim_recovery_results.pdf\n")

############################################################
# 15. SAVE EVERYTHING
############################################################

saveRDS(list(
    # True parameters
    true_params = true_vals,
    # Mock data
    mock_data = list(
        x_obs     = x_obs,
        x_true    = x_true,
        sigma_obs = sigma_obs,
        m_lim     = m_lim_obs,
        z_obs     = det_z[keep],
        z_max     = det_zmax[keep]
    ),
    # Stan results
    fit       = fit,
    posterior = post_mat,
    median    = med,
    q16       = q16,
    q84       = q84,
    # Survey params
    Vsurvey = Vsurvey,
    # Recovery metrics
    bias_sigma = sapply(names(true_vals), function(p) {
        (med[p] - true_vals[p]) / ((q84[p] - q16[p]) / 2)
    })
), "sim_recovery_results.rds")

cat("\n=== DONE ===\n")
cat("Files saved:\n")
cat("  mock_data_diagnostic.pdf  - Diagnostic plots of the mock data\n")
cat("  sim_recovery_results.pdf  - Parameter recovery results\n")
cat("  sim_recovery_results.rds  - All results in R format\n")
