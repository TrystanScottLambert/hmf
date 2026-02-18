############################################################
# STAGE 2 - ADD MASS MEASUREMENT ERRORS
# Each group has sigma_log10M that depends on multiplicity
# We convolve the MRP with the error distribution
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(Cairo)

set.seed(42)

############################################################
# Data preparation
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5

vol_max <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max^3 * 179.92*(pi/180)^2 / (4*pi)

cat("Survey volume:", signif(Vsurvey,4), "Mpc^3\n")

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof >= multi & g3cx$Zfof < zlimit &
             g3cx$Zfof > zmin & g3cx$IterCenDec > -3.5, ]

# Mass calculation (using your GAMA method)
magica <- 13.9
parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30

g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

# Mass errors as function of multiplicity (from your plot)
xx <- seq(3, 22)
yy <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 
        0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506, 
        0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983, 
        0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)

g3c$log10MassErr <- approx(xx, yy, g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)] <- 0.03
g3c$log10MassErr[g3c$log10MassErr < 0.1] <- 0.1

# Mass corrections (from your code)
masscorr <- c(0.0, 0.0, -2.672595e-01, -1.513503e-01, -1.259069e-01, 
              -9.006064e-02, -5.466009e-02, -6.666895e-02, -1.988694e-02,
              -2.439581e-02, -2.067060e-02, -1.812964e-02, -1.556899e-02,
              -1.313664e-02, -1.743112e-02, -7.965513e-03, -1.257178e-02,
              -7.064037e-03, -3.963656e-03, -1.271533e-02, -2.664687e-03,
              -1.691287e-03)

g3c$masscorr <- masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] <- 0.0
g3c$mymasscorr <- g3c$mymass / 10^g3c$masscorr

# Final masses
x_obs <- log10(g3c$mymasscorr[is.finite(g3c$mymasscorr) & g3c$mymasscorr > 0])
sigma_obs <- g3c$log10MassErr[is.finite(g3c$mymasscorr) & g3c$mymasscorr > 0]

# Apply mass cut
mass_limit <- 12.7
keep <- x_obs > mass_limit
x_obs <- x_obs[keep]
sigma_obs <- sigma_obs[keep]

N <- length(x_obs)

cat("N groups:", N, "\n")
cat("Mass range:", round(range(x_obs), 3), "\n")
cat("Error range:", round(range(sigma_obs), 3), "\n\n")

xlo <- mass_limit
xhi <- 15.8
Ng  <- 500   # integration grid for MRP integral
Ne  <- 21    # integration grid for error convolution per group

############################################################
# STAGE 2 STAN MODEL
# Convolves MRP with Gaussian measurement errors
############################################################

stage2_stan <- "
data {
  int<lower=0> N;           // number of groups
  vector[N] x;              // observed log10(mass)
  vector<lower=0>[N] sigma_x;  // mass uncertainty per group
  real V;                   // survey volume
  real xlo;                 // lower mass limit
  real xhi;                 // upper mass limit
  int<lower=1> Ng;          // integration grid for MRP
  int<lower=1> Ne;          // integration grid for errors
}

transformed data {
  // Fixed grid for MRP integral
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

  // ---- MRP integral: Lambda = V * integral[phi(m) dm] ----
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

  // ---- Likelihood with error convolution ----
  // For each observed group, we need:
  //   p(x_obs | theta) = integral[ phi(x_true) * N(x_obs | x_true, sigma) dx_true ]
  //
  // This accounts for the fact that the TRUE mass distribution
  // follows the MRP, but we observe it with Gaussian errors

  for(i in 1:N) {
    // Integration range: x_obs ± 4*sigma
    real x_lo_i = x[i] - 4.0 * sigma_x[i];
    real x_hi_i = x[i] + 4.0 * sigma_x[i];
    
    // Ensure we stay within physical bounds
    x_lo_i = fmax(x_lo_i, xlo);
    x_hi_i = fmin(x_hi_i, xhi);
    
    real dx_i = (x_hi_i - x_lo_i) / (Ne - 1);
    
    // Numerical integration over possible true masses
    vector[Ne] integrand;
    real sqrt_2pi = sqrt(2.0 * pi());
    
    for(e in 1:Ne) {
      real x_true = x_lo_i + (e-1) * dx_i;
      
      // MRP at this true mass
      real u = beta * (x_true - mstar);
      real phi_true = beta * log(10) * pow(10, log_phi)
                      * pow(10, (alpha+1) * (x_true - mstar))
                      * exp(-pow(10, u));
      
      // Gaussian likelihood of observing x[i] given true mass x_true
      real z = (x[i] - x_true) / sigma_x[i];
      real gauss = exp(-0.5 * z * z) / (sigma_x[i] * sqrt_2pi);
      
      integrand[e] = phi_true * gauss;
    }
    
    // Expected rate for observing this group
    real p_obs_i = V * sum(integrand) * dx_i;
    
    // Add to log-likelihood
    target += log(p_obs_i);
  }
}
"

############################################################
# Run Stan
############################################################

stan_data <- list(
    N       = N,
    x       = x_obs,
    sigma_x = sigma_obs,
    V       = Vsurvey,
    xlo     = xlo,
    xhi     = xhi,
    Ng      = Ng,
    Ne      = Ne
)

cat("Compiling Stage 2 model (with mass errors)...\n")
fit2 <- stan(
    model_code = stage2_stan,
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
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\n=== STAGE 2 RESULTS (with mass errors) ===\n")
print(fit2, pars=c("mstar","log_phi","alpha","beta"))

############################################################
# Extract and plot
############################################################

posterior2 <- extract(fit2, pars=c("mstar","log_phi","alpha","beta"))
posterior_matrix <- cbind(
    mstar   = posterior2$mstar,
    log_phi = posterior2$log_phi,
    alpha   = posterior2$alpha,
    beta    = posterior2$beta
)

med <- apply(posterior_matrix, 2, median)
q16 <- apply(posterior_matrix, 2, quantile, 0.16)
q84 <- apply(posterior_matrix, 2, quantile, 0.84)

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))))
}

# Binned data for display
breaks   <- seq(12.5, 16.2, by=0.2)
hist_plt <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_plt  <- hist_plt$counts / (Vsurvey * 0.2)
ok       <- phi_plt > 0
xfit     <- seq(12, 16, length.out=500)

CairoPDF("MRP_STAGE2_WITH_ERRORS.pdf", 10, 7)

plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 2: With Mass Measurement Errors")

grid(col="gray80")

# Posterior samples
idx <- sample(1:nrow(posterior_matrix), 300)
for(i in idx) {
    y <- tryCatch(
        mrp_log10(posterior_matrix[i,"mstar"], posterior_matrix[i,"log_phi"],
                  posterior_matrix[i,"alpha"],  posterior_matrix[i,"beta"], xfit),
        error = function(e) rep(NA, length(xfit))
    )
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
}

# Median fit
lines(xfit, mrp_log10(med["mstar"],med["log_phi"],
                       med["alpha"],med["beta"],xfit),
      col="red", lwd=3)

# Driver et al. reference
lines(xfit, mrp_log10(13.51,-3.19,-1.27,0.47,xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("GAMA (binned)", "With errors (median)",
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

cat("\nPlot: MRP_STAGE2_WITH_ERRORS.pdf\n")

saveRDS(list(
    fit = fit2,
    posterior = posterior_matrix,
    median = med,
    q16 = q16,
    q84 = q84,
    sigma_obs = sigma_obs
), "stage2_results.rds")

cat("Results saved: stage2_results.rds\n")

############################################################
# Compare Stage 1 vs Stage 2
############################################################

if(file.exists("stage1_stan_results.rds")) {
    stage1 <- readRDS("stage1_stan_results.rds")
    
    cat("\n=== COMPARISON: Stage 1 (no errors) vs Stage 2 (with errors) ===\n")
    cat("Parameter | Stage 1 | Stage 2 | Difference\n")
    cat("----------|---------|---------|------------\n")
    
    for(par in c("mstar","log_phi","alpha","beta")) {
        diff <- med[par] - stage1$median[par]
        cat(sprintf("%-9s | %7.3f | %7.3f | %+7.3f\n",
                    par, stage1$median[par], med[par], diff))
    }
    
    cat("\nExpected: Stage 2 should correct for Eddington bias\n")
    cat("- M* should shift higher (scattering moves groups to lower masses)\n")
    cat("- log_phi should adjust to compensate\n")
}
