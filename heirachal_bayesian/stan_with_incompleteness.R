############################################################
# STAGE 3 (ALTERNATIVE) - TRUNCATED MRP
# Each galaxy has a detection limit M_lim[i] based on its
# redshift. We fit the MRP truncated at each galaxy's limit.
#
# For galaxy i at redshift z[i]:
#   p(m_i | detected) = phi(m_i) / integral[phi(m)]_{M_lim[i]}^∞
#
# This accounts for the fact that we can only detect groups
# above some mass threshold at each redshift.
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

vol_max_survey <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max_survey^3 * 179.92*(pi/180)^2 / (4*pi)

cat("Full survey volume:", signif(Vsurvey,4), "Mpc^3\n")

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit & 
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 & 
             g3cx$IterCenDec > -3.5, ]

# Mass calculation (exact from gamahmf.r)
magica <- 13.9
parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30

g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$MassA <- g3c$MassA * 100 / ho
g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

# Mass errors
xx <- seq(3, 22)
yy <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 
        0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506, 
        0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983, 
        0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)

g3c$log10MassErr <- approx(xx, yy, g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)] <- 0.03
g3c$log10MassErr[g3c$log10MassErr < 0.1] <- 0.1

# Mass corrections
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

############################################################
# Calculate mass limit for each galaxy
# This depends on redshift - more distant groups need higher
# mass to be detected
############################################################

cat("Calculating mass limits per galaxy...\n")

# Read galaxy data to get zmax (maximum z at which group could be detected)
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

# Mass limit for each group: 
# If zmax < zlimit, this group is near the detection limit
# We can estimate M_lim from the fact that at z=zmax, this group
# is barely detectable.
#
# Simple model: M_lim ∝ luminosity_limit ∝ D_L^2 ∝ (1+z)^2 for z << 1
# More accurately, use ratio of zmax to zlimit

g3c$mass_limit_factor <- (g3c$zmax / zlimit)^2

# At z=zlimit, assume detection limit is 10^12.5
# At lower z, limit scales downward
base_mass_limit <- 12.5
g3c$log_mass_limit <- base_mass_limit + log10(g3c$mass_limit_factor)

# Groups detected at their limit have observed mass ~ mass limit
# Groups detected well above limit can have lower masses
# Use a conservative estimate: limit is 1 dex below observed mass or the
# calculated limit, whichever is lower

g3c$log_mass_limit <- pmin(g3c$log_mass_limit, 
                            log10(g3c$MassAfunc) - 1.0)

cat("Mass limit range:", range(g3c$log_mass_limit, na.rm=TRUE), "\n")

############################################################
# Prepare data for Stan
############################################################

x_obs     <- log10(g3c$MassAfunc[is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0])
sigma_obs <- g3c$log10MassErr[is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0]
m_lim_obs <- g3c$log_mass_limit[is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0]

# Remove any with NA limits
keep <- !is.na(m_lim_obs) & x_obs > 10
x_obs     <- x_obs[keep]
sigma_obs <- sigma_obs[keep]
m_lim_obs <- m_lim_obs[keep]

N <- length(x_obs)

cat("N groups:", N, "\n")
cat("Mass range:", round(range(x_obs), 3), "\n")
cat("Limit range:", round(range(m_lim_obs), 3), "\n\n")

xlo <- 10.0
xhi <- 16.0
Ng  <- 600
Ne  <- 21

############################################################
# STAGE 3 STAN MODEL - Truncated MRP
############################################################

stage3_truncated <- "
data {
  int<lower=0> N;
  vector[N] x;                 // observed log10(mass)
  vector<lower=0>[N] sigma_x;  // mass uncertainty
  vector[N] m_lim;             // detection limit per galaxy
  real V;                      // total survey volume
  real xlo;
  real xhi;
  int<lower=1> Ng;
  int<lower=1> Ne;
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

  // ---- Full MRP integral (over entire survey) ----
  {
    vector[Ng] phi_grid;
    for(k in 1:Ng) {
      real u = beta * (xgrid[k] - mstar);
      phi_grid[k] = beta * log(10) * pow(10, log_phi)
                    * pow(10, (alpha+1) * (xgrid[k] - mstar))
                    * exp(-pow(10, u));
    }
    real Lambda = V * sum(phi_grid) * dx;
    target += -Lambda;
  }

  // ---- Truncated likelihood for each galaxy ----
  for(i in 1:N) {
    
    // Convolve MRP with measurement error (as before)
    real x_lo_i = x[i] - 4.0 * sigma_x[i];
    real x_hi_i = x[i] + 4.0 * sigma_x[i];
    
    x_lo_i = fmax(x_lo_i, m_lim[i]);  // Can't go below detection limit!
    x_hi_i = fmin(x_hi_i, xhi);
    
    real dx_i = (x_hi_i - x_lo_i) / (Ne - 1);
    
    vector[Ne] integrand;
    real sqrt_2pi = sqrt(2.0 * pi());
    
    for(e in 1:Ne) {
      real x_true = x_lo_i + (e-1) * dx_i;
      real u = beta * (x_true - mstar);
      real phi_true = beta * log(10) * pow(10, log_phi)
                      * pow(10, (alpha+1) * (x_true - mstar))
                      * exp(-pow(10, u));
      
      real z = (x[i] - x_true) / sigma_x[i];
      real gauss = exp(-0.5 * z * z) / (sigma_x[i] * sqrt_2pi);
      
      integrand[e] = phi_true * gauss;
    }
    
    real phi_convolved = sum(integrand) * dx_i;
    
    // Normalization integral: phi integrated from m_lim[i] to infinity
    // This accounts for selection: given that we detected this galaxy,
    // it must be above its mass limit
    vector[Ng] phi_trunc;
    int n_above = 0;
    for(k in 1:Ng) {
      if(xgrid[k] >= m_lim[i]) {
        n_above = n_above + 1;
        real u = beta * (xgrid[k] - mstar);
        phi_trunc[n_above] = beta * log(10) * pow(10, log_phi)
                             * pow(10, (alpha+1) * (xgrid[k] - mstar))
                             * exp(-pow(10, u));
      }
    }
    
    real phi_norm = sum(phi_trunc[1:n_above]) * dx;
    
    // Likelihood: probability of observing x[i] given detection
    target += log(V) + log(phi_convolved) - log(phi_norm);
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
    m_lim   = m_lim_obs,
    V       = Vsurvey,
    xlo     = xlo,
    xhi     = xhi,
    Ng      = Ng,
    Ne      = Ne
)

cat("Compiling Stage 3 model (truncated MRP)...\n")
fit3 <- stan(
    model_code = stage3_truncated,
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

cat("\n=== STAGE 3 RESULTS (truncated MRP) ===\n")
print(fit3, pars=c("mstar","log_phi","alpha","beta"))

############################################################
# Extract and plot
############################################################

posterior3 <- extract(fit3, pars=c("mstar","log_phi","alpha","beta"))
posterior_matrix <- cbind(
    mstar   = posterior3$mstar,
    log_phi = posterior3$log_phi,
    alpha   = posterior3$alpha,
    beta    = posterior3$beta
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
breaks   <- seq(11, 16.5, by=0.2)
hist_plt <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_plt  <- hist_plt$counts / (Vsurvey * 0.2)
ok       <- phi_plt > 0
xfit     <- seq(11, 16, length.out=500)

CairoPDF("MRP_STAGE3_TRUNCATED.pdf", 10, 7)

plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 3: Truncated MRP (Selection Function)")

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
       legend=c("GAMA (binned)", "Truncated fit",
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

cat("\nPlot: MRP_STAGE3_TRUNCATED.pdf\n")

saveRDS(list(
    fit = fit3,
    posterior = posterior_matrix,
    median = med,
    q16 = q16,
    q84 = q84,
    mass_limits = m_lim_obs
), "stage3_truncated_results.rds")

cat("\nNote: Each galaxy has its own truncation mass.\n")
cat("Likelihood accounts for selection effects.\n")
