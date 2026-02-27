############################################################
# Forward-Convolution Eddington Bias (Shuntov+2025 style)
#
# Instead of "correct data then fit" (Driver+22), we:
# 1. Build Vmax-weighted bins (raw, uncorrected)
# 2. Build a mass error kernel from the per-group errors
# 3. In Stan, convolve the MRP model with the kernel and
#    fit the CONVOLVED model to the RAW data
# 4. The fitted MRP parameters are the INTRINSIC (debiased) HMF
#
# This is self-consistent: the Eddington bias correction and
# the model fitting are done simultaneously, not sequentially.
#
# Reference: Shuntov et al. 2025, A&A 695, A20, Section 4.2.4
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
# Vmax-weighted bins (UNCORRECTED - this is what we fit to)
############################################################

massx <- seq(10.3, 16.1, logbin)
gamahmf_uncorrected <- weighted.hist(log10(g3c$MassAfunc), 
                                     w=1/g3c$weightszlimit, 
                                     breaks=massx, 
                                     plot=FALSE)

gamax_all <- gamahmf_uncorrected$mids
gamay_uncorrected_all <- gamahmf_uncorrected$counts / logbin

# Apply mass cut - only bins above completeness limit
# Below 12.7 the survey is incomplete (NFoF≥5 multiplicity limit + flux limit)
# so those bins would bias the fit low. The forward-convolution still extends
# the MODEL below 12.7 to account for scatter-up, but we don't fit to
# incomplete data.
ok <- gamay_uncorrected_all > 0 & gamax_all > mass_limit
gamax <- gamax_all[ok]
gamay_uncorrected <- gamay_uncorrected_all[ok]

N_bins <- length(gamax)
cat("N bins for fitting (M >", mass_limit, "):", N_bins, "\n")

############################################################
# Build per-bin mass error kernel (source-weighted)
#
# For each observed bin at mass M_obs, the Eddington bias comes
# from groups scattering in from nearby masses. The relevant
# kernel width should reflect the errors of groups that COULD
# scatter into this bin — i.e., groups within a few sigma of
# the bin center, weighted by their 1/Vmax contribution.
#
# This is intermediate between:
# - Per-bin median (too narrow at high M, only uses local groups)  
# - Global kernel (too broad at high M, uses all groups equally)
#
# For each bin, we compute a 1/Vmax-weighted mean sigma using
# groups within ±1 dex of the bin center. This naturally
# includes the abundant low-mass groups that can scatter up.
############################################################

cat("\nComputing source-weighted per-bin kernel sigma...\n")

log10mass_all <- log10(g3c$MassAfunc)
vmax_weights <- 1/g3c$weightszlimit

bin_sigma <- numeric(length(gamax))
for(b in seq_along(gamax)) {
    # Include groups within ±1 dex that could scatter into this bin
    # (groups further away contribute negligibly to the convolution)
    nearby <- abs(log10mass_all - gamax[b]) < 1.0
    if(sum(nearby) > 5) {
        # Weight by 1/Vmax AND by how likely they are to scatter to this bin
        # (closer groups contribute more)
        dist_to_bin <- abs(log10mass_all[nearby] - gamax[b])
        scatter_weight <- exp(-0.5 * (dist_to_bin / 0.7)^2)  # ~Gaussian relevance
        w <- vmax_weights[nearby] * scatter_weight
        bin_sigma[b] <- weighted.mean(g3c$log10MassErr[nearby], w)
    } else {
        bin_sigma[b] <- 0.25  # fallback
    }
}

cat("Per-bin source-weighted sigma (dex):", round(bin_sigma, 3), "\n")

# Build a common offset grid
kernel_halfwidth <- 2.0
kernel_dm <- 0.02
kernel_grid <- seq(-kernel_halfwidth, kernel_halfwidth, by=kernel_dm)
cat("Kernel grid:", length(kernel_grid), "points\n")
cat("Kernel range: [", min(kernel_grid), ",", max(kernel_grid), "] dex\n\n")

# Also compute per-bin Poisson errors for the fit
# Count actual number of groups per bin for Poisson stats
# Clip masses to break range before counting (some groups fall outside)
log10mass_clipped <- log10(g3c$MassAfunc)
log10mass_clipped <- log10mass_clipped[log10mass_clipped >= min(massx) & 
                                        log10mass_clipped <= max(massx)]
raw_counts <- hist(log10mass_clipped, breaks=massx, plot=FALSE)$counts
raw_counts_ok <- raw_counts[ok]

# Fractional Poisson error per bin
poisson_frac <- 1.0 / sqrt(pmax(raw_counts_ok, 1))

# Total observational uncertainty per bin in log10(phi)
# Combine Poisson + cosmic variance (7% from Driver+11)
cv_frac <- 0.07
obs_sigma <- sqrt(poisson_frac^2 + cv_frac^2) / log(10)  # convert frac to dex
obs_sigma <- pmax(obs_sigma, 0.10)  # floor at 0.10 dex (accounts for systematic uncertainty)
obs_sigma <- pmin(obs_sigma, 1.0)   # cap at 1.0 dex

cat("Per-bin obs sigma (dex):", round(obs_sigma, 3), "\n")
cat("Per-bin N groups:", raw_counts_ok, "\n\n")

############################################################
# ALSO run the old Monte Carlo method for comparison
############################################################

cat("Running old-style Monte Carlo EdB for comparison (1001 iterations)...\n")

mockcounts <- matrix(0, nrow=1001, ncol=length(gamahmf_uncorrected$mids))

for(iter in 1:1001) {
    if(iter %% 500 == 0) cat("  Iteration", iter, "/", 1001, "\n")
    mockmass <- log10(g3c$MassAfunc) + rnorm(nrow(g3c), 0.0, g3c$log10MassErr)
    mock_hist <- weighted.hist(mockmass, 
                               w=1/g3c$weightszlimit, 
                               breaks=massx, 
                               plot=FALSE)
    mockcounts[iter, ] <- mock_hist$counts
}

meancounts <- apply(mockcounts, 2, mean)
edb <- meancounts / gamahmf_uncorrected$counts
edb[is.infinite(edb)] <- 1.0
edb[is.na(edb)] <- 1.0
gamay_oldcorrected <- (gamahmf_uncorrected$counts / logbin) / edb
gamay_oldcorrected <- gamay_oldcorrected[ok]
edb_ok <- edb[ok]

cat("Old EdB factors:", round(edb_ok, 3), "\n\n")

############################################################
# Stan model: Forward-convolution approach
#
# The MRP is evaluated on a fine grid, convolved with the
# empirical kernel, and the convolved values at bin centers
# are compared to the UNCORRECTED data.
#
# The fitted parameters are the INTRINSIC HMF.
############################################################

stan_model_conv <- "
functions {
    real mrp_phi(real m, real mstar, real log_phi, real alpha, real beta) {
        real dm = m - mstar;
        real bdu = beta * dm;
        if (bdu > 30.0) return 0.0;
        real log_val = log(beta) + log(log(10.0))
                       + log_phi * log(10.0)
                       + (alpha + 1.0) * dm * log(10.0)
                       - pow(10.0, bdu);
        if (log_val < -80.0) return 0.0;
        return exp(log_val);
    }
}

data {
    int<lower=1> N;
    vector[N] x;
    vector<lower=0>[N] y;
    vector<lower=0>[N] obs_err;
    vector<lower=0>[N] bin_sigma;
    int<lower=1> N_kernel;
    vector[N_kernel] kernel_offset;
    real mass_lo;
}

parameters {
    real<lower=13.0, upper=16.0> mstar;
    real<lower=-6.0, upper=-1.0> log_phi;
    real<lower=-2.5, upper=0.0> alpha;
    real<lower=0.1, upper=2.0> beta;
}

model {
    mstar   ~ normal(14.0, 0.8);
    log_phi ~ normal(-3.5, 1.0);
    alpha   ~ normal(-1.5, 0.5);
    beta    ~ normal(0.6, 0.3);

    for (i in 1:N) {
        real phi_conv = 0.0;
        real wt_sum = 0.0;
        real sig = bin_sigma[i];
        
        // Full convolution with NO mass limit truncation
        // The MRP model predicts the HMF below the completeness limit,
        // and halos there genuinely scatter up into observed bins.
        // Truncating would undercount this contribution.
        for (k in 1:N_kernel) {
            real delta = kernel_offset[k];
            real wt = exp(-0.5 * square(delta / sig));
            real m_true = x[i] - delta;
            phi_conv += mrp_phi(m_true, mstar, log_phi, alpha, beta) * wt;
            wt_sum += wt;
        }
        if (wt_sum > 0) phi_conv = phi_conv / wt_sum;
        
        if (phi_conv > 1e-30) {
            log10(y[i]) ~ normal(log10(phi_conv), obs_err[i]);
        }
    }
}

generated quantities {
    vector[N] log_phi_intrinsic;
    vector[N] log_phi_convolved;
    
    for (i in 1:N) {
        real phi_intr = mrp_phi(x[i], mstar, log_phi, alpha, beta);
        log_phi_intrinsic[i] = (phi_intr > 1e-30) ? log10(phi_intr) : -30.0;
        
        real phi_conv = 0.0;
        real wt_sum = 0.0;
        real sig = bin_sigma[i];
        for (k in 1:N_kernel) {
            real delta = kernel_offset[k];
            real wt = exp(-0.5 * square(delta / sig));
            real m_true = x[i] - delta;
            phi_conv += mrp_phi(m_true, mstar, log_phi, alpha, beta) * wt;
            wt_sum += wt;
        }
        if (wt_sum > 0) phi_conv = phi_conv / wt_sum;
        log_phi_convolved[i] = (phi_conv > 1e-30) ? log10(phi_conv) : -30.0;
    }
}
"

# Trim kernel grid
max_sigma <- max(bin_sigma)
kernel_keep <- abs(kernel_grid) <= 3.5 * max_sigma
kernel_grid_trim <- kernel_grid[kernel_keep]

cat("Trimmed kernel grid:", length(kernel_grid_trim), "points (from", length(kernel_grid), ")\n")
cat("  Covers ±", round(max(abs(kernel_grid_trim)), 2), "dex\n\n")

stan_data <- list(
    N             = N_bins,
    x             = gamax,
    y             = gamay_uncorrected,
    obs_err       = obs_sigma,
    bin_sigma     = bin_sigma,
    N_kernel      = length(kernel_grid_trim),
    kernel_offset = kernel_grid_trim,
    mass_lo       = mass_limit
)

cat("Fitting convolved MRP to UNCORRECTED data in Stan...\n")
cat("  N_bins:", N_bins, "\n")
cat("  N_kernel:", length(kernel_grid_trim), "\n\n")

fit <- stan(
    model_code = stan_model_conv,
    data       = stan_data,
    chains     = 4,
    iter       = 6000,
    warmup     = 3000,
    cores      = 4,
    control    = list(adapt_delta = 0.99, max_treedepth = 14)
)

print(fit, pars=c("mstar","log_phi","alpha","beta"))

############################################################
# Extract posteriors
############################################################

posterior <- extract(fit, pars=c("mstar","log_phi","alpha","beta",
                                "log_phi_intrinsic","log_phi_convolved"))

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

# Median intrinsic and convolved phi at bin centers
med_intrinsic <- apply(posterior$log_phi_intrinsic, 2, median)
med_convolved <- apply(posterior$log_phi_convolved, 2, median)
q16_convolved <- apply(posterior$log_phi_convolved, 2, quantile, 0.16)
q84_convolved <- apply(posterior$log_phi_convolved, 2, quantile, 0.84)
q16_intrinsic <- apply(posterior$log_phi_intrinsic, 2, quantile, 0.16)
q84_intrinsic <- apply(posterior$log_phi_intrinsic, 2, quantile, 0.84)

# Implied Eddington bias: difference between convolved and intrinsic
implied_edb <- 10^med_convolved / 10^med_intrinsic
cat("\nImplied Eddington bias factors (convolved/intrinsic):\n")
cat("  Mass bins:", round(gamax, 1), "\n")
cat("  EdB mult: ", round(implied_edb, 3), "\n")
cat("  Old MC EdB:", round(edb_ok, 3), "\n\n")

############################################################
# MRP function for plotting
############################################################

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))))
}

# Convolved MRP on a fine grid for plotting
# Uses per-bin sigma with smooth interpolation
mrp_convolved <- function(mstar, log_phi, alpha, beta, x, 
                          kernel_grid, bin_sigma_vec, gamax_vec,
                          mass_lo = 12.7) {
    result <- numeric(length(x))
    sig_interp <- splinefun(gamax_vec, bin_sigma_vec, method="monoH.FC")
    
    for(i in seq_along(x)) {
        x_clamped <- max(min(x[i], max(gamax_vec)), min(gamax_vec))
        sig <- sig_interp(x_clamped)
        sig <- max(sig, 0.05)
        sig <- min(sig, 0.5)
        
        wt <- exp(-0.5 * (kernel_grid / sig)^2)
        wt <- wt / sum(wt)
        m_true <- x[i] - kernel_grid
        
        phi_vals <- beta * log(10) * 10^log_phi *
                    10^((alpha+1)*(m_true - mstar)) *
                    exp(-10^(beta*(m_true - mstar)))
        phi_vals[!is.finite(phi_vals)] <- 0
        result[i] <- sum(phi_vals * wt)
    }
    result[result <= 0] <- 1e-30
    return(log10(result))
}

############################################################
# LCDM prediction from HMFcalc (as used in Driver+22)
############################################################

betamrp_lcdm  <- 0.7097976
A_lcdm        <- 1.727006E-19
mstarmrp_lcdm <- 14.42947
alphamrp_lcdm <- -1.864908

fitbinwid <- 0.01
rhocrit <- 3*(1000*ho/(1E6*parsec))^2/(8*pi*G)

mrpx_lcdm <- seq(0, 17, 0.001) + log10(100/ho)
mrpy_lcdm <- A_lcdm * betamrp_lcdm * 10^((alphamrp_lcdm+1)*(mrpx_lcdm - mstarmrp_lcdm)) *
             exp(-10^(betamrp_lcdm*(mrpx_lcdm - mstarmrp_lcdm))) * (ho/100)^3
factor_lcdm <- sum(10^mrpx_lcdm * mrpy_lcdm) * 0.001 * msol / (1E6*parsec)^3 / (omegam*rhocrit)
phimrp_lcdm <- A_lcdm / factor_lcdm

# LCDM MRP with the offset Driver uses in his plot (line 352 of gamahmf.r)
# mrpx-0.08, log10(mrpy)-log10(factor)+0.08
lcdm_x <- mrpx_lcdm - 0.08
lcdm_y <- log10(mrpy_lcdm) - log10(factor_lcdm) + 0.08

cat("\nLCDM prediction loaded (HMFcalc, Driver+22 Fig 4)\n")
cat("  MRP params: M*=", mstarmrp_lcdm, " alpha=", alphamrp_lcdm, 
    " beta=", betamrp_lcdm, " phi*=", signif(phimrp_lcdm, 4), "\n\n")

############################################################
# Plot 1: Main HMF plot - shows both methods
############################################################

xfit <- seq(12, 16, length.out=500)
xfit_conv <- seq(min(gamax) - 0.1, max(gamax) + 0.3, length.out=300)  # only plot convolved near data

CairoPDF("MRP_SHUNTOV_CONVOLUTION.pdf", 12, 8)

par(mar=c(5,5,3,2))

# Plot uncorrected data
plot(gamax, log10(gamay_uncorrected),
     pch=17, col="gray40", cex=1.4,
     xlim=c(12.5, 16), ylim=c(-8, -2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Forward-Convolution Eddington Bias (Shuntov+2025 style)")

grid(col="gray85")

# LCDM prediction (HMFcalc)
lcdm_ok <- lcdm_x >= 12 & lcdm_x <= 16 & is.finite(lcdm_y) & lcdm_y > -10
lines(lcdm_x[lcdm_ok], lcdm_y[lcdm_ok], col="darkgreen", lwd=2, lty=2)

# Old corrected data for comparison
points(gamax, log10(gamay_oldcorrected), pch=1, col="orange3", cex=1.2)

# Posterior draws - intrinsic model (debiased) - plot UNDER convolved
idx <- sample(1:nrow(posterior_matrix), 200)
for(i in idx) {
    y_intr <- mrp_log10(posterior_matrix[i,"mstar"], 
                        posterior_matrix[i,"log_phi"],
                        posterior_matrix[i,"alpha"], 
                        posterior_matrix[i,"beta"], xfit)
    if(all(is.finite(y_intr))) lines(xfit, y_intr, col=rgb(1,0,0,0.03))
}

# Posterior draws - convolved model (should match uncorrected data)
for(i in idx) {
    y_conv <- mrp_convolved(posterior_matrix[i,"mstar"], 
                            posterior_matrix[i,"log_phi"],
                            posterior_matrix[i,"alpha"], 
                            posterior_matrix[i,"beta"], 
                            xfit_conv, kernel_grid_trim, bin_sigma, gamax)
    if(all(is.finite(y_conv))) lines(xfit_conv, y_conv, col=rgb(0.2,0.5,1,0.03))
}

# Median intrinsic model (debiased)
lines(xfit, mrp_log10(med["mstar"], med["log_phi"], med["alpha"], med["beta"], xfit),
      col="red", lwd=3)

# Median convolved model (matches data)
lines(xfit_conv, mrp_convolved(med["mstar"], med["log_phi"], med["alpha"], med["beta"],
                          xfit_conv, kernel_grid_trim, bin_sigma, gamax),
      col="blue", lwd=3)

# Driver+22 GAMA5 best fit for reference
lines(xfit, mrp_log10(13.51, -3.19, -1.27, 0.47, xfit),
      col="purple", lwd=2, lty=3)

# Re-plot points on top
points(gamax, log10(gamay_uncorrected), pch=17, col="gray30", cex=1.4)
points(gamax, log10(gamay_oldcorrected), pch=1, col="orange3", cex=1.2)

# Error bars on data
arrows(gamax, log10(gamay_uncorrected) - obs_sigma,
       gamax, log10(gamay_uncorrected) + obs_sigma,
       length=0.03, angle=90, code=3, col="gray40")

legend("bottomleft",
       legend=c("Uncorrected (1/Vmax)", 
                "MC-corrected (Driver method)",
                "Convolved model (fits data)", 
                "Intrinsic model (debiased)",
                "Driver+22 best fit",
                expression(Lambda*"CDM prediction (HMFcalc)")),
       col=c("gray40", "orange3", "blue", "red", "purple", "darkgreen"),
       pch=c(17, 1, NA, NA, NA, NA), 
       lty=c(NA, NA, 1, 1, 3, 2),
       lwd=c(NA, NA, 3, 3, 2, 2), 
       bg="white", cex=0.80)

# Parameter text
text(14.5, -2.3, sprintf(
    "Intrinsic MRP (debiased):\nM* = %.2f (+%.2f/-%.2f)\nlog φ* = %.2f (+%.2f/-%.2f)\nα = %.2f (+%.2f/-%.2f)\nβ = %.2f (+%.2f/-%.2f)",
    med["mstar"], q84["mstar"]-med["mstar"], med["mstar"]-q16["mstar"],
    med["log_phi"], q84["log_phi"]-med["log_phi"], med["log_phi"]-q16["log_phi"],
    med["alpha"], q84["alpha"]-med["alpha"], med["alpha"]-q16["alpha"],
    med["beta"], q84["beta"]-med["beta"], med["beta"]-q16["beta"]),
    col="red", pos=4, cex=0.75)

dev.off()

############################################################
# Plot 2: Kernel visualization
############################################################

CairoPDF("MRP_KERNEL.pdf", 10, 6)

par(mar=c(5,5,3,2))

# Plot Gaussian kernels for representative bins
plot_grid <- seq(-1.5, 1.5, by=0.01)
plot(NULL, xlim=c(-1.5, 1.5), ylim=c(0, 1.05),
     xlab=expression(Delta*" log"[10]*"(M/M"["\u2299"]*") [dex]"),
     ylab="Normalized kernel density",
     main="Source-weighted Per-bin Gaussian Kernels")
grid(col="gray85")

cols <- colorRampPalette(c("blue", "red"))(N_bins)
for(b in seq_along(gamax)) {
    kern <- dnorm(plot_grid, 0, bin_sigma[b])
    kern <- kern / max(kern)
    lines(plot_grid, kern, col=cols[b], lwd=1.5)
}

legend("topright",
       legend=c(sprintf("log M=%.1f (σ=%.2f)", gamax[1], bin_sigma[1]),
                sprintf("log M=%.1f (σ=%.2f)", gamax[round(N_bins/2)], bin_sigma[round(N_bins/2)]),
                sprintf("log M=%.1f (σ=%.2f)", gamax[N_bins], bin_sigma[N_bins])),
       col=c(cols[1], cols[round(N_bins/2)], cols[N_bins]),
       lty=1, lwd=2, bg="white", cex=0.85)

# Inset: sigma vs mass
par(fig=c(0.15, 0.52, 0.45, 0.9), new=TRUE)
plot(gamax, bin_sigma, type="b", pch=19, col="darkblue", cex=0.8,
     xlab="", ylab="", main="Source-weighted σ", cex.main=0.8, cex.axis=0.7)
mtext(expression("log"[10]*"M"), side=1, line=2, cex=0.6)
mtext(expression(sigma*" [dex]"), side=2, line=2, cex=0.6)
grid(col="gray85")

dev.off()

############################################################
# Plot 3: Corner plot
############################################################

CairoPDF("MRP_COVARIANCE_SHUNTOV.pdf", 10, 10)

par(mfrow=c(4,4), mar=c(3,3,1,1), oma=c(0,0,2,0))

params <- c("mstar", "log_phi", "alpha", "beta")
labels <- c("M*", "log φ*", "α", "β")
driver_vals <- c(13.51, -3.19, -1.27, 0.47)

xlims <- list(
    mstar   = c(12, 16),
    log_phi = c(-6, -1),
    alpha   = c(-2.5, 0),
    beta    = c(0, 1.5)
)

for(i in 1:4) {
    for(j in 1:4) {
        if(i == j) {
            hist(posterior_matrix[,params[i]], 
                 breaks=40, col="skyblue", border="white",
                 main="", xlab="", ylab="", xlim=xlims[[params[i]]])
            abline(v=med[params[i]], col="red", lwd=2)
            abline(v=driver_vals[i], col="black", lwd=2, lty=2)
            if(j == 1) mtext(labels[i], side=2, line=2, cex=0.8)
        } else if(i > j) {
            plot(posterior_matrix[,params[j]], posterior_matrix[,params[i]],
                 pch=16, col=rgb(0,0,0,0.05), cex=0.3,
                 xlab="", ylab="", xlim=xlims[[params[j]]], ylim=xlims[[params[i]]])
            points(med[params[j]], med[params[i]], pch=3, col="red", cex=2, lwd=2)
            points(driver_vals[j], driver_vals[i], pch=4, col="black", cex=2, lwd=2)
            if(j == 1) mtext(labels[i], side=2, line=2, cex=0.8)
            if(i == 4) mtext(labels[j], side=1, line=2, cex=0.8)
        } else {
            corr <- cor(posterior_matrix[,params[j]], posterior_matrix[,params[i]])
            plot.new()
            text(0.5, 0.5, sprintf("ρ = %.3f", corr), cex=1.5)
        }
    }
}

mtext("Parameter Covariances - Forward Convolution (Shuntov method)", 
      outer=TRUE, cex=1.0, font=2)

dev.off()

############################################################
# Plot 4: Comparison of implied Eddington bias
############################################################

CairoPDF("MRP_EDB_COMPARISON.pdf", 10, 6)

par(mar=c(5,5,3,2))

plot(gamax, pmin(implied_edb, 15), type="b", pch=19, col="red", lwd=2,
     ylim=c(0.5, min(max(c(pmin(implied_edb, 15), edb_ok), na.rm=TRUE)*1.2, 15)),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab="Eddington Bias Factor (convolved / intrinsic)",
     main="Comparison: Forward-Convolution vs Monte Carlo EdB")
grid(col="gray85")

lines(gamax, edb_ok, type="b", pch=17, col="blue", lwd=2)

abline(h=1.0, lty=2, col="gray50")

legend("topleft",
       legend=c("Forward convolution (Shuntov method)", 
                "Monte Carlo correction (Driver method)"),
       col=c("red", "blue"), pch=c(19, 17), lwd=2, bg="white")

dev.off()

############################################################
# Summary
############################################################

cat("\n============================================\n")
cat("=== FINAL RESULTS (Shuntov-style forward convolution) ===\n")
cat("============================================\n\n")

cat("Method: Fit MRP ⊛ D(M) to UNCORRECTED 1/Vmax data\n")
cat("        where D(M) is empirical kernel from per-group errors\n")
cat("        Fitted parameters = INTRINSIC (debiased) MRP\n\n")

cat(sprintf("           Median    +1sig    -1sig   Driver+22 (GAMA5)\n"))
cat(sprintf("M*:        %6.3f   +%.3f   -%.3f    13.51\n",
            med["mstar"], q84["mstar"]-med["mstar"], med["mstar"]-q16["mstar"]))
cat(sprintf("log_phi:   %6.3f   +%.3f   -%.3f    -3.19\n",
            med["log_phi"], q84["log_phi"]-med["log_phi"], med["log_phi"]-q16["log_phi"]))
cat(sprintf("alpha:     %6.3f   +%.3f   -%.3f    -1.27\n",
            med["alpha"], q84["alpha"]-med["alpha"], med["alpha"]-q16["alpha"]))
cat(sprintf("beta:      %6.3f   +%.3f   -%.3f     0.47\n",
            med["beta"], q84["beta"]-med["beta"], med["beta"]-q16["beta"]))

cat("\nPlots:\n")
cat("  MRP_SHUNTOV_CONVOLUTION.pdf - Main HMF with both methods\n")
cat("  MRP_KERNEL.pdf              - Mass error kernel\n")
cat("  MRP_COVARIANCE_SHUNTOV.pdf  - Parameter covariances\n")
cat("  MRP_EDB_COMPARISON.pdf      - EdB factor comparison\n")

cat("\n✓ Forward-convolution Eddington bias (Shuntov+2025 style)\n")
cat("✓ Empirical kernel from per-group mass errors\n")
cat("✓ Fit to UNCORRECTED data (self-consistent)\n")
cat("✓ Bayesian posteriors via Stan (MCMC)\n")
cat("✓ Intrinsic + convolved model outputs\n")
cat("✓ Old MC method retained for comparison\n")
cat("\nKey advantage: model and bias correction determined simultaneously,\n")
cat("not sequentially. No first-order approximation needed.\n")
