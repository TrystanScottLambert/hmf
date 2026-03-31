############################################################
# PUBLICATION PLOT: GAMA HMF Recovery
#
# Reads the saved RDS results and makes a Driver+22 style
# plot with:
# - LCDM prediction from Murray+21
# - Poisson error bars on binned data
# - MCMC posterior draws
# - Turnover sensitivity (shift mlim by +/- 0.2 dex)
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
library(Cairo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Load saved results
############################################################

cat("Loading results...\n")
res <- readRDS("MRP_GAMA_results.rds")

pmmc   <- res$mcmc$posterior
medmc  <- res$mcmc$median
q16mc  <- res$mcmc$q16
q84mc  <- res$mcmc$q84
x_fit  <- res$data$x_obs
sig_fit <- res$data$sigma
mlim_fit <- res$data$mlim
N_fit  <- res$data$N
Vsurvey <- res$data$Vsurvey
V_sh   <- res$data$V_sh
mlim_sh <- res$data$mlim_sh
driver22 <- res$driver22
map_par <- res$map
cc_to   <- res$mlim$coefs

# Reconstruct mlim_of_z
if(length(cc_to) == 3) {
    mlim_of_z <- function(z) cc_to[1] + cc_to[2]*z + cc_to[3]*z^2
} else {
    mlim_of_z <- function(z) cc_to[1] + cc_to[2]*z
}

############################################################
# 2. Reload GAMA data for full sample info
############################################################

ho <- 67.37; omegam <- 0.3147; zlimit <- 0.25; zmin <- 0.01; multi <- 5
sky_area_deg2 <- 179.92
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit &
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 &
             g3cx$IterCenDec > -3.5, ]

magica <- 13.9; parsec <- 3.0857E16; G <- 6.67408E-11; msol <- 1.988E30
g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

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
g3c$log_mass <- log10(g3c$MassAfunc)

good <- is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0 &
        g3c$log_mass > 10 & g3c$log_mass < 17
x_all <- g3c$log_mass[good]
z_all <- g3c$Zfof[good]

############################################################
# 3. MRP and LCDM prediction
############################################################

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}
mrp_log10 <- function(x, ms, lp, al, be) {
    log10(pmax(mrp_phi(x, ms, lp, al, be), 1e-30))
}

# LCDM prediction: exact computation from Driver+22 code
# Uses Murray+21 MRP with amplitude A normalised to Omega_M
parsec_cgs <- 3.0857E16
G_cgs      <- 6.67408E-11
msol_cgs   <- 1.988E30
rhocrit <- 3*(1000*ho/(1E6*parsec_cgs))^2 / (8*pi*G_cgs)

betamrp  <- 0.7097976
A_raw    <- 1.727006E-19
mstarmrp <- 14.42947
alphamrp <- -1.864908

lcdm_x <- seq(0, 17, 0.001) + log10(100/ho)
lcdm_y <- A_raw * betamrp * 10^((alphamrp+1)*(lcdm_x - mstarmrp)) *
           exp(-10^(betamrp*(lcdm_x - mstarmrp))) * (ho/100)^3
lcdm_factor <- sum(10^lcdm_x * lcdm_y) * 0.001 * msol_cgs / (1E6*parsec_cgs)^3 / (omegam*rhocrit)
lcdm_phi <- lcdm_y / lcdm_factor

# Apply z=0.1 shift (same as Driver+22: -0.08 dex in mass, +0.08 dex in phi)
lcdm_x_plot <- lcdm_x - 0.08
lcdm_y_plot <- log10(lcdm_phi) + 0.08

############################################################
# 4. Binned HMF with Poisson errors
############################################################

bin_width <- 0.2
breaks <- seq(10, 16, by=bin_width)

# All groups
h_all <- hist(x_all, breaks=breaks, plot=FALSE)
phi_all <- h_all$counts / (Vsurvey * bin_width)
phi_all_err <- sqrt(h_all$counts) / (Vsurvey * bin_width)
ok_all <- h_all$counts >= 3

# Above mlim (the fitted sample)
h_fit <- hist(x_fit, breaks=breaks, plot=FALSE)
phi_fit <- h_fit$counts / (Vsurvey * bin_width)
phi_fit_err <- sqrt(h_fit$counts) / (Vsurvey * bin_width)
ok_fit <- h_fit$counts >= 3

############################################################
# 5. Turnover sensitivity: shift mlim by +/- 0.2 dex
############################################################

cat("Computing turnover sensitivity...\n")

# Compile marginalised model for refitting
stan_marg <- "
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;
  vector[N] mlim;
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  real xhi;
  int<lower=2> Ng;
  int<lower=2> Nint;
}
transformed data {
  real ln10 = log(10.0);
  real xlo = min(mlim_sh) - 0.5;
  real dx = (xhi - xlo) / (Ng - 1.0);
  vector[Ng] xg;
  for(k in 1:Ng) xg[k] = xlo + (k-1)*dx;
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
}
model {
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.3, 1.0);
  
  vector[Ng] pg;
  for(k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
  }
  vector[Ng] cum;
  cum[Ng] = 0;
  for(kr in 1:(Ng-1)) {
    int k = Ng - kr;
    cum[k] = cum[k+1] + 0.5*(pg[k]+pg[k+1])*dx;
  }
  
  real Lambda = 0;
  for(j in 1:Nsh) {
    int k0 = 1;
    for(k in 1:Ng) if(xg[k] <= mlim_sh[j]) k0 = k;
    if(k0 >= Ng) k0 = Ng-1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;
  
  for(i in 1:N) {
    real lo_i = mlim[i];
    real hi_i = fmin(xhi, fmax(x_obs[i] + 5*sig[i], mlim[i] + 8*sig[i]));
    real dmt = (hi_i - lo_i) / (Nint - 1.0);
    
    if(dmt < 1e-6) {
      target += -100;
    } else {
      real sum_trap = 0;
      for(g in 1:Nint) {
        real mt_g = lo_i + (g-1) * dmt;
        real u = mt_g - ms;
        real phi_g = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
        real gauss_g = exp(-0.5*((x_obs[i] - mt_g)/sig[i])^2) / (sig[i] * 2.5066283);
        real integrand = phi_g * gauss_g;
        if(g == 1 || g == Nint)
          sum_trap += 0.5 * integrand;
        else
          sum_trap += integrand;
      }
      sum_trap *= dmt;
      
      if(sum_trap > 1e-30)
        target += log(sum_trap);
      else
        target += -100;
    }
  }
}
"

sm <- stan_model(model_code=stan_marg)

# Mass errors for all groups
xx_err <- seq(3, 22)
yy_err <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685,
            0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506,
            0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983,
            0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)
sigma_all <- approx(xx_err, yy_err, g3c$Nfof[good], rule=2)$y
sigma_all[is.na(sigma_all)] <- 0.10
sigma_all[sigma_all < 0.10] <- 0.10

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2

# Function to run MAP for a given mlim shift
run_map_shift <- function(shift_dex) {
    mlim_shifted <- function(z) mlim_of_z(z) + shift_dex
    
    mlim_per <- mlim_shifted(z_all)
    mlim_sh_s <- mlim_shifted(z_mids)
    above <- x_all > mlim_per
    
    if(sum(above) < 50) return(rep(NA, 4))
    
    sd <- list(N=sum(above), x_obs=x_all[above], sig=sigma_all[above],
               mlim=mlim_per[above], Nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh_s,
               xhi=16.5, Ng=300, Nint=100)
    
    init_fn <- function() list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
                               al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
    
    best_opt <- NULL; best_lp <- -Inf
    for(trial in 1:3) {
        this_opt <- tryCatch(
            optimizing(sm, data=sd, init=init_fn(),
                       hessian=FALSE, as_vector=FALSE, iter=20000, algorithm="LBFGS"),
            error=function(e) NULL)
        if(!is.null(this_opt) && this_opt$value > best_lp) {
            best_lp <- this_opt$value; best_opt <- this_opt
        }
    }
    
    if(is.null(best_opt)) return(rep(NA, 4))
    c(best_opt$par$ms, best_opt$par$lp, best_opt$par$al, best_opt$par$be)
}

shifts <- c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)
shift_results <- list()

for(s in shifts) {
    cat(sprintf("  Shift %+.1f... ", s))
    par_s <- run_map_shift(s)
    shift_results[[as.character(s)]] <- par_s
    
    mlim_per_s <- mlim_of_z(z_all) + s
    n_above <- sum(x_all > mlim_per_s)
    cat(sprintf("N=%d  ms=%.3f al=%.3f\n", n_above, par_s[1], par_s[3]))
}

############################################################
# 6. Plot
############################################################

xfit <- seq(10, 16.5, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)

CairoPDF("MRP_GAMA_PUBLICATION.pdf", 12, 14)
layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE), heights=c(1.5, 1))

# ---- Panel 1: Main HMF plot (Driver+22 style) ----
par(mar=c(5, 5, 3, 2))

plot(NA, xlim=c(12, 16), ylim=c(-8, -2),
     xlab=expression("Halo Mass  log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"(number density) [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="MRP Fit to GAMA", cex.lab=1.3, cex.axis=1.1)
grid(col="gray85")

# LCDM prediction (tabulated, from Driver+22 exact computation)
lines(lcdm_x_plot, lcdm_y_plot, col="black", lwd=2, lty=2)

# MCMC posterior draws
idx <- sample(1:nrow(pmmc), min(300, nrow(pmmc)))
for(i in idx) {
    y <- mrp_log10(xfit, pmmc[i,"ms"], pmmc[i,"lp"], pmmc[i,"al"], pmmc[i,"be"])
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0.2,0.4,0.8,0.03))
}

# Best fit median
lines(xfit, mrp_log10(xfit, medmc["ms"], medmc["lp"], medmc["al"], medmc["be"]),
      col="blue", lwd=3)

# Driver+22 GAMA5 fit
lines(xfit, mrp_log10(xfit, driver22[1], driver22[2], driver22[3], driver22[4]),
      col="red", lwd=2, lty=3)

# Binned data: all groups (open, below mlim)
below_mlim <- ok_all & !ok_fit
if(any(ok_all)) {
    # Plot all as open symbols first
    points(h_all$mids[ok_all], log10(phi_all[ok_all]),
           pch=1, col="grey50", cex=1.2)
}

# Binned data: above mlim (filled, with error bars)
if(any(ok_fit)) {
    phi_lo <- pmax(phi_fit[ok_fit] - phi_fit_err[ok_fit], 1e-30)
    phi_hi <- phi_fit[ok_fit] + phi_fit_err[ok_fit]
    
    arrows(h_fit$mids[ok_fit], log10(phi_lo),
           h_fit$mids[ok_fit], log10(phi_hi),
           angle=90, code=3, length=0.04, col="darkgreen", lwd=1.5)
    
    points(h_fit$mids[ok_fit], log10(phi_fit[ok_fit]),
           pch=19, col="darkgreen", cex=1.5)
}

# Turnover sensitivity band
shift_cols <- c("grey70","grey60","grey50","blue","grey50","grey60","grey70")
for(si in 1:length(shifts)) {
    s <- shifts[si]
    par_s <- shift_results[[as.character(s)]]
    if(!any(is.na(par_s)) && s != 0) {
        lines(xfit, mrp_log10(xfit, par_s[1], par_s[2], par_s[3], par_s[4]),
              col=adjustcolor(shift_cols[si], 0.5), lwd=1, lty=3)
    }
}

# Upper inset: histogram of N groups
usr <- par("usr")
rect(12.0, -3.0, 14.5, -2.2, col="white", border="grey80")
text(12.1, -2.35, sprintf("GAMA z<0.25, N>%d\nN groups = %d (fitted: %d)",
                           multi-1, length(x_all), N_fit),
     adj=c(0,0.5), cex=0.7, col="black")

legend("topright",
       legend=c(expression(Lambda*"CDM (Murray+21, z=0.1)"),
                "This work (MCMC median)",
                "Posterior draws",
                "Driver+22 (GAMA5)",
                "GAMA (above mlim)",
                "GAMA (below mlim)",
                "Turnover sensitivity"),
       col=c("black","blue",rgb(0.2,0.4,0.8,0.3),"red","darkgreen","grey50","grey60"),
       pch=c(NA,NA,NA,NA,19,1,NA),
       lty=c(2,1,1,3,NA,NA,3),
       lwd=c(2,3,1,2,NA,NA,1),
       bg="white", cex=0.65)

# Parameter annotation
text(15.5, -2.5,
     sprintf("M* = %.2f (+%.2f/-%.2f)\nlog\u03c6* = %.2f (+%.2f/-%.2f)\n\u03b1 = %.2f (+%.2f/-%.2f)\n\u03b2 = %.2f (+%.2f/-%.2f)",
             medmc["ms"], q84mc["ms"]-medmc["ms"], medmc["ms"]-q16mc["ms"],
             medmc["lp"], q84mc["lp"]-medmc["lp"], medmc["lp"]-q16mc["lp"],
             medmc["al"], q84mc["al"]-medmc["al"], medmc["al"]-q16mc["al"],
             medmc["be"], q84mc["be"]-medmc["be"], medmc["be"]-q16mc["be"]),
     col="blue", cex=0.65, pos=2)

# ---- Panel 2: Mass-redshift ----
par(mar=c(5, 5, 2, 1))

smoothScatter(z_all, x_all,
              xlab="Redshift", ylab=expression("log"[10]*"(M)"),
              xlim=c(zmin, zlimit), ylim=c(10.5, 16),
              colramp=colorRampPalette(c("white","cornflowerblue","blue","darkblue")),
              nbin=200, cex.lab=1.2)

# Turnover lines for different shifts
for(si in 1:length(shifts)) {
    s <- shifts[si]
    if(s == 0) {
        lines(z_plot, mlim_of_z(z_plot), col="red", lwd=3)
    } else {
        lines(z_plot, mlim_of_z(z_plot) + s, col="grey50", lwd=1, lty=3)
    }
}

text(0.02, 15.5, "mlim(z) \u00b1 0.1, 0.2, 0.3 dex", cex=0.7, adj=c(0,1))
grid(col="gray80")

# ---- Panel 3: Turnover sensitivity on parameters ----
par(mar=c(5, 5, 2, 1))

pn <- c("ms","lp","al","be")
plab <- c(expression("M"["*"]), expression("log"*phi["*"]),
          expression(alpha), expression(beta))

# Collect shift results
shift_ms <- sapply(shift_results, function(r) r[1])
shift_lp <- sapply(shift_results, function(r) r[2])
shift_al <- sapply(shift_results, function(r) r[3])
shift_be <- sapply(shift_results, function(r) r[4])

plot(shifts, shift_al, type="b", pch=19, col="blue", lwd=2,
     xlab="mlim shift (dex)", ylab=expression(alpha),
     main=expression("Turnover sensitivity: "*alpha*" vs mlim shift"),
     cex.lab=1.2)
grid(col="gray80")
abline(h=driver22[3], col="red", lwd=2, lty=2)
abline(h=alphamrp, col="black", lwd=1, lty=3)
abline(v=0, col="grey50", lty=3)

legend("topright", c("MAP fit", "Driver+22", expression(Lambda*"CDM")),
       col=c("blue","red","black"), pch=c(19,NA,NA),
       lty=c(1,2,3), lwd=c(2,2,1), cex=0.8)

dev.off()
cat("\nPlot saved: MRP_GAMA_PUBLICATION.pdf\n")

# Print turnover sensitivity table
cat("\n  --- Turnover sensitivity ---\n")
cat(sprintf("  %-8s  %5s  %8s  %8s  %8s  %8s\n", "Shift", "N", "M*", "lp", "al", "be"))
for(si in 1:length(shifts)) {
    s <- shifts[si]
    par_s <- shift_results[[as.character(s)]]
    mlim_per_s <- mlim_of_z(z_all) + s
    n_above <- sum(x_all > mlim_per_s)
    cat(sprintf("  %+.1f dex  %5d  %8.3f  %8.3f  %8.3f  %8.3f\n",
                s, n_above, par_s[1], par_s[2], par_s[3], par_s[4]))
}

cat("\nDone!\n")
