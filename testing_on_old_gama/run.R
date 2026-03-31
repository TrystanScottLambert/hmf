############################################################
# MRP FIT TO GAMA GROUP DATA
#
# Production script combining:
# - GAMA data prep (from Driver+22 / stan_with_incompleteness)
# - Turnover-based mlim(z) (validated on SHARK mocks)
# - Marginalised likelihood (4-param, integrates out errors)
# - Full MCMC for posteriors
#
# The marginalised model computes for each group:
#   L_i = integral[ phi(mt) * N(x_obs|mt,sig) dmt ]
# which accounts for mass measurement errors and
# Eddington bias without latent variables.
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
# 1. Configuration
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5

# GAMA sky area
sky_area_deg2 <- 179.92
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

cat("============================================\n")
cat("  MRP FIT TO GAMA DATA\n")
cat("============================================\n\n")

############################################################
# 2. Read and prepare GAMA group catalogue
############################################################

cat("Reading GAMA group catalogue...\n")

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit &
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 &
             g3cx$IterCenDec > -3.5, ]

# Mass calculation
magica <- 13.9
parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30

g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$MassA <- g3c$MassA * 100 / ho
g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

# Mass errors from multiplicity (vuvuzela)
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

# Observed log mass
g3c$log_mass <- log10(g3c$MassAfunc)

# Filter to good groups
good <- is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0 &
        !is.na(g3c$log10MassErr) &
        g3c$log_mass > 10 & g3c$log_mass < 17

x_obs     <- g3c$log_mass[good]
sigma_obs <- g3c$log10MassErr[good]
z_obs     <- g3c$Zfof[good]
nfof_obs  <- g3c$Nfof[good]
N_all     <- length(x_obs)

cat(sprintf("  N groups: %d\n", N_all))
cat(sprintf("  Mass range: %.2f -- %.2f (median %.2f)\n",
            min(x_obs), max(x_obs), median(x_obs)))
cat(sprintf("  Sigma range: %.2f -- %.2f (median %.2f)\n",
            min(sigma_obs), max(sigma_obs), median(sigma_obs)))

############################################################
# 3. Turnover-based mlim(z)
############################################################

cat("\nComputing turnover mlim(z)...\n")

nbin_z <- 30
z_be <- seq(zmin, zlimit, length.out=nbin_z+1)
z_bm <- (z_be[-1] + z_be[-(nbin_z+1)]) / 2
hist_bw <- 0.3

turnover_bins <- numeric(nbin_z)
for(b in 1:nbin_z) {
    in_bin <- z_obs >= z_be[b] & z_obs < z_be[b+1]
    masses_in_bin <- x_obs[in_bin]
    if(sum(in_bin) > 20) {
        hh <- hist(masses_in_bin, breaks=seq(10, 16, by=hist_bw), plot=FALSE)
        turnover_bins[b] <- hh$mids[which.max(hh$counts)]
    } else {
        turnover_bins[b] <- NA
    }
}

ok_to <- !is.na(turnover_bins)

# Fit polynomial
fit_to_l <- lm(turnover_bins[ok_to] ~ z_bm[ok_to])
fit_to_q <- lm(turnover_bins[ok_to] ~ z_bm[ok_to] + I(z_bm[ok_to]^2))

if(AIC(fit_to_q) < AIC(fit_to_l) - 2) {
    cc_to <- coef(fit_to_q)
    mlim_of_z <- function(z) cc_to[1] + cc_to[2]*z + cc_to[3]*z^2
    cat(sprintf("  mlim(z) = %.2f + %.2f*z + %.2f*z^2\n", cc_to[1], cc_to[2], cc_to[3]))
} else {
    cc_to <- coef(fit_to_l)
    mlim_of_z <- function(z) cc_to[1] + cc_to[2]*z
    cat(sprintf("  mlim(z) = %.2f + %.2f*z\n", cc_to[1], cc_to[2]))
}

cat(sprintf("  mlim(%.2f) = %.2f,  mlim(%.2f) = %.2f\n",
            zmin, mlim_of_z(zmin), zlimit, mlim_of_z(zlimit)))

############################################################
# 4. Redshift shells
############################################################

cat("\nSetting up redshift shells...\n")

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2

d_edges <- sapply(z_edges, function(z) cosdistCoDist(z, OmegaM=omegam, H0=ho))
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
Vsurvey <- sum(V_sh)

mlim_sh <- mlim_of_z(z_mids)

cat(sprintf("  Survey volume: %.4e Mpc^3\n", Vsurvey))
cat(sprintf("  mlim_sh range: %.2f -- %.2f\n", min(mlim_sh), max(mlim_sh)))

############################################################
# 5. Filter groups above mlim
############################################################

mlim_per <- mlim_of_z(z_obs)
above <- x_obs > mlim_per

x_fit     <- x_obs[above]
sig_fit   <- sigma_obs[above]
z_fit     <- z_obs[above]
mlim_fit  <- mlim_per[above]
N_fit     <- length(x_fit)

cat(sprintf("\n  N above mlim: %d / %d (%.1f%%)\n", N_fit, N_all, 100*N_fit/N_all))

############################################################
# 6. Marginalised Stan model
############################################################

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
  // Priors
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.3, 1.0);

  // MRP on grid for Lambda (Poisson normalisation)
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

  // Per-group marginalised likelihood:
  // L_i = integral[ phi(mt) * N(x_obs|mt,sig) dmt ] from mlim to xhi
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

############################################################
# 7. Compile and prepare Stan data
############################################################

cat("\nCompiling Stan model...\n")
sm <- stan_model(model_code=stan_marg)

stan_data <- list(N=N_fit, x_obs=x_fit, sig=sig_fit, mlim=mlim_fit,
                  Nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh,
                  xhi=16.5, Ng=300, Nint=100)

############################################################
# 8. MAP optimisation
############################################################

cat("\n========== MAP ==========\n")

init_fn <- function() {
    list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
         al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
}

best_opt <- NULL; best_lp_val <- -Inf
for(trial in 1:5) {
    cat(sprintf("  Trial %d... ", trial))
    this_opt <- tryCatch(
        optimizing(sm, data=stan_data, init=init_fn(),
                   hessian=FALSE, as_vector=FALSE,
                   iter=20000, algorithm="LBFGS"),
        error=function(e) { cat("failed\n"); NULL })
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f ms=%.3f lp=%.3f al=%.3f be=%.3f\n",
                    this_opt$value, this_opt$par$ms, this_opt$par$lp,
                    this_opt$par$al, this_opt$par$be))
        if(this_opt$value > best_lp_val) {
            best_lp_val <- this_opt$value; best_opt <- this_opt
        }
    }
}

map_par <- c(best_opt$par$ms, best_opt$par$lp, best_opt$par$al, best_opt$par$be)
cat(sprintf("\nMAP: M*=%.3f  lp=%.3f  al=%.3f  be=%.3f\n",
            map_par[1], map_par[2], map_par[3], map_par[4]))

############################################################
# 9. MCMC
############################################################

cat("\n========== MCMC ==========\n")

init_mcmc <- function() {
    list(ms=map_par[1]+rnorm(1,0,.03), lp=map_par[2]+rnorm(1,0,.03),
         al=map_par[3]+rnorm(1,0,.03),
         be=min(1.9, max(0.15, map_par[4]+rnorm(1,0,.03))))
}

t0 <- proc.time()
fit <- sampling(sm, data=stan_data, chains=4, iter=4000, warmup=2000,
                cores=4, init=lapply(1:4, function(i) init_mcmc()),
                control=list(adapt_delta=0.95, max_treedepth=12))
t_mcmc <- (proc.time()-t0)[3]

pmc <- extract(fit)
pmmc <- cbind(ms=pmc$ms, lp=pmc$lp, al=pmc$al, be=pmc$be)
medmc <- apply(pmmc, 2, median)
q16mc <- apply(pmmc, 2, quantile, 0.16)
q84mc <- apply(pmmc, 2, quantile, 0.84)

summ <- summary(fit, pars=c("ms","lp","al","be"))$summary
n_div <- sum(sapply(1:4, function(ch)
    sum(get_sampler_params(fit, inc_warmup=FALSE)[[ch]][,"divergent__"])))

cat(sprintf("MCMC: %.0fs | Rhat %.3f | n_eff %.0f | Div %d\n",
            t_mcmc, max(summ[,"Rhat"]), min(summ[,"n_eff"]), n_div))

############################################################
# 10. Results
############################################################

cat("\n============================================\n")
cat("  RESULTS\n")
cat("============================================\n\n")

driver22 <- c(ms=13.51, lp=-3.19, al=-1.27, be=0.47)
pn <- c("ms","lp","al","be")
plab <- c("M*","log_phi*","alpha","beta")

cat(sprintf("  %-10s  %12s  %10s\n", "Param", "MCMC", "Driver+22"))
for(i in 1:4) {
    unc <- (q84mc[pn[i]] - q16mc[pn[i]]) / 2
    cat(sprintf("  %-10s  %5.3f+/-%.3f  %10.3f\n",
                plab[i], medmc[pn[i]], unc, driver22[i]))
}

############################################################
# 11. Plot
############################################################

mrp_log10 <- function(x, ms, lp, al, be) {
    log10(pmax(be*log(10)*10^lp*10^((al+1)*(x-ms))*exp(-10^(be*(x-ms))), 1e-30))
}

bin_width <- 0.2
breaks <- seq(10, 16, by=bin_width)
h_all <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_all <- h_all$counts / (Vsurvey * bin_width)
ok_all <- h_all$counts >= 5

h_above <- hist(x_fit, breaks=breaks, plot=FALSE)
phi_above <- h_above$counts / (Vsurvey * bin_width)
ok_above <- h_above$counts >= 5

xfit <- seq(10, 16.5, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)

CairoPDF("MRP_GAMA_FIT.pdf", 14, 10)
par(mfrow=c(2,3))

# Panel 1: HMF with fit
plot(h_all$mids[ok_all], log10(phi_all[ok_all]),
     pch=19, col="grey60", cex=1.0,
     xlim=c(11, 16), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")"),
     main="MRP Fit to GAMA")
grid(col="gray80")

points(h_above$mids[ok_above], log10(phi_above[ok_above]),
       pch=19, col="darkgreen", cex=1.3)

# MCMC draws
idx <- sample(1:nrow(pmmc), min(200, nrow(pmmc)))
for(i in idx) {
    y <- mrp_log10(xfit, pmmc[i,"ms"], pmmc[i,"lp"], pmmc[i,"al"], pmmc[i,"be"])
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
}

lines(xfit, mrp_log10(xfit, medmc["ms"], medmc["lp"], medmc["al"], medmc["be"]),
      col="red", lwd=3)
lines(xfit, mrp_log10(xfit, driver22[1], driver22[2], driver22[3], driver22[4]),
      col="black", lwd=2, lty=2)

legend("topright",
       legend=c("All groups", "Above mlim", "MCMC median",
                "Posterior draws", "Driver+22"),
       col=c("grey60","darkgreen","red",rgb(0,0,1,0.3),"black"),
       pch=c(19,19,NA,NA,NA), lty=c(NA,NA,1,1,2),
       lwd=c(NA,NA,3,1,2), bg="white", cex=0.6)

# Panel 2: Mass-redshift with mlim
smoothScatter(z_obs, x_obs,
              xlab="Redshift", ylab=expression("log"[10]*"(M)"),
              main="Mass-redshift + mlim(z)",
              xlim=c(zmin, zlimit), ylim=c(10, 16),
              colramp=colorRampPalette(c("white","cornflowerblue","blue","darkblue")),
              nbin=200)
lines(z_plot, mlim_of_z(z_plot), col="red", lwd=3)
# Turnover points
points(z_bm[ok_to], turnover_bins[ok_to], pch=15, col="cyan", cex=1.2)
grid(col="gray80")

# Panel 3: Posteriors M*
hist(pmmc[,"ms"], breaks=40, col="steelblue", border="white",
     freq=FALSE, main=expression("M"["*"]*" posterior"), xlab="")
abline(v=driver22[1], col="black", lwd=2, lty=2)
abline(v=medmc["ms"], col="red", lwd=2)

# Panel 4: Posteriors log_phi*
hist(pmmc[,"lp"], breaks=40, col="steelblue", border="white",
     freq=FALSE, main=expression("log"*phi["*"]*" posterior"), xlab="")
abline(v=driver22[2], col="black", lwd=2, lty=2)
abline(v=medmc["lp"], col="red", lwd=2)

# Panel 5: Posteriors alpha
hist(pmmc[,"al"], breaks=40, col="steelblue", border="white",
     freq=FALSE, main=expression(alpha*" posterior"), xlab="")
abline(v=driver22[3], col="black", lwd=2, lty=2)
abline(v=medmc["al"], col="red", lwd=2)

# Panel 6: Posteriors beta
hist(pmmc[,"be"], breaks=40, col="steelblue", border="white",
     freq=FALSE, main=expression(beta*" posterior"), xlab="")
abline(v=driver22[4], col="black", lwd=2, lty=2)
abline(v=medmc["be"], col="red", lwd=2)

dev.off()
cat("\nPlot saved: MRP_GAMA_FIT.pdf\n")

############################################################
# 12. Save
############################################################

saveRDS(list(
    mcmc = list(posterior=pmmc, median=medmc, q16=q16mc, q84=q84mc),
    map  = map_par,
    data = list(x_obs=x_fit, sigma=sig_fit, mlim=mlim_fit,
                N=N_fit, Vsurvey=Vsurvey, V_sh=V_sh, mlim_sh=mlim_sh),
    mlim = list(func=mlim_of_z, coefs=cc_to),
    driver22 = driver22,
    diagnostics = list(Rhat=max(summ[,"Rhat"]), n_eff=min(summ[,"n_eff"]),
                       n_div=n_div, time=t_mcmc)
), "MRP_GAMA_results.rds")

cat("Results saved: MRP_GAMA_results.rds\n\nDone!\n")
