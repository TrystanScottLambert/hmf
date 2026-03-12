############################################################
# NESSIE GROUP CATALOGUE: HMF Recovery
#
# Run the MRP model on groups found by the Nessie group
# finder on SHARK mock galaxies. This tests the full
# realistic pipeline: group finding -> mass estimation ->
# HMF recovery.
#
# Uses mass_proxy * 10 as observed mass.
# Runs both Simple and Marginalised models.
# Compares to Driver+22 MRP (abundance-matched truth).
############################################################

library(arrow)
library(Cairo)
library(celestial)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

TRUE_MS <- 13.51; TRUE_LP <- -3.19; TRUE_AL <- -1.27; TRUE_BE <- 0.47

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}

cat("============================================\n")
cat("  NESSIE GROUP CATALOGUE: HMF Recovery\n")
cat("============================================\n\n")

############################################################
# 1. Read Nessie group catalogue
############################################################

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"

cat("Reading Nessie group catalogue...\n")
nessie <- as.data.table(read_parquet("nessie_groups.parquet"))

cat(sprintf("  Total groups: %d\n", nrow(nessie)))
cat(sprintf("  Columns: %s\n", paste(names(nessie), collapse=", ")))
cat(sprintf("  Multiplicity range: %d -- %d\n", min(nessie$multiplicity), max(nessie$multiplicity)))
cat(sprintf("  Redshift range: %.3f -- %.3f\n", min(nessie$iter_redshift), max(nessie$iter_redshift)))

# Observed mass: mass_proxy * 10
nessie$log_mass_obs <- log10(nessie$mass_proxy * 10)

cat(sprintf("  log10(mass_proxy*10) range: %.2f -- %.2f\n",
            min(nessie$log_mass_obs, na.rm=TRUE), max(nessie$log_mass_obs, na.rm=TRUE)))

############################################################
# 2. Selection cuts (GAMA-like)
############################################################

cat("\nApplying selection cuts...\n")

zmin <- 0.01; zlimit <- 0.25; multi <- 5

# Use iter_redshift for group redshift, dec > 0 for SHARK hemisphere
sel <- nessie$iter_redshift > zmin &
       nessie$iter_redshift < zlimit &
       nessie$multiplicity >= multi &
       nessie$iter_dec > 0 &
       is.finite(nessie$log_mass_obs) &
       nessie$log_mass_obs > 8 &   # sanity floor
       nessie$log_mass_obs < 17     # sanity ceiling

groups <- nessie[sel]
cat(sprintf("  After cuts: %d groups (from %d)\n", nrow(groups), nrow(nessie)))

z_obs <- groups$iter_redshift
m_obs <- groups$log_mass_obs
nfof  <- groups$multiplicity

cat(sprintf("  Mass range: %.2f -- %.2f (median %.2f)\n",
            min(m_obs), max(m_obs), median(m_obs)))
cat(sprintf("  Multiplicity range: %d -- %d (median %d)\n",
            min(nfof), max(nfof), median(nfof)))

############################################################
# 3. Mass errors from multiplicity
############################################################

cat("\nAssigning mass errors from multiplicity...\n")

xx_err <- seq(2, 22)
yy_err <- c(0.68, 0.39, 0.40, 0.33, 0.28, 0.24, 0.20, 0.19, 0.17,
            0.14, 0.14, 0.13, 0.14, 0.12, 0.12, 0.10, 0.10, 0.10,
            0.09, 0.08, 0.08)

sigma_obs <- approx(xx_err, yy_err, nfof, rule=2)$y
sigma_obs[sigma_obs < 0.10] <- 0.10

cat(sprintf("  Sigma range: %.2f -- %.2f (median %.2f)\n",
            min(sigma_obs), max(sigma_obs), median(sigma_obs)))

############################################################
# 4. Cosmology setup
############################################################

ho <- 67.37; omegam <- 0.3147

cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

# Sky area (same as SHARK mock)
groups_all <- as.data.table(read_parquet(file.path(data_dir, "groups_shark.parquet")))
groups_all <- groups_all[dec > 0]
ra_range  <- range(groups_all$ra, na.rm=TRUE)
dec_range <- range(groups_all$dec, na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

d_lo <- cosdist(zmin); d_hi <- cosdist(zlimit)
Vsurvey <- (4/3) * pi * (d_hi^3 - d_lo^3) * sky_frac

cat(sprintf("  Sky area: %.1f deg^2, sky fraction: %.5f\n", sky_area_deg2, sky_frac))
cat(sprintf("  Survey volume: %.4e Mpc^3\n", Vsurvey))

############################################################
# 5. Also load the true AM HMF for comparison
############################################################

cat("\nLoading true abundance-matched HMF for reference...\n")
groups_vol <- groups_all[redshift_cosmological > zmin &
                         redshift_cosmological < zlimit &
                         mass_virial > 0]

m_grid <- seq(9, 16.5, by=0.001)
phi_grid <- mrp_phi(m_grid, TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)
cum_counts <- rev(cumsum(rev(phi_grid * 0.001))) * Vsurvey

groups_vol[, rank := frankv(mass_virial, order=-1L)]
groups_vol[, log_mass_am := approx(x=rev(cum_counts), y=rev(m_grid),
                                    xout=rank, rule=2)$y]

bin_width <- 0.2
breaks <- seq(9, 16, by=bin_width)
h_am_true <- hist(groups_vol$log_mass_am, breaks=breaks, plot=FALSE)
phi_am_true <- h_am_true$counts / (Vsurvey * bin_width)
ok_am <- h_am_true$counts >= 5

############################################################
# 6. Turnover-based mlim(z)
############################################################

cat("\nComputing turnover mlim(z)...\n")

nbin_z <- 30
z_be <- seq(zmin, zlimit, length.out=nbin_z+1)
z_bm <- (z_be[-1] + z_be[-(nbin_z+1)]) / 2
hist_bw <- 0.3

turnover_bins <- numeric(nbin_z)
for(b in 1:nbin_z) {
    in_bin <- z_obs >= z_be[b] & z_obs < z_be[b+1]
    masses_in_bin <- m_obs[in_bin]
    if(sum(in_bin) > 20) {
        hh <- hist(masses_in_bin, breaks=seq(10, 16, by=hist_bw), plot=FALSE)
        turnover_bins[b] <- hh$mids[which.max(hh$counts)]
    } else {
        turnover_bins[b] <- NA
    }
}

ok_to <- !is.na(turnover_bins)
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

cat(sprintf("  mlim(0.01) = %.2f,  mlim(0.25) = %.2f\n",
            mlim_of_z(zmin), mlim_of_z(zlimit)))

############################################################
# 7. Setup shells and filter
############################################################

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
mlim_sh <- mlim_of_z(z_mids)

mlim_per <- mlim_of_z(z_obs)
above <- m_obs > mlim_per
x_fit  <- m_obs[above]
sig_fit <- sigma_obs[above]
z_fit  <- z_obs[above]
N_fit <- length(x_fit)

cat(sprintf("\n  N above mlim: %d / %d (%.1f%%)\n", N_fit, length(m_obs), 100*N_fit/length(m_obs)))

############################################################
# 8. Stan models
############################################################

# Simple Poisson (no errors)
stan_simple <- "
data {
  int<lower=1> N;
  vector[N] x_obs;
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  real xhi;
  int<lower=2> Ng;
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
    real u = x_obs[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

# Marginalised (integrate out latent masses)
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

############################################################
# 9. Compile and prepare data
############################################################

cat("\nCompiling Stan models...\n")
sm_simple <- stan_model(model_code=stan_simple)
sm_marg   <- stan_model(model_code=stan_marg)

stan_data_simple <- list(N=N_fit, x_obs=x_fit, Nsh=nsh, V_sh=V_sh,
                         mlim_sh=mlim_sh, xhi=16.5, Ng=300)

stan_data_marg <- list(N=N_fit, x_obs=x_fit, sig=sig_fit,
                       mlim=mlim_per[above],
                       Nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh,
                       xhi=16.5, Ng=300, Nint=100)

############################################################
# 10. Simple model: MAP + MCMC
############################################################

cat("\n========== Simple model ==========\n")

init_simple <- function() {
    list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
         al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
}

best_opt_s <- NULL; best_lp_s <- -Inf
for(trial in 1:5) {
    cat(sprintf("  MAP trial %d... ", trial))
    this_opt <- tryCatch(
        optimizing(sm_simple, data=stan_data_simple, init=init_simple(),
                   hessian=FALSE, as_vector=FALSE, iter=20000, algorithm="LBFGS"),
        error=function(e) { cat("failed\n"); NULL })
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f ms=%.3f al=%.3f\n", this_opt$value, this_opt$par$ms, this_opt$par$al))
        if(this_opt$value > best_lp_s) { best_lp_s <- this_opt$value; best_opt_s <- this_opt }
    }
}

map_s <- c(best_opt_s$par$ms, best_opt_s$par$lp, best_opt_s$par$al, best_opt_s$par$be)
cat(sprintf("  Simple MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n",
            map_s[1], map_s[2], map_s[3], map_s[4]))

cat("  Running MCMC...\n")
init_s_mcmc <- function() {
    list(ms=map_s[1]+rnorm(1,0,.03), lp=map_s[2]+rnorm(1,0,.03),
         al=map_s[3]+rnorm(1,0,.03),
         be=min(1.9, max(0.15, map_s[4]+rnorm(1,0,.03))))
}

t0 <- proc.time()
fit_s <- sampling(sm_simple, data=stan_data_simple, chains=4, iter=4000, warmup=2000,
                  cores=4, init=lapply(1:4, function(i) init_s_mcmc()),
                  control=list(adapt_delta=0.95, max_treedepth=12))
t_s <- (proc.time()-t0)[3]

pmc_s <- extract(fit_s)
pmmc_s <- cbind(ms=pmc_s$ms, lp=pmc_s$lp, al=pmc_s$al, be=pmc_s$be)
medmc_s <- apply(pmmc_s, 2, median)
q16mc_s <- apply(pmmc_s, 2, quantile, 0.16)
q84mc_s <- apply(pmmc_s, 2, quantile, 0.84)

summ_s <- summary(fit_s, pars=c("ms","lp","al","be"))$summary
n_div_s <- sum(sapply(1:4, function(ch)
    sum(get_sampler_params(fit_s, inc_warmup=FALSE)[[ch]][,"divergent__"])))

cat(sprintf("  MCMC: %.0fs | Rhat %.3f | n_eff %.0f | Div %d\n",
            t_s, max(summ_s[,"Rhat"]), min(summ_s[,"n_eff"]), n_div_s))

############################################################
# 11. Marginalised model: MAP + MCMC
############################################################

cat("\n========== Marginalised model ==========\n")

init_marg <- function() {
    list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
         al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
}

best_opt_m <- NULL; best_lp_m <- -Inf
for(trial in 1:5) {
    cat(sprintf("  MAP trial %d... ", trial))
    this_opt <- tryCatch(
        optimizing(sm_marg, data=stan_data_marg, init=init_marg(),
                   hessian=FALSE, as_vector=FALSE, iter=20000, algorithm="LBFGS"),
        error=function(e) { cat("failed\n"); NULL })
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f ms=%.3f al=%.3f\n", this_opt$value, this_opt$par$ms, this_opt$par$al))
        if(this_opt$value > best_lp_m) { best_lp_m <- this_opt$value; best_opt_m <- this_opt }
    }
}

map_m <- c(best_opt_m$par$ms, best_opt_m$par$lp, best_opt_m$par$al, best_opt_m$par$be)
cat(sprintf("  Marg MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n",
            map_m[1], map_m[2], map_m[3], map_m[4]))

cat("  Running MCMC...\n")
init_m_mcmc <- function() {
    list(ms=map_m[1]+rnorm(1,0,.03), lp=map_m[2]+rnorm(1,0,.03),
         al=map_m[3]+rnorm(1,0,.03),
         be=min(1.9, max(0.15, map_m[4]+rnorm(1,0,.03))))
}

t0 <- proc.time()
fit_m <- sampling(sm_marg, data=stan_data_marg, chains=4, iter=4000, warmup=2000,
                  cores=4, init=lapply(1:4, function(i) init_m_mcmc()),
                  control=list(adapt_delta=0.95, max_treedepth=12))
t_m <- (proc.time()-t0)[3]

pmc_m <- extract(fit_m)
pmmc_m <- cbind(ms=pmc_m$ms, lp=pmc_m$lp, al=pmc_m$al, be=pmc_m$be)
medmc_m <- apply(pmmc_m, 2, median)
q16mc_m <- apply(pmmc_m, 2, quantile, 0.16)
q84mc_m <- apply(pmmc_m, 2, quantile, 0.84)

summ_m <- summary(fit_m, pars=c("ms","lp","al","be"))$summary
n_div_m <- sum(sapply(1:4, function(ch)
    sum(get_sampler_params(fit_m, inc_warmup=FALSE)[[ch]][,"divergent__"])))

cat(sprintf("  MCMC: %.0fs | Rhat %.3f | n_eff %.0f | Div %d\n",
            t_m, max(summ_m[,"Rhat"]), min(summ_m[,"n_eff"]), n_div_m))

############################################################
# 12. Results
############################################################

cat("\n============================================\n")
cat("  NESSIE RECOVERY RESULTS\n")
cat("============================================\n\n")

tv <- c(TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)
pn <- c("ms","lp","al","be")
plab <- c("M*","log_phi*","alpha","beta")

cat(sprintf("  TRUE (Driver+22): M*=%.3f lp=%.3f al=%.3f be=%.3f\n\n", tv[1],tv[2],tv[3],tv[4]))

cat("  --- Simple model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc_s[pn[i]] - q16mc_s[pn[i]]) / 2
    bias <- (medmc_s[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc_s[pn[i]], unc, bias))
}

cat("\n  --- Marginalised model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc_m[pn[i]] - q16mc_m[pn[i]]) / 2
    bias <- (medmc_m[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc_m[pn[i]], unc, bias))
}

############################################################
# 13. Plot
############################################################

xfit <- seq(10, 16, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)

# Binned Nessie observed HMF
h_nessie <- hist(m_obs, breaks=breaks, plot=FALSE)
phi_nessie <- h_nessie$counts / (Vsurvey * bin_width)
ok_nes <- h_nessie$counts >= 5

# Binned Nessie above mlim
h_above <- hist(x_fit, breaks=breaks, plot=FALSE)
phi_above <- h_above$counts / (Vsurvey * bin_width)
ok_above <- h_above$counts >= 5

CairoPDF("nessie_recovery.pdf", 14, 14)
par(mfrow=c(3,2))

# ---- Panel 1: HMF recovery ----
plot(h_am_true$mids[ok_am], log10(phi_am_true[ok_am]),
     pch=19, col="grey40", cex=1.3,
     xlim=c(10.5, 15.5), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")"),
     main="HMF Recovery: Nessie groups")
grid(col="gray80")

# Nessie observed (all)
points(h_nessie$mids[ok_nes], log10(phi_nessie[ok_nes]),
       pch=17, col="grey70", cex=0.8)

# Above mlim
points(h_above$mids[ok_above], log10(phi_above[ok_above]),
       pch=2, col="grey50", cex=1.0)

# True MRP
lines(xfit, log10(pmax(mrp_phi(xfit, tv[1], tv[2], tv[3], tv[4]), 1e-30)),
      col="blue", lwd=3)

# Marginalised MCMC draws
idx <- sample(1:nrow(pmmc_m), min(200, nrow(pmmc_m)))
for(i in idx) {
    y <- log10(pmax(mrp_phi(xfit, pmmc_m[i,"ms"], pmmc_m[i,"lp"], pmmc_m[i,"al"], pmmc_m[i,"be"]), 1e-30))
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0.8,0,0,0.03))
}

# Medians
lines(xfit, log10(pmax(mrp_phi(xfit, medmc_m[1], medmc_m[2], medmc_m[3], medmc_m[4]), 1e-30)),
      col="red", lwd=3, lty=2)
lines(xfit, log10(pmax(mrp_phi(xfit, medmc_s[1], medmc_s[2], medmc_s[3], medmc_s[4]), 1e-30)),
      col="orange", lwd=2, lty=3)

legend("topright",
       legend=c("True HMF (AM)", "Nessie (all)", "Nessie (>mlim)",
                "True MRP", "Marg MCMC", "Simple MCMC", "Marg draws"),
       col=c("grey40","grey70","grey50","blue","red","orange",rgb(0.8,0,0,0.3)),
       pch=c(19,17,2,NA,NA,NA,NA), lty=c(NA,NA,NA,1,2,3,1),
       lwd=c(NA,NA,NA,3,3,2,1), bg="white", cex=0.5)

# ---- Panel 2: Mass-redshift ----
smoothScatter(z_obs, m_obs,
              xlab="Redshift", ylab=expression("log"[10]*"(M)"),
              main="Nessie groups + mlim(z)",
              xlim=c(zmin, zlimit), ylim=c(10, 15.5),
              colramp=colorRampPalette(c("white","cornflowerblue","blue","darkblue")),
              nbin=200)
lines(z_plot, mlim_of_z(z_plot), col="red", lwd=3)
grid(col="gray80")

# ---- Panel 3: Multiplicity distribution ----
hist(nfof, breaks=seq(4.5, max(nfof)+0.5, by=1),
     col="steelblue", border="white",
     main="Multiplicity distribution", xlab="N members")
abline(v=multi-0.5, col="red", lwd=2, lty=2)

# ---- Panel 4: Mass error distribution ----
hist(sigma_obs, breaks=30, col="steelblue", border="white",
     main="Mass errors (from multiplicity)", xlab="sigma (dex)")
abline(v=median(sigma_obs), col="red", lwd=2, lty=2)

# ---- Panel 5: Parameter bias ----
bias_s <- bias_m <- numeric(4)
for(i in 1:4) {
    bias_s[i] <- (medmc_s[pn[i]] - tv[i]) / ((q84mc_s[pn[i]] - q16mc_s[pn[i]])/2)
    bias_m[i] <- (medmc_m[pn[i]] - tv[i]) / ((q84mc_m[pn[i]] - q16mc_m[pn[i]])/2)
}

plot(1:4, bias_m, pch=19, col="red", cex=2,
     ylim=range(c(bias_s, bias_m, -3, 3)),
     xlab="", ylab="Bias (sigma)", main="Bias: Marginalised vs Simple", xaxt="n")
points(1:4+0.15, bias_s, pch=17, col="orange", cex=2)
abline(h=0, col="black", lwd=2)
abline(h=c(-2,2), col="grey50", lty=3)
axis(1, at=1:4, labels=plab)
legend("topleft", c("Marginalised","Simple"), pch=c(19,17), col=c("red","orange"), cex=0.8)

# ---- Panel 6: Marginal posteriors (M* and alpha) ----
# Two-panel within one: M* on top half
par(mfrow=c(3,2))  # already set, this is panel 6

xlims_ms <- range(c(pmmc_m[,"ms"], pmmc_s[,"ms"]))
xlims_ms <- xlims_ms + c(-0.2, 0.2) * diff(xlims_ms)

hist(pmmc_m[,"ms"], breaks=40, col=rgb(1,0,0,0.3), border="white",
     freq=FALSE, main=expression("M"["*"]*" posterior"), xlab="", xlim=xlims_ms)
hist(pmmc_s[,"ms"], breaks=40, col=rgb(1,0.6,0,0.3), border="white",
     freq=FALSE, add=TRUE)
abline(v=tv[1], col="blue", lwd=3)
legend("topright", c("Marg","Simple","True"),
       fill=c(rgb(1,0,0,0.3), rgb(1,0.6,0,0.3), NA),
       border=c("white","white",NA),
       col=c(NA,NA,"blue"), lty=c(NA,NA,1), lwd=c(NA,NA,3), cex=0.6)

dev.off()
cat("\nPlot saved: nessie_recovery.pdf\n\nDone!\n")
