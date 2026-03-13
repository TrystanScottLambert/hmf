############################################################
# NESSIE: Test mass_proxy scaling factors
#
# Run the marginalised model with mass_proxy * {2, 5, 10}
# to see how the mass calibration factor affects recovery.
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
tv <- c(TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}

cat("============================================\n")
cat("  NESSIE: Mass scaling comparison\n")
cat("============================================\n\n")

############################################################
# 1. Read data (once)
############################################################

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"

nessie <- as.data.table(read_parquet("nessie_groups.parquet"))

zmin <- 0.01; zlimit <- 0.25; multi <- 5
ho <- 67.37; omegam <- 0.3147

cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

# Sky area
groups_all <- as.data.table(read_parquet(file.path(data_dir, "groups_shark.parquet")))
groups_all <- groups_all[dec > 0]
ra_range  <- range(groups_all$ra, na.rm=TRUE)
dec_range <- range(groups_all$dec, na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

# Shells
nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
Vsurvey <- sum(V_sh)

# True AM HMF for reference
groups_vol <- groups_all[redshift_cosmological > zmin &
                         redshift_cosmological < zlimit & mass_virial > 0]
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

# Base selection (redshift + multiplicity + dec, but NOT mass yet)
sel_base <- nessie$iter_redshift > zmin &
            nessie$iter_redshift < zlimit &
            nessie$multiplicity >= multi &
            nessie$iter_dec > 0 &
            nessie$mass_proxy > 0  # positive mass proxy

groups_base <- nessie[sel_base]
z_base <- groups_base$iter_redshift
nfof_base <- groups_base$multiplicity
mass_proxy_base <- groups_base$mass_proxy

# Error model
xx_err <- seq(2, 22)
yy_err <- c(0.68, 0.39, 0.40, 0.33, 0.28, 0.24, 0.20, 0.19, 0.17,
            0.14, 0.14, 0.13, 0.14, 0.12, 0.12, 0.10, 0.10, 0.10,
            0.09, 0.08, 0.08)
sigma_base <- approx(xx_err, yy_err, nfof_base, rule=2)$y
sigma_base[sigma_base < 0.10] <- 0.10

cat(sprintf("  Base sample: %d groups\n", nrow(groups_base)))

############################################################
# 2. Compile marginalised model (once)
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

cat("Compiling model...\n")
sm_marg <- stan_model(model_code=stan_marg)

############################################################
# 3. Loop over scaling factors
############################################################

scale_factors <- c(2, 5, 10)
scale_cols <- c("darkgreen", "orange", "red")
scale_names <- paste0("x", scale_factors)

results_list <- list()

for(si in 1:length(scale_factors)) {
    sf <- scale_factors[si]
    cat(sprintf("\n========== mass_proxy * %d ==========\n", sf))
    
    # Compute masses with this scaling
    m_obs <- log10(mass_proxy_base * sf)
    
    # Remove non-finite
    good <- is.finite(m_obs) & m_obs > 8 & m_obs < 17
    m_use <- m_obs[good]
    z_use <- z_base[good]
    sig_use <- sigma_base[good]
    
    cat(sprintf("  N groups: %d\n", sum(good)))
    cat(sprintf("  Mass range: %.2f -- %.2f (median %.2f)\n",
                min(m_use), max(m_use), median(m_use)))
    
    # Turnover mlim(z)
    nbin_z <- 30
    z_be <- seq(zmin, zlimit, length.out=nbin_z+1)
    z_bm <- (z_be[-1] + z_be[-(nbin_z+1)]) / 2
    
    turnover_bins <- numeric(nbin_z)
    for(b in 1:nbin_z) {
        in_bin <- z_use >= z_be[b] & z_use < z_be[b+1]
        if(sum(in_bin) > 20) {
            hh <- hist(m_use[in_bin], breaks=seq(8, 17, by=0.3), plot=FALSE)
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
        mlim_func <- function(z) cc_to[1] + cc_to[2]*z + cc_to[3]*z^2
    } else {
        cc_to <- coef(fit_to_l)
        mlim_func <- function(z) cc_to[1] + cc_to[2]*z
    }
    environment(mlim_func) <- list2env(list(cc_to=cc_to))
    
    cat(sprintf("  mlim(0.01)=%.2f, mlim(0.25)=%.2f\n",
                mlim_func(zmin), mlim_func(zlimit)))
    
    # Filter above mlim
    mlim_per <- mlim_func(z_use)
    mlim_sh_this <- mlim_func(z_mids)
    above <- m_use > mlim_per
    
    x_fit <- m_use[above]
    sig_fit <- sig_use[above]
    mlim_fit <- mlim_per[above]
    N_fit <- length(x_fit)
    
    cat(sprintf("  N above mlim: %d (%.1f%%)\n", N_fit, 100*N_fit/length(m_use)))
    
    # Stan data
    stan_data <- list(N=N_fit, x_obs=x_fit, sig=sig_fit, mlim=mlim_fit,
                      Nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh_this,
                      xhi=16.5, Ng=300, Nint=100)
    
    # MAP
    init_fn <- function() {
        list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
             al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
    }
    
    best_opt <- NULL; best_lp_val <- -Inf
    for(trial in 1:5) {
        this_opt <- tryCatch(
            optimizing(sm_marg, data=stan_data, init=init_fn(),
                       hessian=FALSE, as_vector=FALSE, iter=20000, algorithm="LBFGS"),
            error=function(e) NULL)
        if(!is.null(this_opt) && this_opt$value > best_lp_val) {
            best_lp_val <- this_opt$value; best_opt <- this_opt
        }
    }
    
    map_par <- c(best_opt$par$ms, best_opt$par$lp, best_opt$par$al, best_opt$par$be)
    cat(sprintf("  MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n",
                map_par[1], map_par[2], map_par[3], map_par[4]))
    
    # MCMC
    cat("  Running MCMC...\n")
    init_mcmc <- function() {
        list(ms=map_par[1]+rnorm(1,0,.03), lp=map_par[2]+rnorm(1,0,.03),
             al=map_par[3]+rnorm(1,0,.03),
             be=min(1.9, max(0.15, map_par[4]+rnorm(1,0,.03))))
    }
    
    t0 <- proc.time()
    fit <- sampling(sm_marg, data=stan_data, chains=4, iter=4000, warmup=2000,
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
    
    cat(sprintf("  MCMC: %.0fs | Rhat %.3f | n_eff %.0f | Div %d\n",
                t_mcmc, max(summ[,"Rhat"]), min(summ[,"n_eff"]), n_div))
    
    # Store results
    results_list[[scale_names[si]]] <- list(
        sf=sf, map=map_par, pmmc=pmmc, med=medmc, q16=q16mc, q84=q84mc,
        N_fit=N_fit, m_obs=m_use, mlim_func=mlim_func
    )
    
    # Print bias
    pn <- c("ms","lp","al","be")
    plab <- c("M*","log_phi*","alpha","beta")
    cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
    for(i in 1:4) {
        unc <- (q84mc[pn[i]] - q16mc[pn[i]]) / 2
        bias <- (medmc[pn[i]] - tv[i]) / unc
        cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                    plab[i], tv[i], medmc[pn[i]], unc, bias))
    }
}

############################################################
# 4. Summary table
############################################################

cat("\n============================================\n")
cat("  SUMMARY: All scaling factors\n")
cat("============================================\n\n")

pn <- c("ms","lp","al","be")
plab <- c("M*","log_phi*","alpha","beta")

cat(sprintf("  TRUE: M*=%.3f lp=%.3f al=%.3f be=%.3f\n\n", tv[1],tv[2],tv[3],tv[4]))

for(sn in names(results_list)) {
    r <- results_list[[sn]]
    cat(sprintf("  --- %s (N=%d) ---\n", sn, r$N_fit))
    for(i in 1:4) {
        unc <- (r$q84[pn[i]] - r$q16[pn[i]]) / 2
        bias <- (r$med[pn[i]] - tv[i]) / unc
        cat(sprintf("    %-10s  %5.3f +/- %.3f  (bias %+.2f sig)\n",
                    plab[i], r$med[pn[i]], unc, bias))
    }
}

############################################################
# 5. Plot
############################################################

xfit <- seq(10, 16, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)

CairoPDF("nessie_scaling.pdf", 14, 10)
par(mfrow=c(2,3))

# ---- Panel 1: HMF recovery comparison ----
plot(h_am_true$mids[ok_am], log10(phi_am_true[ok_am]),
     pch=19, col="grey40", cex=1.3,
     xlim=c(10.5, 15.5), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M)"),
     ylab=expression("log"[10]*"("*phi*")"),
     main="HMF Recovery: scaling comparison")
grid(col="gray80")

# True MRP
lines(xfit, log10(pmax(mrp_phi(xfit, tv[1], tv[2], tv[3], tv[4]), 1e-30)),
      col="blue", lwd=3)

# Each scaling factor
for(si in 1:length(scale_factors)) {
    sn <- scale_names[si]
    r <- results_list[[sn]]
    
    # MCMC draws
    idx <- sample(1:nrow(r$pmmc), min(100, nrow(r$pmmc)))
    for(i in idx) {
        y <- log10(pmax(mrp_phi(xfit, r$pmmc[i,"ms"], r$pmmc[i,"lp"],
                                r$pmmc[i,"al"], r$pmmc[i,"be"]), 1e-30))
        if(all(is.finite(y))) lines(xfit, y, col=adjustcolor(scale_cols[si], 0.05))
    }
    
    # Median
    lines(xfit, log10(pmax(mrp_phi(xfit, r$med["ms"], r$med["lp"],
                                   r$med["al"], r$med["be"]), 1e-30)),
          col=scale_cols[si], lwd=3, lty=si)
}

legend("topright",
       legend=c("True HMF (AM)", "True MRP",
                paste0("x", scale_factors, " (N=", sapply(results_list, function(r) r$N_fit), ")")),
       col=c("grey40","blue", scale_cols),
       pch=c(19,NA,NA,NA,NA), lty=c(NA,1,1,2,3),
       lwd=c(NA,3,3,3,3), bg="white", cex=0.55)

# ---- Panel 2: Mass-redshift for each scaling ----
plot(NA, xlim=c(zmin, zlimit), ylim=c(10, 16),
     xlab="Redshift", ylab=expression("log"[10]*"(M)"),
     main="Mass-redshift + mlim(z)")
grid(col="gray80")

for(si in 1:length(scale_factors)) {
    sn <- scale_names[si]
    r <- results_list[[sn]]
    points(z_base[is.finite(m_obs)], r$m_obs,
           pch=".", col=adjustcolor(scale_cols[si], 0.1), cex=0.5)
    lines(z_plot, r$mlim_func(z_plot), col=scale_cols[si], lwd=2.5, lty=si)
}
legend("bottomright", paste0("x", scale_factors),
       col=scale_cols, lwd=2.5, lty=1:3, cex=0.7)

# ---- Panel 3: Bias comparison ----
plot(NA, xlim=c(0.5, 4.5), ylim=c(-8, 10),
     xlab="", ylab="Bias (sigma)", main="Parameter bias by scaling", xaxt="n")
grid(col="gray80")
abline(h=0, col="black", lwd=2)
abline(h=c(-2,2), col="grey50", lty=3)

for(si in 1:length(scale_factors)) {
    sn <- scale_names[si]
    r <- results_list[[sn]]
    offset <- (si - 2) * 0.15
    
    bias_vec <- numeric(4)
    for(i in 1:4) {
        unc <- (r$q84[pn[i]] - r$q16[pn[i]]) / 2
        bias_vec[i] <- (r$med[pn[i]] - tv[i]) / unc
    }
    
    points(1:4 + offset, bias_vec, pch=19, col=scale_cols[si], cex=1.5)
    for(j in 1:4) lines(c(j+offset, j+offset), c(0, bias_vec[j]),
                        col=scale_cols[si], lwd=2)
}

axis(1, at=1:4, labels=plab)
legend("topright", paste0("x", scale_factors), col=scale_cols, pch=19, cex=0.7)

# ---- Panels 4-6: Posteriors for M*, log_phi*, alpha ----
for(p in 1:3) {
    pname <- pn[p]
    
    all_vals <- unlist(lapply(results_list, function(r) r$pmmc[,pname]))
    xlims <- range(all_vals, na.rm=TRUE)
    xlims <- xlims + c(-0.1, 0.1) * diff(xlims)
    
    # Empty plot
    plot(NA, xlim=xlims, ylim=c(0, 1), xlab="", ylab="Density",
         main=plab[p], yaxt="n")
    
    for(si in 1:length(scale_factors)) {
        sn <- scale_names[si]
        r <- results_list[[sn]]
        d <- density(r$pmmc[,pname], adjust=1.2)
        d$y <- d$y / max(d$y)  # normalise to peak=1
        lines(d$x, d$y, col=scale_cols[si], lwd=2.5, lty=si)
    }
    
    abline(v=tv[p], col="blue", lwd=3)
    
    if(p == 1) legend("topright", c(paste0("x", scale_factors), "True"),
                      col=c(scale_cols, "blue"),
                      lty=c(1:3, 1), lwd=c(2.5,2.5,2.5,3), cex=0.6)
}

dev.off()
cat("\nPlot saved: nessie_scaling.pdf\n\nDone!\n")
