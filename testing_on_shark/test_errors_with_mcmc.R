############################################################
# SHARK RECOVERY: MCMC with turnover mlim(z) + hierarchical
#
# Best combo from the quantile comparison:
#   - Turnover (Wright+17) for mlim(z)
#   - Hierarchical model with latent masses
#   - Mass errors from multiplicity-dependent model
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
cat("  MCMC: Turnover + Hierarchical model\n")
cat("============================================\n\n")

############################################################
# 1. Data prep (same as comparison script)
############################################################

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"

cat("Reading data...\n")
groups <- as.data.table(read_parquet(file.path(data_dir, "groups_shark.parquet")))
groups <- groups[dec > 0]
galaxies <- as.data.table(read_parquet(file.path(data_dir, "galaxies_shark.parquet")))
galaxies <- galaxies[dec > 0]

zmin <- 0.01; zlimit <- 0.25; multi <- 5; mag_limit <- 19.8

# Abundance matching
cat("Abundance matching...\n")
groups_vol <- groups[redshift_cosmological > zmin & 
                     redshift_cosmological < zlimit & mass_virial > 0]
N_total <- nrow(groups_vol)

ho <- 67.37; omegam <- 0.3147
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

ra_range  <- range(groups_vol$ra, na.rm=TRUE)
dec_range <- range(groups_vol$dec, na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)
d_lo <- cosdist(zmin); d_hi <- cosdist(zlimit)
Vsurvey <- (4/3) * pi * (d_hi^3 - d_lo^3) * sky_frac

m_grid <- seq(9, 16.5, by=0.001)
phi_grid <- mrp_phi(m_grid, TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)
cum_counts <- rev(cumsum(rev(phi_grid * 0.001))) * Vsurvey

groups_vol[, rank := frankv(mass_virial, order=-1L)]
groups_vol[, log_mass_am := approx(x=rev(cum_counts), y=rev(m_grid),
                                    xout=rank, rule=2)$y]

# GAMA selection
cat("GAMA selection...\n")
gal_selected <- galaxies[id_fof != -1 & mass_stellar_total > 1e8 & mag_r_SDSS < mag_limit]
gal_counts <- gal_selected[, .(n_gama = .N), by = id_group_sky]

groups_vol <- merge(groups_vol, gal_counts, by="id_group_sky", all.x=TRUE)
groups_vol[is.na(n_gama), n_gama := 0L]
groups_vol[, detected := n_gama >= multi]

groups_gama <- groups_vol[detected == TRUE]
z_gama <- groups_gama$redshift_cosmological
m_gama_true <- groups_gama$log_mass_am

cat(sprintf("  Detected: %d / %d\n", nrow(groups_gama), N_total))

# Add mass errors
cat("Adding mass errors...\n")
xx_err <- seq(2, 22)
yy_err <- c(0.68, 0.39, 0.40, 0.33, 0.28, 0.24, 0.20, 0.19, 0.17,
            0.14, 0.14, 0.13, 0.14, 0.12, 0.12, 0.10, 0.10, 0.10,
            0.09, 0.08, 0.08)

sigma_obs <- approx(xx_err, yy_err, groups_gama$n_gama, rule=2)$y
sigma_obs[sigma_obs < 0.10] <- 0.10
m_gama <- m_gama_true + rnorm(length(m_gama_true), 0, sigma_obs)

cat(sprintf("  Sigma: median=%.2f  range=%.2f--%.2f\n",
            median(sigma_obs), min(sigma_obs), max(sigma_obs)))

############################################################
# 2. Turnover-based mlim(z)
############################################################

cat("\nComputing turnover mlim(z)...\n")

nbin_z <- 30
z_be <- seq(zmin, zlimit, length.out=nbin_z+1)
z_bm <- (z_be[-1] + z_be[-(nbin_z+1)]) / 2
hist_bw <- 0.3

turnover_bins <- numeric(nbin_z)
for(b in 1:nbin_z) {
    in_bin <- z_gama >= z_be[b] & z_gama < z_be[b+1]
    masses_in_bin <- m_gama[in_bin]
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
} else {
    cc_to <- coef(fit_to_l)
    mlim_of_z <- function(z) cc_to[1] + cc_to[2]*z
}

cat(sprintf("  mlim(0.01) = %.2f,  mlim(0.25) = %.2f\n",
            mlim_of_z(zmin), mlim_of_z(zlimit)))

############################################################
# 3. Setup shells and filter
############################################################

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
mlim_sh <- mlim_of_z(z_mids)

mlim_per <- mlim_of_z(z_gama)
above <- m_gama > mlim_per
x_obs  <- m_gama[above]
sig    <- sigma_obs[above]
mlim_i <- mlim_per[above]
m_true_above <- m_gama_true[above]
N <- length(x_obs)

cat(sprintf("  N above mlim: %d / %d (%.1f%%)\n", N, nrow(groups_gama), 100*N/nrow(groups_gama)))

############################################################
# 4. Hierarchical Stan model
############################################################

stan_hier <- "
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
  vector[N] z_raw;
}
transformed parameters {
  vector[N] mt;
  for(i in 1:N) mt[i] = x_obs[i] + sig[i] * z_raw[i];
}
model {
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.3, 1.0);
  z_raw ~ std_normal();
  
  // MRP on grid for Lambda
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
  
  // Log-likelihood at latent masses
  for(i in 1:N) {
    real u = mt[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(mt[i] > mlim[i] && phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

############################################################
# 5. MAP first
############################################################

cat("\nCompiling Stan model...\n")
sm <- stan_model(model_code=stan_hier)

stan_data <- list(N=N, x_obs=x_obs, sig=sig, mlim=mlim_i,
                  Nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh,
                  xhi=16.5, Ng=300)

init_fn <- function() {
    list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
         al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7),
         z_raw=rnorm(N, 0, 0.3))
}

cat("\n========== MAP ==========\n")
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

map_ms <- best_opt$par$ms; map_lp <- best_opt$par$lp
map_al <- best_opt$par$al; map_be <- best_opt$par$be

cat(sprintf("\nMAP:  M*=%.3f  lp=%.3f  al=%.3f  be=%.3f\n", map_ms, map_lp, map_al, map_be))
cat(sprintf("TRUE: M*=%.3f  lp=%.3f  al=%.3f  be=%.3f\n", TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE))

############################################################
# 6. MCMC
############################################################

cat("\n========== MCMC ==========\n")

init_map <- function() {
    list(ms=map_ms+rnorm(1,0,0.02), lp=map_lp+rnorm(1,0,0.02),
         al=map_al+rnorm(1,0,0.02),
         be=min(1.9, max(0.15, map_be+rnorm(1,0,0.02))),
         z_raw=best_opt$par$z_raw + rnorm(N, 0, 0.1))
}

t0 <- proc.time()
fit_mcmc <- sampling(sm, data=stan_data, chains=4, iter=4000, warmup=2000,
                     thin=1, cores=4,
                     init=lapply(1:4, function(i) init_map()),
                     control=list(adapt_delta=0.95, max_treedepth=12))
t_mcmc <- (proc.time()-t0)[3]
cat(sprintf("MCMC time: %.1f seconds\n", t_mcmc))

pmc <- extract(fit_mcmc)
pmmc <- cbind(ms=pmc$ms, lp=pmc$lp, al=pmc$al, be=pmc$be)
medmc <- apply(pmmc, 2, median)
q16mc <- apply(pmmc, 2, quantile, 0.16)
q84mc <- apply(pmmc, 2, quantile, 0.84)

summ <- summary(fit_mcmc, pars=c("ms","lp","al","be"))$summary
n_div <- sum(sapply(1:4, function(ch)
    sum(get_sampler_params(fit_mcmc, inc_warmup=FALSE)[[ch]][,"divergent__"])))

cat(sprintf("Max Rhat: %.3f | Min n_eff: %.0f | Div: %d\n",
            max(summ[,"Rhat"]), min(summ[,"n_eff"]), n_div))

############################################################
# 7. Also run simple model for comparison
############################################################

cat("\n========== Simple model (for comparison) ==========\n")

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

sm_simple <- stan_model(model_code=stan_simple)

stan_data_simple <- list(N=N, x_obs=x_obs, Nsh=nsh, V_sh=V_sh,
                         mlim_sh=mlim_sh, xhi=16.5, Ng=300)

init_simple <- function() {
    list(ms=map_ms+rnorm(1,0,0.05), lp=map_lp+rnorm(1,0,0.05),
         al=map_al+rnorm(1,0,0.05),
         be=min(1.9, max(0.15, map_be+rnorm(1,0,0.05))))
}

fit_simple <- sampling(sm_simple, data=stan_data_simple, chains=4, iter=4000, warmup=2000,
                       thin=1, cores=4,
                       init=lapply(1:4, function(i) init_simple()),
                       control=list(adapt_delta=0.95, max_treedepth=12))

pmc_s <- extract(fit_simple)
pmmc_s <- cbind(ms=pmc_s$ms, lp=pmc_s$lp, al=pmc_s$al, be=pmc_s$be)
medmc_s <- apply(pmmc_s, 2, median)
q16mc_s <- apply(pmmc_s, 2, quantile, 0.16)
q84mc_s <- apply(pmmc_s, 2, quantile, 0.84)

############################################################
# 8. Results
############################################################

cat("\n============================================\n")
cat("  RECOVERY RESULTS: Turnover + MCMC\n")
cat("============================================\n\n")

pn <- c("ms","lp","al","be")
plab <- c("M*","log_phi*","alpha","beta")
tv <- c(TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)

cat("  --- Hierarchical model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc[pn[i]] - q16mc[pn[i]]) / 2
    bias <- (medmc[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc[pn[i]], unc, bias))
}

cat("\n  --- Simple model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc_s[pn[i]] - q16mc_s[pn[i]]) / 2
    bias <- (medmc_s[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc_s[pn[i]], unc, bias))
}

############################################################
# 9. Plot
############################################################

xfit <- seq(10, 16, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)
bin_width <- 0.2
breaks <- seq(9, 16, by=bin_width)

# Binned HMFs
h_am_true <- hist(groups_vol$log_mass_am, breaks=breaks, plot=FALSE)
phi_am_true <- h_am_true$counts / (Vsurvey * bin_width)
phi_am_err  <- sqrt(h_am_true$counts) / (Vsurvey * bin_width)
ok_am_true  <- h_am_true$counts >= 5

h_gama_plot <- hist(m_gama[above], breaks=breaks, plot=FALSE)
phi_gama_plot <- h_gama_plot$counts / (Vsurvey * bin_width)
ok_gama_plot  <- h_gama_plot$counts >= 5

CairoPDF("shark_recovery.pdf", 14, 10)
par(mfrow=c(2,3))

# ---- Panel 1: HMF recovery ----
plot(h_am_true$mids[ok_am_true], log10(phi_am_true[ok_am_true]),
     pch=19, col="grey40", cex=1.3,
     xlim=c(10.5, 15.5), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")"),
     main="HMF Recovery (turnover mlim)")
grid(col="gray80")

arrows(h_am_true$mids[ok_am_true],
       log10(pmax(phi_am_true[ok_am_true] - phi_am_err[ok_am_true], 1e-30)),
       h_am_true$mids[ok_am_true],
       log10(phi_am_true[ok_am_true] + phi_am_err[ok_am_true]),
       angle=90, code=3, length=0.03, col="grey40")

points(h_gama_plot$mids[ok_gama_plot], log10(phi_gama_plot[ok_gama_plot]),
       pch=17, col="grey60", cex=1.0)

# True MRP
lines(xfit, log10(pmax(mrp_phi(xfit, tv[1], tv[2], tv[3], tv[4]), 1e-30)),
      col="blue", lwd=3)

# Hierarchical MCMC draws
idx <- sample(1:nrow(pmmc), min(200, nrow(pmmc)))
for(i in idx) {
    y <- log10(pmax(mrp_phi(xfit, pmmc[i,"ms"], pmmc[i,"lp"], pmmc[i,"al"], pmmc[i,"be"]), 1e-30))
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0.8,0,0,0.03))
}

# Hierarchical median
lines(xfit, log10(pmax(mrp_phi(xfit, medmc["ms"], medmc["lp"], medmc["al"], medmc["be"]), 1e-30)),
      col="red", lwd=3, lty=2)

# Simple median
lines(xfit, log10(pmax(mrp_phi(xfit, medmc_s["ms"], medmc_s["lp"], medmc_s["al"], medmc_s["be"]), 1e-30)),
      col="orange", lwd=2, lty=3)

legend("topright",
       legend=c("True HMF", "Observed (>mlim)", "True MRP",
                "Hier MCMC draws", "Hier median", "Simple median"),
       col=c("grey40","grey60","blue",rgb(0.8,0,0,0.3),"red","orange"),
       pch=c(19,17,NA,NA,NA,NA), lty=c(NA,NA,1,1,2,3), lwd=c(NA,NA,3,1,3,2),
       bg="white", cex=0.55)

# ---- Panel 2: Mass-redshift ----
smoothScatter(z_gama, m_gama,
              xlab="Redshift", ylab=expression("log"[10]*"(M)"),
              main="Observed masses + mlim(z)",
              xlim=c(zmin, zlimit), ylim=c(11, 15.5),
              colramp=colorRampPalette(c("white","cornflowerblue","blue","darkblue")),
              nbin=200)
lines(z_plot, mlim_of_z(z_plot), col="red", lwd=3)
grid(col="gray80")

# ---- Panel 3: Parameter recovery ----
plot(1:4, tv, pch=4, col="blue", cex=2, lwd=3,
     ylim=range(c(tv, medmc, medmc_s), na.rm=TRUE),
     xlab="", ylab="Value", main="Parameter Recovery", xaxt="n")
grid(col="gray80")

# Hier
points(1:4+0.1, medmc, pch=19, col="red", cex=1.5)
arrows(1:4+0.1, q16mc, 1:4+0.1, q84mc,
       angle=90, code=3, length=0.05, col="red", lwd=2)

# Simple
points(1:4-0.1, medmc_s, pch=17, col="orange", cex=1.5)
arrows(1:4-0.1, q16mc_s, 1:4-0.1, q84mc_s,
       angle=90, code=3, length=0.05, col="orange", lwd=2)

axis(1, at=1:4, labels=plab)
legend("topleft", c("True","Hier MCMC","Simple MCMC"),
       col=c("blue","red","orange"), pch=c(4,19,17),
       pt.lwd=c(3,1,1), cex=0.8)

# ---- Panel 4-6: Marginal posteriors ----
par_names <- c("ms","lp","al","be")
par_labels <- c(expression("M"["*"]), expression("log"*phi["*"]),
                expression(alpha), expression(beta))
par_true <- tv

for(p in 1:3) {
    # Hier posterior
    xlims <- range(c(pmmc[,par_names[p]], pmmc_s[,par_names[p]]))
    xlims <- xlims + c(-0.3, 0.3) * diff(xlims)
    
    hist(pmmc[,par_names[p]], breaks=40, col=rgb(1,0,0,0.3), border="white",
         freq=FALSE, main=par_labels[p], xlab="", xlim=xlims)
    
    # Simple posterior
    hist(pmmc_s[,par_names[p]], breaks=40, col=rgb(1,0.6,0,0.3), border="white",
         freq=FALSE, add=TRUE)
    
    abline(v=par_true[p], col="blue", lwd=3)
    
    if(p == 1) legend("topright", c("Hier","Simple","True"),
                      fill=c(rgb(1,0,0,0.3), rgb(1,0.6,0,0.3), NA),
                      border=c("white","white","blue"),
                      col=c(NA,NA,"blue"), lty=c(NA,NA,1), lwd=c(NA,NA,3),
                      cex=0.7)
}

dev.off()
cat("\nPlot saved: shark_recovery.pdf\n\nDone!\n")
