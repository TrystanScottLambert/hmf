############################################################
# RECOVERY TEST: Marginalised vs Hierarchical model
#
# Validates that marginalising out latent masses analytically
# (via grid integration) gives the same posterior as the
# hierarchical model with z_raw latent variables.
#
# The marginalised model has only 4 parameters (ms, lp, al, be)
# and should sample efficiently with MCMC.
#
# Tests:
# A) Original hierarchical (z_raw) - MAP + short MCMC
# B) Marginalised (grid integral) - MAP + full MCMC
# C) Simple (no errors) - MAP + full MCMC
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(Cairo)

set.seed(42)

############################################################
# 1. Generate data (same as v16 recovery test)
############################################################

TRUE_MSTAR <- 14.13; TRUE_LOG_PHI <- -3.96
TRUE_ALPHA <- -1.68; TRUE_BETA <- 0.63
tv <- c(ms=TRUE_MSTAR, lp=TRUE_LOG_PHI, al=TRUE_ALPHA, be=TRUE_BETA)

mrp_phi <- function(x, mstar, log_phi, alpha, beta) {
    u <- x - mstar
    beta * log(10) * 10^log_phi * 10^((alpha+1)*u) * exp(-10^(beta*u))
}
mrp_log10 <- function(x, mstar, log_phi, alpha, beta)
    log10(pmax(mrp_phi(x, mstar, log_phi, alpha, beta), 1e-30))

ho <- 67.37; omegam <- 0.3147; zlimit <- 0.25; zmin <- 0.01
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

sky_frac <- 60*(pi/180)^2/(4*pi)
nz <- 20
z_edges <- seq(zmin, zlimit, length.out=nz+1)
z_mids <- (z_edges[-1]+z_edges[-(nz+1)])/2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3)*pi*(d_edges[-1]^3 - d_edges[-(nz+1)]^3)*sky_frac
Vsurvey <- sum(V_sh)
mlim_of_z <- function(z) 12.0 + 2.5*(z-zmin)/(zlimit-zmin)
mlim_sh <- mlim_of_z(z_mids)

mass_hi <- 16.0
all_tm <- c(); all_z <- c()
for(j in 1:nz) {
    mg <- seq(mlim_sh[j], mass_hi, by=0.001)
    if(length(mg)<2) next
    pg <- mrp_phi(mg, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
    N_j <- rpois(1, V_sh[j]*sum(pg)*0.001)
    if(N_j==0) next
    pm_j <- max(pg)*1.1; mj <- numeric(0)
    while(length(mj)<N_j) {
        mp <- runif(N_j*10, mlim_sh[j], mass_hi)
        pp <- mrp_phi(mp, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
        mj <- c(mj, mp[runif(length(mp)) < pp/pm_j])
    }
    all_tm <- c(all_tm, mj[1:N_j])
    all_z <- c(all_z, runif(N_j, z_edges[j], z_edges[j+1]))
}
N_true <- length(all_tm)
sigma_all <- pmax(0.10, pmin(0.30, 0.30 - 0.06*(all_tm-13)))
obs_all <- all_tm + rnorm(N_true, 0, sigma_all)
mlim_per <- mlim_of_z(all_z)

cat("=== SETUP ===\n")
cat(sprintf("True: M*=%.2f phi=%.2f a=%.2f b=%.2f\n", tv[1],tv[2],tv[3],tv[4]))
cat(sprintf("N = %d halos, median sigma = %.2f\n\n", N_true, median(sigma_all)))

############################################################
# 2. Stan models
############################################################

# Model A: Hierarchical with latent z_raw (original v16)
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
  al ~ normal(-1.7, 1.0);
  z_raw ~ std_normal();
  
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
    real u = mt[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(mt[i] > mlim[i] && phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

# Model B: Marginalised - integrate out latent masses on a grid
# For each group: log integral(phi(mt) * N(x_obs|mt,sig), mlim_i, xhi) dmt
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
  int<lower=2> Ng;     // grid for Lambda
  int<lower=2> Nint;   // grid for per-group marginalisation
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
  al ~ normal(-1.7, 1.0);
  
  // MRP on grid for Lambda (same as hierarchical)
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
  // L_i = integral from mlim_i to xhi of phi(mt) * N(x_obs_i | mt, sig_i) dmt
  // Compute on adaptive grid centered on x_obs_i
  for(i in 1:N) {
    // Integration range: from mlim[i] to xhi
    // But most weight is near x_obs[i], so use a range that covers
    // mlim[i] to max(x_obs[i] + 4*sig[i], mlim[i] + 6*sig[i])
    real lo_i = mlim[i];
    real hi_i = fmin(xhi, fmax(x_obs[i] + 5*sig[i], mlim[i] + 8*sig[i]));
    real dmt = (hi_i - lo_i) / (Nint - 1.0);
    
    if(dmt < 1e-6) {
      target += -100;  // degenerate case
    } else {
      // Trapezoid integration
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

# Model C: Simple (no errors)
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
  al ~ normal(-1.7, 1.0);
  
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

############################################################
# 3. Compile all models
############################################################

cat("Compiling models...\n")
sm_hier <- stan_model(model_code=stan_hier)
sm_marg <- stan_model(model_code=stan_marg)
sm_simple <- stan_model(model_code=stan_simple)

stan_data_hier <- list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
                       Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200)

stan_data_marg <- list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
                       Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200,
                       Nint=100)

stan_data_simple <- list(N=N_true, x_obs=obs_all, Nsh=nz, V_sh=V_sh,
                         mlim_sh=mlim_sh, xhi=16.5, Ng=200)

############################################################
# 4A. Hierarchical: MAP only (MCMC too slow)
############################################################

cat("\n========== A: Hierarchical (MAP) ==========\n")

init_hier <- function() {
    list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
         al=rnorm(1,-1.7,.1), be=runif(1,.5,.8),
         z_raw=rnorm(N_true, 0, 0.3))
}

best_opt_h <- NULL; best_lp_h <- -Inf
for(trial in 1:5) {
    cat(sprintf("  Trial %d... ", trial))
    this_opt <- tryCatch(
        optimizing(sm_hier, data=stan_data_hier, init=init_hier(),
                   hessian=FALSE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
        error=function(e) { cat("failed\n"); NULL })
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f ms=%.3f al=%.3f\n", this_opt$value, this_opt$par$ms, this_opt$par$al))
        if(this_opt$value > best_lp_h) { best_lp_h <- this_opt$value; best_opt_h <- this_opt }
    }
}

map_hier <- c(best_opt_h$par$ms, best_opt_h$par$lp, best_opt_h$par$al, best_opt_h$par$be)
cat(sprintf("  Hier MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n", map_hier[1],map_hier[2],map_hier[3],map_hier[4]))

############################################################
# 4B. Marginalised: MAP + MCMC
############################################################

cat("\n========== B: Marginalised (MAP + MCMC) ==========\n")

init_marg <- function() {
    list(ms=rnorm(1,14,.3), lp=rnorm(1,-4,.3),
         al=rnorm(1,-1.7,.2), be=runif(1,.3,.7))
}

best_opt_m <- NULL; best_lp_m <- -Inf
for(trial in 1:5) {
    cat(sprintf("  Trial %d... ", trial))
    this_opt <- tryCatch(
        optimizing(sm_marg, data=stan_data_marg, init=init_marg(),
                   hessian=FALSE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
        error=function(e) { cat("failed\n"); NULL })
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f ms=%.3f al=%.3f\n", this_opt$value, this_opt$par$ms, this_opt$par$al))
        if(this_opt$value > best_lp_m) { best_lp_m <- this_opt$value; best_opt_m <- this_opt }
    }
}

map_marg <- c(best_opt_m$par$ms, best_opt_m$par$lp, best_opt_m$par$al, best_opt_m$par$be)
cat(sprintf("  Marg MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n", map_marg[1],map_marg[2],map_marg[3],map_marg[4]))

# MCMC from MAP
cat("  Running MCMC...\n")
init_marg_mcmc <- function() {
    list(ms=map_marg[1]+rnorm(1,0,.02), lp=map_marg[2]+rnorm(1,0,.02),
         al=map_marg[3]+rnorm(1,0,.02),
         be=min(1.9, max(0.15, map_marg[4]+rnorm(1,0,.02))))
}

t0 <- proc.time()
fit_marg <- sampling(sm_marg, data=stan_data_marg, chains=4, iter=4000, warmup=2000,
                     cores=4, init=lapply(1:4, function(i) init_marg_mcmc()),
                     control=list(adapt_delta=0.95, max_treedepth=12))
t_marg <- (proc.time()-t0)[3]

pmc_m <- extract(fit_marg)
pmmc_m <- cbind(ms=pmc_m$ms, lp=pmc_m$lp, al=pmc_m$al, be=pmc_m$be)
medmc_m <- apply(pmmc_m, 2, median)
q16mc_m <- apply(pmmc_m, 2, quantile, 0.16)
q84mc_m <- apply(pmmc_m, 2, quantile, 0.84)

summ_m <- summary(fit_marg, pars=c("ms","lp","al","be"))$summary
n_div_m <- sum(sapply(1:4, function(ch)
    sum(get_sampler_params(fit_marg, inc_warmup=FALSE)[[ch]][,"divergent__"])))

cat(sprintf("  MCMC: %.0f sec | Rhat %.3f | n_eff %.0f | Div %d\n",
            t_marg, max(summ_m[,"Rhat"]), min(summ_m[,"n_eff"]), n_div_m))

############################################################
# 4C. Simple: MAP + MCMC
############################################################

cat("\n========== C: Simple (MAP + MCMC) ==========\n")

init_simple <- function() {
    list(ms=rnorm(1,14,.3), lp=rnorm(1,-4,.3),
         al=rnorm(1,-1.7,.2), be=runif(1,.3,.7))
}

best_opt_s <- NULL; best_lp_s <- -Inf
for(trial in 1:3) {
    this_opt <- tryCatch(
        optimizing(sm_simple, data=stan_data_simple, init=init_simple(),
                   hessian=FALSE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
        error=function(e) NULL)
    if(!is.null(this_opt) && this_opt$value > best_lp_s) {
        best_lp_s <- this_opt$value; best_opt_s <- this_opt
    }
}
map_simple <- c(best_opt_s$par$ms, best_opt_s$par$lp, best_opt_s$par$al, best_opt_s$par$be)
cat(sprintf("  Simple MAP: ms=%.3f lp=%.3f al=%.3f be=%.3f\n", map_simple[1],map_simple[2],map_simple[3],map_simple[4]))

init_simple_mcmc <- function() {
    list(ms=map_simple[1]+rnorm(1,0,.02), lp=map_simple[2]+rnorm(1,0,.02),
         al=map_simple[3]+rnorm(1,0,.02),
         be=min(1.9, max(0.15, map_simple[4]+rnorm(1,0,.02))))
}

fit_simple <- sampling(sm_simple, data=stan_data_simple, chains=4, iter=4000, warmup=2000,
                       cores=4, init=lapply(1:4, function(i) init_simple_mcmc()),
                       control=list(adapt_delta=0.95, max_treedepth=12))

pmc_s <- extract(fit_simple)
pmmc_s <- cbind(ms=pmc_s$ms, lp=pmc_s$lp, al=pmc_s$al, be=pmc_s$be)
medmc_s <- apply(pmmc_s, 2, median)
q16mc_s <- apply(pmmc_s, 2, quantile, 0.16)
q84mc_s <- apply(pmmc_s, 2, quantile, 0.84)

############################################################
# 5. Results comparison
############################################################

cat("\n============================================\n")
cat("  RECOVERY COMPARISON\n")
cat("============================================\n\n")

pn <- c("ms","lp","al","be")
plab <- c("M*","log_phi*","alpha","beta")

cat("  --- MAP comparison ---\n")
cat(sprintf("  %-10s  %8s  %8s  %8s  %8s\n", "Param", "True", "Hier", "Marg", "Simple"))
for(i in 1:4) {
    cat(sprintf("  %-10s  %8.3f  %8.3f  %8.3f  %8.3f\n",
                plab[i], tv[i], map_hier[i], map_marg[i], map_simple[i]))
}

cat("\n  --- MCMC: Marginalised model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc_m[pn[i]] - q16mc_m[pn[i]]) / 2
    bias <- (medmc_m[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc_m[pn[i]], unc, bias))
}

cat("\n  --- MCMC: Simple model ---\n")
cat(sprintf("  %-10s  %8s  %12s  %8s\n", "Param", "True", "MCMC", "Bias(sig)"))
for(i in 1:4) {
    unc <- (q84mc_s[pn[i]] - q16mc_s[pn[i]]) / 2
    bias <- (medmc_s[pn[i]] - tv[i]) / unc
    cat(sprintf("  %-10s  %8.3f  %5.3f+/-%.3f  %+8.2f\n",
                plab[i], tv[i], medmc_s[pn[i]], unc, bias))
}

############################################################
# 6. Plot
############################################################

xfit <- seq(11.5, 16.5, length.out=500)
br <- seq(11,17,by=0.2)
h <- hist(obs_all, breaks=br, plot=FALSE)
ph <- h$counts/(Vsurvey*0.2); ok <- ph>0

CairoPDF("recovery_marginalised.pdf", 14, 10)
par(mfrow=c(2,3))

# Panel 1: HMF recovery
plot(h$mids[ok], log10(ph[ok]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main="HMF recovery: 3 models")
grid(col="gray80")

lines(xfit, mrp_log10(xfit, tv[1], tv[2], tv[3], tv[4]), col="black", lwd=2, lty=2)

# Marg MCMC draws
idx <- sample(1:nrow(pmmc_m), min(200, nrow(pmmc_m)))
for(i in idx) {
    y <- mrp_log10(xfit, pmmc_m[i,"ms"], pmmc_m[i,"lp"], pmmc_m[i,"al"], pmmc_m[i,"be"])
    if(all(is.finite(y))) lines(xfit, y, col=rgb(1,0,0,0.03))
}

lines(xfit, mrp_log10(xfit, map_hier[1], map_hier[2], map_hier[3], map_hier[4]),
      col="blue", lwd=2, lty=1)
lines(xfit, mrp_log10(xfit, medmc_m[1], medmc_m[2], medmc_m[3], medmc_m[4]),
      col="red", lwd=3)
lines(xfit, mrp_log10(xfit, medmc_s[1], medmc_s[2], medmc_s[3], medmc_s[4]),
      col="orange", lwd=2, lty=3)

legend("topright", c("TRUE","Hier MAP","Marg MCMC","Simple MCMC","Marg draws"),
       col=c("black","blue","red","orange",rgb(1,0,0,0.3)),
       lty=c(2,1,1,3,1), lwd=c(2,2,3,2,1), cex=0.7)

# Panels 2-5: Marginal posteriors
for(p in 1:4) {
    pname <- pn[p]
    
    # Get ranges
    all_vals <- c(pmmc_m[,pname], pmmc_s[,pname])
    xlims <- range(all_vals)
    xlims <- xlims + c(-0.2, 0.2) * diff(xlims)
    
    hist(pmmc_m[,pname], breaks=40, col=rgb(1,0,0,0.3), border="white",
         freq=FALSE, main=plab[p], xlab="", xlim=xlims)
    hist(pmmc_s[,pname], breaks=40, col=rgb(1,0.6,0,0.3), border="white",
         freq=FALSE, add=TRUE)
    
    abline(v=tv[p], col="blue", lwd=3)
    abline(v=map_hier[p], col="blue", lwd=2, lty=2)
    
    if(p == 1) legend("topright", c("Marg","Simple","True","Hier MAP"),
                      fill=c(rgb(1,0,0,0.3), rgb(1,0.6,0,0.3), NA, NA),
                      border=c("white","white",NA,NA),
                      col=c(NA,NA,"blue","blue"),
                      lty=c(NA,NA,1,2), lwd=c(NA,NA,3,2), cex=0.6)
}

# Panel 6: Bias comparison
bias_marg <- numeric(4)
bias_simple <- numeric(4)
for(i in 1:4) {
    bias_marg[i] <- (medmc_m[pn[i]] - tv[i]) / ((q84mc_m[pn[i]] - q16mc_m[pn[i]])/2)
    bias_simple[i] <- (medmc_s[pn[i]] - tv[i]) / ((q84mc_s[pn[i]] - q16mc_s[pn[i]])/2)
}

plot(1:4, bias_marg, pch=19, col="red", cex=2,
     ylim=c(-4, 4), xlab="", ylab="Bias (sigma)", main="Bias comparison", xaxt="n")
points(1:4+0.15, bias_simple, pch=17, col="orange", cex=2)
abline(h=0, col="black", lwd=2)
abline(h=c(-2,2), col="grey50", lty=3)
axis(1, at=1:4, labels=plab)
legend("topleft", c("Marginalised","Simple"), pch=c(19,17), col=c("red","orange"), cex=0.8)

dev.off()
cat("\nPlot saved: recovery_marginalised.pdf\n\nDone!\n")
