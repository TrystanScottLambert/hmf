############################################################
# SIMULATION-RECOVERY TEST v16
#
# Same hierarchical model as v14 (proven unbiased), but
# using optimisation and variational inference instead of
# full MCMC to avoid the sampling efficiency problem.
#
# Approach A: optimizing() for MAP + Laplace approx
# Approach B: vb() for variational Bayes (ADVI)
# Approach C: Short MCMC with thin=5 for comparison
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup + generate data (same as all versions)
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
# 2. Stan model (same as v14 hierarchical non-centered)
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

stan_data <- list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
                  Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200)

# Compile once
cat("Compiling Stan model...\n")
sm <- stan_model(model_code=stan_hier)

init_fn <- function() {
    list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
         al=rnorm(1,-1.7,.1), be=runif(1,.5,.8),
         z_raw=rnorm(N_true, 0, 0.3))
}

############################################################
# APPROACH A: MAP via optimizing()
############################################################

cat("\n========== APPROACH A: MAP (optimizing) ==========\n")

t0 <- proc.time()
opt <- optimizing(sm, data=stan_data, init=init_fn(),
                  hessian=TRUE, as_vector=FALSE,
                  iter=5000, algorithm="LBFGS")
t_opt <- (proc.time()-t0)[3]

cat(sprintf("Optimizing time: %.1f seconds\n", t_opt))
cat(sprintf("Converged: %s (return code %d)\n", 
            ifelse(opt$return_code==0, "YES", "NO"), opt$return_code))

# Extract MAP estimates
map_ms <- opt$par$ms
map_lp <- opt$par$lp
map_al <- opt$par$al
map_be <- opt$par$be

cat(sprintf("\nMAP estimates:\n"))
cat(sprintf("  ms  = %.3f (true %.3f, delta = %+.3f)\n", map_ms, tv["ms"], map_ms-tv["ms"]))
cat(sprintf("  lp  = %.3f (true %.3f, delta = %+.3f)\n", map_lp, tv["lp"], map_lp-tv["lp"]))
cat(sprintf("  al  = %.3f (true %.3f, delta = %+.3f)\n", map_al, tv["al"], map_al-tv["al"]))
cat(sprintf("  be  = %.3f (true %.3f, delta = %+.3f)\n", map_be, tv["be"], map_be-tv["be"]))

# Latent mass recovery
map_zraw <- opt$par$z_raw
map_mt <- obs_all + sigma_all * map_zraw
cat(sprintf("\nLatent mass: median |MAP - true| = %.3f dex (vs sigma = %.3f)\n",
            median(abs(map_mt - all_tm)), median(sigma_all)))

# Hessian for uncertainties (just for the 4 MRP params)
# The Hessian is for all parameters including z_raw
# We need the marginal variance for the first 4 params
if(!is.null(opt$hessian)) {
    cat("\nHessian available, computing uncertainties...\n")
    # The parameter ordering in the hessian is: ms, lp, al, be, z_raw[1..N]
    H <- opt$hessian
    # For Laplace approx, Sigma = -H^{-1}
    # But H is huge (784x784). We can try to invert just the 4x4 block
    # after marginalising over z_raw using Schur complement.
    # For now, just use the diagonal of the full inverse if feasible.
    np <- 4  # number of MRP params
    nz_raw <- N_true
    
    # Partition: H = [A B; C D] where A is 4x4, D is NxN
    # Marginal covariance of MRP params = -(A - B D^{-1} C)^{-1}
    # This requires inverting the NxN D matrix.
    # Since z_raw are nearly independent given MRP params,
    # D should be nearly diagonal -> cheap to invert.
    
    tryCatch({
        A <- -H[1:np, 1:np]
        B <- -H[1:np, (np+1):(np+nz_raw)]
        D <- -H[(np+1):(np+nz_raw), (np+1):(np+nz_raw)]
        
        # D is nearly diagonal, use diagonal approximation
        D_diag_inv <- diag(1/diag(D))
        Schur <- A - B %*% D_diag_inv %*% t(B)
        Sigma_mrp <- solve(Schur)
        se_laplace <- sqrt(diag(Sigma_mrp))
        
        cat("Laplace SE (diagonal D approx):\n")
        pnames <- c("ms","lp","al","be")
        for(i in 1:4) {
            bias <- (opt$par[[pnames[i]]] - tv[pnames[i]]) / se_laplace[i]
            cat(sprintf("  %-5s MAP=%.3f SE=%.3f  Bias=%+.2f sig\n",
                        pnames[i], opt$par[[pnames[i]]], se_laplace[i], bias))
        }
    }, error=function(e) {
        cat(sprintf("Hessian analysis failed: %s\n", e$message))
    })
} else {
    cat("Hessian not available.\n")
}

############################################################
# APPROACH B: Variational Bayes
############################################################

cat("\n========== APPROACH B: Variational Bayes (ADVI) ==========\n")

t0 <- proc.time()
vb_fit <- tryCatch(
    vb(sm, data=stan_data, init=init_fn(),
       iter=10000, output_samples=4000,
       algorithm="meanfield"),
    error=function(e) { cat(sprintf("VB failed: %s\n", e$message)); NULL }
)
t_vb <- (proc.time()-t0)[3]

if(!is.null(vb_fit)) {
    cat(sprintf("VB time: %.1f seconds\n", t_vb))
    
    pvb <- extract(vb_fit)
    pmvb <- cbind(ms=pvb$ms, lp=pvb$lp, al=pvb$al, be=pvb$be)
    medvb <- apply(pmvb,2,median)
    q16vb <- apply(pmvb,2,quantile,.16)
    q84vb <- apply(pmvb,2,quantile,.84)
    
    cat("VB Recovery:\n")
    for(p in names(tv)) {
        ea <- ((q84vb[p]-medvb[p])+(medvb[p]-q16vb[p]))/2
        bias <- (medvb[p]-tv[p])/ea
        st <- ifelse(abs(bias)<1,"OK",ifelse(abs(bias)<2,"marginal","FAIL"))
        cat(sprintf("  %-5s True=%.3f Med=%.3f SE~%.3f Bias=%+.2f sig [%s]\n",
                    p, tv[p], medvb[p], ea, bias, st))
    }
}

############################################################
# APPROACH C: Short MCMC with thinning
############################################################

cat("\n========== APPROACH C: Short MCMC (thin=5) ==========\n")

t0 <- proc.time()
fit_mcmc <- sampling(sm, data=stan_data,
    chains=4, iter=2000, warmup=1000, thin=5, cores=4,
    init=lapply(1:4, function(i) init_fn()),
    control=list(adapt_delta=0.85, max_treedepth=12))
t_mcmc <- (proc.time()-t0)[3]

cat(sprintf("MCMC time: %.1f seconds\n", t_mcmc))

pmc <- extract(fit_mcmc)
pmmc <- cbind(ms=pmc$ms, lp=pmc$lp, al=pmc$al, be=pmc$be)
medmc <- apply(pmmc,2,median)
q16mc <- apply(pmmc,2,quantile,.16)
q84mc <- apply(pmmc,2,quantile,.84)

summ_mc <- summary(fit_mcmc, pars=c("ms","lp","al","be"))$summary
n_div_mc <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit_mcmc,inc_warmup=FALSE)[[ch]][,"divergent__"])))
n_mt_mc <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit_mcmc,inc_warmup=FALSE)[[ch]][,"treedepth__"] >= 12)))

cat(sprintf("Max Rhat: %.3f | Min n_eff: %.0f | Div: %d | MaxTree: %d\n",
            max(summ_mc[,"Rhat"]), min(summ_mc[,"n_eff"]), n_div_mc, n_mt_mc))

cat("MCMC Recovery:\n")
for(p in names(tv)) {
    ea <- ((q84mc[p]-medmc[p])+(medmc[p]-q16mc[p]))/2
    bias <- (medmc[p]-tv[p])/ea
    st <- ifelse(abs(bias)<1,"OK",ifelse(abs(bias)<2,"marginal","FAIL"))
    cat(sprintf("  %-5s True=%.3f Med=%.3f SE~%.3f Bias=%+.2f sig [%s]\n",
                p, tv[p], medmc[p], ea, bias, st))
}

############################################################
# PLOTS
############################################################

pdf("sim_recovery_v16.pdf", width=14, height=10)
par(mfrow=c(2,3))

xf <- seq(11.5, 16.5, length=500)
br <- seq(11,17,by=0.2)
h <- hist(obs_all, breaks=br, plot=FALSE)
ph <- h$counts/(Vsurvey*0.2); ok <- ph>0

# HMF with MAP
plot(h$mids[ok], log10(ph[ok]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main="v16: MAP + VB + MCMC")
grid(col="gray80")
lines(xf, mrp_log10(xf, tv["ms"], tv["lp"], tv["al"], tv["be"]), col="black", lwd=2, lty=2)
lines(xf, mrp_log10(xf, map_ms, map_lp, map_al, map_be), col="red", lwd=3)
if(!is.null(vb_fit))
    lines(xf, mrp_log10(xf, medvb["ms"], medvb["lp"], medvb["al"], medvb["be"]), col="blue", lwd=2, lty=3)
lines(xf, mrp_log10(xf, medmc["ms"], medmc["lp"], medmc["al"], medmc["be"]), col="purple", lwd=2, lty=4)
legend("topright", c("TRUE","MAP","VB","MCMC"), col=c("black","red","blue","purple"),
       lty=c(2,1,3,4), lwd=c(2,3,2,2), cex=.8)

# Bias comparison
biases <- data.frame(
    param=c("M*","phi*","alpha","beta"),
    MAP=c(map_ms-tv["ms"], map_lp-tv["lp"], map_al-tv["al"], map_be-tv["be"]),
    MCMC=c(medmc["ms"]-tv["ms"], medmc["lp"]-tv["lp"], medmc["al"]-tv["al"], medmc["be"]-tv["be"])
)
if(!is.null(vb_fit))
    biases$VB <- c(medvb["ms"]-tv["ms"], medvb["lp"]-tv["lp"], medvb["al"]-tv["al"], medvb["be"]-tv["be"])

plot(1:4, biases$MAP, pch=19, col="red", cex=2, ylim=c(-1,1),
     xlab="", ylab="Median - True", main="Absolute bias", xaxt="n")
points(1:4+.1, biases$MCMC, pch=17, col="purple", cex=2)
if(!is.null(vb_fit)) points(1:4-.1, biases$VB, pch=15, col="blue", cex=2)
axis(1, at=1:4, labels=biases$param)
abline(h=0, col="gray", lwd=2)
legend("topleft", c("MAP","MCMC","VB"), pch=c(19,17,15), col=c("red","purple","blue"), cex=.9)

# Latent mass recovery (MAP)
plot(all_tm, map_mt, pch=".", col=rgb(0,0,0,.3),
     xlab="True mass", ylab="MAP latent mass", main="Latent mass (MAP)")
abline(0,1,col="red",lwd=2)

# MCMC posteriors
for(p in c("ms","lp","al")) {
    hist(pmmc[,p], breaks=30, col="steelblue", border="white", freq=FALSE,
         main=c(ms="M* (MCMC)",lp="log_phi (MCMC)",al="alpha (MCMC)")[p], xlab=p)
    abline(v=tv[p], col="red", lwd=3)
    abline(v=medmc[p], col="blue", lwd=2, lty=2)
    if(!is.null(vb_fit)) abline(v=medvb[p], col="darkgreen", lwd=2, lty=3)
}

dev.off()
cat("\nSaved: sim_recovery_v16.pdf\n")

cat("\n=== TIMING SUMMARY ===\n")
cat(sprintf("  MAP:  %.0f seconds\n", t_opt))
cat(sprintf("  VB:   %.0f seconds\n", t_vb))
cat(sprintf("  MCMC: %.0f seconds\n", t_mcmc))

cat("\nDone!\n")
