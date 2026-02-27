############################################################
# SIMULATION-RECOVERY TEST v8
# 
# Back to basics. Pure inhomogeneous Poisson process:
#   log L = sum_i log phi(x_i) - sum_j V_j * int(phi, mlim_j, xhi)
#
# NO per-group denominator. NO error convolution (Stage 1).
# 
# ALSO: verify the likelihood in R before running Stan,
# to catch any implementation bugs.
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. TRUE PARAMETERS
############################################################

TRUE_MSTAR   <- 14.13
TRUE_LOG_PHI <- -3.96
TRUE_ALPHA   <- -1.68
TRUE_BETA    <- 0.63

tv <- c(mstar=TRUE_MSTAR, log_phi=TRUE_LOG_PHI,
        alpha=TRUE_ALPHA, beta=TRUE_BETA)

cat("=== TRUE PARAMETERS ===\n")
cat(sprintf("M*=%.2f  log_phi=%.2f  alpha=%.2f  beta=%.2f\n\n",
            TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA))

############################################################
# 2. MRP
############################################################

mrp_phi <- function(x, mstar, log_phi, alpha, beta) {
    u <- x - mstar
    beta * log(10) * 10^log_phi * 10^((alpha + 1) * u) * exp(-10^(beta * u))
}

mrp_log10 <- function(x, mstar, log_phi, alpha, beta) {
    log10(pmax(mrp_phi(x, mstar, log_phi, alpha, beta), 1e-30))
}

############################################################
# 3. Survey
############################################################

ho <- 67.37; omegam <- 0.3147; zlimit <- 0.25; zmin <- 0.01
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

sky_frac <- 60 * (pi/180)^2 / (4*pi)
nz_shells <- 20
z_edges <- seq(zmin, zlimit, length.out = nz_shells + 1)
z_mids  <- (z_edges[-1] + z_edges[-(nz_shells+1)]) / 2
d_edges <- sapply(z_edges, cosdist)
V_shells <- (4/3)*pi*(d_edges[-1]^3 - d_edges[-(nz_shells+1)]^3) * sky_frac
Vsurvey  <- sum(V_shells)

mlim_of_z <- function(z) 12.0 + 2.5*(z - zmin)/(zlimit - zmin)
mlim_shells <- mlim_of_z(z_mids)

############################################################
# 4. Generate halos (Poisson per shell)
############################################################

mass_hi <- 16.0
all_tm <- c(); all_z <- c(); all_shell <- c()
for(j in 1:nz_shells) {
    mg <- seq(mlim_shells[j], mass_hi, by=0.001)
    if(length(mg) < 2) next
    pg <- mrp_phi(mg, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
    N_j <- rpois(1, V_shells[j] * sum(pg) * 0.001)
    if(N_j == 0) next
    pm_j <- max(pg) * 1.1; mj <- numeric(0)
    while(length(mj) < N_j) {
        mp <- runif(N_j*10, mlim_shells[j], mass_hi)
        pp <- mrp_phi(mp, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
        mj <- c(mj, mp[runif(length(mp)) < pp/pm_j])
    }
    mj <- mj[1:N_j]
    all_tm <- c(all_tm, mj)
    all_z <- c(all_z, runif(N_j, z_edges[j], z_edges[j+1]))
    all_shell <- c(all_shell, rep(j, N_j))
}

# Use true masses directly (no errors, no extra cuts)
x_obs <- all_tm
N <- length(x_obs)

cat(sprintf("Generated N = %d halos\n", N))
cat(sprintf("Mass range: [%.2f, %.2f]\n\n", min(x_obs), max(x_obs)))

############################################################
# 5. VERIFY LIKELIHOOD IN R
#    Compute log L at true parameters and at some wrong ones
############################################################

loglik_R <- function(pars, x, V_shell, mlim_shell, xhi=16.5, Ng=500) {
    mstar <- pars[1]; log_phi <- pars[2]
    alpha <- pars[3]; beta <- pars[4]
    
    # Term 1: sum log phi(x_i)
    phi_vals <- mrp_phi(x, mstar, log_phi, alpha, beta)
    if(any(phi_vals <= 0)) return(-1e10)
    term1 <- sum(log(phi_vals))
    
    # Term 2: Lambda = sum_j V_j * integral(phi, mlim_j, xhi)
    Lambda <- 0
    for(j in seq_along(V_shell)) {
        mg <- seq(mlim_shell[j], xhi, length.out=Ng)
        pg <- mrp_phi(mg, mstar, log_phi, alpha, beta)
        dm <- (xhi - mlim_shell[j]) / (Ng - 1)
        integral_j <- (pg[1]/2 + sum(pg[2:(Ng-1)]) + pg[Ng]/2) * dm
        Lambda <- Lambda + V_shell[j] * integral_j
    }
    
    # Poisson process: sum log phi - Lambda
    return(term1 - Lambda)
}

cat("=== LIKELIHOOD VERIFICATION ===\n")
ll_true <- loglik_R(c(TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA),
                     x_obs, V_shells, mlim_shells)
cat(sprintf("log L at TRUE params:     %.2f\n", ll_true))

# Try some perturbations
perturbs <- list(
    c(TRUE_MSTAR+0.5, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA),
    c(TRUE_MSTAR-0.5, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA),
    c(TRUE_MSTAR, TRUE_LOG_PHI+0.5, TRUE_ALPHA, TRUE_BETA),
    c(TRUE_MSTAR, TRUE_LOG_PHI-0.5, TRUE_ALPHA, TRUE_BETA),
    c(TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA+0.5, TRUE_BETA),
    c(TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA-0.5, TRUE_BETA),
    c(TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA+0.2),
    c(TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA-0.2)
)
labels <- c("M*+0.5","M*-0.5","phi+0.5","phi-0.5",
            "a+0.5","a-0.5","b+0.2","b-0.2")

for(i in seq_along(perturbs)) {
    ll_p <- loglik_R(perturbs[[i]], x_obs, V_shells, mlim_shells)
    cat(sprintf("  %s: %.2f  (delta = %+.2f)\n", 
                labels[i], ll_p, ll_p - ll_true))
}

# Grid search around true values to check for maximum
cat("\nGrid search for log L maximum:\n")
best_ll <- -Inf; best_pars <- tv
for(dm in seq(-1, 1, by=0.25)) {
    for(dp in seq(-1, 1, by=0.25)) {
        for(da in seq(-0.5, 0.5, by=0.25)) {
            for(db in seq(-0.2, 0.2, by=0.1)) {
                pars <- c(TRUE_MSTAR+dm, TRUE_LOG_PHI+dp,
                          TRUE_ALPHA+da, TRUE_BETA+db)
                if(pars[4] < 0.1) next
                ll <- loglik_R(pars, x_obs, V_shells, mlim_shells)
                if(ll > best_ll) {
                    best_ll <- ll; best_pars <- pars
                }
            }
        }
    }
}
cat(sprintf("  True:  M*=%.2f  phi=%.2f  a=%.2f  b=%.2f  => ll=%.2f\n",
            TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA, ll_true))
cat(sprintf("  Best:  M*=%.2f  phi=%.2f  a=%.2f  b=%.2f  => ll=%.2f\n",
            best_pars[1], best_pars[2], best_pars[3], best_pars[4], best_ll))
cat(sprintf("  Delta: %+.2f\n\n", best_ll - ll_true))

############################################################
# 6. STAN MODEL - Pure Poisson process
############################################################

stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  int<lower=1> Nshells;
  vector[Nshells] V_shell;
  vector[Nshells] mlim_shell;
  real xhi;
  int<lower=2> Ng;
}

transformed data {
  real ln10 = log(10.0);
  real xlo_g = min(mlim_shell) - 0.5;
  real dx = (xhi - xlo_g) / (Ng - 1.0);
  vector[Ng] xgrid;
  for(k in 1:Ng) xgrid[k] = xlo_g + (k-1)*dx;
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.1, upper=2.0> beta;
}

model {
  mstar   ~ normal(14.0, 1.5);
  log_phi ~ normal(-4.0, 2.0);
  alpha   ~ normal(-1.7, 1.0);
  
  // MRP on grid
  vector[Ng] phi_g;
  for(k in 1:Ng) {
    real u = xgrid[k] - mstar;
    phi_g[k] = beta * ln10 * pow(10, log_phi)
               * pow(10, (alpha+1)*u) * exp(-pow(10, beta*u));
  }
  
  // Cumulative from right
  vector[Ng] cum;
  cum[Ng] = 0;
  for(kr in 1:(Ng-1)) {
    int k = Ng - kr;
    cum[k] = cum[k+1] + 0.5*(phi_g[k] + phi_g[k+1])*dx;
  }
  
  // Lambda = sum_j V_j * integral(phi above mlim_j)
  real Lambda = 0;
  for(j in 1:Nshells) {
    int k0 = 1;
    for(k in 1:Ng) if(xgrid[k] <= mlim_shell[j]) k0 = k;
    if(k0 >= Ng) k0 = Ng-1;
    Lambda += V_shell[j] * cum[k0];
  }
  
  // Pure Poisson process: sum log phi(x_i) - Lambda
  for(i in 1:N) {
    real u = x[i] - mstar;
    real phi_i = beta * ln10 * pow(10, log_phi)
                 * pow(10, (alpha+1)*u) * exp(-pow(10, beta*u));
    if(phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
  target += -Lambda;
}
"

############################################################
# 7. RUN
############################################################

cat("=== RUNNING STAN ===\n")
cat(sprintf("N=%d, %d shells, Ng=300\n\n", N, nz_shells))

fit <- stan(
    model_code = stan_code,
    data = list(N=N, x=x_obs, Nshells=nz_shells,
                V_shell=V_shells, mlim_shell=mlim_shells,
                xhi=16.5, Ng=300),
    chains=4, iter=2000, warmup=1000, cores=4,
    init = lapply(1:4, function(i) list(
        mstar=rnorm(1,14,0.2), log_phi=rnorm(1,-4,0.2),
        alpha=rnorm(1,-1.7,0.1), beta=runif(1,0.5,0.8))),
    control = list(adapt_delta=0.95, max_treedepth=12)
)

############################################################
# 8. RESULTS
############################################################

cat("\n=== STAN SUMMARY ===\n")
print(fit, pars=c("mstar","log_phi","alpha","beta"))

post <- extract(fit, pars=c("mstar","log_phi","alpha","beta"))
pm <- cbind(mstar=post$mstar, log_phi=post$log_phi,
            alpha=post$alpha, beta=post$beta)
med <- apply(pm, 2, median)
q16 <- apply(pm, 2, quantile, 0.16)
q84 <- apply(pm, 2, quantile, 0.84)

summ <- summary(fit, pars=c("mstar","log_phi","alpha","beta"))$summary
n_div <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit,inc_warmup=FALSE)[[ch]][,"divergent__"])))
cat(sprintf("\nMax Rhat: %.3f | Min n_eff: %.0f | Divergences: %d\n",
            max(summ[,"Rhat"]), min(summ[,"n_eff"]), n_div))

cat("\n=== RECOVERY ===\n")
cat(sprintf("%-10s %8s %8s %8s %8s %10s\n",
            "Param","True","Median","+1sig","-1sig","Bias"))
cat(strrep("-",62),"\n")
for(p in names(tv)) {
    eu <- q84[p]-med[p]; ed <- med[p]-q16[p]; ea <- (eu+ed)/2
    bias <- (med[p]-tv[p])/ea
    st <- ifelse(abs(bias)<1,"OK",ifelse(abs(bias)<2,"marginal","FAIL"))
    cat(sprintf("%-10s %8.3f %8.3f  +%.3f  -%.3f  %+.2f sig  [%s]\n",
                p, tv[p], med[p], eu, ed, bias, st))
}

cat("\nCorrelations:\n")
print(round(cor(pm), 3))

############################################################
# 9. PLOTS
############################################################

pdf("sim_recovery_v8.pdf", width=14, height=10)
par(mfrow=c(2,3))

br <- seq(12, 16, by=0.2)
h <- hist(x_obs, breaks=br, plot=FALSE)
phi_h <- h$counts/(Vsurvey*0.2); ok <- phi_h>0
xf <- seq(11.5, 16.5, length=500)

plot(h$mids[ok], log10(phi_h[ok]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab="log10(M)", ylab="log10(phi)",
     main=sprintf("v8: Poisson process, no errors (N=%d)", N))
grid(col="gray80")
idx <- sample(nrow(pm), min(300, nrow(pm)))
for(i in idx) {
    y <- tryCatch(mrp_log10(xf, pm[i,"mstar"], pm[i,"log_phi"],
                             pm[i,"alpha"], pm[i,"beta"]),
                  error=function(e) rep(NA, length(xf)))
    if(all(is.finite(y))) lines(xf, y, col=rgb(0,0,1,0.02))
}
lines(xf, mrp_log10(xf, med["mstar"], med["log_phi"],
                      med["alpha"], med["beta"]), col="red", lwd=3)
lines(xf, mrp_log10(xf, TRUE_MSTAR, TRUE_LOG_PHI,
                      TRUE_ALPHA, TRUE_BETA), col="black", lwd=2, lty=2)
legend("topright", c("Data","Fit","TRUE"), col=c("darkgreen","red","black"),
       pch=c(19,NA,NA), lty=c(NA,1,2), lwd=c(NA,3,2), cex=0.85)

for(p in c("mstar","log_phi","alpha","beta")) {
    hist(pm[,p], breaks=50, main=p, xlab=p,
         col="steelblue", border="white", freq=FALSE)
    abline(v=tv[p], col="red", lwd=3)
    abline(v=med[p], col="blue", lwd=2, lty=2)
    legend("topright", c(sprintf("True: %.3f",tv[p]),
                          sprintf("Med: %.3f",med[p])),
           col=c("red","blue"), lty=c(1,2), lwd=c(3,2), cex=0.8)
}

plot(pm[,"mstar"], pm[,"log_phi"], pch=".", col=rgb(0,0,0,.1),
     xlab="M*", ylab="log_phi", main="M* vs log_phi")
points(TRUE_MSTAR, TRUE_LOG_PHI, pch=4, col="red", cex=3, lwd=3)
points(med["mstar"], med["log_phi"], pch=4, col="blue", cex=2, lwd=2)

dev.off()
cat("\nSaved: sim_recovery_v8.pdf\n")

saveRDS(list(true=tv, post=pm, med=med, q16=q16, q84=q84, fit=fit,
             data=list(x=x_obs), shells=list(V=V_shells, mlim=mlim_shells)),
        "sim_recovery_v8.rds")
cat("Done!\n")
