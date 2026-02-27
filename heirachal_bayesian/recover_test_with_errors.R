############################################################
# SIMULATION-RECOVERY TEST v11 (optimized + precomputed indices)
#
# Stage 2 optimizations:
# - Reuse phi(mt) on a fixed grid (pg) instead of recomputing per point.
# - Precompute nearest grid indices ke[i,e] in transformed data (data-only),
#   so no real→int casting in the model block.
# - Small local integration (Ne = 5) and ±4σ window.
# - Reduced grid (Ng = 120) for fast testing.
#
# Also:
# - Fix name merging so Stage 1 medians aren’t NA.
# - Robust medians/quantiles with na.rm=TRUE.
############################################################

suppressPackageStartupMessages({
  library(rstan)
})
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup
############################################################

TRUE_MSTAR <- 14.13; TRUE_LOG_PHI <- -3.96
TRUE_ALPHA <- -1.68; TRUE_BETA <- 0.63
tv <- c(mstar=TRUE_MSTAR, log_phi=TRUE_LOG_PHI,
        alpha=TRUE_ALPHA, beta=TRUE_BETA)

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

cat("=== SETUP ===\n")
cat(sprintf("True: M*=%.2f phi=%.2f a=%.2f b=%.2f\n", tv[1],tv[2],tv[3],tv[4]))
cat(sprintf("Survey: %.2e Mpc^3, %d shells\n\n", Vsurvey, nz))

############################################################
# 2. Generate halos + errors
############################################################

mass_hi <- 16.0
all_tm <- c(); all_z <- c()
for (j in 1:nz) {
  mg <- seq(mlim_sh[j], mass_hi, by=0.001)
  if (length(mg) < 2) next
  pg <- mrp_phi(mg, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
  N_j <- rpois(1, V_sh[j]*sum(pg)*0.001)
  if (N_j == 0) next
  pm_j <- max(pg)*1.1; mj <- numeric(0)
  while (length(mj) < N_j) {
    mp <- runif(N_j*10, mlim_sh[j], mass_hi)
    pp <- mrp_phi(mp, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
    mj <- c(mj, mp[runif(length(mp)) < pp/pm_j])
  }
  all_tm <- c(all_tm, mj[1:N_j])
  all_z  <- c(all_z, runif(N_j, z_edges[j], z_edges[j+1]))
}
N_true <- length(all_tm)
sigma_all <- pmax(0.10, pmin(0.30, 0.30 - 0.06*(all_tm-13)))
obs_all <- all_tm + rnorm(N_true, 0, sigma_all)
mlim_per <- mlim_of_z(all_z)

cat(sprintf("N = %d halos, median sigma = %.2f\n\n", N_true, median(sigma_all)))

############################################################
# 3. DIAGNOSTIC PLOT
############################################################

pdf("mock_diagnostic_v11.pdf", width=12, height=8)
par(mfrow=c(2,2))

br <- seq(11, 17, by=0.2)
h_t <- hist(all_tm, breaks=br, plot=FALSE)
h_o <- hist(obs_all, breaks=br, plot=FALSE)
pt <- h_t$counts/(Vsurvey*0.2); po <- h_o$counts/(Vsurvey*0.2)
xf <- seq(12, 16, length=500)
yf <- mrp_phi(xf, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)

plot(h_t$mids[pt>0], log10(pt[pt>0]), pch=19, col="blue", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main=sprintf("HMF: true(blue) obs(red) N=%d", N_true))
points(h_o$mids[po>0], log10(po[po>0]), pch=17, col="red", cex=1)
lines(xf, log10(yf), col="black", lwd=2)

plot(all_tm, obs_all-all_tm, pch=".", col=rgb(0,0,0,.3),
     xlab="True M", ylab="Obs-True", main="Mass errors")
abline(h=0, col="red", lwd=2)

plot(all_tm, all_z, pch=".", col=rgb(0,0,1,.3),
     xlab="True M", ylab="z", main="Detection boundary")
zp <- seq(zmin,zlimit,length=200); lines(mlim_of_z(zp),zp,col="red",lwd=2)

hist(mlim_per, breaks=40, col="steelblue", border="white",
     xlab="m_lim", main="Per-group mass limits")

dev.off()
cat("Saved: mock_diagnostic_v11.pdf\n\n")

############################################################
# STAGE 1: No errors (validation, pure Poisson)
############################################################

cat("========== STAGE 1: No errors ==========\n")

stan_s1 <- "
data {
  int<lower=1> N; vector[N] x;
  int<lower=1> Nsh; vector[Nsh] V_sh; vector[Nsh] mlim_sh;
  real xhi; int<lower=2> Ng;
}
transformed data {
  real ln10=log(10.0);
  real xlo=min(mlim_sh)-0.5;
  real dx=(xhi-xlo)/(Ng-1.0);
  vector[Ng] xg; for(k in 1:Ng) xg[k]=xlo+(k-1)*dx;
}
parameters { real ms; real lp; real al; real<lower=0.1,upper=2.0> be; }
model {
  ms~normal(14,1.5); lp~normal(-4,2); al~normal(-1.7,1);
  vector[Ng] pg;
  for(k in 1:Ng){real u=xg[k]-ms; pg[k]=be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));}
  vector[Ng] cum; cum[Ng]=0;
  for(kr in 1:(Ng-1)){int k=Ng-kr; cum[k]=cum[k+1]+0.5*(pg[k]+pg[k+1])*dx;}
  real Lam=0;
  for(j in 1:Nsh){int k0=1; for(k in 1:Ng) if(xg[k]<=mlim_sh[j]) k0=k; if(k0>=Ng) k0=Ng-1; Lam+=V_sh[j]*cum[k0];}
  for(i in 1:N){real u=x[i]-ms; real p=be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(p>1e-30) target+=log(p); else target+=-100;}
  target+=-Lam;
}
"

fit1 <- stan(
  model_code = stan_s1,
  data = list(N=N_true, x=all_tm, Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=250),
  chains = 4, iter = 1500, warmup = 750, cores = 4,
  init = lapply(1:4, function(i) list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2), al=rnorm(1,-1.7,.1), be=runif(1,.5,.8))),
  control = list(adapt_delta = 0.95)
)

p1 <- extract(fit1)
pm1 <- cbind(ms=p1$ms, lp=p1$lp, al=p1$al, be=p1$be)

# Ensure names align with posterior parameter names
tvs <- setNames(as.numeric(tv[c("mstar","log_phi","alpha","beta")]),
                c("ms","lp","al","be"))

# Robust summaries
med1  <- apply(pm1, 2, median,   na.rm = TRUE)
q16_1 <- apply(pm1, 2, quantile, probs = 0.16, na.rm = TRUE)
q84_1 <- apply(pm1, 2, quantile, probs = 0.84, na.rm = TRUE)

cat("Stage 1:\n")
for (p in names(tvs)) {
  ea <- ((q84_1[p]-med1[p]) + (med1[p]-q16_1[p]))/2
  bias_sig <- if (is.finite(ea) && ea > 0) (med1[p]-tvs[p]) / ea else NA_real_
  cat(sprintf("  %-5s True=%.3f Med=%s Bias=%s sig\n",
              p, tvs[p],
              ifelse(is.finite(med1[p]), sprintf("%.3f", med1[p]), "NA"),
              ifelse(is.finite(bias_sig), sprintf("%+.2f", bias_sig), "NA")))
}

############################################################
# STAGE 2: With errors (FAST)
# Conditional + Poisson(Lambda_conv), precomputed indices
############################################################

cat("\n========== STAGE 2: With errors ==========\n")

# Per-shell mean sigma for Lambda_conv
shell_mean_sigma <- rep(0.25, nz)  # default
for (j in 1:nz) {
  in_shell <- all_z >= z_edges[j] & all_z < z_edges[j+1]
  if (sum(in_shell) > 0) shell_mean_sigma[j] <- mean(sigma_all[in_shell])
}

stan_s2 <- "
functions {
  real normal_cdf_approx(real x) {
    if (x > 5)  return 1.0;
    if (x < -5) return 0.0;
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
  }
}
data {
  int<lower=1> N;
  vector[N] x;
  vector<lower=0>[N] sig;
  vector[N] mlim_i;

  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  vector<lower=0>[Nsh] sig_sh;

  real xhi;
  int<lower=2> Ng;     // grid size
  int<lower=3> Ne;     // local integration pts
}
transformed data {
  real ln10 = log(10.0);
  real s2pi = sqrt(2.0*pi());

  // grid boundaries
  real xlo = min(mlim_sh) - 1.5;
  real dx  = (xhi - xlo) / (Ng - 1);

  vector[Ng] xg;
  for (k in 1:Ng)
    xg[k] = xlo + (k-1) * dx;

  // Precompute nearest grid indices for each halo & eval point
  int ke[N,Ne];
  for (i in 1:N) {
    real si = sig[i];
    real lo_n = fmax(x[i] - 4*si, xlo);
    real hi_n = fmin(x[i] + 4*si, xhi);
    if (hi_n > lo_n) {
      real dxe = (hi_n - lo_n) / (Ne - 1);
      for (e in 1:Ne) {
        real mt = lo_n + (e - 1) * dxe;
        int k = 1 + to_int( (mt - xlo) / dx );
        if (k < 1)  k = 1;
        if (k > Ng) k = Ng;
        ke[i,e] = k;
      }
    } else {
      for (e in 1:Ne) ke[i,e] = 1;
    }
  }
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1,upper=2.0> be;
}
model {
  ms ~ normal(14, 1.5);
  lp ~ normal(-4, 2);
  al ~ normal(-1.7, 1);

  // Precompute phi(mt) on grid
  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
  }

  // === POISSON TERM: Lambda_conv ===
  real Lam = 0;
  for (j in 1:Nsh) {
    real sj = sig_sh[j];
    real sum_j = 0;

    for (k in 1:Ng)
      sum_j += pg[k] * normal_cdf_approx((xg[k] - mlim_sh[j]) / sj);

    // trapezoid endpoints
    sum_j -= 0.5 * pg[1]  * normal_cdf_approx((xg[1]  - mlim_sh[j]) / sj);
    sum_j -= 0.5 * pg[Ng] * normal_cdf_approx((xg[Ng] - mlim_sh[j]) / sj);

    Lam += V_sh[j] * sum_j * dx;
  }
  if (Lam > 0)
    target += N * log(Lam) - Lam;

  // === CONDITIONAL LIKELIHOOD (FAST, uses precomputed ke) ===
  for (i in 1:N) {
    real si = sig[i];

    real lo_n = fmax(x[i] - 4*si, xlo);
    real hi_n = fmin(x[i] + 4*si, xhi);

    real numerator = 0;

    if (hi_n > lo_n) {
      real dxe = (hi_n - lo_n) / (Ne - 1);

      for (e in 1:Ne) {
        int k = ke[i,e];        // precomputed nearest grid index
        real mt = lo_n + (e - 1) * dxe;
        real z  = (x[i] - mt) / si;
        numerator += pg[k] * exp(-0.5 * z*z) / (si * s2pi);
      }
      numerator *= dxe;
    } else {
      numerator = 1e-30;
    }

    // DENOMINATOR: ∫ phi(mt) Phi((mt-mlim_i)/si) dmt on grid
    real denominator = 0;
    for (k in 1:Ng)
      denominator += pg[k] * normal_cdf_approx((xg[k] - mlim_i[i]) / si);

    denominator -= 0.5 * pg[1]  * normal_cdf_approx((xg[1]  - mlim_i[i]) / si);
    denominator -= 0.5 * pg[Ng] * normal_cdf_approx((xg[Ng] - mlim_i[i]) / si);
    denominator *= dx;

    if (numerator > 1e-30 && denominator > 1e-30)
      target += log(numerator) - log(denominator);
    else
      target += -100;
  }
}
"

Ne_fast <- 5
Ng_fast <- 120
cat(sprintf("N=%d, Ne=%d, Ng=%d\n\n", N_true, Ne_fast, Ng_fast))

fit2 <- stan(
  model_code = stan_s2,
  data = list(
    N = N_true,
    x = obs_all,
    sig = sigma_all,
    mlim_i = mlim_per,
    Nsh = nz,
    V_sh = V_sh,
    mlim_sh = mlim_sh,
    sig_sh = shell_mean_sigma,
    xhi = 16.5,
    Ng = Ng_fast,
    Ne = Ne_fast
  ),
  chains = 4, iter = 1500, warmup = 750, cores = 4,
  init = lapply(1:4, function(i) list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
                                      al=rnorm(1,-1.7,.1), be=runif(1,.5,.8))),
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\n=== STAGE 2 RESULTS ===\n")
print(fit2, pars=c("ms","lp","al","be"))

p2 <- extract(fit2); pm2 <- cbind(ms=p2$ms, lp=p2$lp, al=p2$al, be=p2$be)
med2  <- apply(pm2, 2, median,   na.rm = TRUE)
q16_2 <- apply(pm2, 2, quantile, probs = 0.16, na.rm = TRUE)
q84_2 <- apply(pm2, 2, quantile, probs = 0.84, na.rm = TRUE)

summ2 <- summary(fit2, pars=c("ms","lp","al","be"))$summary
n_div2 <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit2,inc_warmup=FALSE)[[ch]][,"divergent__"])))
cat(sprintf("Max Rhat: %.3f | Min n_eff: %.0f | Div: %d\n",
            max(summ2[,"Rhat"]), min(summ2[,"n_eff"]), n_div2))

cat("\nStage 2 Recovery:\n")
cat(sprintf("%-10s %8s %8s %10s\n","Param","True","Median","Bias"))
cat(strrep("-",42),"\n")
for (p in names(tvs)) {
  ea <- ((q84_2[p]-med2[p]) + (med2[p]-q16_2[p]))/2
  bias <- if (is.finite(ea) && ea > 0) (med2[p]-tvs[p]) / ea else NA_real_
  st <- if (is.finite(bias)) ifelse(abs(bias)<1,"OK", ifelse(abs(bias)<2,"marginal","FAIL")) else "NA"
  cat(sprintf("%-10s %8.3f %8s  %s sig [%s]\n",
              p, tvs[p],
              ifelse(is.finite(med2[p]), sprintf("%.3f", med2[p]), "NA"),
              ifelse(is.finite(bias), sprintf("%+.2f", bias), "NA"),
              st))
}

############################################################
# PLOTS
############################################################

pdf("sim_recovery_v11.pdf", width=14, height=10)
par(mfrow=c(2,3))

br <- seq(11,17,by=0.2)
h <- hist(obs_all, breaks=br, plot=FALSE)
ph <- h$counts/(Vsurvey*0.2); ok <- ph>0
xf <- seq(11.5,16.5,length=500)

plot(h$mids[ok], log10(ph[ok]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main=sprintf("v11 (N=%d)", N_true))
grid(col="gray80")
idx <- sample(nrow(pm2), min(200,nrow(pm2)))
for (i in idx) {
  y <- tryCatch(mrp_log10(xf, pm2[i,"ms"], pm2[i,"lp"], pm2[i,"al"], pm2[i,"be"]),
                error=function(e) rep(NA_real_, length(xf)))
  if (all(is.finite(y))) lines(xf, y, col=rgb(0,0,1,.03))
}
lines(xf, mrp_log10(xf, med2["ms"], med2["lp"], med2["al"], med2["be"]), col="red", lwd=3)
lines(xf, mrp_log10(xf, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA), col="black", lwd=2, lty=2)
legend("topright", c("Data","Fit","TRUE"), col=c("darkgreen","red","black"),
       pch=c(19,NA,NA), lty=c(NA,1,2), lwd=c(NA,3,2), cex=.85)

for (p in c("ms","lp","al","be")) {
  hist(pm2[,p], breaks=50, main=c(ms="M*",lp="log_phi",al="alpha",be="beta")[p],
       xlab=p, col="steelblue", border="white", freq=FALSE)
  abline(v=tvs[p], col="red", lwd=3)
  abline(v=med2[p], col="blue", lwd=2, lty=2)
  legend("topright",
         c(sprintf("True: %.3f", tvs[p]),
           sprintf("Med: %.3f", med2[p])),
         col=c("red","blue"), lty=c(1,2), lwd=c(3,2), cex=.8)
}

plot(pm2[,"ms"], pm2[,"lp"], pch=".", col=rgb(0,0,0,.1),
     xlab="M*", ylab="log_phi", main="Degeneracy")
points(TRUE_MSTAR, TRUE_LOG_PHI, pch=4, col="red", cex=3, lwd=3)
points(med2["ms"], med2["lp"], pch=4, col="blue", cex=2, lwd=2)

dev.off()
cat("\nSaved: sim_recovery_v11.pdf, mock_diagnostic_v11.pdf\n")
cat("Done!\n")
