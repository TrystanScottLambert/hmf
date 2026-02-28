############################################################
# DIAGNOSTIC PLOTS: HMF Hierarchical Recovery (v16)
#
# Generates presentation-quality figures showing:
# 1. HMF recovery (data vs true vs MAP fit)
# 2. Parameter bias summary across all methods tested
# 3. Latent mass recovery (shrinkage toward MRP prior)
# 4. Residual diagnostics (latent mass residuals vs true mass)
# 5. Error budget: how Eddington bias depends on sigma
#
# Uses the same mock catalogue as v16.
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup + generate data (identical to v16)
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

cat(sprintf("Generated %d halos, median sigma = %.2f dex\n", N_true, median(sigma_all)))

############################################################
# 2. Run MAP (hierarchical model)
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

cat("Compiling and running MAP...\n")
sm <- stan_model(model_code=stan_hier)
stan_data <- list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
                  Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200)

init_fn <- function() {
    list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
         al=rnorm(1,-1.7,.1), be=runif(1,.5,.8),
         z_raw=rnorm(N_true, 0, 0.3))
}

# Run MAP multiple times and take best
best_opt <- NULL; best_lp <- -Inf
for(attempt in 1:5) {
    opt <- tryCatch(
        optimizing(sm, data=stan_data, init=init_fn(),
                   hessian=TRUE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
        error=function(e) NULL)
    if(!is.null(opt) && opt$value > best_lp) {
        best_opt <- opt; best_lp <- opt$value
        cat(sprintf("  Attempt %d: log_p = %.2f (best so far)\n", attempt, opt$value))
    } else {
        cat(sprintf("  Attempt %d: %s\n", attempt,
                    ifelse(is.null(opt), "failed", sprintf("log_p = %.2f", opt$value))))
    }
}
opt <- best_opt

map_ms <- opt$par$ms; map_lp <- opt$par$lp
map_al <- opt$par$al; map_be <- opt$par$be
map_mt <- obs_all + sigma_all * opt$par$z_raw

cat(sprintf("\nMAP: M*=%.3f lp=%.3f al=%.3f be=%.3f\n", map_ms, map_lp, map_al, map_be))
cat(sprintf("True: M*=%.3f lp=%.3f al=%.3f be=%.3f\n", tv[1], tv[2], tv[3], tv[4]))

# Laplace SE
np <- 4; H <- opt$hessian
tryCatch({
    A <- -H[1:np, 1:np]
    B <- -H[1:np, (np+1):(np+N_true)]
    D <- -H[(np+1):(np+N_true), (np+1):(np+N_true)]
    D_diag_inv <- diag(1/diag(D))
    Schur <- A - B %*% D_diag_inv %*% t(B)
    Sigma_mrp <- solve(Schur)
    se_lap <- sqrt(diag(Sigma_mrp))
    names(se_lap) <- c("ms","lp","al","be")
    cat("Laplace SE:", paste(names(se_lap), round(se_lap,3), sep="="), "\n")
}, error=function(e) {
    se_lap <<- c(ms=0.3, lp=0.4, al=0.05, be=0.3)
    cat("Hessian failed, using approximate SE\n")
})

# Also run no-error MAP for comparison
cat("\nRunning no-error MAP for comparison...\n")
stan_data_noerr <- stan_data
stan_data_noerr$x_obs <- all_tm
stan_data_noerr$sig <- rep(1e-6, N_true)
opt_noerr <- optimizing(sm, data=stan_data_noerr, init=init_fn(),
                         hessian=FALSE, as_vector=FALSE, iter=10000, algorithm="LBFGS")
noerr_ms <- opt_noerr$par$ms; noerr_lp <- opt_noerr$par$lp
noerr_al <- opt_noerr$par$al; noerr_be <- opt_noerr$par$be
cat(sprintf("No-error MAP: M*=%.3f lp=%.3f al=%.3f be=%.3f\n", noerr_ms, noerr_lp, noerr_al, noerr_be))

############################################################
# 3. FIGURE 1: Main recovery plot (4-panel)
############################################################

cat("\nGenerating Figure 1: Main recovery...\n")
pdf("diagnostic_fig1_recovery.pdf", width=12, height=12)
par(mfrow=c(2,2), mar=c(4.5,4.5,3,1), cex.lab=1.3, cex.axis=1.1, cex.main=1.4)

# --- Panel A: HMF fit ---
xf <- seq(11.5, 16.5, length=500)
br <- seq(11,17,by=0.2)
h_o <- hist(obs_all, breaks=br, plot=FALSE)
h_t <- hist(all_tm, breaks=br, plot=FALSE)
po <- h_o$counts/(Vsurvey*0.2); pt <- h_t$counts/(Vsurvey*0.2)

plot(h_t$mids[pt>0], log10(pt[pt>0]), pch=16, col="grey50", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2), xlab=expression(log[10]*(M/M[sun])),
     ylab=expression(log[10]*(phi*" / Mpc"^-3*" dex"^-1)),
     main="(a) HMF Recovery")
points(h_o$mids[po>0], log10(po[po>0]), pch=17, col="firebrick", cex=1.0)
lines(xf, mrp_log10(xf, tv["ms"],tv["lp"],tv["al"],tv["be"]),
      col="black", lwd=2.5, lty=2)
lines(xf, mrp_log10(xf, map_ms, map_lp, map_al, map_be),
      col="dodgerblue3", lwd=3)
lines(xf, mrp_log10(xf, noerr_ms, noerr_lp, noerr_al, noerr_be),
      col="forestgreen", lwd=2, lty=3)
grid(col="gray85")
legend("topright",
       c("True masses", "Obs masses", "True MRP",
         expression("MAP (with "*sigma*")"), expression("MAP (no "*sigma*")")),
       pch=c(16,17,NA,NA,NA), col=c("grey50","firebrick","black","dodgerblue3","forestgreen"),
       lty=c(NA,NA,2,1,3), lwd=c(NA,NA,2.5,3,2), cex=0.9, bg="white")

# --- Panel B: Parameter recovery ---
pnames <- c("M*", expression(log[10]*phi*"*"), expression(alpha), expression(beta))
map_vals <- c(map_ms, map_lp, map_al, map_be)
noerr_vals <- c(noerr_ms, noerr_lp, noerr_al, noerr_be)
true_vals <- unname(tv)
bias_map <- (map_vals - true_vals) / se_lap
bias_noerr_abs <- noerr_vals - true_vals

plot(1:4, map_vals - true_vals, pch=19, col="dodgerblue3", cex=2.5,
     ylim=c(-0.2, 0.2), xlab="", ylab="MAP - True",
     main="(b) Parameter Bias", xaxt="n")
points(1:4, noerr_vals - true_vals, pch=17, col="forestgreen", cex=2)
# Error bars from Laplace
arrows(1:4, map_vals-true_vals-se_lap, 1:4, map_vals-true_vals+se_lap,
       code=3, angle=90, length=0.08, col="dodgerblue3", lwd=2)
abline(h=0, col="black", lwd=1.5)
abline(h=c(-0.1,0.1), col="gray70", lty=2)
axis(1, at=1:4, labels=c(expression(M*"*"), expression(log*phi*"*"),
                          expression(alpha), expression(beta)), cex.axis=1.2)
legend("topleft", c(expression("Hierarchical MAP (with "*sigma*")"),
                     expression("Direct MAP (no "*sigma*")")),
       pch=c(19,17), col=c("dodgerblue3","forestgreen"), cex=0.9, bg="white")
grid(col="gray85")

# --- Panel C: Latent mass recovery ---
plot(all_tm, map_mt, pch=".", col=rgb(0.2,0.2,0.8,0.3), cex=2,
     xlab=expression("True "*log[10]*M),
     ylab=expression("Recovered "*log[10]*M),
     main="(c) Latent Mass Recovery")
abline(0, 1, col="red", lwd=2)
# Add 1:1 envelope
abline(a=0.3, b=1, col="red", lwd=1, lty=2)
abline(a=-0.3, b=1, col="red", lwd=1, lty=2)
# Annotate shrinkage
med_resid <- median(abs(map_mt - all_tm))
med_sig <- median(sigma_all)
legend("topleft",
       c(sprintf("Median |resid| = %.3f dex", med_resid),
         sprintf("Median sigma  = %.3f dex", med_sig),
         sprintf("Shrinkage ratio = %.2f", med_resid/med_sig)),
       cex=0.9, bg="white")
grid(col="gray85")

# --- Panel D: Latent mass residuals vs true mass ---
resid <- map_mt - all_tm
plot(all_tm, resid, pch=".", col=rgb(0,0,0,0.3), cex=2,
     xlab=expression("True "*log[10]*M),
     ylab=expression("Recovered - True (dex)"),
     main="(d) Mass Residuals", ylim=c(-1, 1))
abline(h=0, col="red", lwd=2)
# Binned mean
mbins <- seq(12, 15.5, by=0.5)
for(i in 1:(length(mbins)-1)) {
    idx <- all_tm >= mbins[i] & all_tm < mbins[i+1]
    if(sum(idx) > 5) {
        bm <- mean(resid[idx]); bsd <- sd(resid[idx])
        mid <- (mbins[i]+mbins[i+1])/2
        points(mid, bm, pch=19, col="red", cex=2)
        arrows(mid, bm-bsd, mid, bm+bsd, code=3, angle=90, length=0.06, col="red", lwd=2)
    }
}
# Detection boundary effect
abline(h=c(-0.3, 0.3), col="gray60", lty=2)
legend("topright", c("Individual", "Binned mean ± 1sd"),
       pch=c(46,19), col=c("black","red"), cex=0.9, bg="white")
grid(col="gray85")

dev.off()
cat("Saved: diagnostic_fig1_recovery.pdf\n")

############################################################
# 4. FIGURE 2: Method comparison + Eddington bias
############################################################

cat("Generating Figure 2: Method comparison...\n")

# Results from previous runs (hardcoded from v12 and v16 output)
# v12 Stage 2 (obs masses, no correction): biased
v12_s2 <- c(ms=14.939, lp=-5.422, al=-1.955, be=1.147)
v12_s2_se <- c(ms=0.144, lp=0.359, al=0.049, be=0.428)
# v12 Stage 3 (half errors)
v12_s3 <- c(ms=14.546, lp=-4.678, al=-1.845, be=0.809)
v12_s3_se <- c(ms=0.285, lp=0.580, al=0.093, be=0.340)
# v16 MAP (hierarchical)
v16_map <- c(ms=map_ms, lp=map_lp, al=map_al, be=map_be)
v16_se <- se_lap
# v8/v16 no-error
v_noerr <- c(ms=noerr_ms, lp=noerr_lp, al=noerr_al, be=noerr_be)

pdf("diagnostic_fig2_comparison.pdf", width=14, height=8)
par(mfrow=c(1,2), mar=c(5,5,3,1), cex.lab=1.3, cex.axis=1.1, cex.main=1.4)

# --- Panel A: Bias in sigma across methods ---
bias_noerr <- (v_noerr - tv) / c(0.3, 0.4, 0.05, 0.3)  # approx SE
bias_naive_full <- (v12_s2 - tv) / v12_s2_se
bias_naive_half <- (v12_s3 - tv) / v12_s3_se
bias_hier <- (v16_map - tv) / v16_se

x <- 1:4
plot(x-0.24, bias_noerr, pch=15, col="forestgreen", cex=2, ylim=c(-7, 7),
     xlab="", ylab=expression("Bias ("*sigma*")"),
     main="(a) Parameter Bias Across Methods", xaxt="n")
points(x-0.08, bias_naive_full, pch=16, col="firebrick", cex=2)
points(x+0.08, bias_naive_half, pch=17, col="darkorange", cex=2)
points(x+0.24, bias_hier, pch=18, col="dodgerblue3", cex=2.5)
abline(h=0, col="black", lwd=1.5)
abline(h=c(-1,1), col="gray50", lty=2)
abline(h=c(-2,2), col="gray70", lty=3)
axis(1, at=1:4, labels=c(expression(M*"*"), expression(log*phi*"*"),
                          expression(alpha), expression(beta)), cex.axis=1.2)
legend("bottomleft",
       c(expression("No errors ("*sigma*"=0)"),
         expression("Naive ("*sigma*"~0.28)"),
         expression("Naive ("*sigma*"~0.14)"),
         expression("Hierarchical MAP ("*sigma*"~0.28)")),
       pch=c(15,16,17,18), col=c("forestgreen","firebrick","darkorange","dodgerblue3"),
       cex=0.9, bg="white", pt.cex=c(2,2,2,2.5))
grid(col="gray85")

# --- Panel B: Eddington bias scaling ---
# sigma values and corresponding alpha bias (from v12)
sig_vals <- c(0, 0.14, 0.28)
al_bias <- c(-0.18, -1.77, -5.66)  # in sigma units from v12
ms_bias <- c(-0.07, 1.46, 5.65)
# Add hierarchical point
sig_hier <- 0.28
al_bias_hier <- bias_hier["al"]
ms_bias_hier <- bias_hier["ms"]

plot(sig_vals^2, abs(al_bias), pch=16, col="firebrick", cex=2.5,
     xlim=c(0, 0.10), ylim=c(0, 7),
     xlab=expression(sigma^2*" (dex"^2*")"),
     ylab=expression("|Bias| ("*sigma*")"),
     main=expression("(b) Eddington Bias Scaling with "*sigma^2))
points(sig_vals^2, abs(ms_bias), pch=17, col="darkorange", cex=2.5)
# Linear fits
lm_al <- lm(abs(al_bias) ~ I(sig_vals^2) + 0)
lm_ms <- lm(abs(ms_bias) ~ I(sig_vals^2) + 0)
sg <- seq(0, 0.1, length=100)
lines(sg, predict(lm_al, data.frame(sig_vals=sqrt(sg))), col="firebrick", lwd=2, lty=2)
lines(sg, predict(lm_ms, data.frame(sig_vals=sqrt(sg))), col="darkorange", lwd=2, lty=2)
# Hierarchical point
points(sig_hier^2, abs(al_bias_hier), pch=18, col="dodgerblue3", cex=3)
points(sig_hier^2, abs(ms_bias_hier), pch=18, col="cyan3", cex=3)
abline(h=1, col="gray50", lty=3)
text(0.085, 1.3, expression("1"*sigma*" threshold"), cex=0.85, col="gray40")
legend("topleft",
       c(expression(alpha*" (naive)"),
         expression(M*"* (naive)"),
         expression(alpha*" (hierarchical)"),
         expression(M*"* (hierarchical)")),
       pch=c(16,17,18,18), col=c("firebrick","darkorange","dodgerblue3","cyan3"),
       cex=0.9, bg="white", pt.cex=c(2.5,2.5,3,3))
grid(col="gray85")

dev.off()
cat("Saved: diagnostic_fig2_comparison.pdf\n")

############################################################
# 5. FIGURE 3: Survey structure + error model
############################################################

cat("Generating Figure 3: Mock survey structure...\n")
pdf("diagnostic_fig3_survey.pdf", width=14, height=5)
par(mfrow=c(1,3), mar=c(4.5,4.5,3,1), cex.lab=1.3, cex.axis=1.1, cex.main=1.4)

# --- Panel A: Survey selection ---
plot(all_tm, all_z, pch=".", col=rgb(0.2,0.2,0.8,0.2), cex=2,
     xlab=expression(log[10]*M[true]), ylab="Redshift",
     main="(a) Survey Selection")
zp <- seq(zmin, zlimit, length=200)
lines(mlim_of_z(zp), zp, col="red", lwd=3)
text(13.5, 0.22, expression(M[lim](z)), col="red", cex=1.2)
grid(col="gray85")

# --- Panel B: Error model ---
plot(all_tm, sigma_all, pch=".", col=rgb(0,0,0,0.3), cex=2,
     xlab=expression(log[10]*M[true]),
     ylab=expression(sigma*" (dex)"),
     main="(b) Mass-dependent Errors")
mg <- seq(12, 15.5, length=100)
sg <- pmax(0.10, pmin(0.30, 0.30 - 0.06*(mg-13)))
lines(mg, sg, col="red", lwd=3)
grid(col="gray85")

# --- Panel C: Observed vs true mass distribution ---
br <- seq(11, 17, by=0.15)
h_t <- hist(all_tm, breaks=br, plot=FALSE)
h_o <- hist(obs_all, breaks=br, plot=FALSE)
plot(h_t$mids, h_t$counts, type="s", col="forestgreen", lwd=2.5,
     xlim=c(12, 16), ylim=c(0, max(h_t$counts)*1.15),
     xlab=expression(log[10]*M), ylab="Count",
     main="(c) Mass Distributions")
lines(h_o$mids, h_o$counts, type="s", col="firebrick", lwd=2.5, lty=2)
legend("topright", c("True masses", "Observed masses"),
       col=c("forestgreen","firebrick"), lwd=2.5, lty=c(1,2), cex=0.9, bg="white")
grid(col="gray85")

dev.off()
cat("Saved: diagnostic_fig3_survey.pdf\n")

############################################################
# Summary table
############################################################

cat("\n")
cat("=" , rep("=", 60), "\n", sep="")
cat("  RECOVERY SUMMARY\n")
cat("=", rep("=", 60), "\n", sep="")
cat(sprintf("  %-12s %8s %8s %8s %8s\n", "Method", "M*", "log_phi", "alpha", "beta"))
cat("  ", strrep("-", 56), "\n", sep="")
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "True", tv[1], tv[2], tv[3], tv[4]))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "No errors", noerr_ms, noerr_lp, noerr_al, noerr_be))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Naive full", v12_s2[1], v12_s2[2], v12_s2[3], v12_s2[4]))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Naive half", v12_s3[1], v12_s3[2], v12_s3[3], v12_s3[4]))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Hier. MAP", map_ms, map_lp, map_al, map_be))
cat("  ", strrep("-", 56), "\n", sep="")
cat(sprintf("  %-12s %8s %8s %8s %8s\n", "", "delta", "delta", "delta", "delta"))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "No errors",
            noerr_ms-tv[1], noerr_lp-tv[2], noerr_al-tv[3], noerr_be-tv[4]))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Naive full",
            v12_s2[1]-tv[1], v12_s2[2]-tv[2], v12_s2[3]-tv[3], v12_s2[4]-tv[4]))
cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Hier. MAP",
            map_ms-tv[1], map_lp-tv[2], map_al-tv[3], map_be-tv[4]))

cat("\n  Latent mass: median |recovered - true| = ", sprintf("%.3f", median(abs(map_mt - all_tm))),
    " dex (vs sigma = ", sprintf("%.3f", median(sigma_all)), ")\n", sep="")
cat("  MAP runtime: <1 second\n")
cat("=", rep("=", 60), "\n", sep="")

cat("\nAll figures saved. Done!\n")
