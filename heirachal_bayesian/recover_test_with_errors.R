############################################################
# SIMULATION-RECOVERY TEST v14 (PRODUCTION)
#
# Hierarchical model with latent true masses.
# Non-centered parameterisation for efficient sampling:
#   z_raw[i] ~ Normal(0, 1)
#   mt[i] = x_obs[i] + sig[i] * z_raw[i]
#
# This decorrelates the latent masses from MRP params,
# dramatically improving sampling efficiency.
#
# Also enforces mt[i] > mlim[i] via rejection in the model.
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup
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

cat("=== SETUP ===\n")
cat(sprintf("True: M*=%.2f phi=%.2f a=%.2f b=%.2f\n", tv[1],tv[2],tv[3],tv[4]))

############################################################
# 2. Generate halos + errors
############################################################

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

cat(sprintf("N = %d halos, median sigma = %.2f\n\n", N_true, median(sigma_all)))

############################################################
# 3. DIAGNOSTIC PLOT
############################################################

pdf("mock_diagnostic_v14.pdf", width=12, height=8)
par(mfrow=c(2,2))

br <- seq(11,17,by=0.2)
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
legend("topright",c("True","Observed","MRP"),pch=c(19,17,NA),
       col=c("blue","red","black"),lty=c(NA,NA,1),lwd=c(NA,NA,2))

plot(all_tm, obs_all-all_tm, pch=".", col=rgb(0,0,0,.3),
     xlab="True M", ylab="Obs-True", main="Mass errors")
abline(h=0,col="red",lwd=2); abline(h=c(-.3,.3),col="red",lty=2)

plot(all_tm, all_z, pch=".", col=rgb(0,0,1,.3),
     xlab="True M", ylab="z", main="Detection boundary")
zp <- seq(zmin,zlimit,length=200); lines(mlim_of_z(zp),zp,col="red",lwd=2)

hist(mlim_per, breaks=40, col="steelblue", border="white",
     xlab="m_lim", main=sprintf("Per-group m_lim (spread=%.1f dex)", diff(range(mlim_per))))

dev.off()
cat("Saved: mock_diagnostic_v14.pdf\n\n")

############################################################
# 4. STAGE 1: No errors (validation)
############################################################

cat("========== STAGE 1: No errors ==========\n")

stan_noerr <- "
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

fit1 <- stan(model_code=stan_noerr,
    data=list(N=N_true, x=all_tm, Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=300),
    chains=4, iter=2000, warmup=1000, cores=4,
    init=lapply(1:4, function(i) list(ms=rnorm(1,14,.2),lp=rnorm(1,-4,.2),
                                       al=rnorm(1,-1.7,.1),be=runif(1,.5,.8))),
    control=list(adapt_delta=0.95))

p1 <- extract(fit1)
pm1 <- cbind(ms=p1$ms, lp=p1$lp, al=p1$al, be=p1$be)
med1 <- apply(pm1,2,median)
q16_1 <- apply(pm1,2,quantile,.16); q84_1 <- apply(pm1,2,quantile,.84)

cat("Stage 1:\n")
for(p in names(tv)){ea<-((q84_1[p]-med1[p])+(med1[p]-q16_1[p]))/2
  cat(sprintf("  %-5s True=%.3f Med=%.3f Bias=%+.2f sig\n",p,tv[p],med1[p],(med1[p]-tv[p])/ea))}

############################################################
# 5. STAGE 2: Hierarchical (non-centered)
############################################################

cat("\n========== STAGE 2: Hierarchical (non-centered) ==========\n")

stan_hier_nc <- "
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
  vector[N] z_raw;  // non-centered: mt = x_obs + sig * z_raw
}
transformed parameters {
  vector[N] mt;
  for(i in 1:N)
    mt[i] = x_obs[i] + sig[i] * z_raw[i];
}
model {
  // Priors on MRP params
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.7, 1.0);
  
  // Non-centered prior on z_raw
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
  
  // Lambda (expected true count)
  real Lambda = 0;
  for(j in 1:Nsh) {
    int k0 = 1;
    for(k in 1:Ng) if(xg[k] <= mlim_sh[j]) k0 = k;
    if(k0 >= Ng) k0 = Ng-1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;
  
  // Per-group: MRP evaluated at latent true mass
  for(i in 1:N) {
    real u = mt[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    
    if(mt[i] > mlim[i] && phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
  // Note: the Normal(x_obs | mt, sig) is already accounted for
  // by the z_raw ~ std_normal() prior + the transformation
  // mt = x_obs + sig * z_raw, which implies:
  // x_obs = mt - sig * z_raw => x_obs ~ Normal(mt, sig)
  // The Jacobian is 1 (linear transform), so this is correct.
}
"

cat(sprintf("N=%d + %d latent = %d params\n", N_true, N_true, N_true+4))
cat("Non-centered parameterisation, max_treedepth=14\n\n")

init_fn <- function() {
    list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
         al=rnorm(1,-1.7,.1), be=runif(1,.5,.8),
         z_raw=rnorm(N_true, 0, 0.5))
}

fit2 <- stan(model_code=stan_hier_nc,
    data=list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
              Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=250),
    chains=4, iter=4000, warmup=2000, cores=4,
    init=lapply(1:4, function(i) init_fn()),
    control=list(adapt_delta=0.90, max_treedepth=14))

############################################################
# 6. RESULTS
############################################################

cat("\n=== STAGE 2 RESULTS ===\n")
print(fit2, pars=c("ms","lp","al","be"))

p2 <- extract(fit2)
pm2 <- cbind(ms=p2$ms, lp=p2$lp, al=p2$al, be=p2$be)
med2 <- apply(pm2,2,median)
q16_2 <- apply(pm2,2,quantile,.16); q84_2 <- apply(pm2,2,quantile,.84)

summ2 <- summary(fit2, pars=c("ms","lp","al","be"))$summary
n_div2 <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit2,inc_warmup=FALSE)[[ch]][,"divergent__"])))
n_maxtree <- sum(sapply(1:4, function(ch) sum(get_sampler_params(fit2,inc_warmup=FALSE)[[ch]][,"treedepth__"] >= 14)))

cat(sprintf("Max Rhat: %.3f | Min n_eff: %.0f | Div: %d | MaxTree: %d\n\n",
            max(summ2[,"Rhat"]), min(summ2[,"n_eff"]), n_div2, n_maxtree))

cat("Stage 2 Recovery:\n")
cat(sprintf("%-6s %8s %8s %8s %8s %10s\n","Param","True","Median","+1sig","-1sig","Bias"))
cat(strrep("-",52),"\n")
for(p in names(tv)) {
    eu <- q84_2[p]-med2[p]; ed <- med2[p]-q16_2[p]; ea <- (eu+ed)/2
    bias <- (med2[p]-tv[p])/ea
    st <- ifelse(abs(bias)<1,"OK",ifelse(abs(bias)<2,"marginal","FAIL"))
    cat(sprintf("%-6s %8.3f %8.3f  +%.3f  -%.3f  %+.2f sig [%s]\n",
                p, tv[p], med2[p], eu, ed, bias, st))
}

cat("\nCorrelations:\n")
print(round(cor(pm2), 3))

# Latent mass recovery
mt_med <- apply(p2$mt, 2, median)
cat(sprintf("\nLatent mass: median |recovered - true| = %.3f dex (vs sigma = %.3f)\n",
            median(abs(mt_med - all_tm)), median(sigma_all)))

############################################################
# 7. PLOTS
############################################################

pdf("sim_recovery_v14.pdf", width=14, height=14)
par(mfrow=c(3,3))

xf <- seq(11.5, 16.5, length=500)

# Stage 1 HMF
h1 <- hist(all_tm, breaks=br, plot=FALSE)
ph1 <- h1$counts/(Vsurvey*0.2); ok1 <- ph1>0
plot(h1$mids[ok1], log10(ph1[ok1]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main=sprintf("S1: True masses (N=%d)", N_true))
grid(col="gray80")
idx <- sample(nrow(pm1), min(200,nrow(pm1)))
for(i in idx){y<-tryCatch(mrp_log10(xf,pm1[i,"ms"],pm1[i,"lp"],pm1[i,"al"],pm1[i,"be"]),
    error=function(e) rep(NA,length(xf))); if(all(is.finite(y))) lines(xf,y,col=rgb(0,0,1,.03))}
lines(xf, mrp_log10(xf,tv["ms"],tv["lp"],tv["al"],tv["be"]),col="black",lwd=2,lty=2)
lines(xf, mrp_log10(xf,med1["ms"],med1["lp"],med1["al"],med1["be"]),col="red",lwd=3)

# Stage 2 HMF
h2 <- hist(obs_all, breaks=br, plot=FALSE)
ph2 <- h2$counts/(Vsurvey*0.2); ok2 <- ph2>0
plot(h2$mids[ok2], log10(ph2[ok2]), pch=19, col="darkgreen", cex=1.3,
     xlim=c(12,16), ylim=c(-8,-2), xlab="log10(M)", ylab="log10(phi)",
     main=sprintf("S2: Hierarchical (N=%d)", N_true))
grid(col="gray80")
idx <- sample(nrow(pm2), min(200,nrow(pm2)))
for(i in idx){y<-tryCatch(mrp_log10(xf,pm2[i,"ms"],pm2[i,"lp"],pm2[i,"al"],pm2[i,"be"]),
    error=function(e) rep(NA,length(xf))); if(all(is.finite(y))) lines(xf,y,col=rgb(0,0,1,.03))}
lines(xf, mrp_log10(xf,tv["ms"],tv["lp"],tv["al"],tv["be"]),col="black",lwd=2,lty=2)
lines(xf, mrp_log10(xf,med2["ms"],med2["lp"],med2["al"],med2["be"]),col="red",lwd=3)
legend("topright",c("Data","TRUE","Fit"),col=c("darkgreen","black","red"),
       pch=c(19,NA,NA),lty=c(NA,2,1),lwd=c(NA,2,3),cex=.8)

# Latent mass recovery
plot(all_tm, mt_med, pch=".", col=rgb(0,0,0,.3),
     xlab="True mass", ylab="Recovered mass", main="Latent mass recovery")
abline(0,1,col="red",lwd=2)

# Posteriors for all 4 params
for(p in c("ms","lp","al","be")) {
    hist(pm2[,p], breaks=50,
         main=c(ms="M*",lp="log_phi",al="alpha",be="beta")[p],
         xlab=p, col="steelblue", border="white", freq=FALSE)
    abline(v=tv[p], col="red", lwd=3)
    abline(v=med2[p], col="blue", lwd=2, lty=2)
    legend("topright",c(sprintf("True: %.3f",tv[p]),sprintf("Med: %.3f",med2[p])),
           col=c("red","blue"),lty=c(1,2),lwd=c(3,2),cex=.8)
}

# Bias comparison
bias1 <- sapply(names(tv), function(p) (med1[p]-tv[p])/((q84_1[p]-q16_1[p])/2))
bias2 <- sapply(names(tv), function(p) (med2[p]-tv[p])/((q84_2[p]-q16_2[p])/2))
plot(1:4-.02, bias1, pch=19, col="blue", cex=2, ylim=c(-3,3),
     xlab="", ylab="Bias (sigma)", main="Recovery comparison", xaxt="n")
points(1:4+.02, bias2, pch=17, col="red", cex=2)
axis(1, at=1:4, labels=c("M*","phi*","alpha","beta"))
abline(h=0, col="gray"); abline(h=c(-1,1), col="gray", lty=2)
legend("topleft", c("S1: no err","S2: hierarchical"), pch=c(19,17), col=c("blue","red"))

dev.off()
cat("\nSaved: sim_recovery_v14.pdf, mock_diagnostic_v14.pdf\n")

saveRDS(list(true=tv,
    stage1=list(post=pm1, med=med1, q16=q16_1, q84=q84_1),
    stage2=list(post=pm2, med=med2, q16=q16_2, q84=q84_2,
                mt_med=mt_med, true_mass=all_tm)),
    "sim_recovery_v14.rds")
cat("Done!\n")
