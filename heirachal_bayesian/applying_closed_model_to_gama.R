############################################################
# HIERARCHICAL MRP FIT TO GAMA DATA (v2 - fixed mlim)
#
# Uses the proven hierarchical latent-mass model from
# simulation-recovery testing (v13/v14/v16).
#
# v2 fixes: mass limit handling that was causing the
# optimizer to find unphysical solutions (alpha > 0).
#
# Root cause in v1: the smooth mlim(z) placed detection
# limits ABOVE many observed masses, forcing the model
# into an impossible corner. Now uses:
#   - Per-group mlim from zmax (physically grounded)
#   - Per-shell mlim from the lower envelope of per-group
#     limits (for the Poisson normalisation Lambda)
#   - Consistency check: every group must be plausibly
#     above its limit within measurement error
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
# 1. Data preparation (unchanged)
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5

cat("============================================\n")
cat("  HIERARCHICAL MRP FIT TO GAMA DATA (v2)\n")
cat("============================================\n\n")

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

# Mass errors from multiplicity-error relation
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

############################################################
# 2. Calculate zmax for each group
############################################################

cat("Calculating zmax per group...\n")

gig <- fread("../data/GAMAGalsInGroups.csv")

g3c$zmax <- NA
for (i in 1:nrow(g3c)) {
    if(g3c$Nfof[i] == 2) {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[2]
    } else {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[multi]
    }
}

g3c$zmax <- ifelse(g3c$zmax < g3c$Zfof, g3c$Zfof, g3c$zmax)
g3c$zmax <- ifelse(g3c$zmax > zlimit, zlimit, g3c$zmax)

############################################################
# 3. Per-group mass limits from zmax
#
# CRITICAL: The per-group mlim must satisfy two conditions:
#   (a) Physical: derived from the group's actual detectability
#   (b) Consistent: x_obs must be plausibly above mlim within
#       measurement error (otherwise the latent mass model
#       can't find a valid mt > mlim)
#
# The zmax-based approach from the original script:
#   mlim_i = log10(M_obs_i) - 2*log10(zmax_i / z_obs_i)
# This comes from the idea that if the group were pushed to
# zmax, its flux would drop as ~(d_L)^2, so the minimum
# detectable mass scales similarly.
############################################################

cat("Computing per-group mass limits...\n")

g3c$log_mass <- log10(g3c$MassAfunc)

# zmax-based mass limit (same formula as original)
g3c$log_mass_limit <- g3c$log_mass - 2*log10(g3c$zmax / g3c$Zfof)

# Cap: mlim should not exceed x_obs - 0.5 dex
# (group must be comfortably above its detection limit)
g3c$log_mass_limit <- pmin(g3c$log_mass_limit, g3c$log_mass - 0.5)

# Floor: no detection below 10^10.5
g3c$log_mass_limit <- pmax(g3c$log_mass_limit, 10.5)

cat(sprintf("  Mass limit range: %.2f -- %.2f\n", 
            min(g3c$log_mass_limit, na.rm=TRUE), max(g3c$log_mass_limit, na.rm=TRUE)))
cat(sprintf("  Median (M_obs - mlim): %.2f dex\n",
            median(g3c$log_mass - g3c$log_mass_limit, na.rm=TRUE)))

############################################################
# 4. Redshift shells for Lambda (Poisson normalisation)
#
# For the expected count Lambda, we need a mass limit per
# redshift shell. We derive this from the per-group limits:
# in each shell, use a robust low percentile of per-group
# mlim values as the shell's detection boundary.
############################################################

cat("Setting up redshift shells...\n")

# GAMA sky area
sky_area_deg2 <- 179.92
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2

# Shell volumes using celestial package
d_edges <- sapply(z_edges, function(z) cosdistCoDist(z, OmegaM=omegam, H0=ho))
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
Vsurvey <- sum(V_sh)

# Per-shell mass limit: use 10th percentile of per-group limits
# within each shell for a robust lower envelope.
mlim_sh <- numeric(nsh)
for(j in 1:nsh) {
    in_shell <- !is.na(g3c$Zfof) & !is.na(g3c$log_mass_limit) &
                g3c$Zfof >= z_edges[j] & g3c$Zfof < z_edges[j+1]
    if(sum(in_shell) >= 5) {
        mlim_sh[j] <- quantile(g3c$log_mass_limit[in_shell], 0.10)
    } else if(sum(in_shell) > 0) {
        mlim_sh[j] <- min(g3c$log_mass_limit[in_shell])
    } else {
        mlim_sh[j] <- NA
    }
}

# Fill NAs by interpolation
if(any(is.na(mlim_sh))) {
    ok_sh <- !is.na(mlim_sh)
    mlim_sh <- approx(z_mids[ok_sh], mlim_sh[ok_sh], z_mids, rule=2)$y
}

# Enforce monotonically non-decreasing with redshift
for(j in 2:nsh) {
    mlim_sh[j] <- max(mlim_sh[j], mlim_sh[j-1])
}

cat(sprintf("  Survey volume: %.4e Mpc^3\n", Vsurvey))
cat(sprintf("  N shells: %d\n", nsh))
cat(sprintf("  mlim_sh range: %.2f -- %.2f\n", min(mlim_sh), max(mlim_sh)))

############################################################
# 5. Prepare final data vectors
############################################################

good <- is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0 & 
        !is.na(g3c$log_mass_limit) & !is.na(g3c$log10MassErr) &
        !is.na(g3c$zmax)

x_obs     <- g3c$log_mass[good]
sigma_obs <- g3c$log10MassErr[good]
mlim_per  <- g3c$log_mass_limit[good]

# Quality cuts
keep <- x_obs > 10 & x_obs < 17
x_obs     <- x_obs[keep]
sigma_obs <- sigma_obs[keep]
mlim_per  <- mlim_per[keep]

N <- length(x_obs)

# DIAGNOSTIC: verify every group is above its per-group limit
n_below <- sum(x_obs < mlim_per)
margin  <- x_obs - mlim_per

cat(sprintf("\n  N groups: %d\n", N))
cat(sprintf("  Mass range: %.2f -- %.2f\n", min(x_obs), max(x_obs)))
cat(sprintf("  Sigma range: %.2f -- %.2f (median %.2f)\n", 
            min(sigma_obs), max(sigma_obs), median(sigma_obs)))
cat(sprintf("  mlim_per range: %.2f -- %.2f\n", min(mlim_per), max(mlim_per)))
cat(sprintf("  Groups below own mlim: %d (should be 0)\n", n_below))
cat(sprintf("  Min margin (x_obs - mlim): %.2f dex\n", min(margin)))
cat(sprintf("  Median margin: %.2f dex\n\n", median(margin)))

if(n_below > 0) {
    warning(sprintf("%d groups have x_obs < mlim! Check mass limit calculation.", n_below))
}

############################################################
# 6. Hierarchical Stan model (from recovery test v14/v16)
############################################################

stan_hierarchical <- "
data {
  int<lower=1> N;          // number of groups
  vector[N] x_obs;         // observed log10(M)
  vector<lower=0>[N] sig;  // measurement error per group (dex)
  vector[N] mlim;          // per-group detection limit on true mass
  int<lower=1> Nsh;        // number of redshift shells
  vector[Nsh] V_sh;        // comoving volume per shell (Mpc^3)
  vector[Nsh] mlim_sh;     // mass limit per shell
  real xhi;                // upper mass limit for integration
  int<lower=2> Ng;         // grid points for Lambda integration
}
transformed data {
  real ln10 = log(10.0);
  real xlo = min(mlim_sh) - 0.5;
  real dx = (xhi - xlo) / (Ng - 1.0);
  vector[Ng] xg;
  for(k in 1:Ng) xg[k] = xlo + (k-1)*dx;
}
parameters {
  real ms;                          // M*
  real lp;                          // log10(phi*)
  real al;                          // alpha
  real<lower=0.1, upper=2.0> be;    // beta
  vector[N] z_raw;                  // non-centered latent masses
}
transformed parameters {
  vector[N] mt;
  for(i in 1:N) mt[i] = x_obs[i] + sig[i] * z_raw[i];
}
model {
  // Priors
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.7, 1.0);
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
  
  // Lambda = sum_j V_j * integral(phi, mlim_j, xhi)
  real Lambda = 0;
  for(j in 1:Nsh) {
    int k0 = 1;
    for(k in 1:Ng) if(xg[k] <= mlim_sh[j]) k0 = k;
    if(k0 >= Ng) k0 = Ng-1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;
  
  // MRP at latent true masses
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
# 7. Compile and prepare
############################################################

cat("Compiling hierarchical Stan model...\n")
sm <- stan_model(model_code=stan_hierarchical)

stan_data <- list(
    N       = N,
    x_obs   = x_obs,
    sig     = sigma_obs,
    mlim    = mlim_per,
    Nsh     = nsh,
    V_sh    = V_sh,
    mlim_sh = mlim_sh,
    xhi     = 16.5,
    Ng      = 200
)

init_fn <- function() {
    list(ms    = rnorm(1, 14, 0.2),
         lp    = rnorm(1, -4, 0.2),
         al    = rnorm(1, -1.7, 0.1),
         be    = runif(1, 0.5, 0.8),
         z_raw = rnorm(N, 0, 0.3))
}

############################################################
# 8. MAP via optimizing()
############################################################

cat("\n========== APPROACH A: MAP (optimizing) ==========\n")

# Try multiple initialisations to avoid local optima
best_opt <- NULL
best_lp_val  <- -Inf

for(trial in 1:5) {
    cat(sprintf("  Trial %d... ", trial))
    this_init <- init_fn()
    this_opt <- tryCatch(
        optimizing(sm, data=stan_data, init=this_init,
                   hessian=FALSE, as_vector=FALSE,
                   iter=10000, algorithm="LBFGS"),
        error=function(e) { cat(sprintf("failed: %s\n", e$message)); NULL }
    )
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f  ms=%.2f al=%.2f (code %d)\n", 
                    this_opt$value, this_opt$par$ms, this_opt$par$al, this_opt$return_code))
        if(this_opt$value > best_lp_val) {
            best_lp_val <- this_opt$value
            best_opt <- this_opt
        }
    }
}

# Re-run the best with hessian
if(!is.null(best_opt)) {
    cat("\n  Re-optimizing best solution with hessian...\n")
    rerun_init <- list(
        ms = best_opt$par$ms, lp = best_opt$par$lp,
        al = best_opt$par$al, be = best_opt$par$be,
        z_raw = best_opt$par$z_raw)
    opt <- tryCatch(
        optimizing(sm, data=stan_data, init=rerun_init,
                   hessian=TRUE, as_vector=FALSE,
                   iter=10000, algorithm="LBFGS"),
        error=function(e) { cat("Hessian rerun failed, using best trial\n"); best_opt }
    )
} else {
    stop("All optimization trials failed!")
}

cat(sprintf("\nConverged: %s (return code %d)\n", 
            ifelse(opt$return_code==0, "YES", "NO"), opt$return_code))

map_ms <- opt$par$ms
map_lp <- opt$par$lp
map_al <- opt$par$al
map_be <- opt$par$be

cat(sprintf("\nMAP estimates:\n"))
cat(sprintf("  M*       = %.3f\n", map_ms))
cat(sprintf("  log_phi* = %.3f\n", map_lp))
cat(sprintf("  alpha    = %.3f\n", map_al))
cat(sprintf("  beta     = %.3f\n", map_be))

# Sanity check
if(map_al > 0) {
    warning("alpha > 0 is unphysical! Something is still wrong with mass limits.")
}
if(map_ms < 12) {
    warning("M* < 12 is unphysical! Expected ~13-14.")
}

# Laplace standard errors
se_laplace <- rep(NA, 4)
pnames <- c("ms","lp","al","be")

if(!is.null(opt$hessian)) {
    cat("\nComputing Laplace standard errors...\n")
    H <- opt$hessian
    np <- 4
    nz_h <- N
    
    tryCatch({
        A <- -H[1:np, 1:np]
        B <- -H[1:np, (np+1):(np+nz_h)]
        D <- -H[(np+1):(np+nz_h), (np+1):(np+nz_h)]
        D_diag_inv <- diag(1/diag(D))
        Schur <- A - B %*% D_diag_inv %*% t(B)
        Sigma_mrp <- solve(Schur)
        se_laplace <- sqrt(diag(Sigma_mrp))
        
        cat("Laplace SE:\n")
        for(i in 1:4) {
            cat(sprintf("  %-10s MAP = %7.3f  SE = %.3f\n",
                        c("M*","log_phi*","alpha","beta")[i],
                        opt$par[[pnames[i]]], se_laplace[i]))
        }
    }, error=function(e) {
        cat(sprintf("Hessian analysis failed: %s\n", e$message))
    })
}

# Latent mass summary
map_zraw <- opt$par$z_raw
map_mt <- x_obs + sigma_obs * map_zraw
cat(sprintf("\nLatent mass shift: median |MAP_mt - x_obs| = %.3f dex\n",
            median(abs(map_mt - x_obs))))

n_mt_below <- sum(map_mt < mlim_per)
cat(sprintf("Latent masses below mlim: %d / %d\n", n_mt_below, N))

############################################################
# 9. Short MCMC from MAP
############################################################

cat("\n========== APPROACH B: Short MCMC ==========\n")

init_from_map <- function() {
    list(ms    = map_ms + rnorm(1, 0, 0.05),
         lp    = map_lp + rnorm(1, 0, 0.05),
         al    = map_al + rnorm(1, 0, 0.02),
         be    = min(1.9, max(0.15, map_be + rnorm(1, 0, 0.05))),
         z_raw = map_zraw + rnorm(N, 0, 0.1))
}

t0 <- proc.time()
fit_mcmc <- sampling(sm, data=stan_data,
    chains=4, iter=2000, warmup=1000, thin=5, cores=4,
    init=lapply(1:4, function(i) init_from_map()),
    control=list(adapt_delta=0.90, max_treedepth=14))
t_mcmc <- (proc.time()-t0)[3]

cat(sprintf("MCMC time: %.1f seconds\n", t_mcmc))

pmc <- extract(fit_mcmc)
pmmc <- cbind(ms=pmc$ms, lp=pmc$lp, al=pmc$al, be=pmc$be)
medmc <- apply(pmmc, 2, median)
q16mc <- apply(pmmc, 2, quantile, 0.16)
q84mc <- apply(pmmc, 2, quantile, 0.84)

summ_mc <- summary(fit_mcmc, pars=c("ms","lp","al","be"))$summary
n_div <- sum(sapply(1:4, function(ch) 
    sum(get_sampler_params(fit_mcmc, inc_warmup=FALSE)[[ch]][,"divergent__"])))
n_mt_tree <- sum(sapply(1:4, function(ch) 
    sum(get_sampler_params(fit_mcmc, inc_warmup=FALSE)[[ch]][,"treedepth__"] >= 14)))

cat(sprintf("Max Rhat: %.3f | Min n_eff: %.0f | Div: %d | MaxTree: %d\n",
            max(summ_mc[,"Rhat"]), min(summ_mc[,"n_eff"]), n_div, n_mt_tree))

cat("\nMCMC Results:\n")
cat(sprintf("  %-10s  Median    +1sig    -1sig\n", ""))
for(p in pnames) {
    plab <- c(ms="M*", lp="log_phi*", al="alpha", be="beta")[p]
    cat(sprintf("  %-10s  %7.3f  +%.3f  -%.3f\n",
                plab, medmc[p], q84mc[p]-medmc[p], medmc[p]-q16mc[p]))
}

############################################################
# 10. Results comparison
############################################################

cat("\n============================================\n")
cat("  FINAL RESULTS COMPARISON\n")
cat("============================================\n\n")

driver22 <- c(ms=13.51, lp=-3.19, al=-1.27, be=0.47)

cat(sprintf("  %-10s  %10s  %10s  %10s\n", "Parameter", "MAP", "MCMC med", "Driver+22"))
cat(sprintf("  %-10s  %10s  %10s  %10s\n", "---------", "---", "--------", "---------"))
for(p in pnames) {
    plab <- c(ms="M*", lp="log_phi*", al="alpha", be="beta")[p]
    se_str <- if(!is.na(se_laplace[which(pnames==p)])) 
        sprintf("%.3f+/-%.3f", opt$par[[p]], se_laplace[which(pnames==p)]) 
    else 
        sprintf("%.3f", opt$par[[p]])
    cat(sprintf("  %-10s  %10s  %7.3f+/-%.3f  %10.3f\n",
                plab, se_str,
                medmc[p], (q84mc[p]-q16mc[p])/2,
                driver22[p]))
}

############################################################
# 11. Plotting
############################################################

mrp_log10 <- function(x, mstar, log_phi, alpha, beta) {
    log10(pmax(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))), 1e-30))
}

x_min <- floor(min(x_obs) * 10) / 10
x_max <- ceiling(max(x_obs) * 10) / 10
breaks <- seq(x_min, x_max, by=0.2)
hist_plt <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_plt  <- hist_plt$counts / (Vsurvey * 0.2)
ok       <- phi_plt > 0
xfit     <- seq(x_min, 16.5, length.out=500)

CairoPDF("MRP_HIERARCHICAL_FIT.pdf", 14, 10)
par(mfrow=c(2,3))

# Panel 1: HMF
plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Hierarchical MRP Fit to GAMA (v2)")
grid(col="gray80")

if(nrow(pmmc) > 10) {
    idx <- sample(1:nrow(pmmc), min(300, nrow(pmmc)))
    for(i in idx) {
        y <- tryCatch(
            mrp_log10(xfit, pmmc[i,"ms"], pmmc[i,"lp"], pmmc[i,"al"], pmmc[i,"be"]),
            error = function(e) rep(NA, length(xfit)))
        if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
    }
}

lines(xfit, mrp_log10(xfit, map_ms, map_lp, map_al, map_be), col="red", lwd=3)
lines(xfit, mrp_log10(xfit, medmc["ms"], medmc["lp"], medmc["al"], medmc["be"]),
      col="purple", lwd=2, lty=2)
lines(xfit, mrp_log10(xfit, driver22["ms"], driver22["lp"], driver22["al"], driver22["be"]),
      col="black", lwd=2, lty=3)

legend("bottomleft",
       legend=c("GAMA (binned)", "MAP", "MCMC median", 
                "Posterior draws", "Driver+22"),
       col=c("darkgreen","red","purple",rgb(0,0,1,0.3),"black"),
       pch=c(19,NA,NA,NA,NA), lty=c(NA,1,2,1,3),
       lwd=c(NA,3,2,1,2), bg="white", cex=0.8)

text(14.5, -2.3,
     sprintf("MAP:\nM* = %.2f\nlog\u03c6* = %.2f\n\u03b1 = %.2f\n\u03b2 = %.2f",
             map_ms, map_lp, map_al, map_be),
     col="red", cex=0.8, pos=4)

# Panel 2: Parameter comparison
vals_map  <- c(map_ms, map_lp, map_al, map_be)
vals_mcmc <- medmc
vals_d22  <- driver22

plot(1:4, vals_map, pch=19, col="red", cex=2, 
     ylim=range(c(vals_map, vals_mcmc, vals_d22)),
     xlab="", ylab="Parameter value", main="MAP vs MCMC vs Driver+22", xaxt="n")
points(1:4+0.1, vals_mcmc, pch=17, col="purple", cex=2)
points(1:4-0.1, vals_d22, pch=15, col="black", cex=2)
axis(1, at=1:4, labels=c("M*","log_phi*","alpha","beta"))
legend("topleft", c("MAP","MCMC","Driver+22"), 
       pch=c(19,17,15), col=c("red","purple","black"), cex=0.9)

# Panel 3: Latent mass shift
plot(x_obs, map_mt, pch=".", col=rgb(0,0,0,0.3),
     xlab="Observed mass", ylab="MAP latent mass",
     main="Latent mass recovery")
abline(0, 1, col="red", lwd=2)
points(x_obs[map_mt < mlim_per + 0.1], map_mt[map_mt < mlim_per + 0.1],
       pch=20, col="orange", cex=0.5)

# Panel 4: Mass and limit distributions
hist(x_obs, breaks=50, col=rgb(0.2,0.5,0.2,0.5), border="white",
     main="Observed mass vs detection limits", xlab="log10(M)")
hist(mlim_per, breaks=50, col=rgb(0.8,0.2,0.2,0.3), border="white", add=TRUE)
legend("topright", c("x_obs", "mlim"), 
       fill=c(rgb(0.2,0.5,0.2,0.5), rgb(0.8,0.2,0.2,0.3)), cex=0.9)

# Panel 5: Mass error distribution
hist(sigma_obs, breaks=30, col="steelblue", border="white",
     main="Mass error distribution", xlab="sigma (dex)")
abline(v=median(sigma_obs), col="red", lwd=2, lty=2)

# Panel 6: mlim vs redshift (keep z_groups for plotting)
z_groups <- g3c$Zfof[good][keep]
plot(z_groups, mlim_per, pch=".", col=rgb(0,0,0,0.3),
     xlab="Redshift", ylab="mlim (log10 M)",
     main="Per-group limits & shell limits")
points(z_mids, mlim_sh, col="red", pch=19, cex=1.5)
lines(z_mids, mlim_sh, col="red", lwd=2)
legend("topleft", c("Per-group", "Per-shell"), 
       col=c("black","red"), pch=c(46,19), cex=0.9)

dev.off()
cat("\nPlot saved: MRP_HIERARCHICAL_FIT.pdf\n")

############################################################
# 12. Save results
############################################################

saveRDS(list(
    map = list(ms=map_ms, lp=map_lp, al=map_al, be=map_be,
               se=se_laplace, mt=map_mt, z_raw=map_zraw,
               converged=(opt$return_code==0)),
    mcmc = list(fit=fit_mcmc, posterior=pmmc,
                median=medmc, q16=q16mc, q84=q84mc),
    data = list(x_obs=x_obs, sigma=sigma_obs, mlim=mlim_per,
                N=N, Vsurvey=Vsurvey, nsh=nsh, V_sh=V_sh, mlim_sh=mlim_sh),
    driver22 = driver22
), "hierarchical_mrp_results.rds")

cat("\nResults saved: hierarchical_mrp_results.rds\n")
cat("\nDone!\n")
