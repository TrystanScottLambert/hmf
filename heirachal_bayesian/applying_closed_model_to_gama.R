############################################################
# HIERARCHICAL MRP FIT TO GAMA DATA (v3)
#
# v3 changes:
# - Option A: Uniform mass cut (conservative but clean)
# - Option B: Per-group zmax-based limits (more data, noisier)
# - Both use the proven hierarchical latent-mass Stan model
# - Better diagnostics for understanding mlim behaviour
# - Increased MCMC iterations for better mixing
#
# The key insight: the model needs to know WHERE the
# detection boundary is. A noisy/wrong mlim biases alpha
# more than losing some low-mass groups does.
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
# USER SETTINGS
############################################################

# Choose mass limit strategy:
#   "uniform"  = single mass cut applied everywhere (safest)
#   "pergroup" = zmax-based per-group limits (more data)
MLIM_STRATEGY <- "uniform"

# Uniform mass cut (only used if MLIM_STRATEGY == "uniform")
# This should be above the completeness limit at zlimit.
# From the v2 shell limits, mlim(z=0.25) ~ 12.9, so 12.8
# is a reasonable conservative choice.
UNIFORM_MCUT <- 12.8

############################################################
# 1. Data preparation (unchanged)
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5

cat("============================================\n")
cat("  HIERARCHICAL MRP FIT TO GAMA DATA (v3)\n")
cat(sprintf("  Strategy: %s", MLIM_STRATEGY))
if(MLIM_STRATEGY == "uniform") cat(sprintf(" (cut = %.1f)", UNIFORM_MCUT))
cat("\n============================================\n\n")

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

# Mass errors
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
g3c$log_mass <- log10(g3c$MassAfunc)

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
# 3. Mass limits and shell setup
############################################################

# GAMA sky area
sky_area_deg2 <- 179.92
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2

# Shell volumes
d_edges <- sapply(z_edges, function(z) cosdistCoDist(z, OmegaM=omegam, H0=ho))
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac
Vsurvey <- sum(V_sh)

if(MLIM_STRATEGY == "uniform") {
    ##########################################################
    # OPTION A: Uniform mass cut
    #
    # Every group must have x_obs > UNIFORM_MCUT.
    # Every shell has the same mlim = UNIFORM_MCUT.
    # Per-group mlim = UNIFORM_MCUT for all groups.
    #
    # This is the closest to the recovery test setup:
    # the model knows EXACTLY where the detection boundary
    # is, so alpha can be properly constrained.
    ##########################################################
    
    cat(sprintf("\nApplying uniform mass cut at %.1f...\n", UNIFORM_MCUT))
    
    # Filter groups
    good <- is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0 & 
            !is.na(g3c$log10MassErr) & g3c$log_mass > UNIFORM_MCUT &
            g3c$log_mass < 17
    
    x_obs     <- g3c$log_mass[good]
    sigma_obs <- g3c$log10MassErr[good]
    mlim_per  <- rep(UNIFORM_MCUT, sum(good))  # same for all
    
    # All shells have the same limit
    mlim_sh <- rep(UNIFORM_MCUT, nsh)
    
    N <- length(x_obs)
    cat(sprintf("  Groups above cut: %d / %d\n", N, nrow(g3c)))
    
} else if(MLIM_STRATEGY == "pergroup") {
    ##########################################################
    # OPTION B: Per-group zmax-based limits
    #
    # Same as v2 but with improved limit calculation.
    # The zmax-based limit is:
    #   mlim_i = log10(M_i) - 2*log10(zmax_i / z_i)
    # Capped at x_obs - 0.5 dex, floored at 10.5.
    ##########################################################
    
    cat("\nUsing per-group zmax-based mass limits...\n")
    
    g3c$log_mass_limit <- g3c$log_mass - 2*log10(g3c$zmax / g3c$Zfof)
    g3c$log_mass_limit <- pmin(g3c$log_mass_limit, g3c$log_mass - 0.5)
    g3c$log_mass_limit <- pmax(g3c$log_mass_limit, 10.5)
    
    good <- is.finite(g3c$MassAfunc) & g3c$MassAfunc > 0 & 
            !is.na(g3c$log_mass_limit) & !is.na(g3c$log10MassErr) &
            !is.na(g3c$zmax) & g3c$log_mass > 10 & g3c$log_mass < 17
    
    x_obs     <- g3c$log_mass[good]
    sigma_obs <- g3c$log10MassErr[good]
    mlim_per  <- g3c$log_mass_limit[good]
    
    # Per-shell limits from 10th percentile
    for(j in 1:nsh) {
        in_shell <- g3c$Zfof[good] >= z_edges[j] & g3c$Zfof[good] < z_edges[j+1]
        if(sum(in_shell) >= 5) {
            mlim_sh[j] <- quantile(mlim_per[in_shell], 0.10)
        } else if(sum(in_shell) > 0) {
            mlim_sh[j] <- min(mlim_per[in_shell])
        } else {
            mlim_sh[j] <- NA
        }
    }
    if(any(is.na(mlim_sh))) {
        ok_sh <- !is.na(mlim_sh)
        mlim_sh <- approx(z_mids[ok_sh], mlim_sh[ok_sh], z_mids, rule=2)$y
    }
    for(j in 2:nsh) mlim_sh[j] <- max(mlim_sh[j], mlim_sh[j-1])
    
    N <- length(x_obs)
}

cat(sprintf("\n  N groups: %d\n", N))
cat(sprintf("  Mass range: %.2f -- %.2f\n", min(x_obs), max(x_obs)))
cat(sprintf("  Sigma range: %.2f -- %.2f (median %.2f)\n", 
            min(sigma_obs), max(sigma_obs), median(sigma_obs)))
cat(sprintf("  mlim_per range: %.2f -- %.2f\n", min(mlim_per), max(mlim_per)))
cat(sprintf("  mlim_sh range: %.2f -- %.2f\n", min(mlim_sh), max(mlim_sh)))

margin <- x_obs - mlim_per
n_below <- sum(margin < 0)
cat(sprintf("  Groups below own mlim: %d (should be 0)\n", n_below))
cat(sprintf("  Min margin: %.2f dex | Median margin: %.2f dex\n", 
            min(margin), median(margin)))
cat(sprintf("  Survey volume: %.4e Mpc^3\n\n", Vsurvey))

############################################################
# 4. Stan model (identical to recovery test)
############################################################

stan_hierarchical <- "
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

############################################################
# 5. Compile and prepare
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
# 6. MAP via optimizing() with multi-start
############################################################

cat("\n========== APPROACH A: MAP (optimizing) ==========\n")

best_opt <- NULL
best_lp_val <- -Inf

for(trial in 1:5) {
    cat(sprintf("  Trial %d... ", trial))
    this_init <- init_fn()
    this_opt <- tryCatch(
        optimizing(sm, data=stan_data, init=this_init,
                   hessian=FALSE, as_vector=FALSE,
                   iter=20000, algorithm="LBFGS",
                   tol_rel_grad=1e-12),
        error=function(e) { cat(sprintf("failed: %s\n", e$message)); NULL }
    )
    if(!is.null(this_opt)) {
        cat(sprintf("lp=%.1f  ms=%.2f lp=%.2f al=%.2f be=%.2f (code %d)\n", 
                    this_opt$value, this_opt$par$ms, this_opt$par$lp,
                    this_opt$par$al, this_opt$par$be, this_opt$return_code))
        if(this_opt$value > best_lp_val) {
            best_lp_val <- this_opt$value
            best_opt <- this_opt
        }
    }
}

if(!is.null(best_opt)) {
    cat("\n  Re-optimizing best with hessian...\n")
    rerun_init <- list(
        ms = best_opt$par$ms, lp = best_opt$par$lp,
        al = best_opt$par$al, be = best_opt$par$be,
        z_raw = best_opt$par$z_raw)
    opt <- tryCatch(
        optimizing(sm, data=stan_data, init=rerun_init,
                   hessian=TRUE, as_vector=FALSE,
                   iter=20000, algorithm="LBFGS",
                   tol_rel_grad=1e-12),
        error=function(e) { cat("Hessian rerun failed\n"); best_opt }
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

if(map_al > 0) warning("alpha > 0 is unphysical!")
if(map_ms < 12) warning("M* < 12 is unphysical!")

# Laplace SEs
se_laplace <- rep(NA, 4)
pnames <- c("ms","lp","al","be")

if(!is.null(opt$hessian)) {
    cat("\nComputing Laplace standard errors...\n")
    H <- opt$hessian
    np <- 4; nz_h <- N
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
    }, error=function(e) cat(sprintf("Hessian failed: %s\n", e$message)))
}

map_zraw <- opt$par$z_raw
map_mt <- x_obs + sigma_obs * map_zraw
cat(sprintf("\nLatent mass shift: median |MAP_mt - x_obs| = %.3f dex\n",
            median(abs(map_mt - x_obs))))
cat(sprintf("Latent masses below mlim: %d / %d\n", sum(map_mt < mlim_per), N))

############################################################
# 7. MCMC
############################################################

cat("\n========== APPROACH B: MCMC ==========\n")

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
# 8. Results comparison
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
# 9. Plotting
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
xfit     <- seq(11, 16.5, length.out=500)

CairoPDF("MRP_HIERARCHICAL_FIT.pdf", 14, 10)
par(mfrow=c(2,3))

# Panel 1: HMF
plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main=sprintf("Hierarchical MRP (v3, %s)", MLIM_STRATEGY))
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

# Show mass cut if uniform
if(MLIM_STRATEGY == "uniform") {
    abline(v=UNIFORM_MCUT, col="orange", lwd=2, lty=4)
}

legend("bottomleft",
       legend=c("GAMA (binned)", "MAP", "MCMC median", 
                "Posterior draws", "Driver+22",
                if(MLIM_STRATEGY=="uniform") "Mass cut" else NULL),
       col=c("darkgreen","red","purple",rgb(0,0,1,0.3),"black",
             if(MLIM_STRATEGY=="uniform") "orange" else NULL),
       pch=c(19,NA,NA,NA,NA, if(MLIM_STRATEGY=="uniform") NA else NULL),
       lty=c(NA,1,2,1,3, if(MLIM_STRATEGY=="uniform") 4 else NULL),
       lwd=c(NA,3,2,1,2, if(MLIM_STRATEGY=="uniform") 2 else NULL),
       bg="white", cex=0.7)

text(14.5, -2.3,
     sprintf("MAP:\nM* = %.2f\nlog\u03c6* = %.2f\n\u03b1 = %.2f\n\u03b2 = %.2f",
             map_ms, map_lp, map_al, map_be),
     col="red", cex=0.8, pos=4)

# Panel 2: Parameter comparison
plot(1:4, c(map_ms, map_lp, map_al, map_be), pch=19, col="red", cex=2, 
     ylim=range(c(map_ms, map_lp, map_al, map_be, 
                  medmc, driver22)),
     xlab="", ylab="Parameter value", main="MAP vs MCMC vs Driver+22", xaxt="n")
points(1:4+0.1, medmc, pch=17, col="purple", cex=2)
points(1:4-0.1, driver22, pch=15, col="black", cex=2)
axis(1, at=1:4, labels=c("M*","log_phi*","alpha","beta"))
legend("topleft", c("MAP","MCMC","Driver+22"), 
       pch=c(19,17,15), col=c("red","purple","black"), cex=0.9)

# Panel 3: Latent mass shift
plot(x_obs, map_mt, pch=".", col=rgb(0,0,0,0.3),
     xlab="Observed mass", ylab="MAP latent mass",
     main="Latent mass recovery")
abline(0, 1, col="red", lwd=2)
if(MLIM_STRATEGY == "uniform") abline(h=UNIFORM_MCUT, col="orange", lty=2)

# Panel 4: Mass histogram with cut
hist(g3c$log_mass[is.finite(g3c$log_mass)], breaks=seq(10,17,by=0.2), 
     col=rgb(0.7,0.7,0.7,0.5), border="white",
     main="Sample selection", xlab="log10(M)")
hist(x_obs, breaks=seq(10,17,by=0.2), 
     col=rgb(0.2,0.5,0.2,0.5), border="white", add=TRUE)
if(MLIM_STRATEGY == "uniform") abline(v=UNIFORM_MCUT, col="orange", lwd=3, lty=2)
legend("topright", c("All groups", "Used in fit",
       if(MLIM_STRATEGY=="uniform") sprintf("Cut = %.1f", UNIFORM_MCUT) else NULL),
       fill=c(rgb(0.7,0.7,0.7,0.5), rgb(0.2,0.5,0.2,0.5), NA),
       border=c("white","white",NA),
       col=c(NA, NA, if(MLIM_STRATEGY=="uniform") "orange" else NA),
       lty=c(NA, NA, if(MLIM_STRATEGY=="uniform") 2 else NA),
       lwd=c(NA, NA, if(MLIM_STRATEGY=="uniform") 3 else NA),
       cex=0.8)

# Panel 5: Error distribution
hist(sigma_obs, breaks=30, col="steelblue", border="white",
     main="Mass error distribution", xlab="sigma (dex)")
abline(v=median(sigma_obs), col="red", lwd=2, lty=2)

# Panel 6: MCMC trace for alpha
if(nrow(pmmc) > 10) {
    hist(pmmc[,"al"], breaks=30, col="steelblue", border="white", freq=FALSE,
         main="alpha posterior (MCMC)", xlab="alpha")
    abline(v=map_al, col="red", lwd=3)
    abline(v=driver22["al"], col="black", lwd=2, lty=3)
    legend("topright", c("MAP","Driver+22"), col=c("red","black"), 
           lty=c(1,3), lwd=c(3,2), cex=0.8)
}

dev.off()
cat("\nPlot saved: MRP_HIERARCHICAL_FIT.pdf\n")

############################################################
# 10. Save results
############################################################

saveRDS(list(
    strategy = MLIM_STRATEGY,
    uniform_mcut = if(MLIM_STRATEGY=="uniform") UNIFORM_MCUT else NA,
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
cat(sprintf("\nDone! (Strategy: %s, N=%d groups)\n", MLIM_STRATEGY, N))
