############################################################
# HIERARCHICAL BAYESIAN HMF
# Upgrades the binned chi-squared approach in do_what_simon_did.R
# to a proper Poisson point-process likelihood with exact
# Eddington bias treatment via Gauss-Hermite marginalisation.
#
# Key change vs. old script:
#   OLD: bin halos -> Vmax-weight bins -> chi-squared fit to log(phi)
#        (Eddington correction applied as a separate forward-MC factor)
#   NEW: use every halo individually, marginalise over true mass,
#        Eddington bias falls out automatically from the integral.
#
# Likelihood follows Murray, Robotham & Power (2018) eq. 16:
#   ln L = - sum_i [V_i * INT phi(m) dm]          (global Poisson term)
#           + sum_i ln [V_i * INT phi(m)*N(m_obs|m,sig) dm]  (per-halo)
#
# MRP function: phi(m) = ln10 * phi* * beta * (M/M*)^(alpha+1) * exp(-(M/M*)^beta)
# Parameters: mstar = log10(M*/Msol), log_phi = log10(phi* [Mpc^-3 dex^-1]),
#             alpha (low-mass slope), beta (high-mass softening)
#
# Reference: Driver+22 GSR best fit: M*=14.13, log_phi=-3.96, alpha=-1.68, beta=0.63
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
library(plotrix)
library(Cairo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

############################################################
# DATA PREPARATION - matching gamahmf.r exactly
############################################################

ho         <- 67.37
omegam     <- 0.3147
omegal     <- 1 - omegam
zlimit     <- 0.25
zmin       <- 0.015
multi      <- 5
mass_limit <- 12.7

parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30
magica <- 13.9

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit &
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 &
             g3cx$IterCenDec > -3.5, ]

g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$mymass    <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

# Mass uncertainty as a function of NFoF (Driver+22 Fig 3)
xx <- seq(3, 22)
yy <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685,
        0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506,
        0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983,
        0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)

g3c$log10MassErr <- approx(xx, yy, g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)] <- 0.03
g3c$log10MassErr[g3c$log10MassErr < 0.1]  <- 0.1   # floor

# Low-multiplicity bias correction
masscorr <- c(0.0, 0.0, -2.672595e-01, -1.513503e-01, -1.259069e-01,
              -9.006064e-02, -5.466009e-02, -6.666895e-02, -1.988694e-02,
              -2.439581e-02, -2.067060e-02, -1.812964e-02, -1.556899e-02,
              -1.313664e-02, -1.743112e-02, -7.965513e-03, -1.257178e-02,
              -7.064037e-03, -3.963656e-03, -1.271533e-02, -2.664687e-03,
              -1.691287e-03)

g3c$masscorr    <- masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] <- 0.0
g3c$mymasscorr  <- g3c$mymass / 10^g3c$masscorr
g3c$MassAfunc   <- g3c$mymasscorr

# Flag bad group
g3c$MassAfunc[g3c$GroupID == 100622] <- 1E9

############################################################
# VMAX CALCULATION
############################################################

gig <- fread("../data/GAMAGalsInGroups.csv")

g3c$zmax <- NA
for (i in 1:nrow(g3c)) {
    if (g3c$Nfof[i] == 2) {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8],
                            decreasing=TRUE)[2]
    } else {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8],
                            decreasing=TRUE)[multi]
    }
}

g3c$zmax <- ifelse(g3c$zmax < g3c$Zfof, g3c$Zfof, g3c$zmax)
g3c$zmax <- ifelse(g3c$zmax > zlimit,   zlimit,     g3c$zmax)

vol_zmax <- cosdist(g3c$zmax, OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol
vol_zmin <- cosdist(zmin,     OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol
g3c$vmax <- 179.92/(360^2/pi) * 1E9 * (vol_zmax - vol_zmin)

vlimit    <- 179.92/(360^2/pi) * 1E9 *
             cosdist(zlimit, OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol
vlimitmin <- vlimit / 1000.0

g3c$weightszlimit <- ifelse(g3c$vmax > vlimit,    vlimit,    g3c$vmax)
g3c$weightszlimit <- ifelse(g3c$vmax < vlimitmin, vlimitmin, g3c$vmax)

# Apply mass completeness cut - work in log10 mass from here
g3c$log10m <- log10(g3c$MassAfunc)
g3c        <- g3c[g3c$log10m >= mass_limit, ]

cat("N halos after cuts:", nrow(g3c), "\n")
cat("log10(M) range:", round(range(g3c$log10m), 2), "\n")
cat("sigma_m range:",  round(range(g3c$log10MassErr), 3), "\n\n")

############################################################
# GAUSS-HERMITE QUADRATURE NODES & WEIGHTS (n=7)
# Approximates: integral f(t) exp(-t^2) dt ~ sum_k w_k f(t_k)
# Change of variable: m_true = m_obs + sqrt(2)*sigma*t
# gives: integral phi(m)*Normal(m_obs|m,sig) dm
#       ~ (sqrt(2)*sig/sqrt(pi)) * sum_k w_k * phi(m_obs + sqrt(2)*sig*t_k)
############################################################

gh_nodes   <- c(-2.6519613987, -1.6735516287, -0.8162877999,
                 0.0,
                 0.8162877999,  1.6735516287,  2.6519613987)
gh_weights <- c(9.717812451e-4, 5.45155828191e-2, 0.4256072526,
                0.8102646175,
                0.4256072526,   5.45155828191e-2, 9.717812451e-4)

############################################################
# STAN MODEL
############################################################

hierarchical_stan <- "
functions {

  // Upper incomplete gamma via downward recurrence:
  // Gamma(a, x) = [Gamma(a+1, x) - x^a * exp(-x)] / a
  // Recurse until a > 0, then use Stan's gamma_q * tgamma.
  // Needed because alpha ~ -1.68 -> (alpha+2)/beta ~ 0.5
  real upper_inc_gamma(real a, real x) {
    if (a > 1e-6) {
      return gamma_q(a, x) * tgamma(a);
    } else {
      real g1 = upper_inc_gamma(a + 1.0, x);
      return (g1 - pow(x, a) * exp(-x)) / a;
    }
  }

  // MRP log-density (natural log) in log10-mass space.
  // phi(m) = ln10 * phi* * beta * (M/M*)^(alpha+1) * exp(-(M/M*)^beta)
  // where M = 10^m, M* = 10^mstar.
  real log_mrp(real m, real mstar, real log_phi, real alpha, real beta) {
    real x  = m - mstar;
    return log(log(10.0))
         + log_phi * log(10.0)
         + log(beta)
         + (alpha + 1.0) * log(10.0) * x
         - pow(10.0, beta * x);
  }

  // Integral of MRP from mmin to infinity:
  //   INT phi(m) dm = phi* * Gamma((alpha+2)/beta, (10^(mmin-mstar))^beta)
  real mrp_integral(real mmin, real mstar, real log_phi,
                    real alpha, real beta) {
    real a     = (alpha + 2.0) / beta;
    real x_min = pow(10.0, beta * (mmin - mstar));
    return pow(10.0, log_phi) * upper_inc_gamma(a, x_min);
  }

  // Log of the mass-marginalised likelihood for one halo:
  //   ln INT [ phi(m) * Normal(m_obs | m, sigma) ] dm
  // via 7-point Gauss-Hermite quadrature.
  real log_mrp_marginalised(real m_obs, real sigma,
                             real mstar, real log_phi, real alpha, real beta,
                             vector gh_nodes, vector gh_weights) {
    int K = num_elements(gh_nodes);
    real log_sum = negative_infinity();
    real sqrt2sig = sqrt(2.0) * sigma;
    for (k in 1:K) {
      real m_true = m_obs + sqrt2sig * gh_nodes[k];
      if (m_true > 12.5) {
        real lval = log_mrp(m_true, mstar, log_phi, alpha, beta)
                    + log(gh_weights[k]);
        log_sum = log_sum_exp(log_sum, lval);
      }
    }
    return log_sum + 0.5*log(2.0) + log(sigma) - 0.5*log(pi());
  }

}

data {
  int<lower=1> N;
  vector[N] m_obs;
  vector<lower=0>[N] sig_m;
  vector<lower=0>[N] vmax;

  real mmin;
  real vlimit;
  real dlogm;

  int<lower=1> K;
  vector[K] gh_nodes;
  vector<lower=0>[K] gh_weights;

  int<lower=0> N_penalty;
  vector[N_penalty] m_penalty;
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.05, upper=2.5> beta;
}

model {
  // Weakly informative priors centred on Driver+22 GSR
  mstar   ~ normal(14.13, 1.5);
  log_phi ~ normal(-3.96, 2.0);
  alpha   ~ normal(-1.68, 0.5);
  beta    ~ normal(0.63,  0.3);

  // Global Poisson term: -E[N_total] = -(sum V_i) * INT phi(m) dm
  {
    real phi_int = mrp_integral(mmin, mstar, log_phi, alpha, beta);
    target += -sum(vmax) * phi_int;
  }

  // High-mass empty bin penalty (exact Poisson zero-count probability)
  for (p in 1:N_penalty) {
    target += -exp(log_mrp(m_penalty[p], mstar, log_phi, alpha, beta))
              * vlimit * dlogm;
  }

  // Per-halo marginalised likelihood
  for (i in 1:N) {
    target += log(vmax[i])
            + log_mrp_marginalised(m_obs[i], sig_m[i],
                                   mstar, log_phi, alpha, beta,
                                   gh_nodes, gh_weights);
  }
}

generated quantities {
  vector[42] m_grid;
  vector[42] phi_grid;
  for (g in 1:42) {
    m_grid[g]   = 12.0 + (g-1) * 0.1;
    phi_grid[g] = exp(log_mrp(m_grid[g], mstar, log_phi, alpha, beta));
  }
}
"

############################################################
# PREPARE STAN DATA
############################################################

dlogm     <- 0.2
m_max_obs <- max(g3c$log10m)
m_penalty <- seq(m_max_obs + dlogm, m_max_obs + 10*dlogm, by=dlogm)

stan_data <- list(
    N         = nrow(g3c),
    m_obs     = g3c$log10m,
    sig_m     = g3c$log10MassErr,
    vmax      = g3c$weightszlimit,
    mmin      = mass_limit,
    vlimit    = vlimit,
    dlogm     = dlogm,
    K         = 7L,
    gh_nodes  = gh_nodes,
    gh_weights = gh_weights,
    N_penalty = length(m_penalty),
    m_penalty = m_penalty
)

cat("=== Stan data ===\n")
cat("N halos:", stan_data$N, "\n")
cat("N penalty bins:", stan_data$N_penalty,
    "(from", round(m_penalty[1],2), "to", round(tail(m_penalty,1),2), ")\n\n")

############################################################
# COMPILE AND RUN STAN
############################################################

cat("Compiling Stan model...\n")
sm <- stan_model(model_code = hierarchical_stan,
                 model_name  = "HMF_hierarchical")

cat("Running MCMC (4 chains x 3000 iter, expect 20-60 min)...\n")
fit <- sampling(
    sm,
    data    = stan_data,
    chains  = 4,
    iter    = 3000,
    warmup  = 1500,
    cores   = 4,
    control = list(adapt_delta = 0.92, max_treedepth = 12)
)

############################################################
# DIAGNOSTICS
############################################################

cat("\n=== Posterior summary ===\n")
print(fit, pars = c("mstar","log_phi","alpha","beta"),
      probs = c(0.16, 0.50, 0.84))

sp          <- get_sampler_params(fit, inc_warmup=FALSE)
n_div       <- sum(sapply(sp, function(x) sum(x[,"divergent__"])))
n_treedepth <- sum(sapply(sp, function(x) sum(x[,"treedepth__"] >= 12)))
cat("Divergent transitions:", n_div, "\n")
cat("Treedepth saturations:", n_treedepth, "\n")
if (n_div > 0) cat("WARNING: increase adapt_delta if divergences are numerous.\n")

############################################################
# EXTRACT POSTERIOR
############################################################

post <- extract(fit, pars = c("mstar","log_phi","alpha","beta"))

pnames <- c("mstar","log_phi","alpha","beta")
med    <- sapply(post[pnames], median)
q16    <- sapply(post[pnames], quantile, 0.16)
q84    <- sapply(post[pnames], quantile, 0.84)

cat("\n=== RESULTS vs Driver+22 GSR ===\n")
cat(sprintf("%-10s  %7s  +1sig  -1sig  | Driver+22\n","Parameter","Median"))
cat(sprintf("%-10s  %7.3f  %5.3f  %5.3f  | 14.13\n",
    "M*",      med["mstar"],   q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"]))
cat(sprintf("%-10s  %7.3f  %5.3f  %5.3f  | -3.96\n",
    "log_phi", med["log_phi"], q84["log_phi"]-med["log_phi"], med["log_phi"]-q16["log_phi"]))
cat(sprintf("%-10s  %7.3f  %5.3f  %5.3f  | -1.68\n",
    "alpha",   med["alpha"],   q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"]))
cat(sprintf("%-10s  %7.3f  %5.3f  %5.3f  | 0.63\n",
    "beta",    med["beta"],    q84["beta"]-med["beta"],     med["beta"]-q16["beta"]))

############################################################
# PLOT: HMF with posterior samples
############################################################

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
          10^((alpha+1)*(x-mstar)) *
          exp(-10^(beta*(x-mstar))))
}

# Raw Vmax-weighted bins for context (NOT Eddington-corrected)
logbin   <- 0.2
massx    <- seq(10.3, 16.1, logbin)
gamahmf2 <- weighted.hist(g3c$log10m, w=1/g3c$weightszlimit,
                          breaks=massx, plot=FALSE)
raw_x  <- gamahmf2$mids
raw_y  <- gamahmf2$counts / logbin
ok_raw <- raw_y > 0 & !is.na(raw_y) & raw_x >= mass_limit

xfit       <- seq(12, 16.5, length.out=500)
driver_gsr <- mrp_log10(14.13, -3.96, -1.68, 0.63, xfit)

CairoPDF("MRP_FINAL.pdf", 10, 7)
par(mar=c(4.5, 4.5, 2.5, 1.5))

plot(0, 0, type="n",
     xlim=c(12, 16.2), ylim=c(-8, -1.8),
     xlab=expression("Halo Mass  " * log[10] * "(M / M"["\u2299"] * ")"),
     ylab=expression(log[10] * "(" * phi * ")  [Mpc"^{-3} * " dex"^{-1} * "]"),
     main="Hierarchical Bayesian HMF (Eddington bias via marginalisation)",
     las=1)
grid(col="gray85", lty=1)

idx_samp <- sample(length(post$mstar), 200)
for (i in idx_samp) {
    yy <- mrp_log10(post$mstar[i], post$log_phi[i],
                    post$alpha[i],  post$beta[i], xfit)
    if (all(is.finite(yy)))
        lines(xfit, yy, col=rgb(0.1, 0.4, 0.9, 0.03), lwd=1)
}

points(raw_x[ok_raw], log10(raw_y[ok_raw]),
       pch=5, col="limegreen", cex=0.9)

lines(xfit,
      mrp_log10(med["mstar"], med["log_phi"], med["alpha"], med["beta"], xfit),
      col="red", lwd=2.5)

lines(xfit, driver_gsr, col="black", lwd=2, lty=2)

legend("bottomleft", bg="white", cex=0.88,
       legend=c("Raw Vmax bins (no EdB correction)",
                "Posterior median (hierarchical Bayes)",
                "Posterior samples",
                "Driver+22 GSR (M*=14.13, \u03b1=-1.68)"),
       col=c("limegreen","red",rgb(0.1,0.4,0.9,0.4),"black"),
       pch=c(5,NA,NA,NA), lty=c(NA,1,1,2), lwd=c(NA,2.5,2,2))

text(14.0, -2.0,
     sprintf(paste0("Hierarchical Bayes\n",
                    "M* = %.2f (+%.2f/\u2212%.2f)\n",
                    "log\u03c6* = %.2f (+%.2f/\u2212%.2f)\n",
                    "\u03b1 = %.2f (+%.2f/\u2212%.2f)\n",
                    "\u03b2 = %.2f (+%.2f/\u2212%.2f)"),
             med["mstar"],   q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"],
             med["log_phi"], q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"],
             med["alpha"],   q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"],
             med["beta"],    q84["beta"]-med["beta"],     med["beta"]-q16["beta"]),
     col="red", cex=0.82, pos=4)

dev.off()
cat("Saved: MRP_FINAL.pdf\n")

############################################################
# PLOT: Corner / covariance (mimics gamahmf.r mymagtri style)
############################################################

post_mat <- cbind(
    mstar   = post$mstar,
    log_phi = post$log_phi,
    alpha   = post$alpha,
    beta    = post$beta
)
plab       <- c("log\u2081\u2080(M*) [M\u2299]",
                "log\u2081\u2080(\u03c6*) [Mpc\u207B\u00b3]",
                "\u03b1", "\u03b2")
driver_ref <- c(14.13, -3.96, -1.68, 0.63)

CairoPDF("MRP_covariance.pdf", 12, 12)
par(mfrow=c(4,4), mar=c(2,2,0.5,0.5), oma=c(2,2,2,1))
for (i in 1:4) {
    for (j in 1:4) {
        if (i == j) {
            d <- density(post_mat[,i], n=256)
            plot(d, main="", xlab="", ylab="", axes=FALSE,
                 col="firebrick", lwd=1.5, xlim=range(d$x))
            axis(1, cex.axis=0.7); box()
            abline(v=med[i],        col="firebrick", lty=2, lwd=1.5)
            abline(v=driver_ref[i], col="black",     lty=2, lwd=1.5)
            if (i == 1) mtext(plab[i], side=3, line=0.2, cex=0.72)
        } else if (j < i) {
            idx2 <- sample(nrow(post_mat), min(3000, nrow(post_mat)))
            plot(post_mat[idx2,j], post_mat[idx2,i],
                 pch=".", col=rgb(0,0,0,0.15), axes=FALSE, xlab="", ylab="")
            axis(1, cex.axis=0.7); axis(2, cex.axis=0.7); box()
            points(med[j], med[i], pch=3, col="firebrick", cex=1.5, lwd=2)
            points(driver_ref[j], driver_ref[i], pch=1, col="black", cex=1.5, lwd=2)
        } else {
            plot.new()
            if (i == 1 && j == 2) {
                legend("center", bty="n", cex=0.85,
                       legend=c("Posterior median", "Driver+22 GSR"),
                       col=c("firebrick","black"), lty=c(2,2), lwd=2)
            }
        }
    }
}
mtext("Parameter covariances  --  Hierarchical Bayes HMF",
      outer=TRUE, cex=0.9, line=0.5)
dev.off()
cat("Saved: MRP_covariance.pdf\n")

cat("\n=== DONE ===\n")
cat("Output files:\n")
cat("  MRP_FINAL.pdf       -- HMF with posterior samples vs Driver+22\n")
cat("  MRP_covariance.pdf  -- parameter covariance corner plot\n")
