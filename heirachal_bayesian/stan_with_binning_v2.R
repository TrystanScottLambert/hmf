# ------------------------------------------------------------
# GAMA HMF: Unbinned Stan fit (MRP shape-only) + φ* from 1/Vmax
# Method follows Wright+2017 Sec.3 (unbinned, truncated likelihood)
# and compares to Driver+2022 GAMA5 MRP (z_eff ~ 0.1). 
# Refs: Wright+2017 (MNRAS 470, Sec.3); Driver+2022 (MNRAS 515). 
# ------------------------------------------------------------

library(rstan)
library(Rfits)
library(data.table)
library(celestial)
library(plotrix)
library(Cairo)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(42)

##################################################
# Load GAMA data (same as original script)
##################################################

logbin = 0.2
multi  = 5
zlimit = 0.25
zmin   = 0.015

ho     = 67.37
omegam = 0.3147
omegal = 1 - omegam

g3cx = Rfits_read_table("../data/G3CFoFGroupv10.fits")

g3c = g3cx[
  Nfof > multi-1 &
  Zfof < zlimit &
  Zfof > zmin &
  MassAfunc > 1E1 &
  IterCenDec > -3.5
]

parsec = 3.0857E16
G      = 6.67408E-11
msol   = 1.988E30
magica = 13.9

# Dynamical mass with A=13.9 (MassA) and correction table -> MassAfunc analogue
g3c$mymass =
  magica*(g3c$VelDisp*1000)^2 *
  g3c$Rad50*parsec*1E6/(G*msol) * (100/ho)

masscorr = c(
  0.0,0.0,-0.2672595,-0.1513503,-0.1259069,-0.09006064,
  -0.05466009,-0.06666895,-0.01988694,-0.02439581,
  -0.0206706,-0.01812964,-0.01556899,-0.01313664,
  -0.01743112,-0.007965513,-0.01257178,-0.007064037,
  -0.003963656,-0.01271533,-0.002664687,-0.001691287
)

g3c$masscorr = masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] = 0
g3c$MassAfunc = g3c$mymass / 10^g3c$masscorr

##################################################
# Per-halo mass error (dex) for Eddington bias
# (matches the mapping used in Driver+ scripts)
##################################################
xx = seq(3,22)
yy = c(0.68389355,0.38719116,0.40325591,0.32696735,0.27680685,
       0.24018684,0.20226682,0.18645475,0.17437005,0.14271506,
       0.13922450,0.13482418,0.13741619,0.11715141,0.12134983,
       0.10078830,0.09944761,0.09913166,0.08590223,0.07588408)

g3c$log10MassErr = approx(xx, yy, g3c$Nfof)$y
g3c$log10MassErr[is.na(g3c$log10MassErr)] = 0.03
g3c$log10MassErr[g3c$log10MassErr < 0.10] = 0.10  # floor as in paper/code

##################################################
# Vmax calculation (as in your script)
##################################################
gig = fread("../data/GAMAGalsInGroups.csv")

g3c$zmax = NA_real_
for (i in 1:nrow(g3c)) {
  if (g3c$Nfof[i] == 2) {
    zc = sort(gig[GroupID == g3c$GroupID[i], zmax_19p8], decreasing=TRUE)[2]
  } else {
    zc = sort(gig[GroupID == g3c$GroupID[i], zmax_19p8], decreasing=TRUE)[multi]
  }
  g3c$zmax[i] = zc
}

g3c$zmax  = pmin(g3c$zmax, zlimit)

g3c$vmax  = 179.92/(360^2/pi)*1E9 *
  cosdist(g3c$zmax, OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol

vlimit = 179.92/(360^2/pi)*1E9 *
  cosdist(zlimit, OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol

# Stabilize extremes (as in your earlier code)
vlimitmin = vlimit/1000
v_eff = pmin(pmax(g3c$vmax, vlimitmin), vlimit)

##################################################
# (Only for plotting) Vmax-weighted binned HMF
##################################################
massx   = seq(10.3, 16.1, logbin)
gamahmf2 = weighted.hist(log10(g3c$MassAfunc), w = 1/g3c$vmax,
                         breaks = massx, plot = FALSE)
gamax   = gamahmf2$mids
gamay   = gamahmf2$counts / logbin

##################################################
# Build unbinned data for Stan
##################################################
keep = is.finite(g3c$MassAfunc) & (g3c$Nfof >= multi) &
       is.finite(g3c$log10MassErr) &
       is.finite(v_eff) & (g3c$Zfof > zmin) & (g3c$Zfof < zlimit)

x_obs   = log10(g3c$MassAfunc[keep])
sigma_x = g3c$log10MassErr[keep]
v_eff_k = v_eff[keep]

# Truncation mass floor (GAMA5 completeness)
x_min = 12.7

# Optional shape weights (selection/density). Here: Vmax-based weight.
w_shape = (vlimit / v_eff_k)
# You can clamp if desired to reduce leverage of a few small-Vmax systems:
# w_shape = pmin(w_shape, vlimit/vlimitmin)

stan_data <- list(
  N = length(x_obs),
  x_obs   = x_obs,
  sigma_x = sigma_x,
  w_shape = w_shape,
  x_min   = x_min,
  mstar_prior_mean = 14.1, mstar_prior_sd = 0.7,
  alpha_prior_mean = -1.7, alpha_prior_sd = 0.5,
  beta_prior_mean  = 0.6,  beta_prior_sd  = 0.3
)

##################################################
# Stan model (written to a local .stan file)
##################################################
stan_code <- '
functions {
  real ln10() { return log(10.0); }

  // log of unnormalized MRP density over log10 M
  real mrp_log10_unnorm_lpdf(real x, real mstar, real alpha, real beta) {
    real t = pow(10.0, beta * (x - mstar));        // u = 10^{beta(x-m*)}
    return log(beta) + log(ln10()) + ln10() * (alpha + 1.0) * (x - mstar) - t;
  }

  // log normalization Z(theta; x_min) = Gamma(s, u_min) with s=(alpha+1)/beta
  real log_Z_shape(real mstar, real alpha, real beta, real x_min) {
    real s = (alpha + 1.0) / beta;
    real u_min = pow(10.0, beta * (x_min - mstar));
    // upper incomplete gamma: Gamma(s, u_min) = tgamma(s) * gamma_q(s, u_min)
    real Z = tgamma(s) * gamma_q(s, u_min);
    return log(Z);
  }
}

data {
  int<lower=1> N;
  vector[N] x_obs;                // observed log10(M)
  vector<lower=0>[N] sigma_x;     // mass errors in dex
  vector<lower=0>[N] w_shape;     // per-object shape weights
  real x_min;                     // truncation floor in log10(M)
  // weak priors
  real mstar_prior_mean;
  real mstar_prior_sd;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real beta_prior_mean;
  real beta_prior_sd;
}

parameters {
  real mstar;
  real alpha;
  real<lower=0> beta;
  vector<lower=x_min>[N] x_true;  // latent true log10(M), truncated below x_min
}

model {
  // Priors
  mstar ~ normal(mstar_prior_mean, mstar_prior_sd);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  beta  ~ normal(beta_prior_mean,  beta_prior_sd);

  // Precompute normalization over [x_min, ∞)
  real logZ = log_Z_shape(mstar, alpha, beta, x_min);

  // Likelihood
  for (n in 1:N) {
    target += w_shape[n] * ( mrp_log10_unnorm_lpdf(x_true[n] | mstar, alpha, beta) - logZ );
    x_obs[n] ~ normal(x_true[n], sigma_x[n]);   // Eddington bias generative term
  }
}
'

stan_file <- "mrp_unbinned_trunc_halos.stan"
writeLines(stan_code, con = stan_file)

##################################################
# Fit model
##################################################
fit <- stan(
  file = stan_file,
  data = stan_data,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 200
)

print(fit, pars = c("mstar","alpha","beta"))

post <- rstan::extract(fit)

##################################################
# φ* from 1/Vmax and analytic normalization Z(θ, x_min)
##################################################
logZ_R <- function(mstar, alpha, beta, x_min){
  s <- (alpha + 1)/beta
  u <- 10^( beta*(x_min - mstar) )
  # upper incomplete gamma = Gamma(s)*(1 - P(s,u)), with P the lower regularized gamma
  Z <- gamma(s) * (1 - pgamma(u, shape = s, lower.tail = TRUE))
  log(Z)
}

sum_invV <- sum( 1.0 / v_eff_k )

phi_star_draws <- numeric(length(post$mstar))
for(i in seq_along(phi_star_draws)){
  lZ <- logZ_R(post$mstar[i], post$alpha[i], post$beta[i], x_min)
  phi_star_draws[i] <- sum_invV / exp(lZ)  # Mpc^-3 dex^-1
}

posterior_matrix <- cbind(
  mstar   = post$mstar,
  log_phi = log10(phi_star_draws),
  alpha   = post$alpha,
  beta    = post$beta
)

##################################################
# Posterior summaries
##################################################
med <- apply(posterior_matrix, 2, median)
q16 <- apply(posterior_matrix, 2, quantile, 0.16)
q84 <- apply(posterior_matrix, 2, quantile, 0.84)

##################################################
# Plotting (keep your original style)
##################################################
mrp_log10 <- function(mstar, log_phi, alpha, beta, logM){
  phi <- 10^log_phi
  log10( phi * beta * log(10) *
           10^((alpha+1)*(logM-mstar)) *
           exp(-10^(beta*(logM-mstar))) )
}

xfit <- seq(11, 16, 0.01)

# Driver+2022 GAMA5 reference (z_eff~0.1)
y_ref <- mrp_log10(13.51, -3.19, -1.27, 0.47, xfit)

# Median curve from our unbinned fit
y_med <- mrp_log10(med["mstar"], med["log_phi"], med["alpha"], med["beta"], xfit)

CairoPDF("MRP_STAGE3_VMAX.pdf", 10, 7)

plot(gamax, log10(gamay),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(11,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 3: Vmax Weighting (Incompleteness Corrected)")

grid(col="gray80")

# Posterior spaghetti
set.seed(1)
idx <- sample(1:nrow(posterior_matrix), min(300, nrow(posterior_matrix)))
for(i in idx){
  y <- mrp_log10(posterior_matrix[i,"mstar"],
                 posterior_matrix[i,"log_phi"],
                 posterior_matrix[i,"alpha"],
                 posterior_matrix[i,"beta"],
                 xfit)
  lines(xfit, y, col=rgb(0,0,1,0.03))
}

# Median fit
lines(xfit, y_med, col="red", lwd=3)

# Driver+2022 reference
lines(xfit, y_ref, col="black", lwd=2, lty=2)

# Legend
legend("bottomleft",
       legend=c("GAMA (Vmax weighted)",
                "Stan fit (median)",
                "Posterior samples",
                "Driver+22 GAMA5"),
       col=c("darkgreen","red",rgb(0,0,1,0.3),"black"),
       pch=c(19,NA,NA,NA),
       lty=c(NA,1,1,2),
       lwd=c(NA,3,1,2),
       bg="white", cex=0.9)

# Parameter text box
txt <- sprintf(
  "M* = %.2f (+%.2f/−%.2f)\nlog\u03c6* = %.2f (+%.2f/−%.2f)\n\u03b1 = %.2f (+%.2f/−%.2f)\n\u03b2 = %.2f (+%.2f/−%.2f)",
  med["mstar"],  q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"],
  med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"],
  med["alpha"],  q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"],
  med["beta"],   q84["beta"]-med["beta"],     med["beta"]-q16["beta"]
)

text(14.3,-2.3, txt, col="red", cex=0.85, pos=4)

dev.off()

cat("Saved to MRP_STAGE3_VMAX.pdf\n")

