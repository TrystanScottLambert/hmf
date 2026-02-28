############################################################
# GAMA HALO MASS FUNCTION
# Hierarchical Bayesian Fit using Correct Poisson Process
# Latent true masses with measurement uncertainty
# Shell-based Lambda normalization
# MAP optimization via Stan
############################################################

rm(list=ls())

############################################################
# Libraries
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(rstan)
library(Cairo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# Cosmology and survey parameters
############################################################

H0     <- 67.37
OmegaM <- 0.3147
z_min  <- 0.01
z_max  <- 0.25

sky_area_deg <- 179.92
sky_area_sr  <- sky_area_deg * (pi/180)^2

############################################################
# Load GAMA group catalogue
############################################################

GROUP_FILE <- "../data/G3CFoFGroupv10.fits"

cat("Loading GAMA groups...\n")

g3cx <- Rfits_read_table(GROUP_FILE)

# Selection criteria
multi_min <- 5

mask <-
    g3cx$Nfof >= multi_min &
    g3cx$Zfof > z_min &
    g3cx$Zfof < z_max &
    g3cx$MassAfunc > 1e10


g3c <- g3cx[mask]

cat("Groups after selection:", length(g3c$MassAfunc), "\n")

############################################################
# Observed masses and uncertainties
############################################################

x_obs <- log10(g3c$MassAfunc)

# mass uncertainty model (adjust if better estimate exists)
sigma_obs <- 0.12 + 0.05*(x_obs - 13.5)^2

############################################################
# Limiting mass for each group
############################################################

z_group <- g3c$Zfof

m_lim_obs <- log10(
    g3c$MassAfunc *
    (cosdistLumDist(z_max, OmegaM=OmegaM, H0=H0) /
     cosdistLumDist(z_group, OmegaM=OmegaM, H0=H0))^2
)

############################################################
# Construct redshift shells for Lambda normalization
############################################################

N_shell <- 25

z_edges <- seq(z_min, z_max, length.out=N_shell+1)
z_mid   <- 0.5*(z_edges[-1] + z_edges[-length(z_edges)])

V_shell <- numeric(N_shell)
mlim_shell <- numeric(N_shell)

for(j in 1:N_shell){

    r1 <- cosdistCoDist(z_edges[j], OmegaM=OmegaM, H0=H0)
    r2 <- cosdistCoDist(z_edges[j+1], OmegaM=OmegaM, H0=H0)

    V_shell[j] <-
        (4/3)*pi*(r2^3 - r1^3) * sky_area_sr/(4*pi)

    sel <- z_group >= z_edges[j] & z_group < z_edges[j+1]

    if(any(sel))
        mlim_shell[j] <- median(m_lim_obs[sel])
    else
        mlim_shell[j] <- min(m_lim_obs)
}

cat("Total shell volume:", signif(sum(V_shell),4),"Mpc^3\n")

############################################################
# Stan hierarchical model
############################################################

stan_model_code <- "
data {

  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sigma;
  vector[N] mlim;

  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;

  real xhi;
  int<lower=10> Ng;
}

transformed data {

  real ln10 = log(10);

  real xlo = min(mlim_sh) - 0.5;

  real dx = (xhi - xlo)/(Ng-1);

  vector[Ng] xg;

  for(k in 1:Ng)
    xg[k] = xlo + (k-1)*dx;
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

  for(i in 1:N)
    mt[i] = x_obs[i] + sigma[i]*z_raw[i];
}
model {

  ms ~ normal(14,1.5);
  lp ~ normal(-4,2);
  al ~ normal(-1.5,1);
  z_raw ~ std_normal();

  vector[Ng] phi;

  for(k in 1:Ng){

    real u = xg[k] - ms;

    phi[k] =
      be*ln10 *
      pow(10,lp) *
      pow(10,(al+1)*u) *
      exp(-pow(10,be*u));
  }

  vector[Ng] cum;

  cum[Ng]=0;

  for(kk in 1:(Ng-1)){

    int k = Ng-kk;

    cum[k] =
      cum[k+1] +
      0.5*(phi[k]+phi[k+1])*(xg[k+1]-xg[k]);
  }

  real Lambda=0;

  for(j in 1:Nsh){

    int k0=1;

    for(k in 1:Ng)
      if(xg[k]<=mlim_sh[j])
        k0=k;

    if(k0>=Ng)
      k0=Ng-1;

    Lambda += V_sh[j]*cum[k0];
  }

  target += -Lambda;

  for(i in 1:N){

    real u = mt[i]-ms;

    real phi_i =
      be*ln10 *
      pow(10,lp) *
      pow(10,(al+1)*u) *
      exp(-pow(10,be*u));

    if(mt[i]>mlim[i] && phi_i>1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

############################################################
# Prepare data for Stan
############################################################

stan_data <- list(

  N = length(x_obs),
  x_obs = x_obs,
  sigma = sigma_obs,
  mlim = m_lim_obs,

  Nsh = N_shell,
  V_sh = V_shell,
  mlim_sh = mlim_shell,

  xhi = 16.5,
  Ng = 200
)

############################################################
# Compile model
############################################################

cat("Compiling Stan model...\n")

sm <- stan_model(model_code=stan_model_code)

############################################################
# Run optimization (MAP estimate)
############################################################

cat("Running MAP optimization...\n")

fit <- optimizing(

  sm,
  data=stan_data,

  init=list(

    ms=13.5,
    lp=-3.5,
    al=-1.5,
    be=0.6,
    z_raw=rep(0,length(x_obs))
  ),

  hessian=TRUE,
  as_vector=FALSE
)

############################################################
# Results
############################################################

cat("\nRecovered parameters:\n")
print(fit$par)

############################################################
# Plot
############################################################

ms <- fit$par$ms
lp <- fit$par$lp
al <- fit$par$al
be <- fit$par$be

ln10 <- log(10)

phi_fun <- function(x){

  u <- x-ms

  be*ln10*
  10^lp*
  10^((al+1)*u)*
  exp(-10^(be*u))
}

xgrid <- seq(11,16,length=200)

CairoPDF("HMF_fit.pdf", width=6, height=5)

hist(

  x_obs,
  breaks=30,
  freq=FALSE,
  col="grey",
  border="white",
  main="GAMA Halo Mass Function",
  xlab="log10(Mhalo)"
)

lines(xgrid, phi_fun(xgrid), lwd=3, col="red")

dev.off()

cat("Saved plot to HMF_fit.pdf\n")

############################################################
# END
############################################################

