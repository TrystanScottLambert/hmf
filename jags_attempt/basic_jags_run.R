############################################################
# G3C UNBINNED MRP FIT
############################################################

library(celestial)
library(devtools)
library(mvtnorm)
library(Cairo)
library(sm)
library(Rfits)
library(magicaxis)
library(data.table)
library(plotrix)
library(foreign)
library(MASS)
library(rjags)
library(coda)

############################################################
# 1. Cosmology
############################################################

ho=67.37
omegam=0.3147
omegal=1-omegam
G=6.67408E-11
parsec=3.0857E16

zlimit=0.25
zmin=0.01
multi=3

############################################################
# 2. Survey Volume
############################################################

vol = cosdist(zlimit,
              OmegaM=omegam,
              OmegaL=omegal,
              H0=ho)$CoVol

area_sr = 179.92 * (pi/180)^2
Vsurvey = vol * area_sr / (4*pi)

############################################################
# 3. Read catalogue
############################################################

g3cx=Rfits_read_table("../data/G3CFoFGroupv10.fits")

g3c=g3cx[
  g3cx$Nfof > multi-1 &
  g3cx$Zfof < zlimit &
  g3cx$Zfof > zmin &
  g3cx$MassAfunc > 1E13 &
  g3cx$IterCenDec > -3.5,
]

g3c$MassA = g3c$MassA * 100 / ho

x_obs = log10(g3c$MassA)
N = length(x_obs)

xmin = min(x_obs)
xmax = max(x_obs)

############################################################
# 4. Write unbinned MRP model
############################################################

model_string <- "
model {

  mstar  ~ dnorm(14.4, 1/0.5^2)
  log_phi ~ dnorm(-3, 1/2^2)
  alpha  ~ dnorm(-1.8, 1/0.5^2)
  beta   ~ dnorm(0.7, 1/0.3^2)

  for (i in 1:N) {

    log_intensity[i] <-
        log(V)
      + log(beta)
      + log(log(10))
      + log_phi
      + (alpha+1)*(x[i] - mstar)
      - pow(10, beta*(x[i] - mstar)) / log(10)

    zeros[i] ~ dpois(-log_intensity[i] + C)
  }

  # integral term

  for (k in 1:Ng) {

    log_phi_grid[k] <-
        log(beta)
      + log(log(10))
      + log_phi
      + (alpha+1)*(xgrid[k] - mstar)
      - pow(10, beta*(xgrid[k] - mstar)) / log(10)

    phi_grid[k] <- exp(log_phi_grid[k])
  }

  integral <- V * sum(phi_grid[]) * dx

  zeros_int ~ dpois(integral + C)
}
"

writeLines(model_string, "mrp_unbinned.jags")

############################################################
# 5. Integration grid
############################################################

Ng = 250
xgrid = seq(xmin, xmax, length.out=Ng)
dx = (xmax - xmin)/(Ng-1)

zeros = rep(0, N)
zeros_int = 0
C = 1e6

jags_data <- list(
  x = x_obs,
  N = N,
  V = Vsurvey,
  zeros = zeros,
  zeros_int = zeros_int,
  C = C,
  Ng = Ng,
  xgrid = xgrid,
  dx = dx
)

############################################################
# 6. Run JAGS
############################################################

model <- jags.model(
  "mrp_unbinned.jags",
  data = jags_data,
  n.chains = 4,
  n.adapt = 5000
)

update(model, 10000)

samples <- coda.samples(
  model,
  variable.names = c("mstar","log_phi","alpha","beta"),
  n.iter = 40000,
  thin = 20
)

posterior <- as.matrix(samples)

############################################################
# 7. Plot MRP over empirical HMF
############################################################

# empirical mass function (for plotting only)
logbin=0.2
breaks = seq(floor(xmin), ceiling(xmax), by=logbin)
hist_data = hist(x_obs, breaks=breaks, plot=FALSE)

phi_emp = hist_data$counts / (Vsurvey * logbin)
x_emp = hist_data$mids

# MRP function
mrp_curve <- function(mstar, log_phi, alpha, beta, x) {
  log(beta) +
  log(log(10)) +
  log_phi +
  (alpha+1)*(x - mstar) -
  (10^(beta*(x - mstar))) / log(10)
}

xfit = seq(xmin, xmax, length.out=300)

CairoPDF("MRP_unbinned_fit.pdf", width=7, height=6)

plot(x_emp, log10(phi_emp),
     pch=16,
     xlab=expression(log[10](M)),
     ylab=expression(log[10](phi)),
     main="Unbinned MRP Fit")

# posterior draws
draws = 1000
idx = sample(1:nrow(posterior), draws)

for(i in idx){
  yfit = mrp_curve(
    posterior[i,"mstar"],
    posterior[i,"log_phi"],
    posterior[i,"alpha"],
    posterior[i,"beta"],
    xfit
  )
  lines(xfit, yfit, col=rgb(0,0,1,0.02))
}

# median curve
med = apply(posterior,2,median)
lines(xfit,
      mrp_curve(med["mstar"],
                med["log_phi"],
                med["alpha"],
                med["beta"],
                xfit),
      col="red", lwd=3)

dev.off()

############################################################
# Save posterior
############################################################

saveRDS(posterior,"mrp_unbinned_chain.rds")

