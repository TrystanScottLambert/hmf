############################################################
# G3C UNBINNED MRP FIT — FULLY FIXED NORMALIZATION VERSION
# CORRECTED VERSION - fixes y-axis scaling issue
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(Cairo)
library(rjags)
library(coda)

############################################################
# Cosmology & survey
############################################################

ho      = 67.37
omegam  = 0.3147
omegal  = 1 - omegam

zlimit  = 0.25
zmin    = 0.01
multi   = 3

vol      = cosdist(zlimit, OmegaM=omegam, OmegaL=omegal, H0=ho)$CoVol
area_sr  = 179.92 * (pi/180)^2
Vsurvey  = vol * area_sr / (4*pi)

############################################################
# Data
############################################################

g3cx = Rfits_read_table("../data/G3CFoFGroupv10.fits")

g3c = g3cx[
    g3cx$Nfof > multi-1 &
        g3cx$Zfof < zlimit &
        g3cx$Zfof > zmin &
        g3cx$MassAfunc > 1E9 &
        g3cx$IterCenDec > -3.5,
]

g3c$MassA = g3c$MassA * 100 / ho
g3c = g3c[is.finite(g3c$MassA) & g3c$MassA > 0, ]

x_obs = log10(g3c$MassA)
N     = length(x_obs)

xmin  = min(x_obs)
xmax  = max(x_obs)

############################################################
# JAGS MODEL
############################################################

model_string <- "
model {

  mstar   ~ dnorm(14.4, 1/0.5^2)
  log_phi ~ dnorm(-6,   1/2^2)
  alpha   ~ dnorm(-1.8, 1/0.5^2)
  beta    ~ dnorm(0.7,  1/0.3^2) T(0,)

  # ----- Likelihood -----
  for (i in 1:N) {

    phi[i] <-
        beta * log(10)
      * exp(log_phi)
      * pow(10, (alpha + 1) * (x[i] - mstar))
      * exp(- pow(10, beta * (x[i] - mstar)))

    lambda[i] <- V * phi[i]
  }

  # ----- Integral term -----
  for (k in 1:Ng) {

    phi_grid[k] <-
        beta * log(10)
      * exp(log_phi)
      * pow(10, (alpha + 1) * (xgrid[k] - mstar))
      * exp(- pow(10, beta * (xgrid[k] - mstar)))

    lambda_grid[k] <- V * phi_grid[k]
  }

  Lambda <- sum(lambda_grid[]) * dx

  # ----- Poisson point-process likelihood -----
  logL <- sum(log(lambda[])) - Lambda

  zeros ~ dpois(C - logL)
}
"

writeLines(model_string, "mrp_unbinned.jags")

############################################################
# Integration grid
############################################################

Ng    = 500
xgrid = seq(xmin, xmax, length.out=Ng)
dx    = diff(xgrid)[1]

jags_data <- list(
    x      = x_obs,
    N      = N,
    V      = Vsurvey,
    Ng     = Ng,
    xgrid  = xgrid,
    dx     = dx,
    zeros  = 0,
    C      = 1e6
)

############################################################
# Run JAGS
############################################################

model = jags.model("mrp_unbinned.jags",
                   data=jags_data,
                   n.chains=4,
                   n.adapt=5000)

update(model, 10000)

samples = coda.samples(model,
                       c("mstar","log_phi","alpha","beta"),
                       n.iter=40000,
                       thin=20)

posterior = as.matrix(samples)

############################################################
# Plot - CORRECTED VERSION
############################################################

logbin = 0.2
hist_data = hist(x_obs,
                 breaks=seq(floor(xmin), ceiling(xmax), by=logbin),
                 plot=FALSE)

# KEY FIX: The empirical phi should be counts divided by (Volume * bin_width)
# This gives number density per dex (per log-mass interval)
phi_emp = hist_data$counts / (Vsurvey * logbin)
x_emp   = hist_data$mids

# The MRP function returns phi in units of [Mpc^-3 dex^-1]
# which is exactly what we want to match the empirical data
mrp_log10 <- function(mstar, log_phi, alpha, beta, x){

    phi =
        beta * log(10) *
        exp(log_phi) *
        10^((alpha + 1) * (x - mstar)) *
        exp(- 10^(beta * (x - mstar)))

    log10(phi)
}


xfit = seq(xmin, xmax, length.out=600)

ylims = range(log10(phi_emp), finite=TRUE)

CairoPDF("MRP_unbinned_fit_corrected.pdf", 7, 6)

plot(x_emp, log10(phi_emp),
     pch=16,
     ylim=c(-8, -2),  # Set y-axis limits similar to Driver et al. figure
     xlim=c(12, 16),  # Set x-axis limits similar to Driver et al. figure
     xlab=expression(log[10](M)),
     ylab=expression(log[10](phi)~"[Mpc"^-3~"dex"^-1~"]"),
     main="Unbinned MRP Fit (Corrected)",
     col="green", cex=0.8)

idx = sample(1:nrow(posterior), 800)

for(i in idx){
    lines(xfit,
          mrp_log10(posterior[i,"mstar"],
                    posterior[i,"log_phi"],
                    posterior[i,"alpha"],
                    posterior[i,"beta"],
                    xfit),
          col=rgb(0,0,1,0.02))
}

med = apply(posterior,2,median)

lines(xfit,
      mrp_log10(med["mstar"],
                med["log_phi"],
                med["alpha"],
                med["beta"],
                xfit),
      col="red", lwd=3)

# Add legend
legend("bottomleft",
       legend=c("Binned data", "Best fit", "MCMC samples"),
       col=c("green", "red", rgb(0,0,1,0.5)),
       pch=c(16, NA, NA),
       lty=c(NA, 1, 1),
       lwd=c(NA, 3, 1))

dev.off()

saveRDS(posterior,"mrp_unbinned_chain_corrected.rds")

# Print summary
cat("\n\nPosterior Summary:\n")
cat("==================\n")
print(summary(samples))
