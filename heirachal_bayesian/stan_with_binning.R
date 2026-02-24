library(rstan)
library(Rfits)
library(data.table)
library(celestial)
library(plotrix)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##################################################
# Load GAMA data (same as original script)
##################################################

logbin = 0.2
multi = 5
zlimit = 0.25
zmin = 0.015

ho=67.37
omegam=0.3147
omegal=1-omegam

g3cx=Rfits_read_table("../data/G3CFoFGroupv10.fits")

g3c=g3cx[
  Nfof > multi-1 &
  Zfof < zlimit &
  Zfof > zmin &
  MassAfunc > 1E1 &
  IterCenDec > -3.5
]

parsec=3.0857E16
G=6.67408E-11
msol=1.988E30
magica=13.9

g3c$mymass =
  magica*(g3c$VelDisp*1000)^2*
  g3c$Rad50*parsec*1E6/(G*msol)*(100/ho)

masscorr=c(
0.0,0.0,-0.2672595,-0.1513503,-0.1259069,-0.09006064,
-0.05466009,-0.06666895,-0.01988694,-0.02439581,
-0.0206706,-0.01812964,-0.01556899,-0.01313664,
-0.01743112,-0.007965513,-0.01257178,-0.007064037,
-0.003963656,-0.01271533,-0.002664687,-0.001691287)

g3c$masscorr = masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] = 0

g3c$MassAfunc = g3c$mymass / 10^g3c$masscorr

##################################################
# Vmax calculation
##################################################

gig=fread("../data/GAMAGalsInGroups.csv")

for (i in 1:nrow(g3c)){

  if(g3c$Nfof[i]==2)
    g3c$zmax[i]=sort(
      gig[GroupID==g3c$GroupID[i],zmax_19p8],
      decreasing=TRUE)[2]
  else
    g3c$zmax[i]=sort(
      gig[GroupID==g3c$GroupID[i],zmax_19p8],
      decreasing=TRUE)[multi]
}

g3c$zmax=pmin(g3c$zmax,zlimit)

g3c$vmax =
  179.92/(360^2/pi)*1E9*
  cosdist(g3c$zmax,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoVol

vlimit =
  179.92/(360^2/pi)*1E9*
  cosdist(zlimit,OmegaM=omegam,OmegaL=omegal,H0=ho)$CoVol

##################################################
# Compute HMF
##################################################

massx=seq(10.3,16.1,logbin)

gamahmf2 =
  weighted.hist(
    log10(g3c$MassAfunc),
    w=1/g3c$vmax,
    breaks=massx,
    plot=FALSE
  )

gamax = gamahmf2$mids
gamay = gamahmf2$counts / logbin

valid =
  gamay>0 &
  !is.na(gamay) &
  gamax>12.7

logM_obs = gamax[valid]
log_phi_obs = log10(gamay[valid])

frac_err = rep(0.15, length(logM_obs))

##################################################
# Stan data
##################################################

stan_data = list(

  N = length(logM_obs),

  logM_obs = logM_obs,
  log_phi_obs = log_phi_obs,
  frac_err = frac_err,

  logbin = logbin,
  vlimit = vlimit,

  logM_max = max(logM_obs)
)

##################################################
# Fit model
##################################################

fit = stan(

  file="mrp_binned.stan",

  data=stan_data,

  chains=4,
  iter=4000,
  warmup=2000,
  control=list(adapt_delta=0.95)
)

print(fit)
##################################################
# Extract posterior from Stan
##################################################

post = extract(fit)

posterior_matrix = cbind(
  mstar   = post$logMstar,
  log_phi = post$log_phi,
  alpha   = post$alpha,
  beta    = post$beta
)

##################################################
# Compute posterior summaries
##################################################

med = apply(posterior_matrix, 2, median)
q16 = apply(posterior_matrix, 2, quantile, 0.16)
q84 = apply(posterior_matrix, 2, quantile, 0.84)

##################################################
# Define plotting grid
##################################################

xfit = seq(11, 16, 0.01)

mrp_log10 = function(mstar, log_phi, alpha, beta, logM){

  phi = 10^log_phi

  log10(
    phi *
    beta * log(10) *
    10^((alpha+1)*(logM-mstar)) *
    exp(-10^(beta*(logM-mstar)))
  )
}

##################################################
# Save plot
##################################################

Cairo::CairoPDF("MRP_STAGE3_VMAX.pdf", 10, 7)

plot(gamax, log10(gamay),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(11,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 3: Vmax Weighting (Incompleteness Corrected)")

grid(col="gray80")

##################################################
# Posterior chains
##################################################

set.seed(1)
idx <- sample(1:nrow(posterior_matrix), 300)

for(i in idx){

  y <- mrp_log10(
        posterior_matrix[i,"mstar"],
        posterior_matrix[i,"log_phi"],
        posterior_matrix[i,"alpha"],
        posterior_matrix[i,"beta"],
        xfit
      )

  lines(xfit, y, col=rgb(0,0,1,0.03))
}

##################################################
# Median fit
##################################################

lines(xfit,
      mrp_log10(med["mstar"], med["log_phi"],
                med["alpha"], med["beta"], xfit),
      col="red", lwd=3)

##################################################
# Driver+2022 reference curve
##################################################

lines(xfit,
      mrp_log10(13.51, -3.19, -1.27, 0.47, xfit),
      col="black", lwd=2, lty=2)

##################################################
# Legend
##################################################

legend("bottomleft",
       legend=c("GAMA (Vmax weighted)",
                "Stan fit",
                "Posterior samples",
                "Driver+22"),
       col=c("darkgreen","red",rgb(0,0,1,0.3),"black"),
       pch=c(19,NA,NA,NA),
       lty=c(NA,1,1,2),
       lwd=c(NA,3,1,2),
       bg="white",
       cex=0.9)

##################################################
# Parameter text
##################################################

text(14.3,-2.3,
     sprintf(
       "M* = %.2f (+%.2f/-%.2f)\nlog\u03c6* = %.2f (+%.2f/-%.2f)\n\u03b1 = %.2f (+%.2f/-%.2f)\n\u03b2 = %.2f (+%.2f/-%.2f)",
       med["mstar"],  q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"],
       med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"],
       med["alpha"],  q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"],
       med["beta"],   q84["beta"]-med["beta"],     med["beta"]-q16["beta"]
     ),
     col="red", cex=0.85, pos=4)

dev.off()

cat("Saved to MRP_STAGE3_VMAX.pdf\n")
