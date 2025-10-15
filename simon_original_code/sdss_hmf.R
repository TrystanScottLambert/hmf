#
# Plot segim for individual object
#
library(celestial)
library(devtools)
library(Cairo)
library(ProFound)
library(sm)
library(Rfits)
library("magicaxis")
library("data.table")
library("plotrix")
library(Rfits)
library(foreign)
library(MASS)
#
# Define functions
#
mymagtri <- function(chains, samples = 1000, thin = 1, samptype = "end",
                     grid = FALSE, do.tick = FALSE, refvals = NULL, bestvals = NULL, lab = NULL,
                     ...) {
  chains <- as.data.frame(chains)
  chaincolnames <- colnames(chains)
  Nsamp <- dim(chains)[1]
  Npar <- dim(chains)[2]
  if (!is.null(refvals)) {
    if (length(refvals) != Npar) {
      stop("Length of refvales must be equal to number of parameters!")
    }
  }
  if (Npar <= 1) {
    stop("Need 2+ parameters!")
  }
  if (thin > 1) {
    chains <- chains[seq(1, Nsamp, by = thin), , drop = FALSE]
    Nsamp <- dim(chains)[1]
  }
  if (samples > Nsamp) {
    samples <- Nsamp
  }
  layout(matrix(1:Npar^2, Npar, Npar)[Npar:1, ])
  meanvec <- {}
  sdvec <- {}
  if (samptype == "end") {
    usesamps <- (Nsamp - samples + 1):Nsamp
  }
  if (samptype == "ran") {
    usesamps <- sample(Nsamp, samples)
  }
  if (samptype == "thin") {
    usesamps <- seq(1, Nsamp, length = samples)
  }
  for (i in 1:Npar) {
    meanvec <- c(meanvec, mean(chains[usesamps, i]))
    sdvec <- c(sdvec, sd(chains[usesamps, i]))
  }
  par(oma = c(4.1, 4.1, 1.1, 1.1))
  for (i in 1:Npar) {
    for (j in 1:Npar) {
      par(mar = c(0, 0, 0, 0))
      xrange <- range(chains[usesamps, i])
      yrange <- range(chains[usesamps, j])
      if (xrange[1] == xrange[2]) {
        val <- xrange[1]
        xrange[1] <- val - 0.05
        xrange[2] <- val + 0.05
      }
      if (yrange[1] == yrange[2]) {
        val <- yrange[1]
        yrange[1] <- val - 0.05
        yrange[2] <- val + 0.05
      }
      if (i == 1) {
        xrange <- c(12.25, 14.75)
      } else if (i == 2) {
        xrange <- c(-5.5, -2.5)
      } else if (i == 3) {
        xrange <- c(-2, 0)
      } else if (i == 4) {
        xrange <- c(0, 1)
      }
      if (j == 1) {
        yrange <- c(12.25, 14.75)
      } else if (j == 2) {
        yrange <- c(-5.5, -2.5)
      } else if (j == 3) {
        yrange <- c(-2, 0)
      } else if (j == 4) {
        yrange <- c(0, 1)
      }
      if (i == j) {
        xtemp <- chains[usesamps, i]
        if (sd(xtemp) == 0) {
          xtemp <- xtemp + rnorm(samples, sd = 0.001)
        }
        plot(density(xtemp),
          axes = FALSE, main = "",
          xlim = xrange, col = "darkgrey"
        )
        magaxis(1,
          grid = grid, grid.col = "darkgrey",
          labels = FALSE, do.tick = do.tick
        )
        abline(v = bestvals[i], lty = 1, col = "blue")
        abline(v = meanvec[i], lty = 1, col = "red")
        abline(v = meanvec[i] - sdvec[i], lty = 3, col = "red")
        abline(v = meanvec[i] + sdvec[i], lty = 3, col = "red")
        if (!is.null(refvals)) {
          abline(v = refvals[i], lty = 2, col = "black")
        }
        box()
        if (i == 1) {
          plot.window(xlim = xrange, ylim = yrange)
          if (is.null(lab)) {
            magaxis(1, xlab = chaincolnames[i])
          } else {
            magaxis(1, xlab = lab[[i]])
          }
          if (is.null(lab)) {
            magaxis(2, ylab = chaincolnames[j])
          } else {
            magaxis(2, ylab = lab[[j]])
          }
        }
      } else {
        if (i > j) {
          plot.new()
          plot.window(xlim = xrange, ylim = yrange)
          xtemp <- chains[usesamps, i]
          ytemp <- chains[usesamps, j]
          if (sd(xtemp) == 0) {
            xtemp <- xtemp + rnorm(samples, sd = 0.001)
          }
          if (sd(ytemp) == 0) {
            ytemp <- ytemp + rnorm(samples, sd = 0.001)
          }
          magaxis(1:2,
            grid = grid, grid.col = "darkgrey",
            labels = FALSE, do.tick = do.tick
          )
          magcon(xtemp, ytemp,
            dobar = FALSE, doim = TRUE,
            add = TRUE, lty = c(1, 1, 1), xlim = xrange, imcol = c(NA, terrain.colors(1000, alpha = 1, rev = TRUE)),
            ylim = yrange, h = c(diff(xrange), diff(yrange)) / 50, col = "darkgrey",
            ...
          )
          points(meanvec[i], meanvec[j],
            col = "red",
            pch = 4, cex = 1
          )
          points(refvals[i], refvals[j],
            col = "black",
            pch = 1, cex = 1
          )
          points(bestvals[i], bestvals[j],
            col = "blue",
            pch = 4, cex = 1
          )
          box()
          abline(v = bestvals[i], lty = 1, col = "blue")
          abline(v = meanvec[i], lty = 1, col = "red")
          abline(
            v = meanvec[i] - sdvec[i], lty = 3,
            col = "red"
          )
          abline(
            v = meanvec[i] + sdvec[i], lty = 3,
            col = "red"
          )
          if (!is.null(refvals)) {
            abline(v = refvals[i], lty = 2, col = "black")
          }
          if (j == 1) {
            if (is.null(lab)) {
              magaxis(1, xlab = chaincolnames[i])
            } else {
              magaxis(1, xlab = lab[[i]])
            }
          }
        } else {
          plot.new()
          plot.window(xlim = xrange, ylim = yrange)
          magaxis(1:2,
            grid = grid, grid.col = "darkgrey",
            labels = FALSE, do.tick = do.tick
          )
          points(chains[usesamps, c(i, j)],
            pch = ".",
            col = rgb(169 / 255, 169 / 255, 169 / 255, 0.1)
          )
          points(meanvec[i], meanvec[j],
            col = "red",
            pch = 4, cex = 1
          )
          points(refvals[i], refvals[j],
            col = "black",
            pch = 1, cex = 1
          )
          points(bestvals[i], bestvals[j],
            col = "blue",
            pch = 4, cex = 1
          )
          box()
          if (i == 1) {
            if (is.null(lab)) {
              magaxis(2, ylab = chaincolnames[j])
            } else {
              magaxis(2, ylab = lab[[j]])
            }
          }
        }
      }
    }
  }
  output <- cbind(mean = meanvec, sd = sdvec)
  rownames(output) <- chaincolnames
  return(invisible(output))
}
#
massfn <- function(x) {
  mstar <- x[1]
  phi <- x[2]
  alpha <- x[3]
  beta <- x[4]
  allxxx <- c(max(allx) + seq(1, 10, 1) * logbin)
  penalty <- 2 * volumesdss * sum(log(10) * beta * exp(-10^(beta * (allxxx - mstar))) * (phi * (10^allxxx / 10^mstar)^(alpha + 1))) * logbin
  model <- log10(beta * log(10) * (exp(-10^(beta * (allx - mstar)))) * (phi * ((10^allx / 10^mstar)^(alpha + 1))))
  sum(((ally - model) / (allf / log(10)))^2) + penalty
}
#
# massfn <- function(x){
#  mstar <- x[1]
#  phi <- x[2]
#  alpha <- x[3]
#  sum(((ally-log10(log(10)*(exp(-10^allx/10^mstar))*(phi*((10^allx/10^mstar)^(alpha+1)))))^2)/allw^2)}
#
#
# massfn3 <- function(x){
#  mstar <- x[1]
#  phi <- x[2]
#  alpha <- x[3]
#  beta <- x[4]
#  allxxx=c(max(allx)+logbin,max(allx)+2*logbin,max(allx)+3*logbin)
#  penalty=2*(volumesdss*sum(log(10)*beta*(exp(-10^(beta*(allxxx-mstar))))*(phi*((10^allxxx/10^mstar)^(alpha+1)))))
#  penalty=0
#  sum(((ally-log10(log(10)*beta*(exp(-10^(beta*(allx-mstar))))*(phi*((10^allx/10^mstar)^(alpha+1)))))^2)/allw^2)+penalty}
#
ho <- 67.37
omegam <- 0.3147
omegamerr <- 0.0074
omegal <- 1 - omegam
G <- 6.67408E-11
msol <- 1.988E30
parsec <- 3.0857E16
cosvar <- function(V, N) {
  ((219.7 - 52.4 * log10(V) + 3.21 * (log10(V))^2) / N^0.5) / 100.0
}
logbin <- 0.1
fitbinwid <- 0.01
rhocrit <- 3 * (1000 * ho / (1E6 * parsec))^2 / (8 * pi * G)
zlimitsdss <- 0.08
#
inputargs <- commandArgs(TRUE)
multi <- as.integer(inputargs[1])
mlimit <- as.numeric(inputargs[2])
zmin <- as.numeric(inputargs[3])
# multi=5
# mlimit=12.9
# zmin=0.015
volumesdss <- 7221 / (360^2 / pi) * 1E9 * cosdist(zlimitsdss, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol - 7221.0 / (360^2 / pi) * 1E9 * cosdist(zmin, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol
volumesdssmin <- volumesdss / 1000.0
#
#
betamrp <- 0.7097976
A <- 1.727006E-19
mstarmrp <- 14.42947
alphamrp <- -1.864908
mrpx <- seq(0, 17, 0.001) + log10(100 / ho)
mrpy <- A * betamrp * 10^((alphamrp + 1) * (mrpx - mstarmrp)) * exp(-10^(betamrp * (mrpx - mstarmrp))) * (ho / 100)^3
factor <- sum(10^mrpx * mrpy) * 0.001 * msol / (1E6 * parsec)^3 / (omegam * rhocrit)
phimrp <- A / factor
#
# Prime sdss with zmax and vmax values using multi+1 member
#
sdssgalsx <- Rfits_read_table("/Users/sdriver/Drpbx/active/hmf/sdssdr10table1.fits", ext = 2)
names(sdssgalsx)[6] <- "rank"
names(sdssgalsx)[9] <- "redshift"
names(sdssgalsx)[31] <- "absmag_r"
names(sdssgalsx)[4] <- "idcl"
sdssgals <- sdssgalsx[rank == multi & redshift < zlimitsdss + 0.05, c("idcl", "rank", "redshift", "absmag_r")]
sdssgals$dmax <- 10^(0.2 * (17.77 - sdssgals$absmag_r - 25 - 1.5 * log10(1 + sdssgals$redshift)))
#
#
sdssx <- Rfits_read_table("/Users/sdriver/Drpbx/active/hmf/sdssdr10table2.fits", ext = 2)
names(sdssx)[10] <- "zcl"
names(sdssx)[2] <- "nrich"
names(sdssx)[15] <- "mass"
names(sdssx)[1] <- "idcl"
sdss <- sdssx[zcl < zlimitsdss & zcl > zmin & nrich > multi - 1, c("idcl", "zcl", "nrich", "mass")]
sdss$mass <- sdss$mass * 1E12 * (67.8 / ho)
xx <- seq(3, 22)
yy <- c(0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 0.24018684, 0.20226682, 0.18645475, 0.17437005, 0.14271506, 0.13922450, 0.13482418, 0.13741619, 0.11715141, 0.12134983, 0.10078830, 0.09944761, 0.09913166, 0.08590223, 0.07588408)
sdss$log10MassErr <- approx(xx, yy, sdss$nrich)$y
sdss$log10MassErr[is.na(sdss$log10MassErr)] <- 0.03
# sdss$log10MassErr=(sdss$log10MassErr^2+0.1^2)^0.5
sdss$log10MassErr[sdss$log10MassErr < 0.1] <- 0.1
# sdss$log10MassErr[sdss$mass>10^14.45]=0.15
# sdss$log10MassErr[sdss$mass<10^14.45 & sdss$mass>10^14.25]=0.3
# sdss$log10MassErr[sdss$mass<10^14.25]=0.45
# sdss$log10MassErr=0.01
#
# Build dmax grid and then allocate to galaxies and match to groups
#
try <- 0.0
tryz <- 0.0
for (k in 1:1000) {
  try[k] <- (cosdist(z = as.numeric(k / 1000), OmegaM = omegam, OmegaL = omegal, H0 = ho))$LumDist
  tryz[k] <- k / 1000.0
}
for (i in 1:length(sdssgals$dmax)) {
  sdssgals$zmax[i] <- tryz[which.min(abs(sdssgals$dmax[i] - try))]
}
for (i in 1:length(sdss$idcl)) {
  sdss$zmax[i] <- sdssgals[idcl == sdss$idcl[i], zmax]
}
sdss$zmax <- ifelse(sdss$zmax < sdss$zcl, 1.1 * sdss$zcl, sdss$zmax)
sdss$zmax <- ifelse(sdss$zmax > zlimitsdss, zlimitsdss, sdss$zmax)
sdss$vmax <- 7221 / (360^2 / pi) * 1E9 * cosdist(sdss$zmax, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol
sdss$vmax <- sdss$vmax - 7221.0 / (360^2 / pi) * 1E9 * cosdist(zmin, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol
sdss$weightszlimit <- ifelse(sdss$vmax > volumesdss, volumesdss, sdss$vmax)
sdss$weightszlimit <- ifelse(sdss$vmax < volumesdssmin, volumesdssmin, sdss$vmax)
sdss$weightszlimit[sdss$idcl == 81455] <- volumesdss
sdss$cosvar <- cosvar(sdss$vmax, 1)
#
# SDSS HMF and cosmic variance per galaxy
#
#
massx <- seq(10.3, 16.1, logbin)
sdsshmf <- maghist(log10(sdss$mass), breaks = massx, plot = FALSE, verbose = FALSE)
sdsshmf2 <- weighted.hist(log10(sdss$mass), w = 1 / sdss$weightszlimit, breaks = massx, plot = FALSE)
# cosvariance=((weighted.hist(log10(sdss$mass),w=sdss$cosvar,breaks=massx,plot=FALSE)$counts))/sdsshmf$counts
# cosvariance[is.na(cosvariance)]=0
cosvariance <- cosvar(volumesdss, 1)
#
png(filename = paste0("/Users/sdriver/Drpbx/active/hmf/sdsshmf", multi, ".png"), width = 18.0, height = 12.0, units = "cm", res = 240)
#
par(fig = c(0, 1, 0.8, 1), mar = c(0, 3.35, 0.25, 0.25), oma = c(0, 0, 0, 0))
magplot(sdsshmf$mids - 0.5 * logbin, sdsshmf$counts, xlim = c(12, 16), grid = FALSE, type = "s", ylab = "N", majorn = c(5, 2), labels = c(0, 1))
text(16.2, 0.9 * max(sdsshmf$counts), paste0("SDSS groups at z < ", zlimitsdss), pos = 2, cex = 1.0)
text(16.2, 0.6 * max(sdsshmf$counts), paste0("Multiplicity > ", multi - 1), pos = 2, cex = 1.0)
text(16.2, 0.3 * max(sdsshmf$counts), paste0("N groups = ", sum(sdsshmf$counts)), pos = 2, cex = 1.0)
#
par(fig = c(0, 1, 0, 0.795), mar = c(3.0, 3.35, 0.25, 0.25), oma = c(0, 0, 0, 0), new = TRUE)
magplot(0, 0, xlim = c(12, 16), ylim = c(-8, -2), grid = FALSE, xlab = expression("Halo Mass (M"["\u0298"] ~ ")"), ylab = expression("log"[10] ~ "(umber density) [Mpc"^-3 ~ "dex"^-1 ~ "]"))
# lines(mrpx,log10(mrpy)-log10(factor),lty=2,lwd=2)
# lines(mrpx-0.125,log10(mrpy)-log10(factor)+0.11,lty=2,lwd=2)
lines(mrpx - 0.08, log10(mrpy) - log10(factor) + 0.08, lty = 2, lwd = 2)
#
meancounts <- 0.0
mediancounts <- 0.0
upcounts <- 0.0
docounts <- 0.0
mcerr <- 0.0
mockcounts <- matrix(0, nrow = 1001, ncol = length(sdsshmf$mids))
for (i in 1:1001) {
  mockmass <- log10(sdss$mass) + rnorm(length(sdss$mass), 0.0, sdss$log10MassErr)
  mockcounts[i, 1:length(sdsshmf$mids)] <- weighted.hist(mockmass, w = 1 / sdss$weightszlimit, breaks = massx, plot = FALSE)$counts
  #  points(sdsshmf$mids+rnorm(1,0,0.01),log10(mockcounts[i,1:length(sdsshmf$mids)]/logbin),pch=".",col=rgb(0.5,0.5,0.5,0.1))
}
# meancounts
for (i in 1:length(sdsshmf$mids)) {
  meancounts[i] <- mean(mockcounts[, i])
}
for (i in 1:length(sdsshmf$mids)) {
  mediancounts[i] <- quantile(mockcounts[, i], 0.5)
}
for (i in 1:length(sdsshmf$mids)) {
  upcounts[i] <- quantile(mockcounts[, i], 0.84)
}
for (i in 1:length(sdsshmf$mids)) {
  docounts[i] <- quantile(mockcounts[, i], 0.16)
}
# points(sdsshmf$mids,log10(meancounts/logbin),pch=4,col=rgb(0,0,1,1),cex=0.5)
# points(sdsshmf$mids,log10(mediancounts/logbin),pch=16,col=rgb(0,0,1,1),cex=0.5)
# segments(sdsshmf$mids,log10(upcounts/logbin),sdsshmf$mids,log10(docounts/logbin),col=rgb(0,0,1,1))
# edb
edb <- meancounts / sdsshmf2$counts
edb[is.infinite(edb)] <- 1.0
edb[is.na(edb)] <- 1.0
#
# MC err
#
for (i in 1:length(sdsshmf$mids)) {
  mcerr[i] <- quantile((meancounts[i] - mockcounts[, i])^2, 0.66)^0.5 / (meancounts[i])
}
#
# Poisson Error
#
rootnerr <- (sdsshmf$counts^0.5 / sdsshmf$counts)
#
sdssx <- sdsshmf2$mids
sdssy <- sdsshmf2$counts / (logbin * edb)
sdssf <- (mcerr^2 + rootnerr^2)^0.5
sdssf[is.na(sdssf)] <- 0.9999
sdssf[is.infinite(sdssf)] <- 0.0
sdssf <- ifelse(sdssf >= 1.0, 0.9999, sdssf)
#
elmo14 <- fread("/Users/sdriver/Drpbx/active/hmf/elmo.csv")
elmo14$V1 <- elmo14$V1 + 10.0 + log10(100 / ho)
elmo14$V2 <- 1.0 * elmo14$V2 * (ho / 100)^3
elmo14$V3 <- 1.0 * elmo14$V3 * (ho / 100)^3
elmo14$V4 <- 1.0 * elmo14$V4 * (ho / 100)^3
elmox <- c(elmo14$V1, rev(elmo14$V1))
elmoy <- c(log10(elmo14$V2 - elmo14$V3), rev(log10(elmo14$V2 + elmo14$V4)))
elmoy[is.na(elmoy)] <- -6.5
polygon(elmox, elmoy, col = rgb(0.5, 0.5, 0.5, 0.1))
lines(elmo14$V1, log10(elmo14$V2), pch = 6, col = "pink")
#
allx <- c(sdssx[sdssy > 0 & sdssx > mlimit])
ally <- c(log10(sdssy[sdssy > 0 & sdssx > mlimit]))
allf <- c(sdssf[sdssy > 0 & sdssx > mlimit])
sdssfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfn, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
#
xfit <- seq(0.0, 20, fitbinwid)
yfit <- sdssfit$par[4] * log(10) * (exp(-10^(sdssfit$par[4] * (xfit - sdssfit$par[1])))) * (sdssfit$par[2] * ((10^xfit / 10^sdssfit$par[1])^(sdssfit$par[3] + 1)))
sdssomegamatter <- sum((yfit * 10^xfit) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
#
mstar <- 0
phistar <- 0
alphastar <- 0
betastar <- 0
mass <- 0.0
for (i in 1:10001) {
  cv <- sdssy * rnorm(length(sdssf), 0.0, cosvariance)
  mocksdssy <- sdssy + sdssy * rnorm(length(sdssf), 0.0, sdssf) + cv
  allx <- c(sdssx[mocksdssy > 0 & !is.na(sdssy) & sdssx > mlimit])
  ally <- c(log10(mocksdssy[mocksdssy > 0 & !is.na(sdssy) & sdssx > mlimit]))
  allf <- c(sdssf[mocksdssy > 0 & !is.na(sdssy) & sdssx > mlimit])
  montyfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfn, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
  ymonty <- montyfit$par[4] * log(10) * (exp(-10^(montyfit$par[4] * (xfit - montyfit$par[1])))) * (montyfit$par[2] * ((10^xfit / 10^montyfit$par[1])^(montyfit$par[3] + 1)))
  if (i < 1001) {
    lines(xfit, log10(ymonty), col = rgb((100 / 255), (149 / 255), (237 / 255), 0.01))
  }
  #    points(sdssx,log10(mocksdssy),pch=".",col=rgb(0.5,0.5,0.5,0.1))
  mstar[i] <- montyfit$par[1]
  phistar[i] <- montyfit$par[2]
  alphastar[i] <- montyfit$par[3]
  betastar[i] <- montyfit$par[4]
}
#
points(sdssx, log10(sdsshmf2$counts / logbin), pch = 5, cex = 0.75, col = "limegreen")
# points(sdssx,log10(sdsshmf2$counts/logbin),pch="|",cex=0.75,col="limegreen")
points(sdssx[sdssx > mlimit], log10(sdssy[sdssx > mlimit]), pch = 16, col = "purple")
lines(xfit, log10(yfit), col = rgb((100 / 255), (149 / 255), (237 / 255), 1.0), lwd = 2)
points(sdssx[sdssx < mlimit], log10(sdssy[sdssx < mlimit]), pch = 1, col = "purple", cex = 1)
magerr(sdssx, log10(sdssy), ylo = log10(1 - sdssf), yhi = log10(1 + sdssf), col = "purple")
#
lines(c(14.41, 14.59), c(-4.2, -4.2), col = rgb(0.5, 0.5, 0.5, 0.25), lwd = 5)
lines(c(14.4, 14.6), c(-4.2, -4.2), col = "pink")
text(14.5, -4.2, "  SDSS DR10 (Tempel+ 2014)", col = rgb(0.5, 0.5, 0.5), pos = 4)
lines(c(14.4, 14.6), c(-3.7, -3.7), col = "black", lty = 2, lwd = 2)
text(14.5, -3.7, "  LCDM prediction", col = "black", pos = 4)
lines(c(14.4, 14.6), c(-3.2, -3.2), col = rgb((100 / 255), (149 / 255), (237 / 255)), lwd = 2)
text(14.5, -3.2, "  Best MRP function fit", col = rgb((100 / 255), (149 / 255), (237 / 255)), pos = 4)
points(14.5, -2.7, pch = 5, col = "limegreen", cex = 0.75)
# points(14.5,-2.7,pch="|",col="limegreen",cex=0.75)
text(14.5, -2.7, paste0("  SDSS z<", zlimitsdss, " Raw HMF"), col = "limegreen", pos = 4, cex = 1.0)
points(14.5, -2.2, pch = 16, col = "purple", cex = 1)
text(14.5, -2.2, paste0("  SDSS z<", zlimitsdss, " Corrected HMF"), col = "purple", pos = 4, cex = 1.0)
#
# par(fig=c(0.75,0.9825,0.525,0.775),mar=c(0,0,0.25,0.25),oma=c(0,0,0,0),new=TRUE)
# x=seq(0,50,0.025)
# maghist(mass,breaks=x,xlab=expression(Omega[M]),grid=FALSE,ylab="Freq.",majorn=c(3,2),xlim=c(0,1),col=rgb((100/255),(149/255),(237/255),0.5),verbose=FALSE)
# abline(v=omegam,col="black",lwd=2,lty=2)
tmp <- dev.off()
#
# write
#
sdsstable <- paste(rev(sdssx), rev(round(sdsshmf$counts, 3)), rev(round(log10(sdsshmf2$counts), 3)), rev(round(log10(sdssy), 3)), rev(round(rootnerr, 3)), rev(round(mcerr, 3)), rev(round(cosvariance, 3)), rev(round(sdssf, 3)), sep = " $&$ ")
write(sdsstable, file = paste0("/Users/sdriver/Drpbx/active/hmf/sdsstable", multi, ".txt"), append = FALSE)
sdsstable2 <- as.data.frame(cbind(rev(sdssx), rev(sdsshmf$counts), rev(log10(sdsshmf2$counts)), rev(log10(sdssy)), rev(rootnerr), rev(mcerr), rev(cosvariance), rev(sdssf)))
write.csv(sdsstable2, paste0("/Users/sdriver/Drpbx/active/hmf/sdsshmf", multi, ".csv"), row.names = FALSE)
#
png(filename = paste0("/Users/sdriver/Drpbx/active/hmf/sdsscov", multi, ".png"), width = 16.0, height = 16.0, units = "cm", res = 240)
par(mar = c(3, 3.5, 0.25, 0.25), oma = c(0, 0, 0, 0))
params <- cbind(mstar[!is.na(mstar)], log10(phistar[!is.na(mstar)]), alphastar[!is.na(mstar)], betastar[!is.na(mstar)])
fitvals <- mymagtri(params,
  samples = length(params[, 1]), samptype = "ran",
  lab = c("log\u2081\u2080(M*) [M\u2092]", "log\u2081\u2080(\u03c6*) [Mpc\u207B\u00b3]", "\u03b1", "\u03b2"),
  refvals = c(mstarmrp - 0.075, log10(phimrp) + 0.075, alphamrp, betamrp),
  bestvals = c(sdssfit$par[1], log10(sdssfit$par[2]), sdssfit$par[3], sdssfit$par[4])
)
cat(
  paste0("SDSS", multi),
  sdssfit$par[1], quantile(mstar, 0.16, na.rm = TRUE), quantile(mstar, 0.84, na.rm = TRUE),
  log10(sdssfit$par[2]), log10(quantile(phistar, 0.16, na.rm = TRUE)), log10(quantile(phistar, 0.84, na.rm = TRUE)),
  sdssfit$par[3], quantile(alphastar, 0.16, na.rm = TRUE), quantile(alphastar, 0.84, na.rm = TRUE),
  sdssfit$par[4], quantile(betastar, 0.16, na.rm = TRUE), quantile(betastar, 0.84, na.rm = TRUE), "\n"
)
tmp <- dev.off()
