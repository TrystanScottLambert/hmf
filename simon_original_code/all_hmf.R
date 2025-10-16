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
require(foreign)
require(MASS)
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
  penalty <- 2 * (volume * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumesdss * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumeREFLEXII * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  model <- log10(beta * log(10) * (exp(-10^(beta * (allx - mstar)))) * (phi * ((10^allx / 10^mstar)^(alpha + 1))))
  sum(((ally - model) / (allf / log(10)))^2) + penalty
}
#
massfnfix <- function(x) {
  mstar <- x[1]
  phi <- x[2]
  alpha <- alpha
  beta <- x[3]
  allxxx <- c(max(allx) + seq(1, 10, 1) * logbin)
  penalty <- 2 * (volume * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumesdss * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumeREFLEXII * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  model <- log10(beta * log(10) * (exp(-10^(beta * (allx - mstar)))) * (phi * ((10^allx / 10^mstar)^(alpha + 1))))
  sum(((ally - model) / (allf / log(10)))^2) + penalty
}
#
massfnomegam <- function(x) {
  mstar <- x[1]
  phi <- x[2]
  alpha <- x[3]
  beta <- x[4]
  allxxx <- c(max(allx) + seq(1, 10, 1) * logbin)
  penalty <- 2 * (volume * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumesdss * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  penalty <- penalty + 2 * (volumeREFLEXII * sum(log(10) * beta * (exp(-10^(beta * (allxxx - mstar)))) * (phi * ((10^allxxx / 10^mstar)^(alpha + 1))))) * logbin
  model <- log10(beta * log(10) * (exp(-10^(beta * (allx - mstar)))) * (phi * ((10^allx / 10^mstar)^(alpha + 1))))
  x <- rev(seq(0.0, 18, fitbinwid))
  y <- log(10) * phi * beta * ((10^x / 10^mstar)^(alpha + 1)) * (exp(-10^(beta * (x - mstar))))
  omegamatter <- ((sum((y * 10^x) * fitbinwid * msol / (1E6 * parsec)^3 / rhocrit) - omegam)^2) / omegamerr^2
  omegamatter + sum(((ally - model) / (allf / log(10)))^2) + penalty
}
#
# Define inputs and constants
#
# inputargs=commandArgs(TRUE)
# myoption=as.character(inputargs[1])
myoption <- "Omega"
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
fitbinwid <- 0.001
alpha <- -1.864908
logbin <- 0.2
RAmid <- c(135, 180, 217, 5)
name <- c("G09", "G12", "G15")
zlimit <- 0.25
zlimitsdss <- 0.08
volume <- 175 / (360^2 / pi) * 1E9 * cosdist(zlimit, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol
volumesdss <- 7221 / (360^2 / pi) * 1E9 * cosdist(zlimitsdss, OmegaM = omegam, OmegaL = omegal, H0 = ho)$CoVol
volumeREFLEXII <- 13000000.0
rhocrit <- 3 * (1000 * ho / (1E6 * parsec))^2 / (8 * pi * G)
iters <- 10001
if (myoption == "Omega") {
  iters <- 1001
}
#
# Read in data
#
reflex <- fread("../data/reflex.csv")
reflex2 <- fread("../data/reflex2.csv")
tpigg <- fread("../data/tpigg.dat", sep = " ", data.table = FALSE)
gama <- fread("../gamahmfGAMA5.csv", sep = ",", data.table = FALSE)
gama <- gama[gama$V1 > 12.7 & !is.infinite(gama$V4), ]
sdss <- fread("./sdsshmf5.csv", sep = ",", data.table = FALSE)
sdss <- sdss[sdss$V1 > 12.9 & !is.infinite(sdss$V4), ]
#
# New HMF
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
# Plot combined HMF
#
png(filename = paste0("./allhmf", myoption, ".png"), width = 20.0, height = 12.0, units = "cm", res = 240)
par(fig = c(0, 1, 0, 1), mar = c(3, 3.5, 0.25, 0.25), oma = c(0, 0, 0, 0))
magplot(0, 0, xlim = c(12.75, 16), ylim = c(-8, -2), unlog = "xy", xlab = expression("Halo Mass (M"["\u0298"] ~ ")"), ylab = expression("log"[10] ~ "(number density) [Mpc"^-3 ~ "dex"^-1 ~ "]"))
# lines(mrpx,log10(mrpy)-log10(factor),lty=2,lwd=2)
# lines(mrpx-0.125,log10(mrpy)-log10(factor)+0.11,lty=2,lwd=2)
lines(mrpx - 0.08, log10(mrpy) - log10(factor) + 0.08, lty = 2, lwd = 2)
#
# REFLEXII
#
reflex2$mass <- log10(1E14 * reflex2$mass) + log10(70 / ho)
reflex2$density <- log10(reflex2$density) + 3.0 * log10(ho / 70) - 14.0 + reflex2$mass + 1
reflex2$min <- log10(reflex2$min) + 3.0 * log10(ho / 70) - 14.0 + reflex2$mass + 1
reflex2$max <- log10(reflex2$max) + 3.0 * log10(ho / 70) - 14.0 + reflex2$mass + 1
reflex$x <- reflex$x + log10(70 / ho)
reflex$Curve1 <- reflex$Curve1 + 4.0 * log10(ho / 70) - 14.0 + reflex$x + 1
# lines(reflex2$mass,reflex2$min,col="forestgreen",lty=3)
# lines(reflex2$mass,reflex2$max,col="forestgreen",lty=3)
reflexx <- reflex$x
reflexy <- reflex$Curve1
reflexf <- 1 / 20^0.5 + reflexx * 0.0
reflexf[43] <- 1 / 3^0.5
reflexf[1] <- 1 / 3^0.5
#
# 2PIGG
#
tpigg$mass <- tpigg$V1 + log10(100 / ho)
tpigg$density <- tpigg$V2 + log10((ho / 100)^3)
tpigg$up <- tpigg$V3
tpigg$do <- tpigg$V4
#
# SDSS DR10
#
elmo14 <- fread("./elmo.csv")
elmo14$V1 <- elmo14$V1 + 10.0 + log10(100 / ho)
elmo14$V2 <- elmo14$V2 * (ho / 100)^3
elmo14$V3 <- elmo14$V3 * (ho / 100)^3
elmo14$V4 <- elmo14$V4 * (ho / 100)^3
elmox <- c(elmo14$V1, rev(elmo14$V1))
elmoy <- c(log10(elmo14$V2 - elmo14$V3), rev(log10(elmo14$V2 + elmo14$V4)))
elmoy[is.na(elmoy)] <- -6.5
polygon(elmox, elmoy, col = rgb(0.5, 0.5, 0.5, 0.1), border = rgb(0.5, 0.5, 0.5, 0.1))
lines(elmo14$V1, log10(elmo14$V2), pch = 6, col = "cyan")
#
# Generate fits
#
dm <- 0.0
mstar <- 0
phistar <- 0
alphastar <- 0
betastar <- 0
mass <- 0.0
omegamatter <- 0
omegamatter2 <- 0
xfit <- seq(0, 18, fitbinwid)
#
cosvarreflex <- 0.05
cosvarsdss <- cosvar(volumesdss, 1)
cosvargama <- cosvar(volume / 3, 3)
#
for (i in 1:iters) {
  #
  cvreflex <- reflexy * rnorm(length(gama$V4), 0.0, cosvarreflex)
  cvsdss <- sdss$V4 * rnorm(length(sdss$V4), 0.0, cosvarsdss)
  cvgama <- gama$V4 * rnorm(length(gama$V4), 0.0, cosvargama)
  #
  if (myoption == "GSR" | myoption == "FIX" | myoption == "Omega") {
    allxx <- c(gama$V1, sdss$V1, reflexx)
    allyy <- c(gama$V4 + cvgama, sdss$V4 + cvsdss, reflexy + cvreflex)
    allyy <- 10^allyy
    allff <- c(gama$V8, sdss$V8, reflexf)
  } else if (myoption == "GS") {
    allxx <- c(gama$V1, sdss$V1)
    allyy <- c(gama$V4 + cvgama, sdss$V4 + cvsdss)
    allyy <- 10^allyy
    allff <- c(gama$V8, sdss$V8)
  } else if (myoption == "GR") {
    allxx <- c(gama$V1, reflexx)
    allyy <- c(gama$V4 + cvgama, reflexy + cvreflex)
    allyy <- 10^allyy
    allff <- c(gama$V8, reflexf)
  } else if (myoption == "SR") {
    allxx <- c(sdss$V1, reflexx)
    allyy <- c(sdss$V4 + cvsdss, reflexy + cvreflex)
    allyy <- 10^allyy
    allff <- c(sdss$V8, reflexf)
  } else if (myoption == "R") {
    allxx <- c(reflexx)
    allyy <- c(reflexy + cvreflex)
    allyy <- 10^allyy
    allff <- c(reflexf)
  } else if (myoption == "G") {
    allxx <- c(gama$V1)
    allyy <- c(gama$V4 + cvgama)
    allyy <- 10^allyy
    allff <- c(gama$V8)
  }
  #
  allff[is.na(allff)] <- 0.9999
  allff[is.infinite(allff)] <- 0.0
  allff <- ifelse(allff >= 1.0, 0.9999, allff)
  #
  mockally <- allyy + allyy * rnorm(length(allff), 0.0, allff)
  allx <- c(allxx[mockally > 0 & !is.na(mockally)])
  ally <- c(log10(mockally[mockally > 0 & !is.na(mockally)]))
  allf <- c(allff[mockally > 0 & !is.na(mockally)])
  if (myoption == "FIX") {
    montyfit <- optim(par = c(mstarmrp, phimrp, betamrp), fn = massfnfix, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 0.1)))
    montyfit$par[4] <- montyfit$par[3]
    montyfit$par[3] <- alpha
  } else if (myoption == "Omega") {
    montyfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfnomegam, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
  } else {
    montyfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfn, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
  }
  ymonty <- montyfit$par[4] * log(10) * (exp(-10^(montyfit$par[4] * (xfit - montyfit$par[1])))) * (montyfit$par[2] * ((10^xfit / 10^montyfit$par[1])^(montyfit$par[3] + 1)))
  if (i < 1001) {
    lines(xfit, log10(ymonty), col = rgb((100 / 255), (149 / 255), (237 / 255), 0.01))
  }
  mass[i] <- sum((ymonty * 10^xfit) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
  mstar[i] <- montyfit$par[1]
  phistar[i] <- montyfit$par[2]
  alphastar[i] <- montyfit$par[3]
  betastar[i] <- montyfit$par[4]
  omegamatter[i] <- sum((ymonty * 10^xfit) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
  omegamatter2[i] <- sum((ymonty[xfit > 12.7] * 10^xfit[xfit > 12.7]) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
}
#
if (myoption == "GSR" | myoption == "FIX" | myoption == "Omega") {
  allxx <- c(gama$V1, sdss$V1, reflexx)
  allyy <- c(gama$V4, sdss$V4, reflexy)
  allff <- c(gama$V8, sdss$V8, reflexf)
} else if (myoption == "GS") {
  allxx <- c(gama$V1, sdss$V1)
  allyy <- c(gama$V4, sdss$V4)
  allff <- c(gama$V8, sdss$V8)
} else if (myoption == "GR") {
  allxx <- c(gama$V1, reflexx)
  allyy <- c(gama$V4, reflexy)
  allff <- c(gama$V8, reflexf)
} else if (myoption == "SR") {
  allxx <- c(sdss$V1, reflexx)
  allyy <- c(sdss$V4, reflexy)
  allff <- c(sdss$V8, reflexf)
} else if (myoption == "R") {
  allxx <- c(reflexx)
  allyy <- c(reflexy)
  allff <- c(reflexf)
} else if (myoption == "G") {
  allxx <- c(gama$V1)
  allyy <- c(gama$V4)
  allyy <- 10^allyy
  allff <- c(gama$V8)
}
allx <- allxx[!is.na(allyy)]
allf <- allff[!is.na(allyy)]
allf <- ifelse(allf >= 1.0, 0.9999, allf)
ally <- allyy[!is.na(allyy)]
#
if (myoption == "FIX") {
  allfit <- optim(par = c(mstarmrp, phimrp, betamrp), fn = massfnfix, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 0.1)))
  allfit$par[4] <- allfit$par[3]
  allfit$par[3] <- alpha
} else if (myoption == "Omega") {
  allfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfnomegam, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
} else {
  allfit <- optim(par = c(mstarmrp, phimrp, alphamrp, betamrp), fn = massfn, control = list(maxit = 500, reltol = 1e-8, parscale = c(1, 1, 1, 0.1)))
}
xangus <- seq(0.0, 18, fitbinwid)
yangus <- allfit$par[4] * log(10) * (exp(-10^(allfit$par[4] * (xangus - allfit$par[1])))) * (allfit$par[2] * ((10^xangus / 10^allfit$par[1])^(allfit$par[3] + 1)))
lines(xangus, log10(yangus), col = rgb((100 / 255), (149 / 255), (237 / 255), 1), lwd = 2, lty = 3)
bestomegamatter <- sum((yangus * 10^xangus) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
bestomegamatter2 <- sum((yangus[xangus > 12.7] * 10^xangus[xangus > 12.7]) * fitbinwid) * msol / (1E6 * parsec)^3 / rhocrit
#
# Plot data points
#
points(reflex$x, reflex$Curve1, pch = 18, col = "forestgreen", cex = 1.0)
magerr(reflex$x, reflex$Curve1, yhi = log10(1 + reflexf), ylo = log10(1 - reflexf), col = "forestgreen")
points(tpigg$mass, tpigg$density, pch = 16, col = "grey50", cex = 1.0)
magerr(tpigg$mass, tpigg$density, yhi = tpigg$up, ylo = -tpigg$do, col = "grey50")
points(sdss$V1, sdss$V4, pch = 10, col = "purple", cex = 1.0)
magerr(sdss$V1, sdss$V4, yhi = log10(1 + sdss$V8), ylo = log10(1 - sdss$V8), col = "purple")
points(gama$V1, gama$V4, pch = 16, col = "red", cex = 1.0)
magerr(gama$V1, gama$V4, yhi = log10(1 + gama$V8), ylo = log10(1 - gama$V8), col = "red")
#
# Key
#
points(12.75, -5.6, pch = 16, col = "red", cex = 1)
text(12.75, -5.6, paste0(" GAMA z<", zlimit, " and N>4"), col = "red", pos = 4, cex = 1.0)
points(12.75, -6.0, pch = 10, col = "purple", cex = 1)
text(12.75, -6.0, paste0(" SDSS DR12, z<", zlimitsdss, " and N>4"), col = "purple", pos = 4, cex = 1.0)
points(12.75, -6.4, pch = 18, col = "forestgreen", cex = 1)
text(12.75, -6.4, " REFLEX II, x-ray, z~0.1 (Bohringer et al. 2017)", col = "forestgreen", pos = 4, cex = 1.0)
points(12.75, -6.8, pch = 16, col = "grey50", cex = 1)
text(12.75, -6.8, " 2PIGG z < 0.12 (Eke et al. 2008)", col = "grey50", pos = 4, cex = 1.0)
lines(c(12.71, 12.79), c(-7.2, -7.2), col = rgb(0.5, 0.5, 0.5, 0.1), lwd = 5)
lines(c(12.7, 12.8), c(-7.2, -7.2), col = "cyan")
text(12.75, -7.2, paste0(" SDSS DR10 (Tempel et al 2014)"), col = "cyan", pos = 4, cex = 1.0)
lines(c(12.71, 12.79), c(-7.6, -7.6), col = rgb((100 / 255), (149 / 255), (237 / 255), 0.75), lwd = 2, lty = 3)
# text(13,-7.6,expression(" Best MRP fn fit"),col="black",pos=4)
text(12.75, -7.6, paste0(" Best fit MRP function to ", myoption), col = "black", pos = 4)
lines(c(12.7, 12.8), c(-8.0, -8.0), col = "black", lty = 2, lwd = 2)
text(12.75, -8.0, " LCDM expectation from MRP", col = "black", pos = 4)
#
par(fig = c(0.6, 0.98, 0.6, 0.98), new = TRUE)
x <- seq(-10, 50, 0.01)
info <- maghist(mass, breaks = x, xlab = expression(Omega[M]), grid = FALSE, ylab = "Frequency", majorn = c(3, 2), xlim = c(0.0, 1), col = rgb((100 / 255), (149 / 255), (237 / 255), 1.0), verbose = FALSE)
polygon(c(quantile(mass, 0.84), quantile(mass, 0.84), quantile(mass, 0.16), quantile(mass, 0.16)), c(0, 1E5, 1E5, 0), border = NA, col = rgb(1, 0, 0, 0.25))
abline(v = (omegam), col = "black", lwd = 3, lty = 3)
# abline(v=(omegam+omegamerr),col="black",lwd=1,lty=2)
# abline(v=(omegam-omegamerr),col="black",lwd=1,lty=2)
# abline(v=(quantile(mass,0.16)),col="blue",lwd=1,lty=2)
abline(v = (quantile(mass, 0.50)), col = "red", lwd = 1, lty = 1)
# abline(v=(quantile(mass,0.84)),col="blue",lwd=1,lty=2)
# text(0.6,0.975*max(info$counts),myoption,pos=4,cex=0.75)
text(0.6, 0.75 * max(info$counts), expression(Omega["M"] ~ "="), cex = 1.0, col = "red")
text(0.8, 0.75 * max(info$counts), paste0(signif((quantile(mass, 0.50)), 2), pos = 4), cex = 1, col = "red")
up <- (quantile(mass, 0.84) - quantile(mass, 0.5))
do <- (quantile(mass, 0.5) - quantile(mass, 0.16))
text(0.9, 0.925 * max(info$counts), paste0("+", signif(up, 2), pos = 4), cex = 0.75, col = "red")
text(0.9, 0.575 * max(info$counts), paste0("-", signif(do, 2), pos = 4), cex = 0.75, col = "red")


tmp <- dev.off()
#
png(filename = paste0("./allcov", myoption, ".png"), width = 16.0, height = 16.0, units = "cm", res = 240)
par(mar = c(3, 3.5, 0.25, 0.25), oma = c(0, 0, 0, 0))
params <- cbind(mstar[!is.na(mstar)], log10(phistar[!is.na(mstar)]), alphastar[!is.na(mstar)], betastar[!is.na(mstar)])
fitvals <- mymagtri(params,
  samples = length(params[, 1]), samptype = "ran",
  #               refvals=c(allfit3$par[1],log10(allfit3$par[2]),allfit3$par[3],allfit3$par[4]),
  lab = c("log\u2081\u2080(M*) [M\u2092]", "log\u2081\u2080(\u03c6*) [Mpc\u207B\u00b3]", "\u03b1", "\u03b2"),
  refvals = c(mstarmrp - 0.075, log10(phimrp) + 0.075, alphamrp, betamrp),
  bestvals = c(allfit$par[1], log10(allfit$par[2]), allfit$par[3], allfit$par[4])
)
cat(
  paste0(myoption),
  allfit$par[1], quantile(mstar, 0.16, na.rm = TRUE), quantile(mstar, 0.84, na.rm = TRUE),
  log10(allfit$par[2]), log10(quantile(phistar, 0.16, na.rm = TRUE)), log10(quantile(phistar, 0.84, na.rm = TRUE)),
  allfit$par[3], quantile(alphastar, 0.16, na.rm = TRUE), quantile(alphastar, 0.84, na.rm = TRUE),
  allfit$par[4], quantile(betastar, 0.16, na.rm = TRUE), quantile(betastar, 0.84, na.rm = TRUE),
  bestomegamatter2, bestomegamatter2 - quantile(omegamatter2, 0.16), quantile(omegamatter2, 0.84) - bestomegamatter2,
  bestomegamatter2 / omegam, bestomegamatter2 / omegam - quantile(omegamatter2, 0.16) / omegam, quantile(omegamatter2, 0.84) / omegam - bestomegamatter2 / omegam, "\n"
)
tmp <- dev.off()
#
