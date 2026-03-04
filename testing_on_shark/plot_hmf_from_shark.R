############################################################
# SHARK HMF: Binned mass function + simple MRP fit
############################################################

library(arrow)
library(Cairo)
library(celestial)

set.seed(42)

############################################################
# 1. Read and select
############################################################

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"
groups <- read_parquet(file.path(data_dir, "groups_shark.parquet"))

zmin   <- 0.01
zlimit <- 0.25
multi  <- 5

# All groups (no multiplicity cut) - the TRUE HMF
# dec > 0 selects one footprint only
sel_all <- groups$redshift_cosmological > zmin & 
           groups$redshift_cosmological < zlimit &
           groups$mass_virial > 0 &
           groups$dec > 0
log_mass_all <- log10(groups$mass_virial[sel_all])

# With multiplicity cut (GAMA-like)
sel <- sel_all & groups$number_galaxies_selected >= multi
log_mass <- log10(groups$mass_virial[sel])

cat(sprintf("N groups (all):      %d\n", sum(sel_all)))
cat(sprintf("N groups (Nsel>=%d): %d\n", multi, length(log_mass)))

############################################################
# 2. Survey volume
############################################################

ho <- 67.37; omegam <- 0.3147

cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

ra_range  <- range(groups$ra[sel_all], na.rm=TRUE)
dec_range <- range(groups$dec[sel_all], na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
cat(sprintf("RA range: %.2f -- %.2f\n", ra_range[1], ra_range[2]))
cat(sprintf("Dec range: %.2f -- %.2f\n", dec_range[1], dec_range[2]))
cat(sprintf("Sky area: %.2f deg^2\n", sky_area_deg2))
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)

d_lo <- cosdist(zmin)
d_hi <- cosdist(zlimit)
Vsurvey <- (4/3) * pi * (d_hi^3 - d_lo^3) * sky_frac
cat(sprintf("Volume: %.4e Mpc^3\n", Vsurvey))

############################################################
# 3. Bin to get HMF
############################################################

bin_width <- 0.2
breaks <- seq(9, 16, by=bin_width)

h_all   <- hist(log_mass_all, breaks=breaks, plot=FALSE)
phi_all <- h_all$counts / (Vsurvey * bin_width)
ok_all  <- h_all$counts >= 5

h       <- hist(log_mass, breaks=breaks, plot=FALSE)
phi     <- h$counts / (Vsurvey * bin_width)
ok      <- h$counts >= 5

cat(sprintf("\nBinned HMF (selected bins):\n"))
cat("  log10(M)  N_all    phi_all       N_cut    phi_cut\n")
for(i in 1:length(h$mids)) {
    if(ok_all[i] || ok[i]) {
        cat(sprintf("  %6.2f   %6d  %.3e    %5d  %.3e\n",
                    h$mids[i], h_all$counts[i], phi_all[i], h$counts[i], phi[i]))
    }
}

############################################################
# 4. Fit MRP to the TRUE HMF
############################################################

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}

MASS_FIT_FLOOR <- 12.0
fit_ok <- ok_all & h_all$mids >= MASS_FIT_FLOOR & phi_all > 0

# Data for fitting
fit_x   <- h_all$mids[fit_ok]
fit_y   <- log10(phi_all[fit_ok])
fit_err <- 1 / (sqrt(h_all$counts[fit_ok]) * log(10))

chi2 <- function(par) {
    ms <- par[1]; lp <- par[2]; al <- par[3]; be <- par[4]
    model <- mrp_phi(fit_x, ms, lp, al, be)
    model[model < 1e-30] <- 1e-30
    sum(((fit_y - log10(model)) / fit_err)^2)
}

# L-BFGS-B — M* can go down to 10 now
best_fit <- NULL; best_chi2 <- Inf
starts <- list(
    c(13.5, -3.2, -1.3, 0.5),
    c(14.0, -3.5, -1.5, 0.6),
    c(13.0, -2.8, -1.0, 0.4),
    c(13.8, -4.0, -1.7, 0.7),
    c(12.5, -2.0, -1.0, 0.3),
    c(14.5, -5.0, -1.6, 0.5)
)
for(s in starts) {
    this_fit <- tryCatch(
        optim(s, chi2, method="L-BFGS-B",
              lower=c(10, -8, -3, 0.1), upper=c(16, 0, 0, 2.0),
              control=list(maxit=100000)),
        error=function(e) list(value=Inf))
    cat(sprintf("  start=(%.1f,%.1f,%.1f,%.1f) -> chi2=%.1f par=(%.3f,%.3f,%.3f,%.3f)\n",
                s[1],s[2],s[3],s[4], this_fit$value,
                this_fit$par[1],this_fit$par[2],this_fit$par[3],this_fit$par[4]))
    if(this_fit$value < best_chi2) { best_chi2 <- this_fit$value; best_fit <- this_fit }
}

ms <- best_fit$par[1]; lp <- best_fit$par[2]
al <- best_fit$par[3]; be <- best_fit$par[4]

cat(sprintf("\nMRP fit (M >= %.1f, chi2 = %.1f, %d bins):\n",
            MASS_FIT_FLOOR, best_fit$value, sum(fit_ok)))
cat(sprintf("  M*       = %.3f\n", ms))
cat(sprintf("  log_phi* = %.3f\n", lp))
cat(sprintf("  alpha    = %.3f\n", al))
cat(sprintf("  beta     = %.3f\n", be))

if(ms <= 10.01) cat("  WARNING: M* hit lower bound — true MRP knee may be below fit range\n")

driver22 <- c(ms=13.51, lp=-3.19, al=-1.27, be=0.47)
cat(sprintf("\nDriver+22: M*=%.2f  lp=%.2f  al=%.2f  be=%.2f\n",
            driver22[1], driver22[2], driver22[3], driver22[4]))

############################################################
# 5. Plot
############################################################

xfit <- seq(10, 16, length.out=500)
phi_all_err <- sqrt(h_all$counts) / (Vsurvey * bin_width)

CairoPDF("shark_hmf.pdf", 12, 8)
par(mfrow=c(1,1))

plot(h_all$mids[ok_all], log10(phi_all[ok_all]),
     pch=19, col="grey40", cex=1.3,
     xlim=c(10.5, 16), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M"[vir]*" / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main=sprintf("SHARK HMF (%.2f < z < %.2f, dec > 0)", zmin, zlimit))
grid(col="gray80")

arrows(h_all$mids[ok_all], log10(pmax(phi_all[ok_all] - phi_all_err[ok_all], 1e-30)),
       h_all$mids[ok_all], log10(phi_all[ok_all] + phi_all_err[ok_all]),
       angle=90, code=3, length=0.03, col="grey40")

points(h$mids[ok], log10(phi[ok]), pch=17, col="darkgreen", cex=1.2)

lines(xfit, log10(pmax(mrp_phi(xfit, ms, lp, al, be), 1e-30)),
      col="red", lwd=3)
lines(xfit, log10(pmax(mrp_phi(xfit, driver22[1], driver22[2], driver22[3], driver22[4]), 1e-30)),
      col="blue", lwd=2, lty=2)
abline(v=MASS_FIT_FLOOR, col="orange", lwd=2, lty=4)

legend("topright",
       legend=c("All groups (true HMF)", 
                sprintf("Nsel >= %d", multi),
                sprintf("MRP: M*=%.2f, log\u03c6*=%.2f, \u03b1=%.2f, \u03b2=%.2f", ms, lp, al, be),
                "Driver+22",
                sprintf("Fit range (M > %.1f)", MASS_FIT_FLOOR)),
       col=c("grey40","darkgreen","red","blue","orange"),
       pch=c(19,17,NA,NA,NA), lty=c(NA,NA,1,2,4), lwd=c(NA,NA,3,2,2),
       bg="white", cex=0.7)

dev.off()
cat("\nPlot saved: shark_hmf.pdf\n")
