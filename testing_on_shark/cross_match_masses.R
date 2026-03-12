############################################################
# NESSIE vs TRUE HALO MASS COMPARISON
#
# Cross-match Nessie groups to true SHARK halos by
# position + redshift, then compare mass scales to
# understand the calibration offset.
############################################################

library(arrow)
library(Cairo)
library(celestial)
library(data.table)

set.seed(42)

cat("============================================\n")
cat("  NESSIE vs TRUE MASS CROSS-MATCH\n")
cat("============================================\n\n")

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"

############################################################
# 1. Read catalogues
############################################################

cat("Reading catalogues...\n")
nessie <- as.data.table(read_parquet("nessie_groups.parquet"))
groups_shark <- as.data.table(read_parquet(file.path(data_dir, "groups_shark.parquet")))
groups_shark <- groups_shark[dec > 0]

# Nessie mass columns
nessie$log_mass_proxy <- log10(nessie$mass_proxy * 10)
nessie$log_mass_est   <- log10(nessie$estimated_mass)

cat(sprintf("  Nessie: %d groups\n", nrow(nessie)))
cat(sprintf("  SHARK:  %d halos (dec > 0)\n", nrow(groups_shark)))

# Also compute abundance-matched masses for SHARK
zmin <- 0.01; zlimit <- 0.25
TRUE_MS <- 13.51; TRUE_LP <- -3.19; TRUE_AL <- -1.27; TRUE_BE <- 0.47

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}

ho <- 67.37; omegam <- 0.3147
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

ra_range  <- range(groups_shark$ra, na.rm=TRUE)
dec_range <- range(groups_shark$dec, na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)
d_lo <- cosdist(zmin); d_hi <- cosdist(zlimit)
Vsurvey <- (4/3) * pi * (d_hi^3 - d_lo^3) * sky_frac

groups_vol <- groups_shark[redshift_cosmological > zmin &
                           redshift_cosmological < zlimit &
                           mass_virial > 0]

m_grid <- seq(9, 16.5, by=0.001)
phi_grid <- mrp_phi(m_grid, TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)
cum_counts <- rev(cumsum(rev(phi_grid * 0.001))) * Vsurvey

groups_vol[, rank := frankv(mass_virial, order=-1L)]
groups_vol[, log_mass_am := approx(x=rev(cum_counts), y=rev(m_grid),
                                    xout=rank, rule=2)$y]
groups_vol[, log_mass_virial := log10(mass_virial)]

cat(sprintf("  SHARK in volume: %d halos\n", nrow(groups_vol)))

############################################################
# 2. Cross-match by position + redshift
############################################################

cat("\nCross-matching...\n")

# Filter Nessie to volume
nes_vol <- nessie[iter_redshift > zmin & iter_redshift < zlimit & iter_dec > 0]
cat(sprintf("  Nessie in volume: %d\n", nrow(nes_vol)))

# Match each Nessie group to nearest SHARK halo
# Use angular separation + redshift proximity
# Convert to 3D comoving coordinates for matching

nes_vol[, cos_dec := cos(iter_dec * pi/180)]
nes_vol[, x3d := cosdist(iter_redshift) * cos(iter_dec*pi/180) * cos(iter_ra*pi/180), by=1:nrow(nes_vol)]

# This is slow for large N, use a simpler approach:
# For each Nessie group, find SHARK halos within dz < 0.01 and angular sep < 0.5 deg
# Then pick closest

match_shark_id <- rep(NA_integer_, nrow(nes_vol))
match_sep_deg  <- rep(NA_real_, nrow(nes_vol))
match_dz       <- rep(NA_real_, nrow(nes_vol))

cat("  Matching (this may take a minute)...\n")

for(i in 1:nrow(nes_vol)) {
    # Redshift window
    dz_max <- 0.015
    z_nes <- nes_vol$iter_redshift[i]
    candidates <- which(abs(groups_vol$redshift_cosmological - z_nes) < dz_max)
    
    if(length(candidates) == 0) next
    
    # Angular separation
    ra_nes  <- nes_vol$iter_ra[i]
    dec_nes <- nes_vol$iter_dec[i]
    
    dra  <- (groups_vol$ra[candidates] - ra_nes) * cos(dec_nes * pi/180)
    ddec <- groups_vol$dec[candidates] - dec_nes
    sep  <- sqrt(dra^2 + ddec^2)
    
    best <- which.min(sep)
    if(sep[best] < 1.0) {  # within 1 degree
        match_shark_id[i] <- candidates[best]
        match_sep_deg[i]  <- sep[best]
        match_dz[i]       <- abs(groups_vol$redshift_cosmological[candidates[best]] - z_nes)
    }
}

n_matched <- sum(!is.na(match_shark_id))
cat(sprintf("  Matched: %d / %d (%.1f%%)\n", n_matched, nrow(nes_vol),
            100*n_matched/nrow(nes_vol)))
cat(sprintf("  Median sep: %.3f deg\n", median(match_sep_deg, na.rm=TRUE)))
cat(sprintf("  Median dz: %.4f\n", median(match_dz, na.rm=TRUE)))

############################################################
# 3. Build matched catalogue
############################################################

matched <- !is.na(match_shark_id)

df_match <- data.table(
    log_mass_proxy = nes_vol$log_mass_proxy[matched],
    log_mass_est   = nes_vol$log_mass_est[matched],
    log_mass_am    = groups_vol$log_mass_am[match_shark_id[matched]],
    log_mass_vir   = groups_vol$log_mass_virial[match_shark_id[matched]],
    z_nessie       = nes_vol$iter_redshift[matched],
    z_shark        = groups_vol$redshift_cosmological[match_shark_id[matched]],
    nfof_nessie    = nes_vol$multiplicity[matched],
    sep_deg        = match_sep_deg[matched],
    dz             = match_dz[matched]
)

# Remove bad masses
df_match <- df_match[is.finite(log_mass_proxy) & is.finite(log_mass_am) &
                     is.finite(log_mass_est)]

cat(sprintf("\n  Clean matches: %d\n", nrow(df_match)))

############################################################
# 4. Compare mass scales
############################################################

cat("\n  --- Mass comparison ---\n")

# mass_proxy * 10 vs AM
offset_proxy <- median(df_match$log_mass_proxy - df_match$log_mass_am)
scatter_proxy <- mad(df_match$log_mass_proxy - df_match$log_mass_am)
cat(sprintf("  mass_proxy*10 vs AM:  offset = %+.3f dex, scatter = %.3f dex\n",
            offset_proxy, scatter_proxy))

# estimated_mass vs AM
offset_est <- median(df_match$log_mass_est - df_match$log_mass_am)
scatter_est <- mad(df_match$log_mass_est - df_match$log_mass_am)
cat(sprintf("  estimated_mass vs AM: offset = %+.3f dex, scatter = %.3f dex\n",
            offset_est, scatter_est))

# mass_proxy*10 vs virial
offset_vir <- median(df_match$log_mass_proxy - df_match$log_mass_vir)
scatter_vir <- mad(df_match$log_mass_proxy - df_match$log_mass_vir)
cat(sprintf("  mass_proxy*10 vs virial: offset = %+.3f dex, scatter = %.3f dex\n",
            offset_vir, scatter_vir))

# estimated_mass vs virial
offset_est_vir <- median(df_match$log_mass_est - df_match$log_mass_vir)
scatter_est_vir <- mad(df_match$log_mass_est - df_match$log_mass_vir)
cat(sprintf("  estimated_mass vs virial: offset = %+.3f dex, scatter = %.3f dex\n",
            offset_est_vir, scatter_est_vir))

# Linear fit: proxy vs AM
fit_proxy_am <- lm(log_mass_proxy ~ log_mass_am, data=df_match)
cat(sprintf("\n  Linear fit (proxy*10 = a + b*AM):\n"))
cat(sprintf("    a = %.3f, b = %.3f\n", coef(fit_proxy_am)[1], coef(fit_proxy_am)[2]))

# Linear fit: est vs AM
fit_est_am <- lm(log_mass_est ~ log_mass_am, data=df_match)
cat(sprintf("  Linear fit (est = a + b*AM):\n"))
cat(sprintf("    a = %.3f, b = %.3f\n", coef(fit_est_am)[1], coef(fit_est_am)[2]))

# By multiplicity
cat("\n  --- Offset by multiplicity ---\n")
nfof_bins <- c(2, 5, 10, 20, 50, 200)
for(ib in 1:(length(nfof_bins)-1)) {
    in_bin <- df_match$nfof_nessie >= nfof_bins[ib] & df_match$nfof_nessie < nfof_bins[ib+1]
    if(sum(in_bin) > 10) {
        off <- median(df_match$log_mass_proxy[in_bin] - df_match$log_mass_am[in_bin])
        sc  <- mad(df_match$log_mass_proxy[in_bin] - df_match$log_mass_am[in_bin])
        cat(sprintf("  Nfof %d-%d: offset=%+.3f, scatter=%.3f (N=%d)\n",
                    nfof_bins[ib], nfof_bins[ib+1]-1, off, sc, sum(in_bin)))
    }
}

############################################################
# 5. Plot
############################################################

CairoPDF("nessie_mass_comparison.pdf", 14, 10)
par(mfrow=c(2,3))

# Panel 1: mass_proxy*10 vs AM mass
plot(df_match$log_mass_am, df_match$log_mass_proxy,
     pch=".", col=rgb(0,0,0,0.2), cex=2,
     xlim=c(11, 15.5), ylim=c(11, 15.5),
     xlab="True AM mass (log10)", ylab="mass_proxy*10 (log10)",
     main="mass_proxy*10 vs AM")
abline(0, 1, col="red", lwd=2)
abline(offset_proxy, 1, col="blue", lwd=2, lty=2)
abline(coef(fit_proxy_am), col="darkgreen", lwd=2, lty=3)
grid(col="gray80")
legend("topleft", c("1:1", sprintf("offset %+.2f", offset_proxy), "linear fit"),
       col=c("red","blue","darkgreen"), lty=c(1,2,3), lwd=2, cex=0.7)

# Panel 2: estimated_mass vs AM mass
plot(df_match$log_mass_am, df_match$log_mass_est,
     pch=".", col=rgb(0,0,0,0.2), cex=2,
     xlim=c(11, 15.5), ylim=c(11, 15.5),
     xlab="True AM mass (log10)", ylab="estimated_mass (log10)",
     main="estimated_mass vs AM")
abline(0, 1, col="red", lwd=2)
abline(offset_est, 1, col="blue", lwd=2, lty=2)
abline(coef(fit_est_am), col="darkgreen", lwd=2, lty=3)
grid(col="gray80")
legend("topleft", c("1:1", sprintf("offset %+.2f", offset_est), "linear fit"),
       col=c("red","blue","darkgreen"), lty=c(1,2,3), lwd=2, cex=0.7)

# Panel 3: Residuals vs AM mass
resid_proxy <- df_match$log_mass_proxy - df_match$log_mass_am
plot(df_match$log_mass_am, resid_proxy,
     pch=".", col=rgb(0,0,0,0.2), cex=2,
     xlim=c(11, 15.5), ylim=c(-3, 3),
     xlab="True AM mass (log10)", ylab="mass_proxy*10 - AM",
     main="Residuals vs mass")
abline(h=0, col="red", lwd=2)
abline(h=c(-1,1)*scatter_proxy, col="blue", lwd=1, lty=2)
grid(col="gray80")

# Running median
am_order <- order(df_match$log_mass_am)
am_sorted <- df_match$log_mass_am[am_order]
res_sorted <- resid_proxy[am_order]
nw <- 200
if(length(am_sorted) > nw) {
    run_med <- zoo::rollmedian(res_sorted, nw, fill=NA)
    lines(am_sorted, run_med, col="red", lwd=3)
}

# Panel 4: Histograms of mass distributions
xlims <- c(10, 16)
hist(df_match$log_mass_am, breaks=seq(10,16,by=0.2), col=rgb(0,0,1,0.3),
     border="white", main="Mass distributions (matched)",
     xlab="log10(M)", xlim=xlims)
hist(df_match$log_mass_proxy, breaks=seq(10,16,by=0.2), col=rgb(1,0,0,0.3),
     border="white", add=TRUE)
hist(df_match$log_mass_est, breaks=seq(10,16,by=0.2), col=rgb(0,0.6,0,0.3),
     border="white", add=TRUE)
legend("topright", c("AM (true)", "proxy*10", "estimated"),
       fill=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3), rgb(0,0.6,0,0.3)), cex=0.7)

# Panel 5: Scatter vs multiplicity
boxplot(resid_proxy ~ cut(df_match$nfof_nessie, breaks=c(1,4,7,10,20,50,200)),
        main="Residual vs multiplicity",
        xlab="Nfof bin", ylab="mass_proxy*10 - AM",
        col="steelblue", outline=FALSE)
abline(h=0, col="red", lwd=2)

# Panel 6: Match quality
plot(df_match$sep_deg, df_match$dz,
     pch=".", col=rgb(0,0,0,0.3), cex=2,
     xlab="Angular separation (deg)", ylab="|dz|",
     main="Match quality")
grid(col="gray80")

dev.off()
cat("\nPlot saved: nessie_mass_comparison.pdf\n\nDone!\n")
