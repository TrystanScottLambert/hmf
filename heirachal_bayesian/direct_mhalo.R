############################################################
# ALTERNATIVE: Fit M_halo vs z envelope DIRECTLY
# Skip the M* intermediate step
############################################################

library(Rfits)
library(data.table)
library(Cairo)

############################################################
# Read and prepare data
############################################################

GROUP_FILE <- "../data/G3CFoFGroupv10.fits"
N_CUT <- 5
ho <- 67.37

# Read groups
hdu_groups <- Rfits_read_table(GROUP_FILE)
df_groups <- hdu_groups[hdu_groups$Nfof >= N_CUT & 
                        hdu_groups$GroupID < 400000 &
                        hdu_groups$Zfof > 0.01 &
                        hdu_groups$Zfof < 0.25, ]  # Survey range only!

# Get halo masses and redshifts
group_masses <- log10(df_groups$MassAfunc)
group_redshift <- df_groups$Zfof

# Filter valid
ok <- is.finite(group_masses) & is.finite(group_redshift)
group_masses <- group_masses[ok]
group_redshift <- group_redshift[ok]

cat("N groups:", length(group_masses), "\n")
cat("M_halo range:", round(range(group_masses), 2), "\n")
cat("z range:", round(range(group_redshift), 3), "\n\n")

############################################################
# Fit M_halo vs z envelope (lower boundary)
############################################################

cat("Fitting M_halo-z envelope...\n")

# Bin by redshift and find lower envelope (5th percentile)
z_bins <- seq(0.01, 0.25, by=0.02)
z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
mhalo_envelope <- rep(NA, length(z_centers))

for(i in 1:length(z_centers)) {
    in_bin <- group_redshift >= z_bins[i] & group_redshift < z_bins[i+1]
    if(sum(in_bin) > 10) {
        mhalo_envelope[i] <- quantile(group_masses[in_bin], 0.05, na.rm=TRUE)
    }
}

# Fit linear relation
ok_bins <- !is.na(mhalo_envelope)
fit_mhalo_z <- lm(mhalo_envelope[ok_bins] ~ z_centers[ok_bins])

cat("\n=== M_halo vs z ENVELOPE ===\n")
cat("log M_halo(z) = ", round(coef(fit_mhalo_z)[1], 3),
    " + ", round(coef(fit_mhalo_z)[2], 3), " × z\n", sep="")
cat("R² = ", round(summary(fit_mhalo_z)$r.squared, 3), "\n\n")

############################################################
# Plot
############################################################

CairoPDF("Mhalo_z_envelope.pdf", 10, 7)

plot(group_redshift, group_masses,
     pch=20, cex=0.3, col=rgb(0,0,0,0.3),
     xlab="Redshift",
     ylab=expression("log M"[halo]),
     main="Halo Mass vs Redshift - Direct Detection Envelope",
     xlim=c(0, 0.3), ylim=c(10, 17))

grid(col="gray80")

# Plot envelope points
points(z_centers[ok_bins], mhalo_envelope[ok_bins], 
       pch=19, cex=1.5, col="red")

# Plot envelope fit
z_seq <- seq(0.01, 0.25, length.out=100)
mhalo_lim_pred <- coef(fit_mhalo_z)[1] + coef(fit_mhalo_z)[2] * z_seq
lines(z_seq, mhalo_lim_pred, col="red", lwd=3)

# Shade unobservable region
polygon(c(z_seq, rev(z_seq)), 
        c(mhalo_lim_pred, rep(9, length(z_seq))),
        col=rgb(1,0,0,0.1), border=NA)

legend("topleft",
       legend=c("Data", "5th percentile", "Detection limit", "Unobservable"),
       pch=c(20, 19, NA, 15), lty=c(NA, NA, 1, NA),
       col=c(rgb(0,0,0,0.3), "red", "red", rgb(1,0,0,0.3)),
       pt.cex=c(0.5, 1.5, NA, 2),
       bg="white")

text(0.18, 11, sprintf("log M_halo(z) = %.2f + %.2f × z\nDirect detection limit",
                       coef(fit_mhalo_z)[1],
                       coef(fit_mhalo_z)[2]),
     pos=4, col="red", cex=1.1, font=2)

dev.off()

############################################################
# Save relation
############################################################

direct_relation <- list(
    mhalo_z_intercept = coef(fit_mhalo_z)[1],
    mhalo_z_slope = coef(fit_mhalo_z)[2]
)

saveRDS(direct_relation, "direct_mhalo_z_relation.rds")

cat("Saved: direct_mhalo_z_relation.rds\n")
cat("Plot: Mhalo_z_envelope.pdf\n\n")

cat("=== FOR TRUNCATION ===\n")
cat("M_halo_lim(zmax) = ", 
    round(direct_relation$mhalo_z_intercept, 2), " + ",
    round(direct_relation$mhalo_z_slope, 2), " × zmax\n", sep="")
