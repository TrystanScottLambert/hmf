############################################################
# DIAGNOSTIC: Fit empirical M*-z and M_halo-M* relations
# Converts your Python script to R and fits the relations
# needed for truncation modeling
############################################################

library(Rfits)
library(data.table)
library(Cairo)

############################################################
# Configuration
############################################################

GROUP_FILE <- "../data/G3CFoFGroupv10.fits"
GALAXY_FILE <- "../data/GAMAGalsInGroups.csv"
N_CUT <- 5  # Require at least 5 members (matches multi=5)

############################################################
# Read data (matching your Python script)
############################################################

cat("Reading group and galaxy data...\n")

# Read groups file
hdu_groups <- Rfits_read_table(GROUP_FILE)
df_groups <- hdu_groups[hdu_groups$Nfof >= N_CUT & hdu_groups$GroupID < 400000, ]

# Read galaxies file
df_galaxies <- fread(GALAXY_FILE)

cat("Initial groups:", nrow(df_groups), "\n")
cat("Initial galaxies:", nrow(df_galaxies), "\n")

############################################################
# Get Nth member for each group (matching your Python agg)
############################################################

cat("Extracting Nth brightest member per group...\n")

# Initialize output dataframe
df_members <- data.frame(
    GroupID = integer(),
    logmstar = numeric(),
    Z = numeric()
)

# For each group, find the Nth brightest member
for(group_id in df_groups$GroupID) {
    # Get all members of this group
    members <- df_galaxies[df_galaxies$GroupID == group_id, ]
    
    # Sort by stellar mass (descending)
    members <- members[order(members$logmstar, decreasing=TRUE), ]
    
    # Take Nth member (if it exists)
    if(nrow(members) >= N_CUT) {
        df_members <- rbind(df_members, data.frame(
            GroupID = group_id,
            logmstar = members$logmstar[N_CUT],
            Z = members$Z[N_CUT]
        ))
    }
}

# Match groups to members (only keep groups with Nth member)
df_groups <- df_groups[df_groups$GroupID %in% df_members$GroupID, ]

# Sort both by GroupID to match
df_groups <- df_groups[order(df_groups$GroupID), ]
df_members <- df_members[order(df_members$GroupID), ]

cat("Final matched groups:", nrow(df_groups), "\n\n")

############################################################
# Extract arrays (matching your Python)
############################################################

group_masses <- log10(df_groups$MassAfunc)
galaxy_masses <- df_members$logmstar
galaxy_redshift <- df_members$Z

# Filter valid data (matching your Python 'cut')
# CRITICAL: Also filter to survey redshift range (z < 0.25)
cut <- galaxy_masses > 0 & 
       is.finite(galaxy_masses) & 
       is.finite(group_masses) &
       galaxy_redshift < 0.25 & galaxy_redshift > 0.01  # Survey range

group_masses <- group_masses[cut]
galaxy_masses <- galaxy_masses[cut]
galaxy_redshift <- galaxy_redshift[cut]

cat("Valid data points (within survey z range):", length(group_masses), "\n")
cat("M* range:", round(range(galaxy_masses), 2), "\n")
cat("M_halo range:", round(range(group_masses), 2), "\n")
cat("z range:", round(range(galaxy_redshift), 3), "\n\n")

############################################################
# Plot 1: M_halo vs M* (HSM relation)
############################################################

cat("Creating HSM relation plot...\n")

CairoPDF("HSM_relation_fit.pdf", 10, 7)

plot(group_masses, galaxy_masses, 
     pch=20, cex=0.3, col=rgb(0,0,0,0.3),
     xlab=expression("log M"[halo]),
     ylab=expression("log M"["*"]),
     main="Halo Mass vs Stellar Mass (Nth Member)",
     xlim=c(10, 16), ylim=c(6, 12))

grid(col="gray80")

# Fit linear relation: log M_halo = a + b * log M*
fit_halo_mstar <- lm(group_masses ~ galaxy_masses)

cat("\n=== M_halo vs M* RELATION ===\n")
cat("log M_halo = ", round(coef(fit_halo_mstar)[1], 3), 
    " + ", round(coef(fit_halo_mstar)[2], 3), " × log M*\n", sep="")
cat("R² = ", round(summary(fit_halo_mstar)$r.squared, 3), "\n")
cat("Residual SD = ", round(summary(fit_halo_mstar)$sigma, 3), " dex\n")

# Add fit line
abline(fit_halo_mstar, col="red", lwd=2)

# Add ±1σ lines
sigma_fit <- summary(fit_halo_mstar)$sigma
mstar_seq <- seq(6, 12, length.out=100)
mhalo_pred <- predict(fit_halo_mstar, newdata=data.frame(galaxy_masses=mstar_seq))
lines(mhalo_pred + sigma_fit, mstar_seq, col="red", lty=2)
lines(mhalo_pred - sigma_fit, mstar_seq, col="red", lty=2)

legend("topleft", 
       legend=c("Data", "Linear fit", "±1σ"),
       pch=c(20, NA, NA), lty=c(NA, 1, 2),
       col=c(rgb(0,0,0,0.3), "red", "red"),
       bg="white")

text(15, 7, sprintf("log M_halo = %.2f + %.2f × log M*\nR² = %.3f", 
                    coef(fit_halo_mstar)[1], 
                    coef(fit_halo_mstar)[2],
                    summary(fit_halo_mstar)$r.squared),
     pos=4, col="red")

dev.off()

############################################################
# Plot 2: M* vs z (SMZ relation) - CRITICAL FOR TRUNCATION
############################################################

cat("Creating M* vs z relation plot...\n")

CairoPDF("SMZ_relation_fit.pdf", 10, 7)

plot(galaxy_redshift, galaxy_masses,
     pch=20, cex=0.3, col=rgb(0,0,0,0.3),
     xlab="Redshift",
     ylab=expression("log M"["*"]),
     main="Stellar Mass vs Redshift (Nth Member) - Detection Envelope",
     xlim=c(0, 0.3), ylim=c(6, 12))

grid(col="gray80")

# Fit the LOWER ENVELOPE (detection limit)
# This is the key for truncation modeling!

# Method: Bin by redshift and find 5th percentile
z_bins <- seq(0, 0.25, by=0.02)
z_centers <- (z_bins[-1] + z_bins[-length(z_bins)]) / 2
mstar_envelope <- rep(NA, length(z_centers))

for(i in 1:length(z_centers)) {
    in_bin <- galaxy_redshift >= z_bins[i] & galaxy_redshift < z_bins[i+1]
    if(sum(in_bin) > 10) {
        # 5th percentile gives lower envelope (detection limit)
        mstar_envelope[i] <- quantile(galaxy_masses[in_bin], 0.05, na.rm=TRUE)
    }
}

# Fit linear relation to envelope
ok_bins <- !is.na(mstar_envelope)
fit_envelope <- lm(mstar_envelope[ok_bins] ~ z_centers[ok_bins])

cat("\n=== M* vs z ENVELOPE (DETECTION LIMIT) ===\n")
cat("log M*(z) = ", round(coef(fit_envelope)[1], 3),
    " + ", round(coef(fit_envelope)[2], 3), " × z\n", sep="")
cat("R² = ", round(summary(fit_envelope)$r.squared, 3), "\n")

# Plot envelope points
points(z_centers[ok_bins], mstar_envelope[ok_bins], 
       pch=19, cex=1.2, col="blue")

# Plot envelope fit - use the fitted line directly
z_seq <- seq(0, 0.25, length.out=100)
mstar_lim_pred <- coef(fit_envelope)[1] + coef(fit_envelope)[2] * z_seq
lines(z_seq, mstar_lim_pred, col="blue", lwd=3)

# Add shaded region below detection limit (unobservable)
polygon(c(z_seq, rev(z_seq)), 
        c(mstar_lim_pred, rep(5, length(z_seq))),
        col=rgb(1,0,0,0.1), border=NA)

legend("topleft",
       legend=c("Data", "5th percentile", "Detection limit", "Unobservable"),
       pch=c(20, 19, NA, 15), lty=c(NA, NA, 1, NA),
       col=c(rgb(0,0,0,0.3), "blue", "blue", rgb(1,0,0,0.3)),
       pt.cex=c(0.5, 1.2, NA, 2),
       bg="white")

text(0.18, 7, sprintf("log M*(z) = %.2f + %.2f × z\nDetection limit",
                      coef(fit_envelope)[1],
                      coef(fit_envelope)[2]),
     pos=4, col="blue", cex=1.1, font=2)

dev.off()

############################################################
# Save fitted relations for use in Stan model
############################################################

cat("\nSaving fitted relations...\n")

relations <- list(
    # M_halo vs M* relation
    halo_mstar_intercept = coef(fit_halo_mstar)[1],
    halo_mstar_slope = coef(fit_halo_mstar)[2],
    halo_mstar_sigma = summary(fit_halo_mstar)$sigma,
    
    # M* vs z envelope (detection limit)
    mstar_z_intercept = coef(fit_envelope)[1],
    mstar_z_slope = coef(fit_envelope)[2]
)

saveRDS(relations, "empirical_relations.rds")

cat("\n=== SUMMARY FOR TRUNCATION MODEL ===\n")
cat("Detection limit: log M*(z) = ", 
    round(relations$mstar_z_intercept, 3), " + ",
    round(relations$mstar_z_slope, 3), " × z\n", sep="")
cat("\nFor a group at zmax:\n")
cat("  1. Find M*_lim(zmax) using the envelope relation\n")
cat("  2. Convert to M_halo_lim using: log M_halo = ",
    round(relations$halo_mstar_intercept, 2), " + ",
    round(relations$halo_mstar_slope, 2), " × log M*\n", sep="")
cat("  3. Use M_halo_lim as truncation in Stan model\n\n")

cat("Diagnostic plots saved:\n")
cat("  - HSM_relation_fit.pdf\n")
cat("  - SMZ_relation_fit.pdf\n")
cat("  - empirical_relations.rds\n")
