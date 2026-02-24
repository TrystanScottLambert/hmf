############################################################
# DIAGNOSTIC: Compare our binned data with gamahmf.r
############################################################

library(celestial)
library(Rfits)
library(data.table)
library(plotrix)

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25
zmin   <- 0.01
multi  <- 5
logbin <- 0.2

vol_max_survey <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max_survey^3 * 179.92*(pi/180)^2 / (4*pi)
vlimitmin <- Vsurvey / 1000.0

g3cx <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c  <- g3cx[g3cx$Nfof > multi-1 & g3cx$Zfof < zlimit & 
             g3cx$Zfof > zmin & g3cx$MassAfunc > 1E1 & 
             g3cx$IterCenDec > -3.5, ]

magica <- 13.9
parsec <- 3.0857E16
G      <- 6.67408E-11
msol   <- 1.988E30

g3c$MassAfunc <- g3c$MassAfunc * 100 / ho
g3c$mymass <- magica * (g3c$VelDisp*1000)^2 * g3c$Rad50*parsec*1E6 / (G*msol) * (100/ho)

masscorr <- c(0.0, 0.0, -2.672595e-01, -1.513503e-01, -1.259069e-01, 
              -9.006064e-02, -5.466009e-02, -6.666895e-02, -1.988694e-02,
              -2.439581e-02, -2.067060e-02, -1.812964e-02, -1.556899e-02,
              -1.313664e-02, -1.743112e-02, -7.965513e-03, -1.257178e-02,
              -7.064037e-03, -3.963656e-03, -1.271533e-02, -2.664687e-03,
              -1.691287e-03)

g3c$masscorr <- masscorr[g3c$Nfof]
g3c$masscorr[is.na(g3c$masscorr)] <- 0.0
g3c$mymasscorr <- g3c$mymass / 10^g3c$masscorr
g3c$MassAfunc <- g3c$mymasscorr

# Calculate Vmax EXACTLY as in gamahmf.r
gig <- fread("../data/GAMAGalsInGroups.csv")

g3c$zmax <- NA
for (i in 1:nrow(g3c)) {
    if(g3c$Nfof[i] == 2) {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[2]
    } else {
        g3c$zmax[i] <- sort(gig[GroupID==g3c$GroupID[i], zmax_19p8], 
                           decreasing=TRUE)[multi]
    }
}

g3c$zmax <- ifelse(g3c$zmax < g3c$Zfof, g3c$Zfof, g3c$zmax)
g3c$zmax <- ifelse(g3c$zmax > zlimit, zlimit, g3c$zmax)

vol_zmax <- cosdist(g3c$zmax, OmegaM=omegam, H0=ho)$CoVol
vol_zmin <- cosdist(zmin, OmegaM=omegam, H0=ho)$CoVol

g3c$vmax <- 179.92/(360^2/pi) * 1E9 * (vol_zmax - vol_zmin)

# Clip EXACTLY as in gamahmf.r lines 307-308
g3c$weightszlimit <- ifelse(g3c$vmax > Vsurvey, Vsurvey, g3c$vmax)
g3c$weightszlimit <- ifelse(g3c$vmax < vlimitmin, vlimitmin, g3c$vmax)

# Bin EXACTLY as in gamahmf.r line 316
massx <- seq(10.3, 16.1, logbin)
gamahmf2 <- weighted.hist(log10(g3c$MassAfunc), 
                          w=1/g3c$weightszlimit, 
                          breaks=massx, 
                          plot=FALSE)

gamax <- gamahmf2$mids
gamay <- gamahmf2$counts / logbin

cat("=== OUR BINNED DATA ===\n")
cat("N groups:", nrow(g3c), "\n")
cat("N bins:", length(gamax), "\n")
cat("Bin centers:\n")
for(i in 1:length(gamax)) {
    cat(sprintf("  %.2f: phi=%.4e (counts=%.1f)\n", 
                gamax[i], gamay[i], gamahmf2$counts[i]))
}

cat("\n=== FOR M > 12.7 ===\n")
ok <- gamay > 0 & gamax > 12.7
cat("N bins after cut:", sum(ok), "\n")
cat("Mass range:", range(gamax[ok]), "\n")
cat("Phi range:", range(gamay[ok]), "\n")

cat("\nThese numbers should match gamahmf.r output!\n")
cat("Check: N groups should be ~1733\n")
cat("Check: First bin with data around M~10.8\n")
