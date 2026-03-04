############################################################
# SHARK RECOVERY: Which envelope quantile works best?
#
# Test mlim(z) derived from different quantiles of the
# detected mass distribution (5th, 10th, 25th, 50th)
# and also the true 50% completeness, to see which
# envelope prescription recovers the true MRP best.
#
# This tells us what to use on real data where we
# can't compute true completeness.
############################################################

library(arrow)
library(Cairo)
library(celestial)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup (same as before)
############################################################

TRUE_MS <- 13.51; TRUE_LP <- -3.19; TRUE_AL <- -1.27; TRUE_BE <- 0.47

mrp_phi <- function(x, ms, lp, al, be) {
    be * log(10) * 10^lp * 10^((al+1)*(x-ms)) * exp(-10^(be*(x-ms)))
}

cat("============================================\n")
cat("  QUANTILE COMPARISON FOR mlim(z)\n")
cat("============================================\n\n")

data_dir <- "/Users/00115372/Desktop/masking_mock_cat"

cat("Reading data...\n")
groups <- as.data.table(read_parquet(file.path(data_dir, "groups_shark.parquet")))
groups <- groups[dec > 0]
galaxies <- as.data.table(read_parquet(file.path(data_dir, "galaxies_shark.parquet")))
galaxies <- galaxies[dec > 0]

zmin <- 0.01; zlimit <- 0.25; multi <- 5; mag_limit <- 19.8

############################################################
# 2. Abundance matching
############################################################

cat("Abundance matching...\n")

groups_vol <- groups[redshift_cosmological > zmin & 
                     redshift_cosmological < zlimit & mass_virial > 0]
N_total <- nrow(groups_vol)

ho <- 67.37; omegam <- 0.3147
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

ra_range  <- range(groups_vol$ra, na.rm=TRUE)
dec_range <- range(groups_vol$dec, na.rm=TRUE)
sky_area_deg2 <- skyarea(ra_range, dec_range)["area"]
sky_frac <- sky_area_deg2 * (pi/180)^2 / (4*pi)
d_lo <- cosdist(zmin); d_hi <- cosdist(zlimit)
Vsurvey <- (4/3) * pi * (d_hi^3 - d_lo^3) * sky_frac

m_grid <- seq(9, 16.5, by=0.001)
phi_grid <- mrp_phi(m_grid, TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)
cum_counts <- rev(cumsum(rev(phi_grid * 0.001))) * Vsurvey

groups_vol[, rank := frankv(mass_virial, order=-1L)]
groups_vol[, log_mass_am := approx(x=rev(cum_counts), y=rev(m_grid),
                                    xout=rank, rule=2)$y]

############################################################
# 3. GAMA selection
############################################################

cat("GAMA selection...\n")

gal_selected <- galaxies[id_fof != -1 & mass_stellar_total > 1e8 & mag_r_SDSS < mag_limit]
gal_counts <- gal_selected[, .(n_gama = .N), by = id_group_sky]

groups_vol <- merge(groups_vol, gal_counts, by="id_group_sky", all.x=TRUE)
groups_vol[is.na(n_gama), n_gama := 0L]
groups_vol[, detected := n_gama >= multi]

groups_gama <- groups_vol[detected == TRUE]
z_gama <- groups_gama$redshift_cosmological
m_gama <- groups_gama$log_mass_am

cat(sprintf("  Detected: %d / %d\n", nrow(groups_gama), N_total))

############################################################
# 4. Compute true completeness (for reference)
############################################################

cat("Computing true completeness...\n")

m_bin_edges <- seq(10, 16, by=0.2)
m_bin_mids  <- (m_bin_edges[-1] + m_bin_edges[-length(m_bin_edges)]) / 2
n_m <- length(m_bin_mids)

z_bin_edges_fine <- seq(zmin, zlimit, length.out=21)
z_bin_mids_fine  <- (z_bin_edges_fine[-1] + z_bin_edges_fine[-length(z_bin_edges_fine)]) / 2
n_z_fine <- length(z_bin_mids_fine)

C_mz <- matrix(NA, n_m, n_z_fine)
for(im in 1:n_m) {
    for(iz in 1:n_z_fine) {
        in_bin <- groups_vol$log_mass_am >= m_bin_edges[im] & 
                  groups_vol$log_mass_am <  m_bin_edges[im+1] &
                  groups_vol$redshift_cosmological >= z_bin_edges_fine[iz] &
                  groups_vol$redshift_cosmological <  z_bin_edges_fine[iz+1]
        nt <- sum(in_bin)
        nd <- sum(in_bin & groups_vol$detected)
        if(nt >= 10) C_mz[im, iz] <- nd / nt
    }
}

# True 50% completeness mlim
mlim_true50 <- numeric(n_z_fine)
for(iz in 1:n_z_fine) {
    comp_col <- C_mz[, iz]
    ok <- !is.na(comp_col)
    if(sum(ok) >= 3) {
        m_ok <- m_bin_mids[ok]; c_ok <- comp_col[ok]
        above_50 <- which(c_ok >= 0.5)
        if(length(above_50) > 0 && min(above_50) > 1) {
            ic <- min(above_50)
            mlim_true50[iz] <- m_ok[ic-1] + (0.5 - c_ok[ic-1]) * 
                               (m_ok[ic] - m_ok[ic-1]) / (c_ok[ic] - c_ok[ic-1])
        } else if(length(above_50) > 0) {
            mlim_true50[iz] <- m_ok[1]
        } else {
            mlim_true50[iz] <- NA
        }
    } else {
        mlim_true50[iz] <- NA
    }
}

############################################################
# 5. Compute envelope-based mlim(z) at different quantiles
############################################################

cat("Computing envelope mlim(z) at different quantiles...\n")

quantiles_to_test <- c(0.05, 0.10, 0.25, 0.50)
nbin_z <- 30
z_be <- seq(zmin, zlimit, length.out=nbin_z+1)
z_bm <- (z_be[-1] + z_be[-(nbin_z+1)]) / 2

mlim_envelope <- list()
mlim_funcs    <- list()

for(q in quantiles_to_test) {
    qname <- sprintf("q%02d", round(q*100))
    bins <- numeric(nbin_z)
    for(b in 1:nbin_z) {
        in_bin <- z_gama >= z_be[b] & z_gama < z_be[b+1]
        if(sum(in_bin) > 10) {
            bins[b] <- quantile(m_gama[in_bin], q)
        } else {
            bins[b] <- NA
        }
    }
    mlim_envelope[[qname]] <- bins
    
    # Fit smooth function
    ok_q <- !is.na(bins)
    fit_l <- lm(bins[ok_q] ~ z_bm[ok_q])
    fit_q <- lm(bins[ok_q] ~ z_bm[ok_q] + I(z_bm[ok_q]^2))
    
    if(AIC(fit_q) < AIC(fit_l) - 2) {
        cc <- coef(fit_q)
        mlim_funcs[[qname]] <- function(z) cc[1] + cc[2]*z + cc[3]*z^2
        environment(mlim_funcs[[qname]]) <- list2env(list(cc=cc))
        cat(sprintf("  %s: mlim(z) = %.2f + %.2f*z + %.2f*z^2  [mlim(0.25)=%.2f]\n",
                    qname, cc[1], cc[2], cc[3], cc[1]+cc[2]*0.25+cc[3]*0.25^2))
    } else {
        cc <- coef(fit_l)
        mlim_funcs[[qname]] <- function(z) cc[1] + cc[2]*z
        environment(mlim_funcs[[qname]]) <- list2env(list(cc=cc))
        cat(sprintf("  %s: mlim(z) = %.2f + %.2f*z  [mlim(0.25)=%.2f]\n",
                    qname, cc[1], cc[2], cc[1]+cc[2]*0.25))
    }
}

############################################################
# 6. Shells (shared across all runs)
############################################################

nsh <- 20
z_edges <- seq(zmin, zlimit, length.out=nsh+1)
z_mids  <- (z_edges[-1] + z_edges[-(nsh+1)]) / 2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3) * pi * (d_edges[-1]^3 - d_edges[-(nsh+1)]^3) * sky_frac

############################################################
# 7. Stan model
############################################################

stan_poisson <- "
data {
  int<lower=1> N;
  vector[N] x_obs;
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  real xhi;
  int<lower=2> Ng;
}
transformed data {
  real ln10 = log(10.0);
  real xlo = min(mlim_sh) - 0.5;
  real dx = (xhi - xlo) / (Ng - 1.0);
  vector[Ng] xg;
  for(k in 1:Ng) xg[k] = xlo + (k-1)*dx;
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
}
model {
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.7, 1.0);
  
  vector[Ng] pg;
  for(k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
  }
  vector[Ng] cum;
  cum[Ng] = 0;
  for(kr in 1:(Ng-1)) {
    int k = Ng - kr;
    cum[k] = cum[k+1] + 0.5*(pg[k]+pg[k+1])*dx;
  }
  
  real Lambda = 0;
  for(j in 1:Nsh) {
    int k0 = 1;
    for(k in 1:Ng) if(xg[k] <= mlim_sh[j]) k0 = k;
    if(k0 >= Ng) k0 = Ng-1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;
  
  for(i in 1:N) {
    real u = x_obs[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

cat("\nCompiling Stan model...\n")
sm <- stan_model(model_code=stan_poisson)

############################################################
# 8. Run MAP for each quantile
############################################################

results <- list()
tv <- c(TRUE_MS, TRUE_LP, TRUE_AL, TRUE_BE)

for(qname in names(mlim_funcs)) {
    cat(sprintf("\n========== %s ==========\n", qname))
    
    this_mlim <- mlim_funcs[[qname]]
    mlim_sh <- this_mlim(z_mids)
    mlim_per <- this_mlim(z_gama)
    above <- m_gama > mlim_per
    x_obs <- m_gama[above]
    N <- length(x_obs)
    
    cat(sprintf("  mlim_sh range: %.2f -- %.2f\n", min(mlim_sh), max(mlim_sh)))
    cat(sprintf("  N above mlim: %d / %d (%.1f%%)\n", N, length(m_gama), 100*N/length(m_gama)))
    
    stan_data <- list(N=N, x_obs=x_obs, Nsh=nsh, V_sh=V_sh,
                      mlim_sh=mlim_sh, xhi=16.5, Ng=300)
    
    init_fn <- function() {
        list(ms=rnorm(1,14,0.3), lp=rnorm(1,-3.5,0.3),
             al=rnorm(1,-1.3,0.2), be=runif(1,0.3,0.7))
    }
    
    best_opt <- NULL; best_lp_val <- -Inf
    for(trial in 1:5) {
        this_opt <- tryCatch(
            optimizing(sm, data=stan_data, init=init_fn(),
                       hessian=FALSE, as_vector=FALSE,
                       iter=20000, algorithm="LBFGS"),
            error=function(e) NULL)
        if(!is.null(this_opt) && this_opt$value > best_lp_val) {
            best_lp_val <- this_opt$value; best_opt <- this_opt
        }
    }
    
    # Hessian run
    opt <- tryCatch(
        optimizing(sm, data=stan_data,
                   init=list(ms=best_opt$par$ms, lp=best_opt$par$lp,
                             al=best_opt$par$al, be=best_opt$par$be),
                   hessian=TRUE, as_vector=FALSE, iter=20000, algorithm="LBFGS"),
        error=function(e) best_opt)
    
    map_par <- c(opt$par$ms, opt$par$lp, opt$par$al, opt$par$be)
    
    se <- rep(NA, 4)
    if(!is.null(opt$hessian)) {
        tryCatch({
            se <- sqrt(diag(solve(-opt$hessian)))
        }, error=function(e) {})
    }
    
    bias <- (map_par - tv) / se
    
    results[[qname]] <- list(map=map_par, se=se, bias=bias, N=N,
                              mlim_func=this_mlim)
    
    cat(sprintf("  MAP: M*=%.3f lp=%.3f al=%.3f be=%.3f\n",
                map_par[1], map_par[2], map_par[3], map_par[4]))
    cat(sprintf("  Bias(sig): ms=%+.2f lp=%+.2f al=%+.2f be=%+.2f\n",
                bias[1], bias[2], bias[3], bias[4]))
}

############################################################
# 9. Summary table
############################################################

cat("\n============================================\n")
cat("  COMPARISON: Which quantile works best?\n")
cat("============================================\n\n")

cat(sprintf("  TRUE: M*=%.3f  lp=%.3f  al=%.3f  be=%.3f\n\n", tv[1], tv[2], tv[3], tv[4]))

plab <- c("M*","log_phi*","alpha","beta")
cat(sprintf("  %-8s  %5s  %8s  %8s  %8s  %8s  %6s  %6s  %6s  %6s\n",
            "Quantile", "N", "M*", "lp", "al", "be", 
            "dM*", "dlp", "dal", "dbe"))

for(qname in names(results)) {
    r <- results[[qname]]
    cat(sprintf("  %-8s  %5d  %8.3f  %8.3f  %8.3f  %8.3f  %+6.2f  %+6.2f  %+6.2f  %+6.2f\n",
                qname, r$N,
                r$map[1], r$map[2], r$map[3], r$map[4],
                r$bias[1], r$bias[2], r$bias[3], r$bias[4]))
}

############################################################
# 10. Plot
############################################################

xfit <- seq(10, 16, length.out=500)
z_plot <- seq(zmin, zlimit, length.out=200)

bin_width <- 0.2
breaks <- seq(9, 16, by=bin_width)
h_am_true <- hist(groups_vol$log_mass_am, breaks=breaks, plot=FALSE)
phi_am_true <- h_am_true$counts / (Vsurvey * bin_width)
ok_am_true  <- h_am_true$counts >= 5

h_gama_plot <- hist(m_gama, breaks=breaks, plot=FALSE)
phi_gama_plot <- h_gama_plot$counts / (Vsurvey * bin_width)
ok_gama_plot  <- h_gama_plot$counts >= 5

qcols <- c(q05="red", q10="orange", q25="darkgreen", q50="purple")

CairoPDF("shark_recovery.pdf", 14, 14)
par(mfrow=c(3,2))

# ---- Panel 1: Completeness heatmap with all mlim(z) curves ----
image(z_bin_mids_fine, m_bin_mids, t(C_mz),
      col=colorRampPalette(c("white","yellow","orange","red","darkred"))(100),
      xlab="Redshift", ylab=expression("log"[10]*"(M)"),
      main="True completeness + envelope mlim(z)",
      xlim=c(zmin, zlimit), ylim=c(11, 15.5), zlim=c(0, 1))
contour(z_bin_mids_fine, m_bin_mids, t(C_mz),
        levels=c(0.1, 0.25, 0.5, 0.75, 0.9),
        add=TRUE, col="black", labcex=0.7)

for(qname in names(mlim_funcs)) {
    lines(z_plot, mlim_funcs[[qname]](z_plot), col=qcols[qname], lwd=2.5)
}

# True 50% line
ok50 <- !is.na(mlim_true50)
fit50_l <- lm(mlim_true50[ok50] ~ z_bin_mids_fine[ok50])
fit50_q <- lm(mlim_true50[ok50] ~ z_bin_mids_fine[ok50] + I(z_bin_mids_fine[ok50]^2))
if(AIC(fit50_q) < AIC(fit50_l) - 2) {
    cc50 <- coef(fit50_q)
    lines(z_plot, cc50[1]+cc50[2]*z_plot+cc50[3]*z_plot^2, col="blue", lwd=3, lty=2)
} else {
    cc50 <- coef(fit50_l)
    lines(z_plot, cc50[1]+cc50[2]*z_plot, col="blue", lwd=3, lty=2)
}

legend("topleft",
       legend=c("5th pctile", "10th pctile", "25th pctile", 
                "50th pctile (median)", "True 50% completeness"),
       col=c("red","orange","darkgreen","purple","blue"),
       lty=c(1,1,1,1,2), lwd=c(2.5,2.5,2.5,2.5,3),
       bg="white", cex=0.6)

# ---- Panel 2: Completeness vs mass slices ----
z_slices <- c(0.05, 0.10, 0.15, 0.20)
cols_z <- c("blue","darkgreen","orange","red")

plot(NA, xlim=c(11, 15.5), ylim=c(0, 1.05),
     xlab=expression("log"[10]*"(M)"),
     ylab="Completeness", main="Completeness vs mass")
grid(col="gray80")
abline(h=0.5, col="grey50", lty=3)

for(s in 1:length(z_slices)) {
    iz <- which.min(abs(z_bin_mids_fine - z_slices[s]))
    comp_slice <- C_mz[, iz]
    ok_s <- !is.na(comp_slice)
    lines(m_bin_mids[ok_s], comp_slice[ok_s], col=cols_z[s], lwd=2, type="b", pch=19, cex=0.8)
}
legend("bottomright", sprintf("z = %.2f", z_slices), col=cols_z, lwd=2, cex=0.7)

# ---- Panel 3: HMF with all recovered curves ----
plot(h_am_true$mids[ok_am_true], log10(phi_am_true[ok_am_true]),
     pch=19, col="grey40", cex=1.3,
     xlim=c(10.5, 15.5), ylim=c(-8, -1),
     xlab=expression("log"[10]*"(M / M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")"),
     main="HMF recovery by quantile")
grid(col="gray80")

points(h_gama_plot$mids[ok_gama_plot], log10(phi_gama_plot[ok_gama_plot]),
       pch=17, col="grey60", cex=1.0)

# True MRP
lines(xfit, log10(pmax(mrp_phi(xfit, tv[1], tv[2], tv[3], tv[4]), 1e-30)),
      col="blue", lwd=3)

# Each quantile's MAP recovery
for(qname in names(results)) {
    r <- results[[qname]]$map
    lines(xfit, log10(pmax(mrp_phi(xfit, r[1], r[2], r[3], r[4]), 1e-30)),
          col=qcols[qname], lwd=2, lty=2)
}

legend("topright",
       legend=c("True HMF (AM)", "GAMA-selected", "True MRP",
                "q05 MAP", "q10 MAP", "q25 MAP", "q50 MAP"),
       col=c("grey40","grey60","blue","red","orange","darkgreen","purple"),
       pch=c(19,17,NA,NA,NA,NA,NA), lty=c(NA,NA,1,2,2,2,2),
       lwd=c(NA,NA,3,2,2,2,2), bg="white", cex=0.55)

# ---- Panel 4: Bias summary ----
bias_mat <- sapply(results, function(r) r$bias)

plot(NA, xlim=c(0.5, 4.5), ylim=c(-8, 8),
     xlab="", ylab="Bias (sigma)", main="Parameter bias by quantile", xaxt="n")
grid(col="gray80")
abline(h=0, col="black", lwd=2)
abline(h=c(-2,2), col="grey50", lty=3)

for(iq in 1:length(results)) {
    qname <- names(results)[iq]
    offset <- (iq - 2.5) * 0.15
    points(1:4 + offset, results[[qname]]$bias, pch=19, col=qcols[qname], cex=1.5)
    for(j in 1:4) {
        lines(c(j+offset, j+offset), c(0, results[[qname]]$bias[j]),
              col=qcols[qname], lwd=2)
    }
}

axis(1, at=1:4, labels=plab)
legend("topright", names(results), col=unname(qcols[names(results)]),
       pch=19, cex=0.7)

# ---- Panel 5: N groups used per quantile ----
n_used <- sapply(results, function(r) r$N)
barplot(n_used, col=unname(qcols[names(results)]),
        main="N groups used", ylab="N", names.arg=names(results))

# ---- Panel 6: Mass-redshift with mlim curves ----
smoothScatter(z_gama, m_gama,
              xlab="Redshift", ylab=expression("log"[10]*"(M)"),
              main="GAMA-selected + mlim(z) curves",
              xlim=c(zmin, zlimit), ylim=c(11, 15.5),
              colramp=colorRampPalette(c("white","cornflowerblue","blue","darkblue")),
              nbin=200)
for(qname in names(mlim_funcs)) {
    lines(z_plot, mlim_funcs[[qname]](z_plot), col=qcols[qname], lwd=2.5)
}
grid(col="gray80")
legend("bottomright", names(mlim_funcs),
       col=unname(qcols[names(mlim_funcs)]), lwd=2.5, cex=0.7)

dev.off()
cat("\nPlot saved: shark_recovery.pdf\n\nDone!\n")
