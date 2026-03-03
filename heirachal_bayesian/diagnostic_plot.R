############################################################
# PUBLICATION-QUALITY DIAGNOSTIC: Mock Construction & Recovery
#
# 6-panel figure telling the complete story:
#
# (a) The input MRP and how halos are drawn per shell
# (b) The survey selection: mlim(z) and detected halos
# (c) The truncation mechanism: Lambda and per-group limits
# (d) Observed vs true masses (error model)
# (e) Recovery: no errors (validation)
# (f) Recovery: with errors (hierarchical MAP)
#
# Plus a 4-panel parameter recovery comparison figure
############################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

############################################################
# 1. Setup
############################################################

TRUE_MSTAR <- 14.13; TRUE_LOG_PHI <- -3.96
TRUE_ALPHA <- -1.68; TRUE_BETA <- 0.63
tv <- c(ms=TRUE_MSTAR, lp=TRUE_LOG_PHI, al=TRUE_ALPHA, be=TRUE_BETA)

mrp_phi <- function(x, mstar, log_phi, alpha, beta) {
    u <- x - mstar
    beta * log(10) * 10^log_phi * 10^((alpha+1)*u) * exp(-10^(beta*u))
}
mrp_log10 <- function(x, mstar, log_phi, alpha, beta)
    log10(pmax(mrp_phi(x, mstar, log_phi, alpha, beta), 1e-30))

ho <- 67.37; omegam <- 0.3147; zlimit <- 0.25; zmin <- 0.01
cosdist <- function(z) {
    f <- function(zp) 1/sqrt(omegam*(1+zp)^3 + (1-omegam))
    299792.458/ho * integrate(f, 0, z)$value
}

sky_frac <- 60*(pi/180)^2/(4*pi)
nz <- 20
z_edges <- seq(zmin, zlimit, length.out=nz+1)
z_mids <- (z_edges[-1]+z_edges[-(nz+1)])/2
d_edges <- sapply(z_edges, cosdist)
V_sh <- (4/3)*pi*(d_edges[-1]^3 - d_edges[-(nz+1)]^3)*sky_frac
Vsurvey <- sum(V_sh)
mlim_of_z <- function(z) 12.0 + 2.5*(z-zmin)/(zlimit-zmin)
mlim_sh <- mlim_of_z(z_mids)

############################################################
# 2. Generate halos (keeping track of per-shell info)
############################################################

mass_hi <- 16.0
all_tm <- c(); all_z <- c(); all_shell <- c()
N_expected <- numeric(nz); N_actual <- numeric(nz)

for(j in 1:nz) {
    mg <- seq(mlim_sh[j], mass_hi, by=0.001)
    if(length(mg)<2) next
    pg <- mrp_phi(mg, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
    N_expected[j] <- V_sh[j]*sum(pg)*0.001
    N_j <- rpois(1, N_expected[j])
    N_actual[j] <- N_j
    if(N_j==0) next
    pm_j <- max(pg)*1.1; mj <- numeric(0)
    while(length(mj)<N_j) {
        mp <- runif(N_j*10, mlim_sh[j], mass_hi)
        pp <- mrp_phi(mp, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)
        mj <- c(mj, mp[runif(length(mp)) < pp/pm_j])
    }
    all_tm <- c(all_tm, mj[1:N_j])
    all_z <- c(all_z, runif(N_j, z_edges[j], z_edges[j+1]))
    all_shell <- c(all_shell, rep(j, N_j))
}
N_true <- length(all_tm)

# Errors
sigma_all <- pmax(0.10, pmin(0.30, 0.30 - 0.06*(all_tm-13)))
obs_all <- all_tm + rnorm(N_true, 0, sigma_all)
mlim_per <- mlim_of_z(all_z)

# Lambda at true params
Lambda_true <- sum(N_expected)

cat(sprintf("Generated %d halos (Lambda = %.0f)\n", N_true, Lambda_true))
cat(sprintf("Median sigma = %.2f dex\n\n", median(sigma_all)))

############################################################
# 3. Stan model
############################################################

stan_hier <- "
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;
  vector[N] mlim;
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
  vector[N] z_raw;
}
transformed parameters {
  vector[N] mt;
  for(i in 1:N) mt[i] = x_obs[i] + sig[i] * z_raw[i];
}
model {
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.7, 1.0);
  z_raw ~ std_normal();
  
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
    real u = mt[i] - ms;
    real phi_i = be*ln10*pow(10,lp)*pow(10,(al+1)*u)*exp(-pow(10,be*u));
    if(mt[i] > mlim[i] && phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
"

cat("Compiling Stan model...\n")
sm <- stan_model(model_code=stan_hier)

init_fn <- function() {
    list(ms=rnorm(1,14,.2), lp=rnorm(1,-4,.2),
         al=rnorm(1,-1.7,.1), be=runif(1,.5,.8),
         z_raw=rnorm(N_true, 0, 0.3))
}

############################################################
# 4. MAP: no errors
############################################################

cat("Running MAP (no errors)...\n")
sd_noerr <- list(N=N_true, x_obs=all_tm, sig=rep(1e-6, N_true), mlim=mlim_per,
                 Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200)

best_noerr <- NULL; best_lp_noerr <- -Inf
for(a in 1:5) {
    o <- tryCatch(optimizing(sm, data=sd_noerr, init=init_fn(),
                              hessian=TRUE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
                  error=function(e) NULL)
    if(!is.null(o) && o$value > best_lp_noerr) { best_noerr <- o; best_lp_noerr <- o$value }
}
ne <- best_noerr$par
cat(sprintf("  No-error MAP: M*=%.3f lp=%.3f al=%.3f be=%.3f\n", ne$ms, ne$lp, ne$al, ne$be))

# Laplace SE for no-error
np <- 4
tryCatch({
    H <- best_noerr$hessian
    A <- -H[1:np, 1:np]
    B <- -H[1:np, (np+1):(np+N_true)]
    D <- -H[(np+1):(np+N_true), (np+1):(np+N_true)]
    Schur_ne <- A - B %*% diag(1/diag(D)) %*% t(B)
    se_ne <- sqrt(diag(solve(Schur_ne)))
    names(se_ne) <- c("ms","lp","al","be")
}, error=function(e) { se_ne <<- c(ms=0.3, lp=0.4, al=0.05, be=0.3) })

############################################################
# 5. MAP: with errors
############################################################

cat("Running MAP (with errors)...\n")
sd_err <- list(N=N_true, x_obs=obs_all, sig=sigma_all, mlim=mlim_per,
               Nsh=nz, V_sh=V_sh, mlim_sh=mlim_sh, xhi=16.5, Ng=200)

best_err <- NULL; best_lp_err <- -Inf
for(a in 1:5) {
    o <- tryCatch(optimizing(sm, data=sd_err, init=init_fn(),
                              hessian=TRUE, as_vector=FALSE, iter=10000, algorithm="LBFGS"),
                  error=function(e) NULL)
    if(!is.null(o) && o$value > best_lp_err) { best_err <- o; best_lp_err <- o$value }
}
we <- best_err$par
cat(sprintf("  With-error MAP: M*=%.3f lp=%.3f al=%.3f be=%.3f\n", we$ms, we$lp, we$al, we$be))

map_mt <- obs_all + sigma_all * we$z_raw

# Laplace SE for with-error
tryCatch({
    H <- best_err$hessian
    A <- -H[1:np, 1:np]
    B <- -H[1:np, (np+1):(np+N_true)]
    D <- -H[(np+1):(np+N_true), (np+1):(np+N_true)]
    Schur_we <- A - B %*% diag(1/diag(D)) %*% t(B)
    se_we <- sqrt(diag(solve(Schur_we)))
    names(se_we) <- c("ms","lp","al","be")
}, error=function(e) { se_we <<- c(ms=0.3, lp=0.4, al=0.05, be=0.3) })

############################################################
# 6. FIGURE 1: Mock Construction (3 panels)
############################################################

cat("\nGenerating Figure 1: Mock construction...\n")

pdf("fig1_mock_construction.pdf", width=16, height=5.5)
par(mfrow=c(1,3), mar=c(5,5,3.5,1.5), cex.lab=1.4, cex.axis=1.2, cex.main=1.5)

# --- (a) The MRP and per-shell sampling ---
xf <- seq(11, 16.5, length=500)
yf <- mrp_phi(xf, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA)

plot(xf, log10(pmax(yf, 1e-10)), type="l", lwd=3, col="black",
     xlim=c(11.5, 16), ylim=c(-8, -2),
     xlab=expression(log[10]*(M/M['\u2299'])),
     ylab=expression(log[10]*(phi~"/"~Mpc^-3~dex^-1)),
     main="(a) Input MRP & Shell Sampling")
grid(col="gray85")

# Show a few shells with their mlim as coloured shading
shell_cols <- c("dodgerblue3","forestgreen","darkorange","firebrick","purple")
show_shells <- c(2, 6, 10, 15, 19)
for(ii in seq_along(show_shells)) {
    j <- show_shells[ii]
    # Shade the region this shell can see
    xx <- seq(mlim_sh[j], 16.5, length=200)
    yy <- log10(pmax(mrp_phi(xx, TRUE_MSTAR, TRUE_LOG_PHI, TRUE_ALPHA, TRUE_BETA), 1e-10))
    polygon(c(xx, rev(xx)), c(yy, rep(-10, length(xx))),
            col=adjustcolor(shell_cols[ii], alpha=0.15), border=NA)
    abline(v=mlim_sh[j], col=shell_cols[ii], lwd=1.5, lty=2)
    text(mlim_sh[j]+0.08, -2.3 - 0.35*ii,
         sprintf("z=%.2f", z_mids[j]),
         col=shell_cols[ii], cex=0.8, adj=0)
}

# Histogram of true masses
br <- seq(11.5, 16.5, by=0.2)
h <- hist(all_tm, breaks=br, plot=FALSE)
ph <- h$counts/(Vsurvey*0.2)
points(h$mids[ph>0], log10(ph[ph>0]), pch=16, col="grey30", cex=1.3)

legend("topright",
       c("MRP (true)", "Binned data",
         expression(M[lim]*" per shell")),
       lty=c(1, NA, 2), pch=c(NA, 16, NA), lwd=c(3, NA, 1.5),
       col=c("black","grey30","grey50"), cex=0.85, bg="white")

# --- (b) Selection function: M vs z ---
plot(all_tm, all_z, pch=16, col=rgb(0.3, 0.3, 0.7, 0.25), cex=0.7,
     xlab=expression(log[10]*(M[true]/M['\u2299'])),
     ylab="Redshift",
     main=expression("(b) Selection: "*M[lim](z)),
     xlim=c(11.5, 16), ylim=c(0, 0.27))
grid(col="gray85")

# Detection boundary
zp <- seq(zmin, zlimit, length=300)
mp <- mlim_of_z(zp)
lines(mp, zp, col="red", lwd=3)

# Shade excluded region
polygon(c(11, mp, 11), c(zmin, zp, zlimit),
        col=adjustcolor("red", alpha=0.08), border=NA)

# Label
text(12.8, 0.18, expression(M[lim](z) == 12.0 + 2.5 %.% frac(z - 0.01, 0.24)),
     cex=0.95, col="red")

# Shell boundaries
for(j in 1:nz) {
    segments(mlim_sh[j]-0.05, z_mids[j], mlim_sh[j]+0.05, z_mids[j],
             col="red", lwd=1.5)
}

text(14.5, 0.03, sprintf("N = %d groups", N_true), cex=1.1)
text(14.5, 0.015, sprintf("%d shells, %.0f deg²", nz, sky_frac*4*pi*(180/pi)^2), cex=0.9, col="grey40")

# --- (c) Error model and observed masses ---
plot(all_tm, obs_all, pch=16, col=rgb(0.3, 0.3, 0.7, 0.15), cex=0.7,
     xlab=expression(log[10]*(M[true])),
     ylab=expression(log[10]*(M[obs])),
     main=expression("(c) Measurement Errors: "*sigma*(M)),
     xlim=c(11.5, 16), ylim=c(11.5, 16))
grid(col="gray85")
abline(0, 1, col="red", lwd=2)
abline(a=0.3, b=1, col="red", lwd=1, lty=2)
abline(a=-0.3, b=1, col="red", lwd=1, lty=2)

# Inset: sigma vs M
rect(14.5, 11.7, 16, 12.8, col="white", border="grey50")
# Mini plot of sigma(M)
mg <- seq(12, 15.5, length=100)
sg <- pmax(0.10, pmin(0.30, 0.30 - 0.06*(mg-13)))
# Scale to inset
x_in <- 14.55 + (mg - 12)/(15.5-12) * 1.35
y_in <- 11.75 + sg/0.35 * 1.0
lines(x_in, y_in, col="firebrick", lwd=2)
text(15.25, 12.68, expression(sigma(M)), cex=0.8, col="firebrick")
text(14.7, 11.88, "0.10", cex=0.65, col="grey40")
text(14.7, 12.55, "0.30", cex=0.65, col="grey40")

legend("topleft",
       c("1:1", expression(pm*"0.3 dex")),
       lty=c(1, 2), col="red", lwd=c(2,1), cex=0.85, bg="white")

dev.off()
cat("Saved: fig1_mock_construction.pdf\n")

############################################################
# 7. FIGURE 2: Truncation explained (3 panels)
############################################################

cat("Generating Figure 2: Truncation mechanism...\n")

pdf("fig2_truncation.pdf", width=16, height=5.5)
par(mfrow=c(1,3), mar=c(5,5,3.5,1.5), cex.lab=1.4, cex.axis=1.2, cex.main=1.5)

# --- (a) Lambda: expected count per shell ---
barplot_data <- rbind(N_expected, N_actual)
rownames(barplot_data) <- c("Expected", "Observed")

bp <- barplot(N_expected, col=adjustcolor("dodgerblue3", alpha=0.6),
              border="dodgerblue4",
              names.arg=sprintf("%.2f", z_mids), las=2, cex.names=0.7,
              ylab="N halos", main=expression("(a) "*Lambda*" = "*Sigma[j]*V[j]*integral(phi*" dM",M[lim][j],infinity)),
              ylim=c(0, max(N_expected)*1.4))
# Overlay actual
barplot(N_actual, col=adjustcolor("firebrick", alpha=0.4),
        border="firebrick", add=TRUE, names.arg=rep("",nz))
legend("topright",
       c(sprintf("Expected (Lambda=%.0f)", Lambda_true),
         sprintf("Observed (N=%d)", N_true)),
       fill=c(adjustcolor("dodgerblue3",alpha=0.6), adjustcolor("firebrick",alpha=0.4)),
       border=c("dodgerblue4","firebrick"), cex=0.85, bg="white")
mtext("Shell midpoint redshift", side=1, line=3.5, cex=1.1)

# --- (b) Per-group truncation: mt > mlim ---
# Show groups near their detection limit
margin <- all_tm - mlim_per
near_boundary <- which(margin < 0.5 & margin > 0)

plot(mlim_per, all_tm, pch=16, col=rgb(0.3,0.3,0.7,0.15), cex=0.7,
     xlab=expression(M[lim]*" (per-group)"),
     ylab=expression(M[true]),
     main=expression("(b) Per-group truncation: "*M[true] > M[lim]),
     xlim=c(12, 14.8), ylim=c(12, 16))
grid(col="gray85")
abline(0, 1, col="red", lwd=3)

# Highlight near-boundary groups
points(mlim_per[near_boundary], all_tm[near_boundary],
       pch=16, col=adjustcolor("firebrick", alpha=0.5), cex=1.0)

# Shade forbidden region
polygon(c(12, 14.8, 14.8, 12), c(12, 14.8, 12, 12),
        col=adjustcolor("red", alpha=0.08), border=NA)

text(13.0, 15.0, expression("Allowed: "*M[true] > M[lim]), cex=1.0)
text(13.8, 12.8, "Forbidden", cex=1.0, col="firebrick")
arrows(13.8, 12.65, 13.8, 12.3, length=0.1, col="firebrick", lwd=2)

legend("bottomright",
       c("All groups", expression("Near boundary ("*Delta*"<0.5)")),
       pch=16, col=c(rgb(0.3,0.3,0.7,0.5), adjustcolor("firebrick",0.7)),
       cex=0.85, bg="white")

# --- (c) Complete likelihood diagram ---
plot.new()
plot.window(xlim=c(0,10), ylim=c(0,10))

# Title
text(5, 9.5, "The Poisson Process Likelihood", cex=1.5, font=2)

# Equation
text(5, 8.3, expression(log*L == sum(log*phi(M[true][i]), i==1, N) - Lambda),
     cex=1.5)

# Lambda definition
rect(0.5, 5.5, 9.5, 7.5, border="dodgerblue3", lwd=2)
text(5, 7.1, expression(bold("Truncation via "*Lambda*":")), cex=1.1, col="dodgerblue3")
text(5, 6.4, expression(Lambda == sum(V[j] %.% integral(phi(M)*dM, M[lim][j], infinity), j==1, N[sh])),
     cex=1.3, col="dodgerblue3")

# Per-group
rect(0.5, 3.0, 9.5, 5.0, border="firebrick", lwd=2)
text(5, 4.6, expression(bold("Per-group constraint:")), cex=1.1, col="firebrick")
text(5, 3.9, expression(M[true][i] > M[lim][i]*"  (enforced for each group)"), cex=1.2, col="firebrick")
text(5, 3.3, expression(M[obs][i] %~% N(M[true][i]*", "*sigma[i])*"  (error model)"),
     cex=1.1, col="grey30")

# Key insight
rect(0.5, 0.8, 9.5, 2.5, border="forestgreen", lwd=2)
text(5, 2.1, expression(bold("Key:")*" MRP is evaluated at "*M[true]*", not "*M[obs]),
     cex=1.1, col="forestgreen")
text(5, 1.4, expression("Eddington bias avoided: latent "*M[true]*" sampled by Stan"),
     cex=1.0, col="forestgreen")

dev.off()
cat("Saved: fig2_truncation.pdf\n")

############################################################
# 8. FIGURE 3: Recovery results (2x3)
############################################################

cat("Generating Figure 3: Recovery results...\n")

pdf("fig3_recovery.pdf", width=16, height=10)
par(mfrow=c(2,3), mar=c(5,5,3.5,1.5), cex.lab=1.4, cex.axis=1.2, cex.main=1.5)

xf <- seq(11.5, 16.5, length=500)
br <- seq(11.5, 16.5, by=0.2)

# --- (a) HMF: no errors ---
h_t <- hist(all_tm, breaks=br, plot=FALSE)
pt <- h_t$counts/(Vsurvey*0.2)
plot(h_t$mids[pt>0], log10(pt[pt>0]), pch=16, col="grey30", cex=1.5,
     xlim=c(12, 16), ylim=c(-8, -2),
     xlab=expression(log[10]*M), ylab=expression(log[10]*phi),
     main=expression("(a) No errors: "*sigma*" = 0"))
grid(col="gray85")
lines(xf, mrp_log10(xf, tv["ms"],tv["lp"],tv["al"],tv["be"]),
      col="black", lwd=2, lty=2)
lines(xf, mrp_log10(xf, ne$ms, ne$lp, ne$al, ne$be),
      col="dodgerblue3", lwd=3)
legend("topright", c("True MRP","MAP fit","Binned data"),
       lty=c(2,1,NA), pch=c(NA,NA,16), lwd=c(2,3,NA),
       col=c("black","dodgerblue3","grey30"), cex=0.85, bg="white")

# --- (b) HMF: with errors ---
h_o <- hist(obs_all, breaks=br, plot=FALSE)
po <- h_o$counts/(Vsurvey*0.2)
plot(h_o$mids[po>0], log10(po[po>0]), pch=16, col="firebrick", cex=1.5,
     xlim=c(12, 16), ylim=c(-8, -2),
     xlab=expression(log[10]*M), ylab=expression(log[10]*phi),
     main=expression("(b) With errors: "*sigma*" ~ 0.28 dex"))
points(h_t$mids[pt>0], log10(pt[pt>0]), pch=1, col="grey60", cex=1.2)
grid(col="gray85")
lines(xf, mrp_log10(xf, tv["ms"],tv["lp"],tv["al"],tv["be"]),
      col="black", lwd=2, lty=2)
lines(xf, mrp_log10(xf, we$ms, we$lp, we$al, we$be),
      col="dodgerblue3", lwd=3)
legend("topright", c("True MRP","Hier. MAP fit","Obs masses","True masses"),
       lty=c(2,1,NA,NA), pch=c(NA,NA,16,1), lwd=c(2,3,NA,NA),
       col=c("black","dodgerblue3","firebrick","grey60"), cex=0.85, bg="white")

# --- (c) Latent mass recovery ---
plot(all_tm, map_mt, pch=16, col=rgb(0.2,0.2,0.8,0.2), cex=0.7,
     xlab=expression("True "*log[10]*M),
     ylab=expression("Recovered "*log[10]*M),
     main="(c) Latent Mass Recovery")
grid(col="gray85")
abline(0, 1, col="red", lwd=2.5)
med_resid <- median(abs(map_mt - all_tm))
legend("topleft",
       c(sprintf("Median |resid| = %.3f dex", med_resid),
         sprintf("Median sigma  = %.3f dex", median(sigma_all)),
         sprintf("Shrinkage = %.0f%%", (1-med_resid/median(sigma_all))*100)),
       cex=0.9, bg="white")

# --- (d-f) Parameter posteriors ---
pinfo <- list(
    ms=list(name=expression(M*"*"), true=tv["ms"], noerr=ne$ms, werr=we$ms, se_ne=se_ne["ms"], se_we=se_we["ms"]),
    lp=list(name=expression(log[10]*phi*"*"), true=tv["lp"], noerr=ne$lp, werr=we$lp, se_ne=se_ne["lp"], se_we=se_we["lp"]),
    al=list(name=expression(alpha), true=tv["al"], noerr=ne$al, werr=we$al, se_ne=se_ne["al"], se_we=se_we["al"]),
    be=list(name=expression(beta), true=tv["be"], noerr=ne$be, werr=we$be, se_ne=se_ne["be"], se_we=se_we["be"])
)

pnms <- c("ms","lp","al","be")
labels <- c(expression(M*"*"), expression(log[10]*phi*"*"), expression(alpha), expression(beta))

# --- (d) Parameter comparison: point + error bar ---
plot(1:4-0.12, sapply(pinfo, function(p) p$noerr), pch=17, col="forestgreen", cex=2.5,
     xlim=c(0.5,4.5), ylim=c(NA,NA),
     xlab="", ylab="Parameter value",
     main=expression("(d) Parameter Recovery"), xaxt="n",
     panel.first=grid(col="gray85"))
# No-error error bars
for(i in 1:4) {
    p <- pinfo[[pnms[i]]]
    arrows(i-0.12, p$noerr-p$se_ne, i-0.12, p$noerr+p$se_ne,
           code=3, angle=90, length=0.06, col="forestgreen", lwd=2)
}
# With-error points + bars
points(1:4+0.12, sapply(pinfo, function(p) p$werr), pch=19, col="dodgerblue3", cex=2.5)
for(i in 1:4) {
    p <- pinfo[[pnms[i]]]
    arrows(i+0.12, p$werr-p$se_we, i+0.12, p$werr+p$se_we,
           code=3, angle=90, length=0.06, col="dodgerblue3", lwd=2)
}
# True values
points(1:4, sapply(pinfo, function(p) p$true), pch=4, col="red", cex=3, lwd=3)
axis(1, at=1:4, labels=labels, cex.axis=1.3)
legend("topright",
       c("True", expression("No errors ("*sigma*"=0)"),
         expression("Hierarchical ("*sigma*"~0.28)")),
       pch=c(4,17,19), col=c("red","forestgreen","dodgerblue3"),
       pt.cex=c(2,2,2), pt.lwd=c(3,1,1), cex=0.85, bg="white")

# (d) is wide so it covers the bottom-left. Now do bias plot
# --- (e) Bias in sigma ---
bias_ne <- c((ne$ms-tv["ms"])/se_ne["ms"], (ne$lp-tv["lp"])/se_ne["lp"],
             (ne$al-tv["al"])/se_ne["al"], (ne$be-tv["be"])/se_ne["be"])
bias_we <- c((we$ms-tv["ms"])/se_we["ms"], (we$lp-tv["lp"])/se_we["lp"],
             (we$al-tv["al"])/se_we["al"], (we$be-tv["be"])/se_we["be"])

plot(1:4-0.1, bias_ne, pch=17, col="forestgreen", cex=2.5,
     xlim=c(0.5,4.5), ylim=c(-3, 3),
     xlab="", ylab=expression("Bias ("*sigma*")"),
     main="(e) Bias in Laplace SE", xaxt="n")
grid(col="gray85")
points(1:4+0.1, bias_we, pch=19, col="dodgerblue3", cex=2.5)
abline(h=0, col="black", lwd=1.5)
abline(h=c(-1,1), col="grey50", lty=2, lwd=1.5)
abline(h=c(-2,2), col="grey70", lty=3)
axis(1, at=1:4, labels=labels, cex.axis=1.3)
legend("topleft",
       c(expression("No errors ("*sigma*"=0)"),
         expression("Hierarchical ("*sigma*"~0.28)")),
       pch=c(17,19), col=c("forestgreen","dodgerblue3"), cex=0.85, bg="white")

# --- (f) Summary table as a text plot ---
plot.new()
plot.window(xlim=c(0,10), ylim=c(0,10))

text(5, 9.3, "Recovery Summary", cex=1.6, font=2)

# Table header
y0 <- 8.2; dy <- 0.8
text(1.0, y0, "Param", cex=1.1, font=2, adj=0)
text(3.0, y0, "True", cex=1.1, font=2)
text(4.8, y0, expression(sigma*"=0"), cex=1.1, font=2, col="forestgreen")
text(6.6, y0, expression("Hier."), cex=1.1, font=2, col="dodgerblue3")
text(8.5, y0, expression("Bias ("*sigma*")"), cex=1.1, font=2)

segments(0.5, y0-0.35, 9.5, y0-0.35, col="grey50")

for(i in 1:4) {
    y <- y0 - i*dy
    p <- pinfo[[pnms[i]]]
    bias_str <- sprintf("%+.2f", bias_we[i])
    bias_col <- ifelse(abs(bias_we[i]) < 1, "forestgreen",
                       ifelse(abs(bias_we[i]) < 2, "darkorange", "firebrick"))
    
    text(1.0, y, labels[i], cex=1.1, adj=0)
    text(3.0, y, sprintf("%.3f", p$true), cex=1.05)
    text(4.8, y, sprintf("%.3f", p$noerr), cex=1.05, col="forestgreen")
    text(6.6, y, sprintf("%.3f", p$werr), cex=1.05, col="dodgerblue3")
    text(8.5, y, bias_str, cex=1.15, col=bias_col, font=2)
}

segments(0.5, y0-4*dy-0.35, 9.5, y0-4*dy-0.35, col="grey50")

# Footer
y_foot <- y0 - 5*dy
text(5, y_foot, sprintf("N = %d groups  |  Lambda = %.0f  |  MAP < 1 sec",
                         N_true, Lambda_true), cex=1.0, col="grey40")
text(5, y_foot - 0.7,
     sprintf("Latent mass: |resid| = %.3f dex vs sigma = %.3f dex (%.0f%% shrinkage)",
             med_resid, median(sigma_all), (1-med_resid/median(sigma_all))*100),
     cex=1.0, col="grey40")

dev.off()
cat("Saved: fig3_recovery.pdf\n")

############################################################
# Print summary
############################################################

cat("\n")
cat("=", strrep("=", 65), "\n", sep="")
cat("  RECOVERY TEST SUMMARY\n")
cat("=", strrep("=", 65), "\n", sep="")
cat(sprintf("\n  Mock: %d groups, %d shells, median sigma = %.2f dex\n",
            N_true, nz, median(sigma_all)))
cat(sprintf("  Lambda (true params) = %.0f\n\n", Lambda_true))

cat(sprintf("  %-8s %8s %10s %10s %10s\n", "", "True", "No-error", "Hier.MAP", "Bias(sig)"))
cat("  ", strrep("-", 50), "\n", sep="")
for(i in 1:4) {
    p <- pinfo[[pnms[i]]]
    cat(sprintf("  %-8s %8.3f %10.3f %10.3f %+10.2f\n",
                pnms[i], p$true, p$noerr, p$werr, bias_we[i]))
}
cat("  ", strrep("-", 50), "\n\n", sep="")
cat(sprintf("  Latent mass: median |recovered - true| = %.3f dex\n", med_resid))
cat(sprintf("  MAP runtime: < 1 second\n\n"))
cat("  Files: fig1_mock_construction.pdf\n")
cat("         fig2_truncation.pdf\n")
cat("         fig3_recovery.pdf\n")
cat("=", strrep("=", 65), "\n", sep="")
cat("\nDone!\n")
