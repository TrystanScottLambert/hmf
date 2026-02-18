############################################################
# STAGE 1 - UNBINNED POISSON POINT PROCESS
# Using Metropolis-Hastings in R directly
#
# Abandoning the JAGS zeros trick - it's fundamentally
# unstable for this problem because logL ~ -3000 to -5000
# which makes the Poisson rate in dpois() numerically
# impossible to handle.
#
# Instead: pure R MCMC with full numerical control.
# This is the foundation for the hierarchical model.
############################################################

library(celestial)
library(Rfits)
library(Cairo)
library(coda)

set.seed(42)

############################################################
# Data
############################################################

ho     <- 67.37
omegam <- 0.3147
zlimit <- 0.25

vol_max <- cosdistCoDist(zlimit, OmegaM=omegam, H0=ho)
Vsurvey <- (4/3)*pi*vol_max^3 * 179.92*(pi/180)^2 / (4*pi)

cat("Survey volume:", signif(Vsurvey,4), "Mpc^3\n")

g3cx  <- Rfits_read_table("../data/G3CFoFGroupv10.fits")
g3c   <- g3cx[g3cx$Nfof >= 5 & g3cx$Zfof < zlimit &
                  g3cx$Zfof > 0.01 & g3cx$IterCenDec > -3.5, ]
g3c$MassA <- g3c$MassA * 100 / ho
x_obs <- log10(g3c$MassA[is.finite(g3c$MassA) & g3c$MassA > 0])
x_obs <- sort(x_obs[x_obs > 12.7])
N     <- length(x_obs)

xlo <- 12.7
xhi <- 15.8

cat("N groups:", N, "\n")
cat("Mass range:", round(range(x_obs),3), "\n\n")

############################################################
# Integration grid (fixed)
############################################################

Ng    <- 1000
xgrid <- seq(xlo, xhi, length.out=Ng)
dx    <- xgrid[2] - xgrid[1]

############################################################
# MRP function and log-likelihood
############################################################

# MRP: phi(x) = beta * ln(10) * phi* * 10^((a+1)(x-m*)) * exp(-10^(b(x-m*)))
# where x = log10(M)

log_mrp <- function(x, mstar, log_phi, alpha, beta) {
    # Returns log(phi(x)) in natural log
    log(beta) + log(log(10)) +
        log_phi * log(10) +
        (alpha + 1) * (x - mstar) * log(10) -
        10^(beta * (x - mstar))
}

mrp <- function(x, mstar, log_phi, alpha, beta) {
    exp(log_mrp(x, mstar, log_phi, alpha, beta))
}

log_likelihood <- function(mstar, log_phi, alpha, beta) {

    # Guard against degenerate parameters
    if(beta <= 0) return(-Inf)
    if(is.nan(beta) || is.nan(mstar) || is.nan(log_phi) || is.nan(alpha)) return(-Inf)

    # Sum of log phi over observed groups
    lp <- log_mrp(x_obs, mstar, log_phi, alpha, beta)

    # Check for numerical issues
    if(any(!is.finite(lp))) return(-Inf)

    sum_log_lam <- sum(lp) + N * log(Vsurvey)

    # Integral: Lambda = V * integral[phi(x) dx]
    phi_g <- mrp(xgrid, mstar, log_phi, alpha, beta)
    if(any(!is.finite(phi_g))) return(-Inf)

    Lambda <- Vsurvey * sum(phi_g) * dx

    # Poisson point process log-likelihood
    sum_log_lam - Lambda
}

############################################################
# Log prior
############################################################

log_prior <- function(mstar, log_phi, alpha, beta) {

    if(beta < 0.1 || beta > 2.0) return(-Inf)

    # dnorm in log space
    lp_mstar   <- dnorm(mstar,   13.5, 1.0, log=TRUE)
    lp_log_phi <- dnorm(log_phi, -3.5, 1.5, log=TRUE)
    lp_alpha   <- dnorm(alpha,   -1.5, 0.8, log=TRUE)
    lp_beta    <- 0  # uniform on [0.1, 2.0]

    lp_mstar + lp_log_phi + lp_alpha + lp_beta
}

############################################################
# Log posterior
############################################################

log_posterior <- function(params) {
    mstar   <- params[1]
    log_phi <- params[2]
    alpha   <- params[3]
    beta    <- params[4]

    lprior <- log_prior(mstar, log_phi, alpha, beta)
    if(!is.finite(lprior)) return(-Inf)

    llik <- log_likelihood(mstar, log_phi, alpha, beta)
    if(!is.finite(llik)) return(-Inf)

    lprior + llik
}

############################################################
# Test log-posterior at expected values
############################################################

test_params <- c(13.5, -3.5, -1.5, 0.5)
test_lp     <- log_posterior(test_params)
cat("Log-posterior at test params:", round(test_lp, 2), "\n")

if(!is.finite(test_lp)) {
    stop("Log-posterior is not finite at test parameters! Check data and volume.")
}

############################################################
# Adaptive Metropolis-Hastings MCMC
############################################################

run_chain <- function(init, n_iter, n_burnin, proposal_sd,
                      adapt_interval=200, target_rate=0.25,
                      chain_id=1) {

    cat("\nChain", chain_id, "starting...\n")

    params      <- init
    lp_current  <- log_posterior(params)
    n_params    <- length(params)

    # Storage
    chain       <- matrix(NA, nrow=n_iter, ncol=n_params)
    colnames(chain) <- c("mstar","log_phi","alpha","beta")
    accept      <- rep(0, n_params)
    prop_sd     <- proposal_sd

    total_iter  <- n_burnin + n_iter

    for(iter in 1:total_iter) {

        # Update each parameter individually (Metropolis-within-Gibbs)
        for(j in 1:n_params) {

            params_prop      <- params
            params_prop[j]   <- params[j] + rnorm(1, 0, prop_sd[j])

            lp_prop          <- log_posterior(params_prop)
            log_accept_ratio <- lp_prop - lp_current

            if(log(runif(1)) < log_accept_ratio) {
                params     <- params_prop
                lp_current <- lp_prop
                if(iter > n_burnin) accept[j] <- accept[j] + 1
            }
        }

        # Adaptive proposal during burn-in
        if(iter <= n_burnin && iter %% adapt_interval == 0) {
            rate <- accept / pmax(iter - n_burnin, adapt_interval)
            prop_sd <- prop_sd * ifelse(rate < target_rate,
                                        0.8, 1.2)
            prop_sd <- pmax(prop_sd, 1e-4)
            accept  <- rep(0, n_params)  # reset acceptance count
        }

        # Store after burn-in
        if(iter > n_burnin) {
            chain[iter - n_burnin, ] <- params
        }

        # Progress
        if(iter %% 5000 == 0) {
            pct <- round(100 * iter / total_iter)
            cat(sprintf("  Chain %d: %d%% | lp=%.1f | params=[%.2f, %.2f, %.2f, %.2f]\n",
                        chain_id, pct, lp_current,
                        params[1], params[2], params[3], params[4]))
        }
    }

    accept_rate <- accept / n_iter
    cat(sprintf("  Chain %d done. Acceptance rates: mstar=%.2f, log_phi=%.2f, alpha=%.2f, beta=%.2f\n",
                chain_id, accept_rate[1], accept_rate[2],
                accept_rate[3], accept_rate[4]))

    list(chain=chain, prop_sd=prop_sd, accept_rate=accept_rate)
}

############################################################
# Run 3 chains
############################################################

n_burnin <- 10000
n_iter   <- 20000
thin     <- 10

# Starting points spread around expected values
starts <- list(
    c(13.5, -3.2, -1.3, 0.50),
    c(13.8, -3.6, -1.6, 0.45),
    c(13.2, -3.8, -1.1, 0.55)
)

# Initial proposal SDs (tuned roughly to expected posterior widths)
prop_sd_init <- c(0.15, 0.15, 0.10, 0.05)

chains <- list()
for(k in 1:3) {
    result    <- run_chain(starts[[k]], n_iter, n_burnin,
                           prop_sd_init, chain_id=k)
    chains[[k]] <- result$chain
}

############################################################
# Thin and combine
############################################################

thin_idx  <- seq(1, n_iter, by=thin)
posterior <- do.call(rbind, lapply(chains, function(ch) ch[thin_idx,]))

cat("\n=== POSTERIOR SUMMARY ===\n")
cat("Total samples:", nrow(posterior), "\n\n")

for(par in c("mstar","log_phi","alpha","beta")) {
    q <- quantile(posterior[,par], c(0.16, 0.5, 0.84))
    cat(sprintf("%-8s: %.3f  (+%.3f / -%.3f)\n",
                par, q[2], q[3]-q[2], q[2]-q[1]))
}

############################################################
# Gelman-Rubin convergence diagnostic
############################################################

cat("\n=== GELMAN-RUBIN ===\n")
mc_list <- lapply(chains, function(ch) mcmc(ch[thin_idx,]))
mc_coda <- as.mcmc.list(mc_list)
gr      <- gelman.diag(mc_coda)
print(gr)

if(max(gr$psrf[,1]) < 1.1) {
    cat("CONVERGED!\n")
} else {
    cat("Not fully converged - consider running longer.\n")
}

cat("\nEffective sample size:\n")
print(effectiveSize(mc_coda))

############################################################
# Plot
############################################################

mrp_log10 <- function(mstar, log_phi, alpha, beta, x) {
    log10(beta * log(10) * 10^log_phi *
              10^((alpha+1)*(x-mstar)) *
              exp(-10^(beta*(x-mstar))))
}

breaks   <- seq(12.5, 16.2, by=0.2)
hist_plt <- hist(x_obs, breaks=breaks, plot=FALSE)
phi_plt  <- hist_plt$counts / (Vsurvey * 0.2)
ok       <- phi_plt > 0

xfit <- seq(12, 16, length.out=500)
med  <- apply(posterior, 2, median)
q16  <- apply(posterior, 2, quantile, 0.16)
q84  <- apply(posterior, 2, quantile, 0.84)

CairoPDF("MRP_STAGE1_FINAL.pdf", 10, 7)

plot(hist_plt$mids[ok], log10(phi_plt[ok]),
     pch=19, col="darkgreen", cex=1.2,
     xlim=c(12,16), ylim=c(-8,-2),
     xlab=expression("Halo Mass  log"[10]*"(M/M"["\u2299"]*")"),
     ylab=expression("log"[10]*"("*phi*")  [Mpc"^{-3}*" dex"^{-1}*"]"),
     main="Stage 1: Unbinned Poisson Point Process (MH-MCMC)")

grid(col="gray80")

idx <- sample(1:nrow(posterior), 300)
for(i in idx) {
    y <- tryCatch(mrp_log10(posterior[i,"mstar"], posterior[i,"log_phi"],
                            posterior[i,"alpha"],  posterior[i,"beta"], xfit),
                  error=function(e) rep(NA, length(xfit)))
    if(all(is.finite(y))) lines(xfit, y, col=rgb(0,0,1,0.03))
}

lines(xfit, mrp_log10(med["mstar"],med["log_phi"],med["alpha"],med["beta"],xfit),
      col="red", lwd=3)

lines(xfit, mrp_log10(13.51,-3.19,-1.27,0.47,xfit),
      col="black", lwd=2, lty=2)

legend("bottomleft",
       legend=c("GAMA (binned for display)",
                "Posterior median", "Posterior samples", "Driver+22"),
       col=c("darkgreen","red",rgb(0,0,1,0.3),"black"),
       pch=c(19,NA,NA,NA), lty=c(NA,1,1,2),
       lwd=c(NA,3,1,2), bg="white", cex=0.9)

text(14.3,-2.3,
     sprintf("M* = %.2f (+%.2f/-%.2f)\nlog\u03c6* = %.2f (+%.2f/-%.2f)\n\u03b1 = %.2f (+%.2f/-%.2f)\n\u03b2 = %.2f (+%.2f/-%.2f)",
             med["mstar"],  q84["mstar"]-med["mstar"],   med["mstar"]-q16["mstar"],
             med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"],
             med["alpha"],  q84["alpha"]-med["alpha"],   med["alpha"]-q16["alpha"],
             med["beta"],   q84["beta"]-med["beta"],     med["beta"]-q16["beta"]),
     col="red", cex=0.85, pos=4)

dev.off()

# Trace plots
CairoPDF("MCMC_traces_stage1.pdf", 12, 9)
par(mfrow=c(4,1), mar=c(2,4,1,1))
for(par in c("mstar","log_phi","alpha","beta")) {
    plot(posterior[,par], type="l", col=rgb(0,0,0,0.3),
         ylab=par, xlab="", main="")
    abline(h=median(posterior[,par]), col="red", lwd=2)
}
dev.off()

cat("\n=== FINAL RESULTS ===\n")
cat(sprintf("           Median    +1sig    -1sig   Driver+22\n"))
cat(sprintf("M*:        %6.3f   +%.3f   -%.3f    13.51\n",
            med["mstar"],  q84["mstar"]-med["mstar"],  med["mstar"]-q16["mstar"]))
cat(sprintf("log_phi:   %6.3f   +%.3f   -%.3f    -3.19\n",
            med["log_phi"],q84["log_phi"]-med["log_phi"],med["log_phi"]-q16["log_phi"]))
cat(sprintf("alpha:     %6.3f   +%.3f   -%.3f    -1.27\n",
            med["alpha"],  q84["alpha"]-med["alpha"],  med["alpha"]-q16["alpha"]))
cat(sprintf("beta:      %6.3f   +%.3f   -%.3f     0.47\n",
            med["beta"],   q84["beta"]-med["beta"],    med["beta"]-q16["beta"]))

saveRDS(list(posterior=posterior, median=med, q16=q16, q84=q84,
             chains=chains, N=N, V=Vsurvey),
        "stage1_results.rds")

cat("\nSaved: stage1_results.rds\n")
cat("Plot:  MRP_STAGE1_FINAL.pdf\n")
cat("Trace: MCMC_traces_stage1.pdf\n")
