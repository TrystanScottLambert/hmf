
functions {
  real ln10() { return log(10.0); }

  // log of unnormalized MRP density over log10 M
  real mrp_log10_unnorm_lpdf(real x, real mstar, real alpha, real beta) {
    real t = pow(10.0, beta * (x - mstar));        // u = 10^{beta(x-m*)}
    return log(beta) + log(ln10()) + ln10() * (alpha + 1.0) * (x - mstar) - t;
  }

  // log normalization Z(theta; x_min) = Gamma(s, u_min) with s=(alpha+1)/beta
  real log_Z_shape(real mstar, real alpha, real beta, real x_min) {
    real s = (alpha + 1.0) / beta;
    real u_min = pow(10.0, beta * (x_min - mstar));
    // upper incomplete gamma: Gamma(s, u_min) = tgamma(s) * gamma_q(s, u_min)
    real Z = tgamma(s) * gamma_q(s, u_min);
    return log(Z);
  }
}

data {
  int<lower=1> N;
  vector[N] x_obs;                // observed log10(M)
  vector<lower=0>[N] sigma_x;     // mass errors in dex
  vector<lower=0>[N] w_shape;     // per-object shape weights
  real x_min;                     // truncation floor in log10(M)
  // weak priors
  real mstar_prior_mean;
  real mstar_prior_sd;
  real alpha_prior_mean;
  real alpha_prior_sd;
  real beta_prior_mean;
  real beta_prior_sd;
}

parameters {
  real mstar;
  real alpha;
  real<lower=0> beta;
  vector<lower=x_min>[N] x_true;  // latent true log10(M), truncated below x_min
}

model {
  // Priors
  mstar ~ normal(mstar_prior_mean, mstar_prior_sd);
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  beta  ~ normal(beta_prior_mean,  beta_prior_sd);

  // Precompute normalization over [x_min, ∞)
  real logZ = log_Z_shape(mstar, alpha, beta, x_min);

  // Likelihood
  for (n in 1:N) {
    target += w_shape[n] * ( mrp_log10_unnorm_lpdf(x_true[n] | mstar, alpha, beta) - logZ );
    x_obs[n] ~ normal(x_true[n], sigma_x[n]);   // Eddington bias generative term
  }
}

