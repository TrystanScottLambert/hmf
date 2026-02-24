
functions {
  real ln10() { return log(10.0); }

  // log of unnormalized MRP over log10 M
  real mrp_log10_unnorm_lpdf(real x, real mstar, real alpha, real beta) {
    real t = pow(10.0, beta * (x - mstar));        // u = 10^{beta(x-m*)}
    return log(beta) + log(ln10()) + ln10() * (alpha + 1.0) * (x - mstar) - t;
  }

  // positive unnormalized density in log space (for integration)
  real mrp_unnorm(real x, real[] theta, real[] x_r, int[] x_i) {
    real mstar = theta[1];
    real alpha = theta[2];
    real beta  = theta[3];
    return exp( mrp_log10_unnorm_lpdf(x, mstar, alpha, beta) );
  }

  // log normalization Z(theta; x_min) = ∫_{x_min}^{XUP} f(x|θ) dx  (XUP large)
  real log_Z_numeric(real mstar, real alpha, real beta, real x_min) {
    real theta[3];
    theta[1] = mstar; theta[2] = alpha; theta[3] = beta;
    real x_r[0];
    int x_i[0];
    // Fixed large upper limit in log10(M/Msun) safely beyond cluster masses:
    real XUP = 18.5; // 10^18.5 Msun >> any realistic halo; tail negligible
    real Z = integrate_1d(mrp_unnorm, x_min, XUP, theta, x_r, x_i, 1e-6);
    return log(Z);
  }
}

data {
  int<lower=1> N;
  vector[N] x_obs;                // observed log10(M)
  vector<lower=0>[N] sigma_x;     // mass errors in dex
  vector<lower=0>[N] w_shape;     // per-object shape weights (usually 1)
  real x_min;                     // truncation floor in log10(M)
  // priors
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

  // Numeric normalization over [x_min, ∞)
  real logZ = log_Z_numeric(mstar, alpha, beta, x_min);

  // Likelihood: shape-only + measurement
  for (n in 1:N) {
    target += w_shape[n] * ( mrp_log10_unnorm_lpdf(x_true[n] | mstar, alpha, beta) - logZ );
    x_obs[n] ~ normal(x_true[n], sigma_x[n]);   // Eddington bias generative term
  }
}

