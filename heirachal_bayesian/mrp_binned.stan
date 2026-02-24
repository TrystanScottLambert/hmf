functions {

  real mrp_phi(real logM, real logMstar, real log_phi, real alpha, real beta) {

    real phi = pow(10, log_phi);
    real term1 = beta * log(10);
    real term2 = pow(10, (alpha+1) * (logM - logMstar));
    real term3 = exp(-pow(10, beta * (logM - logMstar)));

    return phi * term1 * term2 * term3;
  }

}

data {

  int<lower=1> N;
  vector[N] logM_obs;
  vector[N] log_phi_obs;
  vector[N] frac_err;

  real logbin;
  real vlimit;

  real logM_max;
}

parameters {

  real logMstar;
  real log_phi;
  real alpha;
  real<lower=0> beta;

}

model {

  // Priors tuned to realistic ranges
  logMstar ~ normal(14.4, 0.5);
  log_phi  ~ normal(-3, 1);
  alpha    ~ normal(-1.8, 0.5);
  beta     ~ normal(0.7, 0.3);

  // Likelihood
  for (n in 1:N) {

    real model_phi =
      mrp_phi(logM_obs[n], logMstar, log_phi, alpha, beta);

    real model_log = log10(model_phi);

    real sigma = frac_err[n] / log(10);

    target += normal_lpdf(log_phi_obs[n] | model_log, sigma);
  }

  // Correct penalty term (matching original code exactly)

  real integral = 0;

  for (i in 1:10) {

    real logM = logM_max + i * logbin;

    integral += mrp_phi(logM, logMstar, log_phi, alpha, beta);

  }

  integral *= logbin;

  target += -2 * vlimit * integral;
}

generated quantities {

  vector[N] model_logphi;

  for (n in 1:N)
    model_logphi[n] =
      log10(mrp_phi(logM_obs[n], logMstar, log_phi, alpha, beta));

}
