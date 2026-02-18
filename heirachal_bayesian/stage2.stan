
data {
  int<lower=0> N;
  vector[N] x;             // observed log10(mass)
  vector[N] sigma_x;       // measurement uncertainty per group
  real V;
  real xlo;
  real xhi;
  int<lower=1> Ng;
  int<lower=1> Ne;         // error integration grid points
}

transformed data {
  vector[Ng] xgrid;
  real dx;
  dx = (xhi - xlo) / (Ng - 1);
  for(k in 1:Ng)
    xgrid[k] = xlo + (k-1) * dx;
}

parameters {
  real mstar;
  real log_phi;
  real alpha;
  real<lower=0.1, upper=2> beta;
}

model {
  mstar   ~ normal(13.5, 1.0);
  log_phi ~ normal(-3.5, 1.5);
  alpha   ~ normal(-1.5, 0.8);

  // ---- Integral term (same as Stage 1) ----
  {
    vector[Ng] phi_grid;
    real Lambda;
    for(k in 1:Ng) {
      real u = beta * (xgrid[k] - mstar);
      phi_grid[k] = beta * log(10) * pow(10, log_phi)
                    * pow(10, (alpha+1) * (xgrid[k] - mstar))
                    * exp(-pow(10, u));
    }
    Lambda = V * sum(phi_grid) * dx;
    target += -Lambda;
  }

  // ---- Likelihood: each group's observed mass is
  //      drawn from the true MRP convolved with its
  //      measurement error Gaussian ----
  //
  //  p(x_obs | theta) = integral[ phi(x_true) * N(x_obs | x_true, sigma) dx_true ]
  //                   / integral[ phi(x_true) dx_true ]   <- normalisation
  //
  //  We approximate this integral numerically for each group

  for(i in 1:N) {
    // Fine grid around this observation (+-3 sigma)
    real x_lo_i = x[i] - 4 * sigma_x[i];
    real x_hi_i = x[i] + 4 * sigma_x[i];
    real dx_i   = (x_hi_i - x_lo_i) / (Ne - 1);

    vector[Ne] integrand;
    for(e in 1:Ne) {
      real x_true = x_lo_i + (e-1) * dx_i;
      real u      = beta * (x_true - mstar);
      real phi_true = beta * log(10) * pow(10, log_phi)
                      * pow(10, (alpha+1) * (x_true - mstar))
                      * exp(-pow(10, u));
      real gauss  = exp(-0.5 * square((x[i] - x_true) / sigma_x[i]))
                    / (sigma_x[i] * sqrt(2 * pi()));
      integrand[e] = phi_true * gauss;
    }

    real p_obs_i = sum(integrand) * dx_i * V;
    target += log(p_obs_i);
  }
}

