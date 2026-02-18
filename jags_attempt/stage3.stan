
data {
  int<lower=0> N;
  vector[N] x;             // observed log10(mass)
  vector[N] sigma_x;       // per-group mass uncertainty
  real V;
  real xlo;
  real xhi;
  int<lower=1> Ng;
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

  // Latent TRUE masses for each group
  // This is what enables proper Eddington bias correction
  vector<lower=xlo, upper=xhi>[N] x_true;
}

model {
  mstar   ~ normal(13.5, 1.0);
  log_phi ~ normal(-3.5, 1.5);
  alpha   ~ normal(-1.5, 0.8);

  // ---- Integral term ----
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

  // ---- True mass distribution: MRP ----
  for(i in 1:N) {
    real u = beta * (x_true[i] - mstar);
    real log_phi_i = log(beta) + log(log(10))
                     + log_phi * log(10)
                     + (alpha+1) * (x_true[i] - mstar) * log(10)
                     - pow(10, u);
    target += log(V) + log_phi_i;
  }

  // ---- Observed mass | true mass: Gaussian error ----
  x ~ normal(x_true, sigma_x);
}

