
data {
  int<lower=0> N;          // number of observed groups
  vector[N] x;             // observed log10(mass)
  real V;                  // survey volume [Mpc^3]
  real xlo;                // lower mass limit
  real xhi;                // upper mass limit
  int<lower=1> Ng;         // integration grid size
}

transformed data {
  // Fixed integration grid
  vector[Ng] xgrid;
  real dx;
  dx = (xhi - xlo) / (Ng - 1);
  for(k in 1:Ng)
    xgrid[k] = xlo + (k-1) * dx;
}

parameters {
  real mstar;              // characteristic mass [log10 Msun]
  real log_phi;            // log10(phi*) [Mpc^-3 dex^-1]
  real alpha;              // low-mass slope
  real<lower=0.1, upper=2> beta;  // high-mass cutoff shape
}

model {
  // ---- Priors ----
  mstar   ~ normal(13.5, 1.0);
  log_phi ~ normal(-3.5, 1.5);
  alpha   ~ normal(-1.5, 0.8);
  // beta has implicit uniform prior from bounds

  // ---- MRP helper: log phi(x) ----
  // log phi(x) = log(beta) + log(log10) + log_phi*log(10)
  //            + (alpha+1)*(x-mstar)*log(10) - 10^(beta*(x-mstar))

  // ---- Integral: expected total number of groups ----
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

    // Integral term of Poisson point process
    target += -Lambda;
  }

  // ---- Sum of log phi over observed groups ----
  for(i in 1:N) {
    real u = beta * (x[i] - mstar);
    real log_phi_i = log(beta) + log(log(10))
                     + log_phi * log(10)
                     + (alpha+1) * (x[i] - mstar) * log(10)
                     - pow(10, u);
    target += log(V) + log_phi_i;
  }
}

