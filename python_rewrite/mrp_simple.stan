
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
  for (k in 1:Ng) xg[k] = xlo + (k - 1) * dx;
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

  // MRP on the grid
  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
  }
  // upper-cumulative integral: cum[k] = int_{xg[k]}^{xhi} phi dm
  vector[Ng] cum;
  cum[Ng] = 0;
  for (kr in 1:(Ng - 1)) {
    int k = Ng - kr;
    cum[k] = cum[k + 1] + 0.5 * (pg[k] + pg[k + 1]) * dx;
  }

  // Poisson normalisation
  real Lambda = 0;
  for (j in 1:Nsh) {
    int k0 = 1;
    for (k in 1:Ng) if (xg[k] <= mlim_sh[j]) k0 = k;
    if (k0 >= Ng) k0 = Ng - 1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;

  // Poisson point term (observed mass treated as truth)
  for (i in 1:N) {
    real u = x_obs[i] - ms;
    real phi_i = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
    if (phi_i > 1e-30)
      target += log(phi_i);
    else
      target += -100;
  }
}
