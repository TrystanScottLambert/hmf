
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;          // per-object log-mass error
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  vector<lower=0>[Nsh] sig_sh;     // representative boundary sigma per shell
  real xhi;
  int<lower=2> Ng;
}
transformed data {
  real ln10 = log(10.0);
  real xlo = min(mlim_sh) - 2.0;   // go well below mlim to capture down-scatter
  real dx = (xhi - xlo) / (Ng - 1.0);
  real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
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

  // phi on the global grid (computed once)
  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
  }

  // ---- Lambda: boundary-consistent normalisation (soft Phi cut) ----
  real Lambda = 0;
  for (j in 1:Nsh) {
    real acc = 0;
    for (k in 1:Ng) {
      real w = Phi((xg[k] - mlim_sh[j]) / sig_sh[j]);  // P(scatter above limit)
      real term = pg[k] * w;
      acc += (k == 1 || k == Ng) ? 0.5 * term : term;  // trapezoid
    }
    Lambda += V_sh[j] * acc * dx;
  }
  target += -Lambda;

  // ---- per-object: marginalise the latent true mass over the full grid ----
  for (i in 1:N) {
    real inv_s = 1.0 / sig[i];
    real sm = 0;
    for (k in 1:Ng) {
      real z = (x_obs[i] - xg[k]) * inv_s;
      real gauss = inv_sqrt2pi * inv_s * exp(-0.5 * z * z);
      real term = pg[k] * gauss;
      sm += (k == 1 || k == Ng) ? 0.5 * term : term;  // trapezoid
    }
    sm *= dx;
    if (sm > 1e-300)
      target += log(sm);
    else
      target += -300;
  }
}
