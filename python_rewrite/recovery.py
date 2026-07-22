"""
============================================================
SHARK SIMPLE RECOVERY TEST  (cmdstanpy version)
============================================================

Closed-loop recovery test for the MRP halo mass function.

Pipeline:
  1. Read SHARK/WAVES groups + galaxies.
  2. Abundance-match group masses onto the Driver+22 MRP  -> TRUTH.
  3. Apply a GAMA-like detection rule (>= `multi` members brighter
     than `mag_limit`).
  4. Add multiplicity-dependent Gaussian mass errors (observed masses).
  5. Build mlim(z) from the *observed* masses (turnover / histogram mode).
  6. Fit the SIMPLE POISSON model (document-1 `stan_simple`, run via
     cmdstanpy) and check we recover the injected M*, log_phi*, alpha, beta.

This is the BASELINE bug-check only: simple model (no error de-biasing),
single mock realisation, turnover mlim. The coverage loop, the marginalised
model with the boundary fix, and the calibrated mass proxy come later.

Requires: numpy, pandas, pyarrow, cmdstanpy, matplotlib, and a CmdStan
toolchain. If CmdStan is not yet installed:
    python -c "import cmdstanpy; cmdstanpy.install_cmdstan()"

Run:
    python recovery.py            # real pipeline (needs the parquet files)
    python recovery.py --selftest # synthetic recovery, no data files needed
============================================================
"""

import argparse
import os

import re
import numpy as np
import pandas as pd
from scipy.stats import rankdata

rng_global = np.random.default_rng(42)

# numpy>=2.0 renamed trapz -> trapezoid; support both
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ----------------------------------------------------------
# 0. Configuration  (edit paths / numbers here)
# ----------------------------------------------------------
DATA_DIR = "/Users/00115372/Desktop/mock_catalogs/offical_waves_mocks/v0.5.0"

# Driver+22 MRP = injected truth
TRUE = dict(
    ms=14.13, lp=-3.96, al=-1.68, be=0.63
)  # Driver+22 (2022) MRP abstract values
PARAMS = ["ms", "lp", "al", "be"]

# Cosmology
H0, OMEGA_M = 67.37, 0.3147

# Survey / selection
ZMIN, ZLIMIT = 0.01, 0.25
MULTI = 5  # min members for a detection
# Selection band + limit. GAMA: r_SDSS < 19.8 ; WAVES: Z_VISTA < 21.1 (deeper NIR).
SEL_COL = "total_ap_dust_r_SDSS"
MAG_LIMIT = 19.8
# Dynamical-mass calibration prefactor M = A * sigma^2 R / G. Driver's fiducial
# is 13.9; A=10 is his variant (Fig. A2). Applied to BOTH surveys for a common
# mass scale.
A_SCALE = 10.0

# Completeness ramp C(Delta) = 0.5(1+erf((Delta-D50)/(sqrt2 w))), Delta=m-mlim(z),
# measured from the mock (measure_completeness.py). z-dependent: interpolate
# (D50, w) from these per-z-bin values. CMIN = floor (drop groups below 20%
# completeness; corrections below that are unreliable, cf. Driver+22).
COMP_Z_PTS = [0.045, 0.115, 0.20]
COMP_D50_PTS = [-0.148, -0.193, -0.232]
COMP_W_PTS = [0.326, 0.256, 0.227]
CMIN = 0.2

# Mock only: inject mass errors of SIGMA_INJECT_SCALE * sigma_reported while the
# fit is still given sigma_reported. 1.0 -> the reported errors are correct
# (s_scale truth = 1). Set to e.g. 1.5 to check marg_comp_serr recovers a scale
# it was not handed, which is the real test of the hyperparameter.
SIGMA_INJECT_SCALE = 1.0

# Use Driver+22's GSR chains as a multivariate-normal prior (preserves the
# rho ~ -0.97 parameter covariances that independent Gaussians discard).
# Set via --driver-prior. DRIVER_PRIOR_INFLATE widens the covariance to soften
# the double-counting (Driver's GSR posterior used the same GAMA+SDSS data).
USE_DRIVER_PRIOR = False
DRIVER_PRIOR_PATH = "../data/hmfparams_gsr.csv"
DRIVER_PRIOR_INFLATE = 1.0

# ---------------------------------------------------------------------------
# alpha bias calibration (diagnose_pinned_mstar.py, mock with known truth).
# The completeness model returns an alpha that is biased by an amount that
# depends on where M* sits: fitting the mock (truth alpha = -1.68) gives
#     M* = 14.130 (pinned)  ->  alpha = -1.516   bias = +0.164
#     M* = 14.175 (free)    ->  alpha = -1.563   bias = +0.117
#     M* = 14.600 (pinned)  ->  alpha = -1.894   bias = -0.214
# i.e. bias falls ~linearly with M*. Correcting real fits by this removes the
# apparent 0.39 dex spread between the REFLEX-anchored and unanchored results
# (both -> alpha ~ -1.65).
# CAVEAT: ONE mock realisation, so there is no uncertainty on the correction
# yet. Run calibrate_alpha_bias() over many realisations before quoting it.
ALPHA_BIAS_MS = [14.130, 14.175, 14.600]
ALPHA_BIAS_DA = [0.164, 0.117, -0.214]


def alpha_bias(ms):
    """Fitted-minus-true alpha as a function of where M* sits (linear fit to
    the mock calibration points; extrapolates linearly outside them)."""
    ms = np.asarray(ms, float)
    b, a = np.polyfit(ALPHA_BIAS_MS, ALPHA_BIAS_DA, 1)  # slope, intercept
    return a + b * ms


def corrected_alpha(flat):
    """Per-draw bias-corrected alpha: alpha_i - bias(ms_i). Applying it draw by
    draw propagates the M*-alpha covariance instead of correcting the median."""
    return flat[:, 2] - alpha_bias(flat[:, 0])


ADD_ERRORS = True

# Likelihood grid / integration (passed to Stan as data)
XHI = 16.5  # upper mass bound for the Lambda integral
NG = 300  # grid points for the phi integral
NINT = 30  # local per-object integration points (marginalised model; exact to <1e-4)
NSH = 20  # redshift shells for the Poisson normalisation
SIG_SH_FLOOR = 0.25  # min boundary sigma for the Phi soft cut (keeps HMC stable)


# ----------------------------------------------------------
# 1. Physics helpers
# ----------------------------------------------------------
def mrp_phi(x, ms, lp, al, be):
    """MRP number density per dex.  x = log10(M)."""
    u = x - ms
    return be * np.log(10) * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def comoving_distance(z, H0=H0, Om=OMEGA_M, ngrid=4000):
    """Line-of-sight comoving distance [Mpc].  Vectorised via a fine grid."""
    c = 299792.458
    zg = np.linspace(0.0, float(np.max(z)) * 1.001 + 1e-6, ngrid)
    Ez = np.sqrt(Om * (1 + zg) ** 3 + (1 - Om))
    integ = np.concatenate(
        [[0.0], np.cumsum(0.5 * (1 / Ez[1:] + 1 / Ez[:-1]) * np.diff(zg))]
    )
    return (c / H0) * np.interp(z, zg, integ)


def sky_area_deg2(ra, dec):
    """Spherical area of the RA/Dec bounding box [deg^2].
    NOTE: bounding-box overestimate if the footprint is not a filled
    rectangle -- identical to celestial::skyarea as used in doc 1."""
    ra_min, ra_max = np.min(ra), np.max(ra)
    d_lo, d_hi = np.radians(np.min(dec)), np.radians(np.max(dec))
    return (ra_max - ra_min) * (180 / np.pi) * (np.sin(d_hi) - np.sin(d_lo))


def survey_volume(sky_frac, zmin=ZMIN, zmax=ZLIMIT):
    d = comoving_distance(np.array([zmin, zmax]))
    return (4 / 3) * np.pi * (d[1] ** 3 - d[0] ** 3) * sky_frac


def shell_volumes(sky_frac, zmin=ZMIN, zmax=ZLIMIT, nsh=NSH):
    z_edges = np.linspace(zmin, zmax, nsh + 1)
    z_mids = 0.5 * (z_edges[1:] + z_edges[:-1])
    d = comoving_distance(z_edges)
    V_sh = (4 / 3) * np.pi * (d[1:] ** 3 - d[:-1] ** 3) * sky_frac
    return z_mids, V_sh


# ----------------------------------------------------------
# 2. Data layer
# ----------------------------------------------------------
def load_catalogues(data_dir):
    groups = pd.read_parquet(f"{data_dir}/waves_wide_groups.parquet")
    galaxies = pd.read_parquet(f"{data_dir}/waves_wide_gals.parquet")
    groups = groups[groups["dec"] > 0].copy()
    galaxies = galaxies[galaxies["dec"] > 0].copy()
    return groups, galaxies


def abundance_match(groups, sky_frac):
    """Assign each group a mass by matching its rank to the cumulative MRP
    count in the survey volume.  Deterministic -> truth is exactly MRP."""
    gv = groups[
        (groups["zcos"] > ZMIN) & (groups["zcos"] < ZLIMIT) & (groups["mvir"] > 0)
    ].copy()

    Vsurvey = survey_volume(sky_frac)
    m_grid = np.arange(9.0, 16.5 + 1e-9, 0.001)
    phi = mrp_phi(m_grid, **TRUE)
    cum_counts = np.cumsum((phi * 0.001)[::-1])[::-1] * Vsurvey  # counts ABOVE each m

    rank = rankdata(-gv["mvir"].values, method="average")  # 1 = most massive
    gv["log_mass_am"] = np.interp(rank, cum_counts[::-1], m_grid[::-1])
    return gv, Vsurvey


def gama_select(gv, galaxies):
    """Detection = >= MULTI galaxies brighter than MAG_LIMIT in SEL_COL
    (and log Mstar > 8).

    NOTE: log_mstar_total is log10(Mstar/Msun) (~8-12), so the stellar-mass
    floor is `> 8`, i.e. log10(1e8) -- NOT `> 1e8`."""
    if SEL_COL not in galaxies.columns:
        raise KeyError(
            f"selection column {SEL_COL!r} not in galaxy catalogue; "
            f"available magnitude-like columns: "
            f"{[c for c in galaxies.columns if 'ap_dust' in c or 'mag' in c.lower()]}"
        )
    gal = galaxies[
        (galaxies["id_fof"] != -1)
        & (galaxies["log_mstar_total"] > 8)
        & (galaxies[SEL_COL] < MAG_LIMIT)
    ]
    counts = gal.groupby("id_group_sky").size().rename("n_gama")
    gv = gv.merge(counts, left_on="id_group_sky", right_index=True, how="left")
    gv["n_gama"] = gv["n_gama"].fillna(0).astype(int)
    gv["detected"] = gv["n_gama"] >= MULTI
    return gv


def sigma_from_nfof(n_gama):
    """Multiplicity-dependent log-mass error model (GAMA-like).
    Single source of truth used by both the single-run and coverage paths."""
    xx = np.arange(2, 23)
    yy = np.array(
        [
            0.68,
            0.39,
            0.40,
            0.33,
            0.28,
            0.24,
            0.20,
            0.19,
            0.17,
            0.14,
            0.14,
            0.13,
            0.14,
            0.12,
            0.12,
            0.10,
            0.10,
            0.10,
            0.09,
            0.08,
            0.08,
        ]
    )
    return np.maximum(np.interp(np.asarray(n_gama, float), xx, yy), 0.10)


def add_mass_errors(gv, rng):
    """Add multiplicity-dependent Gaussian errors to the detected groups."""
    det = gv[gv["detected"]].copy()
    sigma = sigma_from_nfof(det["n_gama"].values)
    if ADD_ERRORS:
        m_obs = det["log_mass_am"].values + rng.normal(0, sigma * SIGMA_INJECT_SCALE)
    else:
        m_obs = det["log_mass_am"].values.copy()
        sigma = np.zeros_like(sigma)
    return det["zcos"].values, m_obs, sigma, det["n_gama"].values


# ----------------------------------------------------------
# 3. mlim(z) -- turnover / histogram mode (validated winner)
# ----------------------------------------------------------
def _aic_ols(y, X):
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    rss = np.sum((y - X @ beta) ** 2)
    n, k = len(y), X.shape[1]
    return n * np.log(rss / n) + 2 * k, beta


def turnover_mlim(
    z_obs, m_obs, nbin_z=30, hist_bw=0.3, min_in_bin=20, zmin=ZMIN, zmax=ZLIMIT
):
    z_edges = np.linspace(zmin, zmax, nbin_z + 1)
    z_mids = 0.5 * (z_edges[1:] + z_edges[:-1])
    mass_edges = np.arange(10, 16 + 1e-9, hist_bw)
    mass_mids = 0.5 * (mass_edges[1:] + mass_edges[:-1])

    turn = np.full(nbin_z, np.nan)
    for b in range(nbin_z):
        m = m_obs[(z_obs >= z_edges[b]) & (z_obs < z_edges[b + 1])]
        if m.size > min_in_bin:
            counts, _ = np.histogram(m, bins=mass_edges)
            turn[b] = mass_mids[np.argmax(counts)]

    ok = np.isfinite(turn)
    if ok.sum() < 3:
        raise RuntimeError(
            f"turnover_mlim: only {ok.sum()} usable z-bins -- too few detected "
            "groups to build mlim(z). Check the selection cuts."
        )
    zb, tb = z_mids[ok], turn[ok]
    aic_l, beta_l = _aic_ols(tb, np.column_stack([np.ones_like(zb), zb]))
    aic_q, beta_q = _aic_ols(tb, np.column_stack([np.ones_like(zb), zb, zb**2]))

    if aic_q < aic_l - 2:
        c = beta_q
        func = lambda z: c[0] + c[1] * z + c[2] * z**2
        kind = "quad"
    else:
        c = beta_l
        func = lambda z: c[0] + c[1] * z
        kind = "linear"
    return func, c, kind, (z_mids, turn)


# ----------------------------------------------------------
# 4. Stan models -- simple (doc-1) and marginalised+boundary
# ----------------------------------------------------------
# SIMPLE: observed mass treated as truth, sharp cut at mlim. Baseline.
SIMPLE_CODE = r"""
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
  ms ~ normal(14.13, 0.42);   // Driver+22 M* (broad -> data-driven)
  lp ~ normal(-3.96, 0.69);   // Driver+22 logphi* (broad -> data-driven)
  al ~ normal(-1.68, 0.22);   // Driver+22 alpha (informative)
  be ~ normal(0.63, 0.02);    // beta pinned at Driver+22 value

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
"""

# MARGINALISED + BOUNDARY: integrate out the latent true mass per object, and
# replace the sharp mlim cut in Lambda with a Phi-weighted soft boundary.
#
#   per-object:  L_i = int phi(m_t) * Normal(x_obs_i | m_t, sig_i) dm_t   (all m_t)
#   Lambda    :  sum_j V_sh[j] * int phi(m_t) * Phi((m_t - mlim_sh[j]) / sig_sh[j]) dm_t
#
# Both integrals use the same global grid xg (phi computed once). xlo extends
# well below min(mlim) so down-scattered groups (m_t < mlim) are captured.
MARG_CODE = r"""
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;          // per-object log-mass error
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  vector<lower=0>[Nsh] sig_sh;     // representative boundary sigma per shell
  real xhi;
  int<lower=2> Ng;        // global grid for the Lambda integral
  int<lower=2> Nint;      // local grid points per object (per-object integral)
}
transformed data {
  real ln10 = log(10.0);
  real xlo = min(mlim_sh) - 1.0;   // Lambda grid; Phi kills the integrand below mlim
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
  ms ~ normal(14.13, 0.42);   // Driver+22 M* (broad -> data-driven)
  lp ~ normal(-3.96, 0.69);   // Driver+22 logphi* (broad -> data-driven)
  al ~ normal(-1.68, 0.22);   // Driver+22 alpha (informative)
  be ~ normal(0.63, 0.02);    // beta pinned at Driver+22 value

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

  // ---- per-object: numerical marginalisation on a local grid (+-5 sigma) ----
  // L_i = int phi(m) Normal(x_obs|m,sig) dm, over a narrow window around x_obs
  // (still reaches below mlim -> down-scatter captured). Nint=30 / +-5 sigma is
  // exact to <1e-4 vs a fine integral. Numerical (not the Laplace analytic form)
  // because the analytic version, while accurate near truth, develops spurious
  // structure at extreme parameters and lets the sampler wander off.
  for (i in 1:N) {
    real lo_i = x_obs[i] - 5 * sig[i];
    real hi_i = x_obs[i] + 5 * sig[i];
    real dmt = (hi_i - lo_i) / (Nint - 1.0);
    real inv_s = 1.0 / sig[i];
    real sm = 0;
    for (g in 1:Nint) {
      real mt = lo_i + (g - 1) * dmt;
      real u = mt - ms;
      real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
      real zsc = (x_obs[i] - mt) * inv_s;
      real term = phi_g * exp(-0.5 * zsc * zsc);
      sm += (g == 1 || g == Nint) ? 0.5 * term : term;   // trapezoid
    }
    sm *= dmt * inv_s * inv_sqrt2pi;                       // Gaussian 1/(sig sqrt(2pi))
    if (sm > 1e-300)
      target += log(sm);
    else
      target += -300;
  }
}
"""

# Verbatim port of the R production model (run.R stan_marg): wide priors
# (al ~ N(-1.3,1.0), be unconstrained within [0.1,2]), per-object integral
# from mlim upward, Lambda via a sharp cut on a cumulative grid. This exists
# ONLY to reproduce the R fit on the same GAMA data (a port-consistency
# check), NOT for the mock recovery -- do not confuse it with MARG_CODE.
GAMA_CODE = r"""
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
  int<lower=2> Nint;
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
  al ~ normal(-1.3, 1.0);

  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
  }
  vector[Ng] cum;
  cum[Ng] = 0;
  for (kr in 1:(Ng - 1)) {
    int k = Ng - kr;
    cum[k] = cum[k + 1] + 0.5 * (pg[k] + pg[k + 1]) * dx;
  }
  real Lambda = 0;
  for (j in 1:Nsh) {
    int k0 = 1;
    for (k in 1:Ng) if (xg[k] <= mlim_sh[j]) k0 = k;
    if (k0 >= Ng) k0 = Ng - 1;
    Lambda += V_sh[j] * cum[k0];
  }
  target += -Lambda;

  for (i in 1:N) {
    real lo_i = mlim[i];
    real hi_i = fmin(xhi, fmax(x_obs[i] + 5 * sig[i], mlim[i] + 8 * sig[i]));
    real dmt = (hi_i - lo_i) / (Nint - 1.0);
    if (dmt < 1e-6) {
      target += -100;
    } else {
      real sum_trap = 0;
      for (g in 1:Nint) {
        real mt_g = lo_i + (g - 1) * dmt;
        real u = mt_g - ms;
        real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
        real gauss_g = exp(-0.5 * square((x_obs[i] - mt_g) / sig[i])) / (sig[i] * 2.5066283);
        real integrand = phi_g * gauss_g;
        if (g == 1 || g == Nint) sum_trap += 0.5 * integrand;
        else sum_trap += integrand;
      }
      sum_trap *= dmt;
      if (sum_trap > 1e-30) target += log(sum_trap);
      else target += -100;
    }
  }
}
"""

# Combined two-survey model: one shared MRP, each survey contributes its own
# -Lambda + Sum(per-object) with its own volumes/mlim/sig_sh. The per-survey
# likelihood is written once as a Stan function and called for each survey.
MARG_COMBINED_CODE = r"""
functions {
  real survey_contrib(vector x_obs, vector sig, vector V_sh, vector mlim_sh,
                      vector sig_sh, real xhi, int Ng, int Nint,
                      real ms, real lp, real al, real be) {
    real ln10 = log(10.0);
    real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
    int N = num_elements(x_obs);
    int Nsh = num_elements(V_sh);
    real xlo = min(mlim_sh) - 1.0;
    real dx = (xhi - xlo) / (Ng - 1.0);
    vector[Ng] xg;
    vector[Ng] pg;
    real Lambda = 0;
    real out;
    for (k in 1:Ng) xg[k] = xlo + (k - 1) * dx;
    for (k in 1:Ng) {
      real u = xg[k] - ms;
      pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
    }
    for (j in 1:Nsh) {
      real acc = 0;
      for (k in 1:Ng) {
        real w = Phi((xg[k] - mlim_sh[j]) / sig_sh[j]);
        real term = pg[k] * w;
        acc += (k == 1 || k == Ng) ? 0.5 * term : term;
      }
      Lambda += V_sh[j] * acc * dx;
    }
    out = -Lambda;
    for (i in 1:N) {
      real lo_i = x_obs[i] - 5 * sig[i];
      real hi_i = x_obs[i] + 5 * sig[i];
      real dmt = (hi_i - lo_i) / (Nint - 1.0);
      real inv_s = 1.0 / sig[i];
      real sm = 0;
      for (g in 1:Nint) {
        real mt = lo_i + (g - 1) * dmt;
        real u = mt - ms;
        real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
        real zsc = (x_obs[i] - mt) * inv_s;
        real term = phi_g * exp(-0.5 * zsc * zsc);
        sm += (g == 1 || g == Nint) ? 0.5 * term : term;
      }
      sm *= dmt * inv_s * inv_sqrt2pi;
      out += (sm > 1e-300) ? log(sm) : -300;
    }
    return out;
  }
}
data {
  int<lower=1> N_a; vector[N_a] x_obs_a; vector<lower=0>[N_a] sig_a;
  int<lower=1> Nsh_a; vector[Nsh_a] V_sh_a; vector[Nsh_a] mlim_sh_a;
  vector<lower=0>[Nsh_a] sig_sh_a;
  int<lower=1> N_b; vector[N_b] x_obs_b; vector<lower=0>[N_b] sig_b;
  int<lower=1> Nsh_b; vector[Nsh_b] V_sh_b; vector[Nsh_b] mlim_sh_b;
  vector<lower=0>[Nsh_b] sig_sh_b;
  real xhi; int<lower=2> Ng; int<lower=2> Nint;
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
}
model {
  ms ~ normal(14.13, 0.42);
  lp ~ normal(-3.96, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.02);    // beta pinned at Driver+22 value
  target += survey_contrib(x_obs_a, sig_a, V_sh_a, mlim_sh_a, sig_sh_a,
                           xhi, Ng, Nint, ms, lp, al, be);
  target += survey_contrib(x_obs_b, sig_b, V_sh_b, mlim_sh_b, sig_sh_b,
                           xhi, Ng, Nint, ms, lp, al, be);
}
"""

# Completeness forward-model: replaces the sharp mlim cut with the measured
# erf ramp C(m,z), applied CONSISTENTLY in both the per-object term and Lambda.
# z-dependent (d50/w passed per object and per shell), floored at cmin. Keeps
# all detected groups above the floor (no mlim cut). This is the boundary fix.
MARG_COMP_CODE = r"""
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;
  vector[N] mlim_obj;            // mlim(z_i) per object
  vector[N] d50_obj;             // completeness D50(z_i)
  vector<lower=0>[N] w_obj;      // completeness width w(z_i)
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  vector[Nsh] d50_sh;
  vector<lower=0>[Nsh] w_sh;
  real xhi;
  int<lower=2> Ng;
  int<lower=2> Nint;
  real cmin;
}
transformed data {
  real ln10 = log(10.0);
  real sqrt2 = sqrt(2.0);
  real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
  real xlo = min(mlim_sh) - 2.5;
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
  ms ~ normal(14.13, 0.42);
  lp ~ normal(-3.96, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.02);    // beta pinned at Driver+22 value

  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
  }

  // Lambda = sum_j V_sh int phi(m) C_j(m) dm, C floored at cmin
  real Lambda = 0;
  for (j in 1:Nsh) {
    real acc = 0;
    for (k in 1:Ng) {
      real C = 0.5 * (1 + erf((xg[k] - mlim_sh[j] - d50_sh[j]) / (sqrt2 * w_sh[j])));
      real Cf = C > cmin ? C : 0.0;
      real term = pg[k] * Cf;
      acc += (k == 1 || k == Ng) ? 0.5 * term : term;
    }
    Lambda += V_sh[j] * acc * dx;
  }
  target += -Lambda;

  // per-object: int phi(m) C_i(m) N(x|m,sig) dm on a local +-6 sigma grid
  for (i in 1:N) {
    real lo_i = x_obs[i] - 6 * sig[i];
    real hi_i = x_obs[i] + 6 * sig[i];
    real dmt = (hi_i - lo_i) / (Nint - 1.0);
    real inv_s = 1.0 / sig[i];
    real sm = 0;
    for (g in 1:Nint) {
      real mt = lo_i + (g - 1) * dmt;
      real u = mt - ms;
      real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
      real C = 0.5 * (1 + erf((mt - mlim_obj[i] - d50_obj[i]) / (sqrt2 * w_obj[i])));
      real zsc = (x_obs[i] - mt) * inv_s;
      real term = phi_g * C * exp(-0.5 * zsc * zsc);
      sm += (g == 1 || g == Nint) ? 0.5 * term : term;
    }
    sm *= dmt * inv_s * inv_sqrt2pi;
    if (sm > 1e-300)
      target += log(sm);
    else
      target += -300;
  }
}
"""

# Combined GAMA+SDSS with completeness. One shared MRP. GAMA's ramp is FIXED
# (measured from the GAMA-selected mock). SDSS's ramp (D50, w) is FITTED, with
# priors informed by GAMA's measurement -- the WAVES lightcone is too small in
# area to build an SDSS-like mock (71 groups), so we marginalise over the SDSS
# selection rather than assume it. The ramp shape is imprinted on the observed
# counts near SDSS's limit, and 4894 groups constrain it.
COMBINED_COMP_CODE = r"""
functions {
  // Both surveys use FIXED, precomputed completeness (Cobj per-object grid,
  // Csh per-shell grid) passed as data -- no completeness parameters, so no
  // erf inside the sampler and no ramp/MRP degeneracy.
  real survey_ll_fixC(vector x_obs, vector sig, matrix Cobj, matrix mt_obj,
                      vector V_sh, matrix Csh, vector xg, real dx,
                      int Nint, real ms, real lp, real al, real be) {
    real ln10 = log(10.0);
    real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
    int N = num_elements(x_obs);
    int Nsh = num_elements(V_sh);
    int Ng = num_elements(xg);
    vector[Ng] pg;
    real Lambda = 0;
    real out;
    for (k in 1:Ng) {
      real u = xg[k] - ms;
      pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
    }
    for (j in 1:Nsh) {
      real acc = 0;
      for (k in 1:Ng) {
        real term = pg[k] * Csh[j, k];
        acc += (k == 1 || k == Ng) ? 0.5 * term : term;
      }
      Lambda += V_sh[j] * acc * dx;
    }
    out = -Lambda;
    for (i in 1:N) {
      real inv_s = 1.0 / sig[i];
      real dmt = mt_obj[i, 2] - mt_obj[i, 1];
      real sm = 0;
      for (g in 1:Nint) {
        real mt = mt_obj[i, g];
        real u = mt - ms;
        real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
        real zsc = (x_obs[i] - mt) * inv_s;
        real term = phi_g * Cobj[i, g] * exp(-0.5 * zsc * zsc);
        sm += (g == 1 || g == Nint) ? 0.5 * term : term;
      }
      sm *= dmt * inv_s * inv_sqrt2pi;
      out += (sm > 1e-300) ? log(sm) : -300;
    }
    return out;
  }
}
data {
  // GAMA block (ramp FIXED, measured)
  int<lower=1> N_a; vector[N_a] x_obs_a; vector<lower=0>[N_a] sig_a;
  vector[N_a] mlim_obj_a; vector[N_a] d50_obj_a; vector<lower=0>[N_a] w_obj_a;
  int<lower=1> Nsh_a; vector[Nsh_a] V_sh_a; vector[Nsh_a] mlim_sh_a;
  vector[Nsh_a] d50_sh_a; vector<lower=0>[Nsh_a] w_sh_a;
  // SDSS block (ramp FIXED too: d50/w passed as data, not fitted)
  int<lower=1> N_b; vector[N_b] x_obs_b; vector<lower=0>[N_b] sig_b;
  vector[N_b] mlim_obj_b; vector[N_b] d50_obj_b; vector<lower=0>[N_b] w_obj_b;
  int<lower=1> Nsh_b; vector[Nsh_b] V_sh_b; vector[Nsh_b] mlim_sh_b;
  vector[Nsh_b] d50_sh_b; vector<lower=0>[Nsh_b] w_sh_b;
  real xhi; int<lower=2> Ng; int<lower=2> Nint; real cmin;
}
transformed data {
  real sqrt2 = sqrt(2.0);
  // ---- GAMA precompute ----
  real xlo_a = min(mlim_sh_a) - 2.5;
  real dx_a = (xhi - xlo_a) / (Ng - 1.0);
  vector[Ng] xg_a;
  matrix[Nsh_a, Ng] Csh_a;
  matrix[N_a, Nint] mt_a;
  matrix[N_a, Nint] Cobj_a;
  // ---- SDSS precompute ----
  real xlo_b = min(mlim_sh_b) - 2.5;
  real dx_b = (xhi - xlo_b) / (Ng - 1.0);
  vector[Ng] xg_b;
  matrix[Nsh_b, Ng] Csh_b;
  matrix[N_b, Nint] mt_b;
  matrix[N_b, Nint] Cobj_b;
  for (k in 1:Ng) xg_a[k] = xlo_a + (k - 1) * dx_a;
  for (j in 1:Nsh_a) for (k in 1:Ng) {
    real C = 0.5 * (1 + erf((xg_a[k] - mlim_sh_a[j] - d50_sh_a[j]) / (sqrt2 * w_sh_a[j])));
    Csh_a[j, k] = C > cmin ? C : 0.0;
  }
  for (i in 1:N_a) {
    real lo_i = x_obs_a[i] - 5 * sig_a[i];
    real dmt = (10 * sig_a[i]) / (Nint - 1.0);
    for (g in 1:Nint) {
      mt_a[i, g] = lo_i + (g - 1) * dmt;
      Cobj_a[i, g] = 0.5 * (1 + erf((mt_a[i, g] - mlim_obj_a[i] - d50_obj_a[i])
                                    / (sqrt2 * w_obj_a[i])));
    }
  }
  for (k in 1:Ng) xg_b[k] = xlo_b + (k - 1) * dx_b;
  for (j in 1:Nsh_b) for (k in 1:Ng) {
    real C = 0.5 * (1 + erf((xg_b[k] - mlim_sh_b[j] - d50_sh_b[j]) / (sqrt2 * w_sh_b[j])));
    Csh_b[j, k] = C > cmin ? C : 0.0;
  }
  for (i in 1:N_b) {
    real lo_i = x_obs_b[i] - 5 * sig_b[i];
    real dmt = (10 * sig_b[i]) / (Nint - 1.0);
    for (g in 1:Nint) {
      mt_b[i, g] = lo_i + (g - 1) * dmt;
      Cobj_b[i, g] = 0.5 * (1 + erf((mt_b[i, g] - mlim_obj_b[i] - d50_obj_b[i])
                                    / (sqrt2 * w_obj_b[i])));
    }
  }
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
}
model {
  ms ~ normal(14.13, 0.42);
  lp ~ normal(-3.96, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.02);    // beta pinned at Driver+22 value
  target += survey_ll_fixC(x_obs_a, sig_a, Cobj_a, mt_a, V_sh_a, Csh_a,
                           xg_a, dx_a, Nint, ms, lp, al, be);
  target += survey_ll_fixC(x_obs_b, sig_b, Cobj_b, mt_b, V_sh_b, Csh_b,
                           xg_b, dx_b, Nint, ms, lp, al, be);
}
"""

_STAN = {
    "simple": SIMPLE_CODE,
    "marg": MARG_CODE,
    "gama": GAMA_CODE,
    "combined": MARG_COMBINED_CODE,
    "marg_comp": MARG_COMP_CODE,
    "combined_comp": COMBINED_COMP_CODE,
}

# GAMA+SDSS (per-object completeness) PLUS REFLEX II as binned chi^2 anchoring the
# high-mass cutoff. Built from COMBINED_COMP_CODE by adding a REFLEX data block and
# a Gaussian term: -0.5 * sum(((log10 phi_MRP(m_r + dX) - y_r)/sig_r)^2). dX is a
# FIXED X-ray->dynamical mass offset (data, default 0) -- fitting it is degenerate
# with M*, so it is asserted, not fitted (set it explicitly if you have a value).
COMBINED_COMP_REFLEX_CODE = (
    COMBINED_COMP_CODE.replace(
        "  real xhi; int<lower=2> Ng; int<lower=2> Nint; real cmin;\n}",
        "  real xhi; int<lower=2> Ng; int<lower=2> Nint; real cmin;\n"
        "  int<lower=1> N_r; vector[N_r] m_r; vector[N_r] y_r;\n"
        "  vector<lower=0>[N_r] sig_r;\n"
        "  real dXa_mu; real dXa_sd; real dXb_mu; real dXb_sd;\n}",
    )
    .replace(
        "  real<lower=0.1, upper=2.0> be;\n}",
        "  real<lower=0.1, upper=2.0> be;\n"
        "  real dXa;   // REFLEX mass offset at logM=14 (dex)\n"
        "  real dXb;   // REFLEX offset slope d(offset)/d(logM) -- shallow M-sigma\n}",
    )
    .replace(
        "                           xg_b, dx_b, Nint, ms, lp, al, be);\n}",
        "                           xg_b, dx_b, Nint, ms, lp, al, be);\n\n"
        "  // REFLEX II binned chi^2, with a MASS-DEPENDENT X-ray->dynamical offset\n"
        "  // dX(M) = dXa + dXb*(m_r - 14): grounded in the shallow dynamical M-sigma\n"
        "  // relation (Han+15, Viola+15), so massive clusters shift more than groups.\n"
        "  dXa ~ normal(dXa_mu, dXa_sd);\n"
        "  dXb ~ normal(dXb_mu, dXb_sd);\n"
        "  {\n"
        "    real ln10r = log(10.0);\n"
        "    for (r in 1:N_r) {\n"
        "      real dXr = dXa + dXb * (m_r[r] - 14.0);\n"
        "      real u = (m_r[r] + dXr) - ms;\n"
        "      real phir = be * ln10r * pow(10, lp) * pow(10, (al + 1) * u)\n"
        "                  * exp(-pow(10, be * u));\n"
        "      real logphi = log10(phir > 1e-300 ? phir : 1e-300);\n"
        "      target += -0.5 * square((logphi - y_r[r]) / sig_r[r]);\n"
        "    }\n"
        "  }\n}",
    )
)

_STAN["combined_comp_reflex"] = COMBINED_COMP_REFLEX_CODE

# Same, but the mass-dependent offset is FIXED DATA (dXa_fix, dXb_fix) rather than
# a parameter. Pinning a parameter with a spike prior (sd~1e-6) destroys HMC's step
# size; removing it from the parameter space entirely samples at normal speed.
COMBINED_COMP_REFLEX_FIXED_CODE = COMBINED_COMP_CODE.replace(
    "  real xhi; int<lower=2> Ng; int<lower=2> Nint; real cmin;\n}",
    "  real xhi; int<lower=2> Ng; int<lower=2> Nint; real cmin;\n"
    "  int<lower=1> N_r; vector[N_r] m_r; vector[N_r] y_r;\n"
    "  vector<lower=0>[N_r] sig_r; real dXa_fix; real dXb_fix;\n}",
).replace(
    "                           xg_b, dx_b, Nint, ms, lp, al, be);\n}",
    "                           xg_b, dx_b, Nint, ms, lp, al, be);\n\n"
    "  // REFLEX II binned chi^2 with a FIXED mass offset dX(M)\n"
    "  {\n"
    "    real ln10r = log(10.0);\n"
    "    for (r in 1:N_r) {\n"
    "      real dXr = dXa_fix + dXb_fix * (m_r[r] - 14.0);\n"
    "      real u = (m_r[r] + dXr) - ms;\n"
    "      real phir = be * ln10r * pow(10, lp) * pow(10, (al + 1) * u)\n"
    "                  * exp(-pow(10, be * u));\n"
    "      real logphi = log10(phir > 1e-300 ? phir : 1e-300);\n"
    "      target += -0.5 * square((logphi - y_r[r]) / sig_r[r]);\n"
    "    }\n"
    "  }\n}",
)
_STAN["combined_comp_reflex_fixed"] = COMBINED_COMP_REFLEX_FIXED_CODE

# GAMA-only completeness model PLUS REFLEX II binned chi^2 with a mass-dependent
# X-ray->dynamical offset dX(M)=dXa+dXb*(M-14). Same REFLEX treatment as the
# combined model, but on the single-survey (marg_comp) likelihood -- for testing
# GAMA + REFLEX without SDSS.
MARG_COMP_REFLEX_CODE = MARG_COMP_CODE.replace(
    "  real cmin;\n}",
    "  real cmin;\n"
    "  int<lower=1> N_r; vector[N_r] m_r; vector[N_r] y_r;\n"
    "  vector<lower=0>[N_r] sig_r;\n"
    "  real dXa_mu; real dXa_sd; real dXb_mu; real dXb_sd;\n}",
).replace(
    "  real<lower=0.1, upper=2.0> be;\n}",
    "  real<lower=0.1, upper=2.0> be;\n"
    "  real dXa;   // REFLEX mass offset at logM=14 (dex)\n"
    "  real dXb;   // REFLEX offset slope d(offset)/d(logM)\n}",
)
# append the REFLEX chi^2 before the final closing brace of the model block
_i = MARG_COMP_REFLEX_CODE.rfind("}")
MARG_COMP_REFLEX_CODE = (
    MARG_COMP_REFLEX_CODE[:_i]
    + "\n  // REFLEX II binned chi^2 with mass-dependent X-ray->dynamical offset\n"
    "  dXa ~ normal(dXa_mu, dXa_sd);\n"
    "  dXb ~ normal(dXb_mu, dXb_sd);\n"
    "  {\n"
    "    real ln10r = log(10.0);\n"
    "    for (r in 1:N_r) {\n"
    "      real dXr = dXa + dXb * (m_r[r] - 14.0);\n"
    "      real u = (m_r[r] + dXr) - ms;\n"
    "      real phir = be * ln10r * pow(10, lp) * pow(10, (al + 1) * u)\n"
    "                  * exp(-pow(10, be * u));\n"
    "      real logphi = log10(phir > 1e-300 ? phir : 1e-300);\n"
    "      target += -0.5 * square((logphi - y_r[r]) / sig_r[r]);\n"
    "    }\n"
    "  }\n" + MARG_COMP_REFLEX_CODE[_i:]
)
_STAN["marg_comp_reflex"] = MARG_COMP_REFLEX_CODE

# marg_comp with the per-object mass errors scaled by a fitted hyperparameter
# s_scale: sigma_i -> s_scale * sigma_i. The reported sigmas come from run.R's
# multiplicity lookup and are asserted, not measured; since the Eddington
# deconvolution is driven entirely by sigma, this tests whether the data
# themselves prefer larger/smaller errors. Cobj/mt must be rebuilt inside the
# model (they depend on sigma), so this is slower than marg_comp.
MARG_COMP_SERR_CODE = r"""
data {
  int<lower=1> N;
  vector[N] x_obs;
  vector<lower=0>[N] sig;
  vector[N] mlim_obj;
  vector[N] d50_obj;
  vector<lower=0>[N] w_obj;
  int<lower=1> Nsh;
  vector[Nsh] V_sh;
  vector[Nsh] mlim_sh;
  vector[Nsh] d50_sh;
  vector<lower=0>[Nsh] w_sh;
  real xhi;
  int<lower=2> Ng;
  int<lower=2> Nint;
  real cmin;
}
transformed data {
  real ln10 = log(10.0);
  real sqrt2 = sqrt(2.0);
  real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
  real xlo = min(mlim_sh) - 2.5;
  real dx = (xhi - xlo) / (Ng - 1.0);
  vector[Ng] xg;
  matrix[Nsh, Ng] Csh;
  for (k in 1:Ng) xg[k] = xlo + (k - 1) * dx;
  for (j in 1:Nsh) {
    for (k in 1:Ng) {
      real C = 0.5 * (1 + erf((xg[k] - mlim_sh[j] - d50_sh[j]) / (sqrt2 * w_sh[j])));
      Csh[j, k] = C > cmin ? C : 0.0;
    }
  }
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
  real<lower=0.2, upper=4.0> s_scale;   // multiplies every reported sigma
}
model {
  ms ~ normal(14.13, 0.42);
  lp ~ normal(-3.96, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.02);
  s_scale ~ lognormal(0, 0.30);         // prior median 1, ~ +/-35 per cent

  vector[Ng] pg;
  for (k in 1:Ng) {
    real u = xg[k] - ms;
    pg[k] = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
  }
  real Lambda = 0;
  for (j in 1:Nsh) {
    real acc = 0;
    for (k in 1:Ng) {
      real term = pg[k] * Csh[j, k];
      acc += (k == 1 || k == Ng) ? 0.5 * term : term;
    }
    Lambda += V_sh[j] * acc * dx;
  }
  target += -Lambda;

  for (i in 1:N) {
    real si = s_scale * sig[i];
    real lo_i = x_obs[i] - 5 * si;
    real dmt = (10 * si) / (Nint - 1.0);
    real inv_s = 1.0 / si;
    real sm = 0;
    for (g in 1:Nint) {
      real mt = lo_i + (g - 1) * dmt;
      real u = mt - ms;
      real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
      real C = 0.5 * (1 + erf((mt - mlim_obj[i] - d50_obj[i]) / (sqrt2 * w_obj[i])));
      real zsc = (x_obs[i] - mt) * inv_s;
      real term = phi_g * C * exp(-0.5 * zsc * zsc);
      sm += (g == 1 || g == Nint) ? 0.5 * term : term;
    }
    sm *= dmt * inv_s * inv_sqrt2pi;
    target += (sm > 1e-300) ? log(sm) : -300;
  }
}
"""
_STAN["marg_comp_serr"] = MARG_COMP_SERR_CODE

# GAMA-only completeness + REFLEX with the offset as FIXED DATA (4 parameters).
MARG_COMP_REFLEX_FIXED_CODE = MARG_COMP_CODE.replace(
    "  real cmin;\n}",
    "  real cmin;\n"
    "  int<lower=1> N_r; vector[N_r] m_r; vector[N_r] y_r;\n"
    "  vector<lower=0>[N_r] sig_r; real dXa_fix; real dXb_fix;\n}",
)
_i2 = MARG_COMP_REFLEX_FIXED_CODE.rfind("}")
MARG_COMP_REFLEX_FIXED_CODE = (
    MARG_COMP_REFLEX_FIXED_CODE[:_i2]
    + "\n  // REFLEX II binned chi^2, fixed mass offset\n"
    "  {\n"
    "    real ln10r = log(10.0);\n"
    "    for (r in 1:N_r) {\n"
    "      real dXr = dXa_fix + dXb_fix * (m_r[r] - 14.0);\n"
    "      real u = (m_r[r] + dXr) - ms;\n"
    "      real phir = be * ln10r * pow(10, lp) * pow(10, (al + 1) * u)\n"
    "                  * exp(-pow(10, be * u));\n"
    "      real logphi = log10(phir > 1e-300 ? phir : 1e-300);\n"
    "      target += -0.5 * square((logphi - y_r[r]) / sig_r[r]);\n"
    "    }\n"
    "  }\n" + MARG_COMP_REFLEX_FIXED_CODE[_i2:]
)
_STAN["marg_comp_reflex_fixed"] = MARG_COMP_REFLEX_FIXED_CODE
_MODELS = {}

# Compile OUTSIDE any iCloud-synced tree (e.g. ~/Desktop, ~/Documents).
# A compiled Stan binary placed under a synced folder can be evicted mid-run,
# which shows up as "No such file or directory: .../mrp_marg" and retcode -1.
# Override with the HMF_STAN_DIR environment variable if you like.
STAN_BUILD_DIR = os.environ.get(
    "HMF_STAN_DIR", os.path.join(os.path.expanduser("~"), ".cache", "hmf_mrp_stan")
)
os.makedirs(STAN_BUILD_DIR, exist_ok=True)


def driver_prior(path="../data/hmfparams_gsr.csv", inflate=1.0):
    """Mean vector and covariance of Driver+22's GSR chains, for use as a
    multivariate-normal prior that preserves the strong parameter covariances
    (rho(M*,phi*) = -0.97) that independent Gaussians discard.

    NOTE this prior is NOT independent of the data: Driver's GSR posterior was
    derived from GAMA+SDSS+REFLEX. `inflate` widens the covariance (x inflate^2)
    to partially offset double-counting."""
    d = np.genfromtxt(path, delimiter=",", skip_header=1)
    d = d[np.isfinite(d).all(axis=1)]
    mu = d.mean(axis=0)
    cov = np.cov(d.T) * float(inflate) ** 2
    return mu, cov, d


def _with_driver_prior(code):
    """Swap a model's independent-Gaussian prior block for a multivariate normal."""
    old = (
        "  ms ~ normal(14.13, 0.42);\n"
        "  lp ~ normal(-3.96, 0.69);\n"
        "  al ~ normal(-1.68, 0.22);\n"
    )
    if old not in code:
        raise ValueError("prior block not found")
    code = code.replace(
        old, "  [ms, lp, al, be]' ~ multi_normal(prior_mu, prior_Sigma);\n", 1
    )
    # drop the now-duplicated beta prior line that followed
    code = re.sub(r"\n\s*be ~ normal\(0\.63, 0\.\d+\);[^\n]*", "", code, count=1)
    # declare the prior data just before the close of the data block
    i = code.index("data {")
    j, depth = i + len("data {"), 1
    while depth:
        if code[j] == "{":
            depth += 1
        elif code[j] == "}":
            depth -= 1
        j += 1
    code = (
        code[: j - 1]
        + "  vector[4] prior_mu;\n  matrix[4, 4] prior_Sigma;\n"
        + code[j - 1 :]
    )
    return code


def apply_driver_prior(kind, data):
    """If enabled, switch to the multivariate-prior variant and attach its data."""
    if not USE_DRIVER_PRIOR or kind.endswith("_dp"):
        return kind, data
    mu, cov, _ = driver_prior(DRIVER_PRIOR_PATH, DRIVER_PRIOR_INFLATE)
    data = dict(
        data,
        prior_mu=[float(v) for v in mu],
        prior_Sigma=[[float(v) for v in row] for row in cov],
    )
    print(
        f"  [prior] Driver+22 GSR chains as multivariate normal "
        f"(inflate={DRIVER_PRIOR_INFLATE}); mu={np.round(mu, 3).tolist()}"
    )
    return kind + "_dp", data


def get_model(kind):
    """Compile a Stan model once, in a non-synced build dir (cmdstanpy caches)."""
    if kind.endswith("_dp") and kind not in _STAN:
        _STAN[kind] = _with_driver_prior(_STAN[kind[:-3]])
    if kind not in _MODELS:
        from cmdstanpy import CmdStanModel

        path = os.path.join(STAN_BUILD_DIR, f"mrp_{kind}.stan")
        with open(path, "w") as f:
            f.write(_STAN[kind])
        model = CmdStanModel(stan_file=path)
        exe = getattr(model, "exe_file", None)
        if exe and not os.path.exists(exe):
            raise RuntimeError(f"Stan exe missing right after compile: {exe}")
        print(f"  [{kind}] compiled -> {STAN_BUILD_DIR}")
        _MODELS[kind] = model
    return _MODELS[kind]


def build_data(kind, x_fit, sig_fit, mlim_sh, V_sh, sig_sh=None, mlim_obj=None):
    """Assemble the Stan data dict for the chosen model."""
    data = dict(
        N=int(np.asarray(x_fit).size),
        x_obs=np.asarray(x_fit, float),
        Nsh=int(np.asarray(V_sh).size),
        V_sh=np.asarray(V_sh, float),
        mlim_sh=np.asarray(mlim_sh, float),
        xhi=float(XHI),
        Ng=int(NG),
    )
    if kind == "marg":
        if sig_fit is None or sig_sh is None:
            raise ValueError(
                "marg model needs sig_fit (per object) and sig_sh (per shell)"
            )
        data["sig"] = np.asarray(sig_fit, float)
        data["sig_sh"] = np.asarray(sig_sh, float)
        data["Nint"] = int(NINT)
    elif kind == "gama":
        # verbatim R port: needs per-object sigma AND per-object mlim
        if sig_fit is None or mlim_obj is None:
            raise ValueError("gama model needs sig_fit and mlim_obj (per-object mlim)")
        data["sig"] = np.asarray(sig_fit, float)
        data["mlim"] = np.asarray(mlim_obj, float)
        data["Nint"] = 100  # match run.R
    return data


def run_stan(
    kind,
    data,
    chains=4,
    warmup=1500,
    sampling=1500,
    adapt_delta=0.95,
    max_treedepth=10,
    seed=42,
    show_progress=True,
    output_dir=None,
):
    """MAP (optimize) + posterior (sample) for the chosen model.
    Returns (map_par, flat) with flat an (Ndraws, 4) array in PARAMS order."""
    kind, data = apply_driver_prior(kind, data)
    model = get_model(kind)
    if kind == "marg":
        ss = np.asarray(data["sig_sh"])
        print(
            f"  [marg] sig_sh: min={ss.min():.2f} max={ss.max():.2f} med={np.median(ss):.2f}"
        )

    # ---- MAP: multi-start LBFGS (Stan = exact gradients, unconstrained scale) ----
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(5):
        init = dict(
            ms=float(rng.normal(14, 0.3)),
            lp=float(rng.normal(-3.5, 0.3)),
            al=float(rng.normal(-1.3, 0.2)),
            be=float(rng.uniform(0.3, 0.7)),
        )
        try:
            opt = model.optimize(
                data=data, inits=init, algorithm="lbfgs", iter=20000, seed=seed
            )
            val = float(opt.optimized_params_dict["lp__"])
            if best is None or val > best[0]:
                best = (
                    val,
                    np.array([float(opt.optimized_params_dict[p]) for p in PARAMS]),
                )
        except Exception as e:
            print("  optimize trial failed:", e)
    map_par = best[1] if best is not None else np.array([14.0, -3.5, -1.3, 0.5])
    print(f"  [{kind}] MAP:", dict(zip(PARAMS, np.round(map_par, 3))))

    # ---- MCMC: NUTS, chains dispersed around the MAP ----
    inits = [
        dict(
            ms=float(map_par[0] + rng.normal(0, 0.03)),
            lp=float(map_par[1] + rng.normal(0, 0.03)),
            al=float(map_par[2] + rng.normal(0, 0.03)),
            be=float(np.clip(map_par[3] + rng.normal(0, 0.03), 0.12, 1.9)),
        )
        for _ in range(chains)
    ]

    fit = model.sample(
        data=data,
        chains=chains,
        parallel_chains=chains,
        iter_warmup=warmup,
        iter_sampling=sampling,
        adapt_delta=adapt_delta,
        max_treedepth=max_treedepth,
        inits=inits,
        seed=seed,
        show_progress=show_progress,
        output_dir=output_dir,
    )

    flat = np.column_stack([fit.stan_variable(p) for p in PARAMS])

    # diagnostics (column names vary slightly across cmdstanpy versions)
    try:
        summ = fit.summary()
        rhat = float(summ.loc[PARAMS, "R_hat"].max())
        ess_col = next((c for c in ("ESS_bulk", "N_Eff") if c in summ.columns), None)
        ess = float(summ.loc[PARAMS, ess_col].min()) if ess_col else float("nan")
    except Exception:
        rhat, ess = float("nan"), float("nan")
    try:
        ndiv = int(np.sum(fit.method_variables()["divergent__"]))
    except Exception:
        ndiv = -1
    # treedepth saturation: high fraction => NUTS fighting the geometry (slow)
    try:
        td = fit.method_variables()["treedepth__"]
        td_frac = float(np.mean(td >= max_treedepth))
        td_max = int(np.max(td))
    except Exception:
        td_frac, td_max = float("nan"), -1
    print(
        f"  [{kind}] Rhat={rhat:.3f}  min ESS={ess:.0f}  divergences={ndiv}  "
        f"treedepth>={max_treedepth}: {td_frac:.0%} (max {td_max})"
    )
    if "_serr" in kind:
        try:
            ss = fit.stan_variable("s_scale")
            print(
                f"  [{kind}] s_scale = {np.median(ss):.3f} "
                f"+{np.percentile(ss, 84) - np.median(ss):.3f}/"
                f"-{np.median(ss) - np.percentile(ss, 16):.3f}  "
                f"(1.0 = reported sigmas correct; prior lognormal(0, 0.30))"
            )
            np.savetxt(
                "s_scale_draws.csv", ss, delimiter=",", header="s_scale", comments=""
            )
        except Exception as e:
            print(f"  [s_scale unavailable: {e}]")
    return map_par, flat


def sigma_eff_per_shell(z, m_obs, sigma, mlim_sh, nsh=NSH, zmin=ZMIN, zmax=ZLIMIT):
    """Representative boundary sigma per shell for the marginalised Lambda.
    Uses the median sigma of detected groups NEAR the limit in each shell
    (that is the scatter scale that smooths the cut), with fallbacks."""
    z_edges = np.linspace(zmin, zmax, nsh + 1)
    sig_global = float(np.median(sigma))
    out = np.full(nsh, sig_global)
    for j in range(nsh):
        in_sh = (z >= z_edges[j]) & (z < z_edges[j + 1])
        near = in_sh & (m_obs >= mlim_sh[j] - 0.25) & (m_obs <= mlim_sh[j] + 0.75)
        if near.sum() >= 10:
            out[j] = float(np.median(sigma[near]))
        elif in_sh.sum() >= 5:
            out[j] = float(np.median(sigma[in_sh]))
    # Floor the boundary width: a very small sig_sh makes Phi((m-mlim)/sig_sh)
    # a near-vertical step, which wrecks the HMC geometry (treedepth blow-up).
    # The floor is well below the near-limit group errors (~0.28 dex) so it
    # rarely binds and does not over-smooth the boundary.
    return np.maximum(out, SIG_SH_FLOOR)


def summarise(flat):
    med = np.median(flat, axis=0)
    sd = np.std(flat, axis=0)
    q16, q84 = np.percentile(flat, [16, 84], axis=0)
    tv = np.array([TRUE[p] for p in PARAMS])
    print(
        f"\n  {'param':9s} {'true':>8s} {'median':>9s} {'sd':>7s} "
        f"{'q16':>8s} {'q84':>8s} {'bias(sd)':>9s}"
    )
    for i, p in enumerate(PARAMS):
        bias = (med[i] - tv[i]) / sd[i]
        print(
            f"  {p:9s} {tv[i]:8.3f} {med[i]:9.3f} {sd[i]:7.3f} "
            f"{q16[i]:8.3f} {q84[i]:8.3f} {bias:+9.2f}"
        )
    ac = corrected_alpha(flat)
    acm = float(np.median(ac))
    print(
        f"  {'al_corr':9s} {TRUE['al']:8.3f} {acm:9.3f} {np.std(ac):7.3f} "
        f"{np.percentile(ac, 16):8.3f} {np.percentile(ac, 84):8.3f} "
        f"{(acm - TRUE['al']) / np.std(ac):+9.2f}   <- mock bias-corrected"
    )
    print(
        f"  [alpha bias at M*={med[0]:.2f} is {float(alpha_bias(med[0])):+.3f} dex; "
        f"single-realisation calibration, no error bar yet]"
    )
    return dict(
        median=med, sd=sd, q16=q16, q84=q84, al_corr=acm, al_corr_sd=float(np.std(ac))
    )


# ----------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------
def plot_recovery(
    flat,
    z_obs,
    m_obs,
    x_fit,
    mlim_func,
    Vsurvey,
    turn_pts=None,
    fname="recovery_simple.pdf",
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    med = np.median(flat, axis=0)
    tv = np.array([TRUE[p] for p in PARAMS])
    xfit = np.linspace(10, 16.5, 500)
    z_plot = np.linspace(ZMIN, ZLIMIT, 200)

    bw = 0.2
    edges = np.arange(10, 16 + 1e-9, bw)
    mids = 0.5 * (edges[1:] + edges[:-1])
    c_all, _ = np.histogram(m_obs, bins=edges)
    c_fit, _ = np.histogram(x_fit, bins=edges)
    phi_all = c_all / (Vsurvey * bw)
    phi_fit = c_fit / (Vsurvey * bw)

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))

    a = ax[0, 0]
    ok = c_all >= 5
    a.plot(
        mids[ok], np.log10(phi_all[ok]), "o", color="grey", ms=4, label="all detected"
    )
    ok = c_fit >= 5
    a.plot(
        mids[ok],
        np.log10(phi_fit[ok]),
        "o",
        color="darkgreen",
        ms=6,
        label="above mlim",
    )
    idx = np.random.default_rng(0).choice(
        len(flat), size=min(200, len(flat)), replace=False
    )
    for i in idx:
        y = np.log10(np.maximum(mrp_phi(xfit, *flat[i]), 1e-30))
        a.plot(xfit, y, color="red", alpha=0.02)
    a.plot(
        xfit,
        np.log10(np.maximum(mrp_phi(xfit, *med), 1e-30)),
        "r-",
        lw=2,
        label="MCMC median",
    )
    a.plot(
        xfit,
        np.log10(np.maximum(mrp_phi(xfit, *tv), 1e-30)),
        "b--",
        lw=2,
        label="truth",
    )
    a.set(
        xlim=(11, 16),
        ylim=(-8, -1),
        xlabel=r"$\log_{10}(M)$",
        ylabel=r"$\log_{10}\phi$",
        title="HMF recovery",
    )
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    a = ax[0, 1]
    a.scatter(z_obs, m_obs, s=2, alpha=0.2, color="steelblue")
    a.plot(z_plot, mlim_func(z_plot), "r-", lw=2, label="mlim(z)")
    if turn_pts is not None:
        zb, tb = turn_pts
        ok = np.isfinite(tb)
        a.plot(zb[ok], tb[ok], "cs", ms=5, label="turnover")
    a.set(
        xlim=(ZMIN, ZLIMIT),
        ylim=(10, 15.5),
        xlabel="z",
        ylabel=r"$\log_{10}(M)$",
        title="mass-redshift + mlim(z)",
    )
    a.legend(fontsize=7)
    a.grid(alpha=0.3)

    labels = [r"$M_*$", r"$\log\phi_*$", r"$\alpha$", r"$\beta$"]
    prior_mu = [14.13, -3.96, -1.68, 0.63]  # must match the Stan model priors
    prior_sd = [0.42, 0.69, 0.22, 0.02]
    for k, (i, j) in enumerate([(1, 0), (1, 1), (1, 2), (0, 2)]):
        a = ax[i, j]
        a.hist(flat[:, k], bins=40, color="steelblue", density=True)
        xs = np.linspace(*a.get_xlim(), 200)
        prior = np.exp(-0.5 * ((xs - prior_mu[k]) / prior_sd[k]) ** 2) / (
            prior_sd[k] * np.sqrt(2 * np.pi)
        )
        a.plot(xs, prior, color="green", lw=1.5, ls=":", label="prior")
        a.axvline(tv[k], color="blue", lw=2, ls="--", label="truth")
        a.axvline(med[k], color="red", lw=2, label="median")
        a.set(title=labels[k])
        a.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")
    return fname


# ----------------------------------------------------------
# 7. Coverage loop  (step 1 validation)
# ----------------------------------------------------------
# Wraps the GAMA-mock in N realisations to ask: are the credible intervals
# honest?  The abundance-matched TRUTH and the detected-group set are fixed
# (deterministic selection on one mock volume); what is re-drawn each
# realisation is the Gaussian MASS ERROR, which then propagates through the
# turnover mlim(z) and the above-mlim cut -- i.e. the model-relevant noise.
#
# LIMITATION (be honest about it): because there is a single mock volume,
# the Poisson count of groups is NOT re-sampled, so logphi* coverage here is
# CONDITIONAL on that fixed count.  Fuller coverage needs independent volumes
# or an injected count/photometric-scatter term -- a later refinement.
def _prepare_fixed():
    """Deterministic part of the GAMA-mock, run once."""
    groups, galaxies = load_catalogues(DATA_DIR)
    sky_frac = (
        sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = abundance_match(groups, sky_frac)
    gv = gama_select(gv, galaxies)
    det = gv[gv["detected"]].copy()
    z_mids, V_sh = shell_volumes(sky_frac)
    print(f"  fixed detected-group set: {len(det)} groups")
    return dict(
        z=det["zcos"].values,
        m_true=det["log_mass_am"].values,
        n_gama=det["n_gama"].values,
        sky_frac=sky_frac,
        Vsurvey=Vsurvey,
        z_mids=z_mids,
        V_sh=V_sh,
    )


def run_coverage(
    model_kind="marg",
    n_real=20,
    seed0=1000,
    chains=4,
    warmup=500,
    sampling=500,
    checkpoint=None,
):
    # Fewer iterations than a single headline fit: coverage only needs a stable
    # median + 16/84 interval per realisation (ESS ~150 suffices), and we run
    # many of them. The marginalised model's per-object integral is the cost, so
    # this keeps the full loop to ~hours rather than ~overnight.
    if checkpoint is None:
        checkpoint = f"coverage_{model_kind}.csv"
    print("=" * 60)
    print(f"  COVERAGE LOOP [{model_kind}]: {n_real} realisations")
    print("=" * 60)
    fixed = _prepare_fixed()
    z, m_true, n_gama = fixed["z"], fixed["m_true"], fixed["n_gama"]
    z_mids, V_sh = fixed["z_mids"], fixed["V_sh"]
    sigma = sigma_from_nfof(n_gama)  # error SIZE fixed per group; the DRAW varies

    # resume from checkpoint if present
    rows, done = [], set()
    if os.path.exists(checkpoint):
        prev = pd.read_csv(checkpoint)
        rows = prev.to_dict("records")
        done = set(prev["real"].astype(int))
        print(f"  resuming: {len(done)} realisations already on disk")

    for r in range(n_real):
        if r in done:
            continue
        rng = np.random.default_rng(seed0 + r)
        m_obs = m_true + rng.normal(0, sigma)

        try:
            mlim_func, _, tkind, _ = turnover_mlim(z, m_obs)
        except RuntimeError as e:
            print(f"  [real {r:02d}] mlim failed ({e}); skipped")
            continue
        mlim_sh = mlim_func(z_mids)
        if model_kind in ("marg_comp", "marg_comp_serr"):
            data, keep = prep_comp(z, m_obs, sigma, mlim_func, z_mids, mlim_sh, V_sh)
            x_fit = m_obs[keep]
            print(
                f"\n[real {r:02d}] N_kept={x_fit.size} (C>{CMIN})  mlim[{tkind}] "
                f"{mlim_func(ZMIN):.2f}->{mlim_func(ZLIMIT):.2f}"
            )
        else:
            above = m_obs > mlim_func(z)
            x_fit = m_obs[above]
            sig_fit = sigma[above]
            sig_sh = (
                sigma_eff_per_shell(z, m_obs, sigma, mlim_sh)
                if model_kind == "marg"
                else None
            )
            print(
                f"\n[real {r:02d}] N_fit={x_fit.size}  mlim[{tkind}] "
                f"{mlim_func(ZMIN):.2f}->{mlim_func(ZLIMIT):.2f}"
            )
            data = build_data(model_kind, x_fit, sig_fit, mlim_sh, V_sh, sig_sh=sig_sh)

        _, flat = run_stan(
            model_kind,
            data,
            chains=chains,
            warmup=warmup,
            sampling=sampling,
            seed=seed0 + r,
            show_progress=False,
        )

        med = np.median(flat, axis=0)
        sd = np.std(flat, axis=0)
        q025, q16, q84, q975 = np.percentile(flat, [2.5, 16, 84, 97.5], axis=0)
        row = dict(real=r, N_fit=int(x_fit.size))
        for i, p in enumerate(PARAMS):
            row.update(
                {
                    f"{p}_med": med[i],
                    f"{p}_sd": sd[i],
                    f"{p}_q025": q025[i],
                    f"{p}_q16": q16[i],
                    f"{p}_q84": q84[i],
                    f"{p}_q975": q975[i],
                }
            )
        rows.append(row)
        pd.DataFrame(rows).to_csv(checkpoint, index=False)  # checkpoint every iter

    df = pd.DataFrame(rows).sort_values("real").reset_index(drop=True)
    report_coverage(df)
    plot_coverage(df, checkpoint.replace(".csv", ".pdf"))
    return df


def report_coverage(df):
    # prior widths (must match the Stan model priors) -- to flag which
    # parameters are data-constrained vs prior-driven on this sample.
    prior_sd = {"ms": 0.42, "lp": 0.69, "al": 0.22, "be": 0.02}
    print("\n" + "=" * 78)
    print(f"  COVERAGE SUMMARY over {len(df)} realisations")
    print("=" * 78)
    print(
        f"  {'param':6s} {'true':>8s} {'mean_med':>9s} {'bias_dex':>9s} "
        f"{'bias_sd':>8s} {'cov68':>7s} {'cov95':>7s} {'post/prior':>10s} {'constrained':>12s}"
    )
    for p in PARAMS:
        t = TRUE[p]
        med = df[f"{p}_med"].values
        sd = df[f"{p}_sd"].values
        cov68 = np.mean((df[f"{p}_q16"].values <= t) & (t <= df[f"{p}_q84"].values))
        cov95 = np.mean((df[f"{p}_q025"].values <= t) & (t <= df[f"{p}_q975"].values))
        ratio = np.mean(sd) / prior_sd[p]  # posterior width / prior width
        tag = (
            "data-driven"
            if ratio < 0.6
            else ("prior-driven" if ratio > 0.9 else "mixed")
        )
        print(
            f"  {p:6s} {t:8.3f} {np.mean(med):9.3f} {np.mean(med - t):+9.3f} "
            f"{np.mean((med - t) / sd):+8.2f} {cov68:7.0%} {cov95:7.0%} "
            f"{ratio:10.2f} {tag:>12s}"
        )
    print("\n  cov68 ~ 68%, cov95 ~ 95%, |bias_sd| small  ->  calibrated & usable.")
    print("  'post/prior' ~ 1 (prior-driven): the data did NOT constrain that")
    print("  parameter -- report it as prior-informed, not measured (expected for")
    print("  M* and the knee, which this sample does not reach).  'data-driven'")
    print("  with good coverage -> a genuine measurement (expect this for the")
    print("  high-mass/cutoff behaviour).")


def plot_coverage(df, fname="coverage_results.pdf"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = dict(ms=r"$M_*$", lp=r"$\log\phi_*$", al=r"$\alpha$", be=r"$\beta$")
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    for a, p in zip(ax.flat, PARAMS):
        t = TRUE[p]
        x = df["real"].values
        med = df[f"{p}_med"].values
        lo = med - df[f"{p}_q16"].values
        hi = df[f"{p}_q84"].values - med
        inside = (df[f"{p}_q16"].values <= t) & (t <= df[f"{p}_q84"].values)
        for xi, mi, loi, hii, ins in zip(x, med, lo, hi, inside):
            a.errorbar(
                xi,
                mi,
                yerr=[[loi], [hii]],
                fmt="o",
                ms=3,
                capsize=2,
                color="steelblue" if ins else "crimson",
            )
        a.axhline(t, color="k", ls="--", lw=1.5)
        a.set(
            title=f"{labels[p]}   68% coverage = {inside.mean():.0%}",
            xlabel="realisation",
            ylabel=labels[p],
        )
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")


# ----------------------------------------------------------
# 9. Drivers
# ----------------------------------------------------------
def build_comp_data(
    x_fit, sig_fit, mlim_obj, d50_obj, w_obj, mlim_sh, V_sh, d50_sh, w_sh, Nint=61
):
    """Stan data dict for the completeness (marg_comp) model."""
    f = lambda v: np.asarray(v, float)
    return dict(
        N=int(f(x_fit).size),
        x_obs=f(x_fit),
        sig=f(sig_fit),
        mlim_obj=f(mlim_obj),
        d50_obj=f(d50_obj),
        w_obj=f(w_obj),
        Nsh=int(f(V_sh).size),
        V_sh=f(V_sh),
        mlim_sh=f(mlim_sh),
        d50_sh=f(d50_sh),
        w_sh=f(w_sh),
        xhi=float(XHI),
        Ng=int(NG),
        Nint=int(Nint),
        cmin=float(CMIN),
    )


def prep_comp(z_obs, m_obs, sigma, mlim_func, z_mids, mlim_sh, V_sh):
    """Build the completeness-model data: z-dependent ramp (D50(z), w(z)
    interpolated from the mock measurement) per object and per shell, keep all
    detected groups whose own completeness exceeds CMIN (no mlim cut). Returns
    (data, keep_mask)."""
    from scipy.stats import norm

    mlim_obj = mlim_func(z_obs)
    d50_obj = np.interp(z_obs, COMP_Z_PTS, COMP_D50_PTS)
    w_obj = np.interp(z_obs, COMP_Z_PTS, COMP_W_PTS)
    C_obj = norm.cdf((m_obs - mlim_obj - d50_obj) / w_obj)  # = 0.5(1+erf(./sqrt2 w))
    keep = C_obj > CMIN
    d50_sh = np.interp(z_mids, COMP_Z_PTS, COMP_D50_PTS)
    w_sh = np.interp(z_mids, COMP_Z_PTS, COMP_W_PTS)
    data = build_comp_data(
        m_obs[keep],
        sigma[keep],
        mlim_obj[keep],
        d50_obj[keep],
        w_obj[keep],
        mlim_sh,
        V_sh,
        d50_sh,
        w_sh,
    )
    return data, keep


def run_real_pipeline(model_kind="marg"):
    print("Reading catalogues ...")
    groups, galaxies = load_catalogues(DATA_DIR)
    sky_frac = (
        sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )

    print("Abundance matching ...")
    gv, Vsurvey = abundance_match(groups, sky_frac)

    print(f"Selection ... [{SEL_COL} < {MAG_LIMIT}, >= {MULTI} members]")
    gv = gama_select(gv, galaxies)
    n_det = int(gv["detected"].sum())
    print(f"  detected: {n_det} / {len(gv)}")
    if n_det == 0:
        raise RuntimeError(
            "No groups detected -- check the selection cuts in gama_select()."
        )

    z_obs, m_obs, sigma_obs, nfof = add_mass_errors(gv, rng_global)

    print("Turnover mlim(z) ...")
    mlim_func, coefs, kind, turn_pts = turnover_mlim(z_obs, m_obs)
    print(
        f"  mlim(z) [{kind}]: mlim({ZMIN})={mlim_func(ZMIN):.2f} "
        f"mlim({ZLIMIT})={mlim_func(ZLIMIT):.2f}"
    )

    z_mids, V_sh = shell_volumes(sky_frac)
    mlim_sh = mlim_func(z_mids)

    if model_kind in ("marg_comp", "marg_comp_serr"):
        data, keep = prep_comp(
            z_obs, m_obs, sigma_obs, mlim_func, z_mids, mlim_sh, V_sh
        )
        x_fit = m_obs[keep]
        print(
            f"  completeness model: kept {x_fit.size} / {m_obs.size} "
            f"(C > {CMIN}); no mlim cut"
        )
    else:
        above = m_obs > mlim_func(z_obs)
        x_fit = m_obs[above]
        sig_fit = sigma_obs[above]
        print(
            f"  N above mlim: {x_fit.size} / {m_obs.size} ({100 * x_fit.size / m_obs.size:.1f}%)"
        )
        sig_sh = (
            sigma_eff_per_shell(z_obs, m_obs, sigma_obs, mlim_sh)
            if model_kind == "marg"
            else None
        )
        data = build_data(model_kind, x_fit, sig_fit, mlim_sh, V_sh, sig_sh=sig_sh)

    print(f"\nFitting [{model_kind}] (cmdstanpy) ...")
    map_par, flat = run_stan(model_kind, data)
    res = summarise(flat)
    plot_recovery(
        flat,
        z_obs,
        m_obs,
        x_fit,
        mlim_func,
        Vsurvey,
        turn_pts=turn_pts,
        fname=f"recovery_{model_kind}.pdf",
    )
    emit_publication(
        flat,
        {"mock": dict(x_fit=x_fit, Vsurvey=Vsurvey)},
        tag=f"mock_{model_kind}",
        title="Mock HMF",
    )
    return res


def load_real_gama(fits_path):
    """Read the GAMA G3C group catalogue and build observed log-masses exactly
    as run.R does: A=13.9 dynamical mass, multiplicity error model, and the
    empirical masscorr(Nfof) calibration. Returns (log_mass, sigma, z, Nfof).
    Column names/units follow run.R; adjust here if the FITS schema differs."""
    from astropy.io import fits as afits

    parsec, Gnewton, msol = 3.0857e16, 6.67408e-11, 1.988e30
    with afits.open(fits_path) as hdul:
        t = hdul[1].data
    Nfof = np.asarray(t["Nfof"], float)
    Zfof = np.asarray(t["Zfof"], float)
    MassAfunc = np.asarray(t["MassAfunc"], float)
    VelDisp = np.asarray(t["VelDisp"], float)
    Rad50 = np.asarray(t["Rad50"], float)
    IterCenDec = np.asarray(t["IterCenDec"], float)

    sel = (
        (Nfof > MULTI - 1)
        & (Zfof < ZLIMIT)
        & (Zfof > ZMIN)
        & (MassAfunc > 1e1)
        & (IterCenDec > -3.5)
    )
    Nfof, Zfof, VelDisp, Rad50 = Nfof[sel], Zfof[sel], VelDisp[sel], Rad50[sel]

    # A=13.9 dynamical mass (Msun), h-scaled as in run.R
    mymass = (
        A_SCALE
        * (VelDisp * 1000) ** 2
        * Rad50
        * parsec
        * 1e6
        / (Gnewton * msol)
        * (100 / H0)
    )

    # multiplicity error model (run.R table, xx = 3..22; NA -> 0.03; floor 0.1)
    xx = np.arange(3, 23)
    yy = np.array(
        [
            0.68389355,
            0.38719116,
            0.40325591,
            0.32696735,
            0.27680685,
            0.24018684,
            0.20226682,
            0.18645475,
            0.17437005,
            0.14271506,
            0.13922450,
            0.13482418,
            0.13741619,
            0.11715141,
            0.12134983,
            0.10078830,
            0.09944761,
            0.09913166,
            0.08590223,
            0.07588408,
        ]
    )
    err = np.interp(Nfof, xx, yy, left=np.nan, right=np.nan)
    err = np.where(np.isfinite(err), err, 0.03)
    err = np.where(err < 0.1, 0.1, err)

    # empirical mass corrections indexed by Nfof (run.R, 1-based; out-of-range -> 0)
    masscorr = np.array(
        [
            0.0,
            0.0,
            -2.672595e-01,
            -1.513503e-01,
            -1.259069e-01,
            -9.006064e-02,
            -5.466009e-02,
            -6.666895e-02,
            -1.988694e-02,
            -2.439581e-02,
            -2.067060e-02,
            -1.812964e-02,
            -1.556899e-02,
            -1.313664e-02,
            -1.743112e-02,
            -7.965513e-03,
            -1.257178e-02,
            -7.064037e-03,
            -3.963656e-03,
            -1.271533e-02,
            -2.664687e-03,
            -1.691287e-03,
        ]
    )
    idx = Nfof.astype(int) - 1
    mc = np.where(
        (idx >= 0) & (idx < masscorr.size),
        masscorr[np.clip(idx, 0, masscorr.size - 1)],
        0.0,
    )
    log_mass = np.log10(mymass / 10**mc)

    good = np.isfinite(log_mass) & (log_mass > 10) & (log_mass < 17) & np.isfinite(err)
    return log_mass[good], err[good], Zfof[good], Nfof[good].astype(int)


def run_real_gama(
    fits_path,
    sky_area_deg2_val=179.92,
    model_kind="marg",
    reflex=False,
    reflex_dX=0.0,
    fix_offset=False,
    reflex_mmin=None,
):
    """Fit the MRP to the REAL GAMA group catalogue using OUR developed model
    ('marg' by default, with the current tight cutoff/slope priors; 'gama' =
    verbatim R port, kept only for a pure port-check). load_real_gama does the
    run.R data prep (A=13.9 masses, multiplicity error model); the FIT is our
    marginalised model. 'bias vs truth' and the plot 'truth' lines are Driver+22,
    so read them as 'offset vs Driver+22' -- with tight priors, be/al are
    prior-set near Driver by construction; M* and logphi* are the measurements."""
    print(f"Reading GAMA catalogue: {fits_path}")
    log_mass, sigma, z, nfof = load_real_gama(fits_path)
    sky_frac = sky_area_deg2_val * (np.pi / 180) ** 2 / (4 * np.pi)
    print(
        f"  N groups: {log_mass.size}   mass {log_mass.min():.2f}..{log_mass.max():.2f} "
        f"(med {np.median(log_mass):.2f})   sigma med {np.median(sigma):.2f}"
    )

    print("Turnover mlim(z) ...")
    mlim_func, coefs, tkind, turn_pts = turnover_mlim(z, log_mass)
    print(
        f"  mlim(z) [{tkind}]: mlim({ZMIN})={mlim_func(ZMIN):.2f} "
        f"mlim({ZLIMIT})={mlim_func(ZLIMIT):.2f}"
    )

    z_mids, V_sh = shell_volumes(sky_frac)
    Vsurvey = float(V_sh.sum())
    mlim_sh = mlim_func(z_mids)

    mlim_per = mlim_func(z)

    if model_kind in ("marg_comp", "marg_comp_serr"):
        # Completeness ramp from the GAMA-selected mock (same mag limit, same
        # >=MULTI members, same group finder) -- i.e. injection-recovery applied
        # to real GAMA. mlim(z) is derived from the REAL data above.
        data, keep = prep_comp(z, log_mass, sigma, mlim_func, z_mids, mlim_sh, V_sh)
        x_fit = log_mass[keep]
        print(
            f"  completeness model: kept {x_fit.size} / {log_mass.size} "
            f"(C > {CMIN}); no mlim cut; ramp from GAMA mock"
        )
        if reflex:
            comp = _load_comparison("../data")
            rf = comp.get("REFLEX II (Böhringer+17)")
            if rf is None:
                print("  [reflex requested but reflex.csv not loaded -> without it]")
            else:
                sig_r = np.clip(
                    0.5 * (np.abs(rf["elo"]) + np.abs(rf["ehi"])), 0.03, None
                )
                fa = lambda v: np.asarray(v, float)
                rx, ry = fa(rf["x"]), fa(rf["y"])
                if reflex_mmin is not None:
                    kr = rx > float(reflex_mmin)
                    print(
                        f"  [REFLEX] mass cut logM > {reflex_mmin}: {int(kr.sum())}/{rx.size}"
                    )
                    rx, ry, sig_r = rx[kr], ry[kr], sig_r[kr]
                data.update(N_r=int(rx.size), m_r=rx, y_r=ry, sig_r=fa(sig_r))
                if fix_offset:
                    data.update(dXa_fix=float(reflex_dX), dXb_fix=0.0)
                    model_kind = "marg_comp_reflex_fixed"
                    print(
                        f"  [REFLEX] {rx.size} points, logM {rx.min():.2f}-{rx.max():.2f}, "
                        f"dX fixed at {reflex_dX:+.2f}"
                    )
                else:
                    data.update(
                        dXa_mu=float(reflex_dX), dXa_sd=0.25, dXb_mu=0.0, dXb_sd=0.3
                    )
                    model_kind = "marg_comp_reflex"
                    print(f"  [REFLEX] {rx.size} points, dX(M)=a+b(M-14) fitted")
    else:
        above = log_mass > mlim_per
        x_fit, sig_fit = log_mass[above], sigma[above]
        print(
            f"  N above mlim: {x_fit.size} / {log_mass.size} "
            f"({100 * x_fit.size / log_mass.size:.1f}%)"
        )
        if model_kind == "marg":
            sig_sh = sigma_eff_per_shell(z, log_mass, sigma, mlim_sh)
            data = build_data("marg", x_fit, sig_fit, mlim_sh, V_sh, sig_sh=sig_sh)
        elif model_kind == "gama":
            data = build_data(
                "gama", x_fit, sig_fit, mlim_sh, V_sh, mlim_obj=mlim_per[above]
            )
        else:
            data = build_data(model_kind, x_fit, sig_fit, mlim_sh, V_sh)

    print(f"\nFitting [{model_kind}] on real GAMA (cmdstanpy) ...")
    map_par, flat = run_stan(model_kind, data)
    res = summarise(flat)  # 'bias vs truth' = offset vs Driver+22
    if model_kind == "marg_comp":
        MSTAR_LCDM = 14.13  # anchor: Driver's M* (sits on the Murray+21 LCDM curve)
        A_draws = A_SCALE * 10 ** (MSTAR_LCDM - flat[:, 0])  # flat[:,0] = ms draws
        amed = np.median(A_draws)
        print(f"\n  === fitted mass calibration A (M* fixed to LCDM {MSTAR_LCDM}) ===")
        print(
            f"  A = {amed:.2f}  +{np.percentile(A_draws, 84) - amed:.2f}"
            f"/-{amed - np.percentile(A_draws, 16):.2f}   (A_assumed={A_SCALE})"
        )
        print(f"  vs Robotham+11 sim-calibrated 13.9, Driver variant 10, Zwicky 1.667")
        print(f"  [caveat: M* has a mild informative prior N(14.13,0.42); data-driven")
        print(f"   posterior ~0.1 dominates, so prior pull on A is small (~7%)]")
    plot_recovery(
        flat,
        z,
        log_mass,
        x_fit,
        mlim_func,
        Vsurvey,
        turn_pts=turn_pts,
        fname=f"recovery_gama_{model_kind}.pdf",
    )
    emit_publication(
        flat,
        {"GAMA": dict(x_fit=x_fit, Vsurvey=Vsurvey)},
        tag=f"gama_{model_kind}",
        title="GAMA HMF",
    )
    return res


def load_sdss_groups(
    parquet_path, zmin=ZMIN, zmax=ZLIMIT, mass_col="mass_proxy", A=None
):
    """Read a per-object SDSS group catalogue (sdss_groups.parquet).

    mass_col='mass_proxy' (default): apply the SAME A-scaling as GAMA, i.e.
      mass = A * mass_proxy, so both surveys share one mass definition. Assumes
      mass_proxy is the sigma^2 R / G combination in Msun (the same quantity
      GAMA scales). A scale-check is printed -- if log10 lands outside ~13-15,
      the units differ and mass_proxy needs GAMA's unit conversion.
    mass_col='estimated_mass': use the catalogue's own linear mass as-is.

    Per-object sigma reuses the GAMA multiplicity error model for now.
    Returns (log_mass, sigma, z, multiplicity)."""
    import pandas as pd

    df = pd.read_parquet(parquet_path)
    z = df["median_redshift"].values.astype(float)
    mult = df["multiplicity"].values.astype(float)

    if mass_col == "mass_proxy":
        if A is None:
            A = A_SCALE
        m_lin = A * df["mass_proxy"].values.astype(float)
    else:
        m_lin = df[mass_col].values.astype(float)

    good = (
        np.isfinite(m_lin)
        & (m_lin > 0)
        & np.isfinite(z)
        & (z > zmin)
        & (z < zmax)
        & (mult >= MULTI)
    )
    log_mass = np.log10(m_lin[good])
    z, mult = z[good], mult[good]
    sigma = sigma_from_nfof(mult)
    keep = np.isfinite(log_mass) & (log_mass > 10) & (log_mass < 17)

    # scale sanity check vs the catalogue's own estimated_mass
    if mass_col == "mass_proxy" and "estimated_mass" in df.columns:
        em = df["estimated_mass"].values.astype(float)[good][keep]
        em = em[np.isfinite(em) & (em > 0)]
        print(
            f"  [mass] {A} x mass_proxy -> log10 med={np.median(log_mass[keep]):.2f} "
            f"range {log_mass[keep].min():.2f}..{log_mass[keep].max():.2f}; "
            f"catalogue estimated_mass log10 med={np.median(np.log10(em)):.2f}"
        )
    return log_mass[keep], sigma[keep], z[keep], mult[keep].astype(int)


def run_real_sdss(
    parquet_path,
    sky_area_deg2_val=None,
    sky_frac=None,
    sdss_zmin=0.01,
    sdss_zmax=0.08,
    model_kind="marg",
):
    """Fit the MRP to a REAL per-object SDSS group catalogue with OUR marginalised
    model. z-range defaults to Driver's SDSS cut (0.01-0.08): beyond that the
    turnover mlim runs away (only massive groups survive at high z). Pass either
    sky_frac (fractional) or sky_area_deg2_val. Per-object sigma reuses the GAMA
    error model for now. 'truth' lines are Driver+22."""
    print(f"Reading SDSS group catalogue: {parquet_path}")
    log_mass, sigma, z, mult = load_sdss_groups(
        parquet_path, zmin=sdss_zmin, zmax=sdss_zmax
    )
    if sky_frac is None:
        sky_frac = sky_area_deg2_val * (np.pi / 180) ** 2 / (4 * np.pi)
    print(
        f"  N groups: {log_mass.size}   mass {log_mass.min():.2f}..{log_mass.max():.2f} "
        f"(med {np.median(log_mass):.2f})   sigma med {np.median(sigma):.2f}"
    )
    print(
        f"  z range {z.min():.3f}..{z.max():.3f}  (fit z {sdss_zmin}-{sdss_zmax}, "
        f"frac={sky_frac:.5f})"
    )

    print("Turnover mlim(z) ...")
    mlim_func, coefs, tkind, turn_pts = turnover_mlim(
        z, log_mass, zmin=sdss_zmin, zmax=sdss_zmax
    )
    print(
        f"  mlim(z) [{tkind}]: mlim({sdss_zmin})={mlim_func(sdss_zmin):.2f} "
        f"mlim({sdss_zmax})={mlim_func(sdss_zmax):.2f}"
    )

    z_mids, V_sh = shell_volumes(sky_frac, zmin=sdss_zmin, zmax=sdss_zmax)
    Vsurvey = float(V_sh.sum())
    mlim_sh = mlim_func(z_mids)

    mlim_per = mlim_func(z)
    above = log_mass > mlim_per
    x_fit, sig_fit = log_mass[above], sigma[above]
    print(
        f"  N above mlim: {x_fit.size} / {log_mass.size} "
        f"({100 * x_fit.size / log_mass.size:.1f}%)"
    )

    if model_kind == "marg":
        sig_sh = sigma_eff_per_shell(
            z, log_mass, sigma, mlim_sh, zmin=sdss_zmin, zmax=sdss_zmax
        )
        data = build_data("marg", x_fit, sig_fit, mlim_sh, V_sh, sig_sh=sig_sh)
    elif model_kind == "gama":
        data = build_data(
            "gama", x_fit, sig_fit, mlim_sh, V_sh, mlim_obj=mlim_per[above]
        )
    else:
        data = build_data(model_kind, x_fit, sig_fit, mlim_sh, V_sh)

    print(f"\nFitting [{model_kind}] on real SDSS (cmdstanpy) ...")
    map_par, flat = run_stan(model_kind, data)
    res = summarise(flat)  # offset vs Driver+22
    plot_recovery(
        flat,
        z,
        log_mass,
        x_fit,
        mlim_func,
        Vsurvey,
        turn_pts=turn_pts,
        fname=f"recovery_sdss_{model_kind}.pdf",
    )
    emit_publication(
        flat,
        {"SDSS": dict(x_fit=x_fit, Vsurvey=Vsurvey)},
        tag=f"sdss_{model_kind}",
        title="SDSS HMF",
    )
    return res


def _prep_survey(z, log_mass, sigma, sky_frac, zmin, zmax, label=""):
    """Turnover mlim, shells, sig_sh, and the above-mlim cut for one survey,
    over its own z-range and area. Returns everything the combined data needs."""
    mlim_func, coefs, tkind, turn_pts = turnover_mlim(z, log_mass, zmin=zmin, zmax=zmax)
    z_mids, V_sh = shell_volumes(sky_frac, zmin=zmin, zmax=zmax)
    mlim_sh = mlim_func(z_mids)
    sig_sh = sigma_eff_per_shell(z, log_mass, sigma, mlim_sh, zmin=zmin, zmax=zmax)
    mlim_per = mlim_func(z)
    above = log_mass > mlim_per
    x_fit, sig_fit = log_mass[above], sigma[above]
    print(
        f"  [{label}] mlim(z)[{tkind}] {mlim_func(zmin):.2f}->{mlim_func(zmax):.2f}"
        f"   N above mlim: {x_fit.size}/{log_mass.size}  Vsurvey={V_sh.sum():.3e}"
    )
    return dict(
        x_fit=x_fit,
        sig_fit=sig_fit,
        mlim_sh=mlim_sh,
        V_sh=V_sh,
        sig_sh=sig_sh,
        mlim_func=mlim_func,
        turn_pts=turn_pts,
        Vsurvey=float(V_sh.sum()),
    )


def build_data_combined(a, b):
    """Assemble the two-survey Stan data dict (survey a and b share the MRP)."""
    f = lambda v: np.asarray(v, float)
    return dict(
        N_a=int(a["x_fit"].size),
        x_obs_a=f(a["x_fit"]),
        sig_a=f(a["sig_fit"]),
        Nsh_a=int(a["V_sh"].size),
        V_sh_a=f(a["V_sh"]),
        mlim_sh_a=f(a["mlim_sh"]),
        sig_sh_a=f(a["sig_sh"]),
        N_b=int(b["x_fit"].size),
        x_obs_b=f(b["x_fit"]),
        sig_b=f(b["sig_fit"]),
        Nsh_b=int(b["V_sh"].size),
        V_sh_b=f(b["V_sh"]),
        mlim_sh_b=f(b["mlim_sh"]),
        sig_sh_b=f(b["sig_sh"]),
        xhi=float(XHI),
        Ng=int(NG),
        Nint=int(NINT),
    )


def plot_combined(flat, surveys, fname="recovery_combined.pdf"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    med = np.median(flat, axis=0)
    tv = np.array([TRUE[p] for p in PARAMS])
    mgrid = np.linspace(11, 16, 300)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))

    a = ax[0, 0]
    idx = np.random.default_rng(0).choice(
        flat.shape[0], size=min(300, flat.shape[0]), replace=False
    )
    for k in idx:
        a.plot(mgrid, np.log10(mrp_phi(mgrid, *flat[k])), color="red", alpha=0.02)
    a.plot(
        mgrid, np.log10(mrp_phi(mgrid, *med)), color="red", lw=2, label="MCMC median"
    )
    a.plot(mgrid, np.log10(mrp_phi(mgrid, *tv)), "b--", lw=2, label="Driver+22 GSR")
    for name, s in surveys.items():
        edges = np.arange(12, 16, 0.2)
        cen = 0.5 * (edges[:-1] + edges[1:])
        cnt, _ = np.histogram(s["x_fit"], bins=edges)
        with np.errstate(divide="ignore"):
            phi = np.log10(cnt / (s["Vsurvey"] * 0.2))
        a.scatter(cen, phi, s=25, label=name, zorder=5)
    a.set(
        xlim=(11, 16),
        ylim=(-8, -1),
        xlabel=r"$\log_{10} M$",
        ylabel=r"$\log_{10}\phi$",
        title="Combined HMF (GAMA + SDSS)",
    )
    a.legend(fontsize=8)

    labels = [r"$M_*$", r"$\log\phi_*$", r"$\alpha$", r"$\beta$"]
    prior_mu, prior_sd = [14.13, -3.96, -1.68, 0.63], [0.42, 0.69, 0.22, 0.02]
    for i, (r, c) in enumerate([(0, 1), (0, 2), (1, 0), (1, 1)]):
        aa = ax[r, c]
        aa.hist(flat[:, i], bins=40, color="steelblue", density=True)
        xs = np.linspace(*aa.get_xlim(), 200)
        aa.plot(
            xs,
            np.exp(-0.5 * ((xs - prior_mu[i]) / prior_sd[i]) ** 2)
            / (prior_sd[i] * np.sqrt(2 * np.pi)),
            "g:",
            label="prior",
        )
        aa.axvline(tv[i], color="blue", ls="--", label="Driver+22")
        aa.axvline(med[i], color="red", label="median")
        aa.set(title=labels[i])
        aa.legend(fontsize=7)
    ax[1, 2].axis("off")
    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")
    return fname


# ---------------------------------------------------------------------------
# Publication plotting: LCDM curve (exact Murray+21, from Driver's allhmf.r),
# a Driver-style multi-survey HMF figure, and a corner plot.
# ---------------------------------------------------------------------------
def lcdm_curve():
    """Murray+21 LCDM MRP, Omega_M-normalised and z=0.1-shifted, exactly as
    Driver+22 allhmf.r. Returns (x, log10 phi) for the curve plus the (M*, logphi*,
    alpha, beta) reference point for the corner plot."""
    parsec, Gn, msol = 3.0857e16, 6.67408e-11, 1.988e30
    rhocrit = 3 * (1000 * H0 / (1e6 * parsec)) ** 2 / (8 * np.pi * Gn)
    be, A, ms, al = 0.7097976, 1.727006e-19, 14.42947, -1.864908
    mrpx = np.arange(0, 17, 0.001) + np.log10(100 / H0)
    mrpy = (
        A
        * be
        * 10 ** ((al + 1) * (mrpx - ms))
        * np.exp(-(10 ** (be * (mrpx - ms))))
        * (H0 / 100) ** 3
    )
    factor = (
        np.sum(10**mrpx * mrpy)
        * 0.001
        * msol
        / (1e6 * parsec) ** 3
        / (OMEGA_M * rhocrit)
    )
    x = mrpx - 0.08
    y = np.log10(mrpy) - np.log10(factor) + 0.08
    ref = (ms - 0.075, np.log10(A / factor) + 0.075, al, be)  # Driver's corner refvals
    return x, y, ref


def _load_comparison(data_dir):
    """Load REFLEX II, 2PIGG, and ELMO/Tempel binned HMFs from data_dir, with the
    exact h-unit conversions from allhmf.r. Missing files are skipped."""
    import os

    out = {}
    try:
        import pandas as pd

        rf = pd.read_csv(os.path.join(data_dir, "reflex.csv"))
        rx = rf["x"].values + np.log10(70 / H0)
        ry = rf["Curve1"].values + 4.0 * np.log10(H0 / 70) - 14.0 + rx + 1
        f = np.full(rx.size, 1 / np.sqrt(20))
        f[0] = f[-1] = 1 / np.sqrt(3)
        f = np.clip(f, 0, 0.99)
        elo = -np.log10(1 - f)  # positive distance below
        ehi = np.log10(1 + f)  # positive distance above
        out["REFLEX II (Böhringer+17)"] = dict(
            x=rx, y=ry, elo=elo, ehi=ehi, c="forestgreen", m="D"
        )
    except Exception as e:
        print(f"  [reflex.csv skipped: {e}]")
    try:
        tp = np.loadtxt(os.path.join(data_dir, "tpigg.dat"))
        tx = tp[:, 0] + np.log10(100 / H0)
        ty = tp[:, 1] + np.log10((H0 / 100) ** 3)
        out["2PIGG (Eke+08)"] = dict(
            x=tx, y=ty, elo=np.abs(tp[:, 3]), ehi=np.abs(tp[:, 2]), c="0.5", m="o"
        )
    except Exception as e:
        print(f"  [tpigg.dat skipped: {e}]")
    try:
        import pandas as pd

        el = pd.read_csv(os.path.join(data_dir, "elmo.csv"), header=None)
        ex = el[0].values + 10.0 + np.log10(100 / H0)
        v2, v3, v4 = (
            el[1].values * (H0 / 100) ** 3,
            el[2].values * (H0 / 100) ** 3,
            el[3].values * (H0 / 100) ** 3,
        )
        out["_elmo"] = dict(
            x=ex,
            ycen=np.log10(v2),
            ylo=np.log10(np.clip(v2 - v3, 1e-30, None)),
            yhi=np.log10(v2 + v4),
        )
    except Exception as e:
        print(f"  [elmo.csv skipped: {e}]")
    return out


def plot_publication(
    flat,
    my_surveys=None,
    driver=(14.13, -3.96, -1.68, 0.63),
    data_dir="../data",
    fname="hmf_publication.pdf",
    title="HMF",
    show_corrected=True,
    corrected_band=False,
):
    """Driver-style HMF: our MCMC band (highlighted), the LCDM curve, the Driver+22
    MRP curve, our own survey binned points (grey, 'not fitted'), and external
    comparison data (REFLEX/2PIGG/ELMO) from data_dir."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    med = np.median(flat, axis=0)
    mgrid = np.linspace(12.5, 16, 400)
    fig, ax = plt.subplots(figsize=(8, 6))

    # our posterior band
    idx = np.random.default_rng(0).choice(
        flat.shape[0], size=min(400, flat.shape[0]), replace=False
    )
    for k in idx:
        ax.plot(
            mgrid,
            np.log10(mrp_phi(mgrid, *flat[k])),
            color=(100 / 255, 149 / 255, 237 / 255),
            alpha=0.02,
        )
    ax.plot(
        mgrid,
        np.log10(mrp_phi(mgrid, *med)),
        color=(100 / 255, 149 / 255, 237 / 255),
        lw=2,
        ls="--",
        label="Best-fit MRP (this work)",
    )

    # alpha bias-corrected curve: alpha -> alpha - bias(M*), other params as fitted.
    # NOTE this curve is NOT a fit to the data -- only alpha is corrected, while
    # M*/phi*/beta keep their (also biased) fitted values, so it will sit off the
    # points. It is shown to indicate the size of the alpha systematic.
    if show_corrected:
        if corrected_band:
            for k in idx:
                pc = flat[k].copy()
                pc[2] = pc[2] - float(alpha_bias(pc[0]))
                ax.plot(
                    mgrid, np.log10(mrp_phi(mgrid, *pc)), color="darkorange", alpha=0.02
                )
        medc = med.copy()
        medc[2] = float(np.median(corrected_alpha(flat)))
        ax.plot(
            mgrid,
            np.log10(mrp_phi(mgrid, *medc)),
            color="darkorange",
            lw=2,
            ls="-.",
            label=r"$\alpha$ bias-corrected (not a fit)",
        )

    # LCDM + Driver MRP curves
    lx, ly, _ = lcdm_curve()
    ax.plot(lx, ly, "k--", lw=2, label=r"$\Lambda$CDM (Murray+21, z=0.1)")
    ax.plot(
        mgrid,
        np.log10(mrp_phi(mgrid, *driver)),
        color="red",
        lw=1.5,
        ls=":",
        label="Driver+22 MRP (GSR)",
    )

    # our own survey binned points (grey, not fitted)
    if my_surveys:
        edges = np.arange(12.5, 16, 0.2)
        cen = 0.5 * (edges[:-1] + edges[1:])
        for name, s in my_surveys.items():
            cnt, _ = np.histogram(s["x_fit"], bins=edges)
            ok = cnt >= 3
            with np.errstate(divide="ignore"):
                phi = np.log10(cnt / (s["Vsurvey"] * 0.2))
                err = 0.4343 / np.sqrt(np.maximum(cnt, 1))
            ax.errorbar(
                cen[ok],
                phi[ok],
                yerr=err[ok],
                fmt="s",
                ms=4,
                color="0.35",
                alpha=0.7,
                capsize=2,
                label=f"{name} (this work, not fitted)",
            )

    # external comparison data
    comp = _load_comparison(data_dir)
    for name, d in comp.items():
        if name == "_elmo":
            ax.fill_between(d["x"], d["ylo"], d["yhi"], color="0.5", alpha=0.12, lw=0)
            ax.plot(
                d["x"], d["ycen"], color="cyan", lw=1, label="SDSS DR10 (Tempel+14)"
            )
        else:
            elo = np.nan_to_num(np.abs(d["elo"]), nan=0.0)
            ehi = np.nan_to_num(np.abs(d["ehi"]), nan=0.0)
            ax.errorbar(
                d["x"],
                d["y"],
                yerr=[elo, ehi],
                fmt=d["m"],
                ms=5,
                color=d["c"],
                capsize=2,
                lw=1,
                label=name,
            )

    ax.set(
        xlim=(12.75, 16),
        ylim=(-8, -2),
        xlabel=r"$\log_{10}(M_{\rm halo}/M_\odot)$",
        ylabel=r"$\log_{10}$ number density [Mpc$^{-3}$ dex$^{-1}$]",
        title=title,
    )
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")
    return fname


def _cred_levels(H):
    """68/95% contour levels for a 2D histogram."""
    Hs = np.sort(H.ravel())[::-1]
    cs = np.cumsum(Hs)
    cs = cs / cs[-1]
    lv = []
    for frac in (0.95, 0.68):
        i = np.searchsorted(cs, frac)
        lv.append(Hs[min(i, Hs.size - 1)])
    return sorted(set(lv))


def plot_corner(
    flat,
    driver=(14.13, -3.96, -1.68, 0.63),
    fname="corner.pdf",
    use_corrected_alpha=False,
):
    """Corner plot for (M*, logphi*, alpha, beta): 2D credible contours off the
    diagonal, marginals on it, with the median (red) and Driver+22 (black)
    reference points. (LCDM single-parameter markers are intentionally omitted:
    comparing a theoretical curve's MRP parameters to a data fit is misleading --
    LCDM consistency belongs on the HMF plot, not here.)"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from scipy.ndimage import gaussian_filter

        smooth = lambda H: gaussian_filter(H, 1.0)
    except Exception:
        smooth = lambda H: H

    labels = [r"$\log_{10} M_*$", r"$\log_{10}\phi_*$", r"$\alpha$", r"$\beta$"]
    if use_corrected_alpha:
        flat = flat.copy()
        flat[:, 2] = corrected_alpha(flat)
        labels[2] = r"$\alpha_{\rm corr}$"
    med = np.median(flat, axis=0)
    P = flat.shape[1]
    fig, axes = plt.subplots(P, P, figsize=(9, 9))
    for i in range(P):
        for j in range(P):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
                continue
            if i == j:
                ax.hist(flat[:, i], bins=40, color="steelblue", density=True)
                ax.axvline(med[i], color="red")
                if driver[i] is not None:
                    ax.axvline(driver[i], color="k", ls="--")
                ax.set_yticks([])
            else:
                H, xe, ye = np.histogram2d(flat[:, j], flat[:, i], bins=40)
                Hs = smooth(H.T)
                xc, yc = 0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:])
                lv = _cred_levels(Hs)
                ax.contourf(
                    xc,
                    yc,
                    Hs,
                    levels=[lv[0], lv[-1], Hs.max()],
                    colors=["#c6dbef", "#6baed6"],
                )
                ax.contour(xc, yc, Hs, levels=lv, colors="#2171b5", linewidths=0.6)
                ax.plot(med[j], med[i], "r+", ms=8)
                if driver[j] is not None and driver[i] is not None:
                    ax.plot(driver[j], driver[i], "k*", ms=9)
            if i == P - 1:
                ax.set_xlabel(labels[j])
            else:
                ax.set_xticklabels([])
            if j == 0 and i > 0:
                ax.set_ylabel(labels[i])
            elif j > 0:
                ax.set_yticklabels([])
    axes[0, 0].set_yticklabels([])
    from matplotlib.lines import Line2D

    proxies = [
        Line2D([0], [0], color="red", marker="+", ls="", label="median (this work)"),
        Line2D([0], [0], color="k", marker="*", ls="", label="Driver+22 GSR"),
    ]
    fig.legend(handles=proxies, loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")
    return fname


def emit_publication(flat, my_surveys, tag, title="HMF"):
    """Always emit the corner + publication HMF plots for a fit."""
    try:
        plot_corner(flat, fname=f"corner_{tag}.pdf")
    except Exception as e:
        print(f"  [corner_{tag} failed: {e}]")
    try:
        plot_publication(
            flat, my_surveys=my_surveys, fname=f"hmf_{tag}.pdf", title=title
        )
    except Exception as e:
        print(f"  [hmf_{tag} failed: {e}]")


def run_combined(
    gama_fits,
    sdss_parquet,
    gama_area=179.92,
    sdss_frac=0.2126803,
    sdss_zmin=0.01,
    sdss_zmax=0.20,
):
    """Joint GAMA+SDSS fit: one shared MRP, each survey with its own volume,
    z-range, mlim(z), and per-object errors. The valid multi-survey version of
    the marginalised method (per-object, not binned)."""
    print(f"Reading GAMA: {gama_fits}")
    g_lm, g_sig, g_z, _ = load_real_gama(gama_fits)
    print(f"Reading SDSS: {sdss_parquet}")
    s_lm, s_sig, s_z, _ = load_sdss_groups(sdss_parquet, zmin=sdss_zmin, zmax=sdss_zmax)

    gama_frac = gama_area * (np.pi / 180) ** 2 / (4 * np.pi)
    print(f"  GAMA: {g_lm.size} groups (frac={gama_frac:.5f}, z {ZMIN}-{ZLIMIT})")
    print(
        f"  SDSS: {s_lm.size} groups (frac={sdss_frac:.5f}, z {sdss_zmin}-{sdss_zmax})"
    )

    A = _prep_survey(g_z, g_lm, g_sig, gama_frac, ZMIN, ZLIMIT, "GAMA")
    B = _prep_survey(s_z, s_lm, s_sig, sdss_frac, sdss_zmin, sdss_zmax, "SDSS")

    data = build_data_combined(A, B)
    print("\nFitting [combined] GAMA + SDSS (cmdstanpy) ...")
    map_par, flat = run_stan("combined", data)
    res = summarise(flat)  # offset vs Driver+22 GSR
    plot_combined(flat, {"GAMA": A, "SDSS": B})
    return res


def run_combined_comp(
    gama_fits,
    sdss_parquet,
    gama_area=179.92,
    sdss_frac=0.2126803,
    sdss_zmin=0.01,
    sdss_zmax=0.08,
    full_sample=False,
    reflex=False,
    reflex_dX=0.0,
    fix_offset=False,
    reflex_mmin=None,
    data_dir="../data",
):
    """Joint GAMA+SDSS with completeness. GAMA's ramp is fixed (measured from the
    GAMA-selected mock); SDSS's ramp (D50, w) is FITTED with GAMA-informed priors,
    because the WAVES lightcone is too small in area to build an SDSS-like mock.
    One shared MRP; per-survey volumes, z-ranges, mlim(z), and completeness."""
    from scipy.stats import norm

    print(f"Reading GAMA: {gama_fits}")
    g_lm, g_sig, g_z, _ = load_real_gama(gama_fits)
    print(f"Reading SDSS: {sdss_parquet}")
    s_lm, s_sig, s_z, _ = load_sdss_groups(sdss_parquet, zmin=sdss_zmin, zmax=sdss_zmax)

    gama_frac = gama_area * (np.pi / 180) ** 2 / (4 * np.pi)

    # --- GAMA block: measured ramp, C>CMIN floor, no mlim cut ---
    g_mlim_func, _, gk, _ = turnover_mlim(g_z, g_lm, zmin=ZMIN, zmax=ZLIMIT)
    gz_mids, gV_sh = shell_volumes(gama_frac, zmin=ZMIN, zmax=ZLIMIT)
    g_mlim_sh = g_mlim_func(gz_mids)
    g_mlim_obj = g_mlim_func(g_z)
    g_d50 = np.interp(g_z, COMP_Z_PTS, COMP_D50_PTS)
    g_w = np.interp(g_z, COMP_Z_PTS, COMP_W_PTS)
    g_keep = norm.cdf((g_lm - g_mlim_obj - g_d50) / g_w) > CMIN
    print(
        f"  [GAMA] mlim[{gk}] {g_mlim_func(ZMIN):.2f}->{g_mlim_func(ZLIMIT):.2f}  "
        f"kept {int(g_keep.sum())}/{g_lm.size} (C>{CMIN})  V={gV_sh.sum():.3e}"
    )

    # --- SDSS block: ramp FIXED (adopted = GAMA-measured; the WAVES lightcone
    #     is too small to measure an SDSS-specific ramp). C precomputed, C>CMIN. ---
    s_mlim_func, _, sk, _ = turnover_mlim(s_z, s_lm, zmin=sdss_zmin, zmax=sdss_zmax)
    sz_mids, sV_sh = shell_volumes(sdss_frac, zmin=sdss_zmin, zmax=sdss_zmax)
    s_mlim_sh = s_mlim_func(sz_mids)
    s_mlim_obj = s_mlim_func(s_z)
    s_d50 = np.interp(s_z, COMP_Z_PTS, COMP_D50_PTS)  # adopted GAMA ramp
    s_w = np.interp(s_z, COMP_Z_PTS, COMP_W_PTS)
    s_keep = norm.cdf((s_lm - s_mlim_obj - s_d50) / s_w) > CMIN
    print(
        f"  [SDSS] mlim[{sk}] {s_mlim_func(sdss_zmin):.2f}->{s_mlim_func(sdss_zmax):.2f}  "
        f"kept {int(s_keep.sum())}/{s_lm.size} (C>{CMIN}, ramp adopted from GAMA)  "
        f"V={sV_sh.sum():.3e}"
    )

    f = lambda v: np.asarray(v, float)
    data = dict(
        N_a=int(g_keep.sum()),
        x_obs_a=f(g_lm[g_keep]),
        sig_a=f(g_sig[g_keep]),
        mlim_obj_a=f(g_mlim_obj[g_keep]),
        d50_obj_a=f(g_d50[g_keep]),
        w_obj_a=f(g_w[g_keep]),
        Nsh_a=int(gV_sh.size),
        V_sh_a=f(gV_sh),
        mlim_sh_a=f(g_mlim_sh),
        d50_sh_a=f(np.interp(gz_mids, COMP_Z_PTS, COMP_D50_PTS)),
        w_sh_a=f(np.interp(gz_mids, COMP_Z_PTS, COMP_W_PTS)),
        N_b=int(s_keep.sum()),
        x_obs_b=f(s_lm[s_keep]),
        sig_b=f(s_sig[s_keep]),
        mlim_obj_b=f(s_mlim_obj[s_keep]),
        d50_obj_b=f(s_d50[s_keep]),
        w_obj_b=f(s_w[s_keep]),
        Nsh_b=int(sV_sh.size),
        V_sh_b=f(sV_sh),
        mlim_sh_b=f(s_mlim_sh),
        d50_sh_b=f(np.interp(sz_mids, COMP_Z_PTS, COMP_D50_PTS)),
        w_sh_b=f(np.interp(sz_mids, COMP_Z_PTS, COMP_W_PTS)),
        xhi=float(XHI),
        Ng=int(NG),
        Nint=31,
        cmin=float(CMIN),
    )

    model_kind = "combined_comp"
    if reflex:
        comp = _load_comparison(data_dir)
        rf = comp.get("REFLEX II (Böhringer+17)")
        if rf is None:
            print(
                "  [reflex requested but reflex.csv not loaded -> running without it]"
            )
        else:
            # symmetric log-error from the fractional errors (mean of elo/ehi)
            sig_r = 0.5 * (np.abs(rf["elo"]) + np.abs(rf["ehi"]))
            sig_r = np.clip(sig_r, 0.03, None)
            rx, ry = np.asarray(rf["x"], float), np.asarray(rf["y"], float)
            if reflex_mmin is not None:
                keepr = rx > float(reflex_mmin)
                print(
                    f"  [REFLEX] mass cut logM > {reflex_mmin}: "
                    f"{int(keepr.sum())}/{rx.size} points kept"
                )
                rx, ry, sig_r = rx[keepr], ry[keepr], sig_r[keepr]
            base = dict(N_r=int(rx.size), m_r=f(rx), y_r=f(ry), sig_r=f(sig_r))
            if fix_offset:
                # offset as DATA -> parameter space stays 4D, samples normally
                base.update(dXa_fix=float(reflex_dX), dXb_fix=0.0)
                model_kind = "combined_comp_reflex_fixed"
            else:
                base.update(
                    dXa_mu=float(reflex_dX), dXa_sd=0.25, dXb_mu=0.0, dXb_sd=0.3
                )
                model_kind = "combined_comp_reflex"
            data.update(base)
            print(
                f"  [REFLEX] {rx.size} binned points, "
                f"logM {rx.min():.2f}-{rx.max():.2f}, "
                f"dX(M) fixed at {reflex_dX:+.2f}"
                if fix_offset
                else f"dX(M)=a+b(M-14): a~N({reflex_dX:+.2f},0.25), b~N(0,0.3) (fitted)"
            )

    model_kind, data = apply_driver_prior(model_kind, data)
    print(
        f"\nFitting [{model_kind}] GAMA + SDSS"
        f"{' + REFLEX' if 'reflex' in model_kind else ''} (cmdstanpy) ..."
    )
    model = get_model(model_kind)

    # --- Stan MAP first (same model, L-BFGS): a fast point result + a sane init ---
    print("  MAP (Stan optimize, L-BFGS) ...")
    init0 = dict(ms=14.13, lp=-3.96, al=-1.68, be=0.63)
    if model_kind.endswith("reflex"):
        init0["dXa"] = float(reflex_dX)
        init0["dXb"] = 0.0
    opt = model.optimize(
        data=data, inits=init0, algorithm="lbfgs", iter=20000, show_console=False
    )
    mp = {p: float(opt.optimized_params_dict[p]) for p in PARAMS}
    print(
        f"  MAP:  ms={mp['ms']:.3f}  lp={mp['lp']:.3f}  al={mp['al']:.3f}  be={mp['be']:.3f}"
    )
    if model_kind.endswith("reflex"):
        da = float(opt.optimized_params_dict["dXa"])
        db = float(opt.optimized_params_dict["dXb"])
        print(
            f"        REFLEX offset dX(M) = {da:+.3f} {db:+.3f}*(logM-14)  "
            f"[at 10^14: {da:+.2f}, at 10^15: {da + db:+.2f} dex]"
        )
    print(f"        vs Driver GSR 14.13 / -3.96 / -1.68 / 0.63")
    print(
        f"  SDSS ramp adopted from GAMA: D50 ~ {np.mean(COMP_D50_PTS):+.3f}, "
        f"w ~ {np.mean(COMP_W_PTS):.3f} (fixed, not fitted)"
    )

    if not full_sample:
        print("  (MAP only; pass full_sample=True to also run MCMC)")
        return mp

    fit = model.sample(
        data=data,
        chains=4,
        iter_warmup=1500,
        iter_sampling=1500,
        adapt_delta=0.95,
        max_treedepth=10,
        seed=42,
        show_progress=True,
        inits=init0,
    )
    flat = np.column_stack([fit.stan_variable(p) for p in PARAMS])
    otag = "combined_comp_reflex" if "reflex" in model_kind else "combined_comp"
    save_arr, save_hdr = flat, "ms,lp,al,be"
    if model_kind.endswith("reflex"):
        da, db = fit.stan_variable("dXa"), fit.stan_variable("dXb")
        print(
            f"  fitted REFLEX offset: a(10^14) = {np.median(da):+.3f}"
            f"±{0.5 * (np.percentile(da, 84) - np.percentile(da, 16)):.3f}, "
            f"slope b = {np.median(db):+.3f}"
            f"±{0.5 * (np.percentile(db, 84) - np.percentile(db, 16)):.3f} dex/dex"
        )
        print(
            f"  [prior was a~N({reflex_dX:+.2f},0.25), b~N(0,0.30) -- compare widths:"
            f" posterior/prior = {0.5 * (np.percentile(db, 84) - np.percentile(db, 16)) / 0.30:.2f}"
            f" (near 1 => prior-driven, << 1 => measured)]"
        )
        save_arr = np.column_stack([flat, da, db])
        save_hdr = "ms,lp,al,be,dXa,dXb"
    np.savetxt(
        f"{otag}_draws.csv", save_arr, delimiter=",", header=save_hdr, comments=""
    )
    print(f"  saved draws -> {otag}_draws.csv")
    try:
        print(
            f"  Rhat/ESS: {fit.diagnose().splitlines()[0] if hasattr(fit, 'diagnose') else 'n/a'}"
        )
    except Exception:
        pass
    res = summarise(flat)
    A = dict(x_fit=g_lm[g_keep], Vsurvey=float(gV_sh.sum()))
    B = dict(x_fit=s_lm[s_keep], Vsurvey=float(sV_sh.sum()))
    ttl = "GAMA + SDSS + REFLEX HMF" if "reflex" in model_kind else "GAMA + SDSS HMF"
    emit_publication(flat, {"GAMA": A, "SDSS": B}, tag=otag, title=ttl)
    try:
        plot_combined(flat, {"GAMA": A, "SDSS": B}, fname=f"recovery_{otag}.pdf")
    except Exception as e:
        print(f"  [plot failed: {e}] draws are saved; re-plot from the CSV.")
    return res


def calibrate_alpha_bias(n_real=20, ms_pins=(14.13, 14.35, 14.60), seed0=500):
    """Measure the alpha bias vs M* over many mock realisations (numpy, no Stan).
    For each realisation and each pinned M*, fit (logphi*, alpha) with the
    production completeness likelihood and record alpha_fit - alpha_true.
    Prints mean bias +/- scatter at each pin -> the error bar that
    ALPHA_BIAS_DA currently lacks. Writes alpha_bias_calibration.csv."""
    from scipy import optimize
    from scipy.stats import norm

    groups, galaxies = load_catalogues(DATA_DIR)
    sky_frac = (
        sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, _ = abundance_match(groups, sky_frac)
    gv = gama_select(gv, galaxies)
    z_mids, V_sh = shell_volumes(sky_frac, zmin=ZMIN, zmax=ZLIMIT)
    ln10 = np.log(10.0)
    tv_al, be = TRUE["al"], TRUE["be"]

    def phi(m, ms, lp, al):
        u = m - ms
        return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))

    out = np.full((n_real, len(ms_pins)), np.nan)
    for r in range(n_real):
        rng = np.random.default_rng(seed0 + r)
        z, m_obs, sig, _ = add_mass_errors(gv, rng)
        try:
            mlim_func, _, _, _ = turnover_mlim(z, m_obs, zmin=ZMIN, zmax=ZLIMIT)
        except Exception:
            continue
        mlim_sh, mlim_o = mlim_func(z_mids), mlim_func(z)
        d50_o = np.interp(z, COMP_Z_PTS, COMP_D50_PTS)
        w_o = np.interp(z, COMP_Z_PTS, COMP_W_PTS)
        d50_s = np.interp(z_mids, COMP_Z_PTS, COMP_D50_PTS)
        w_s = np.interp(z_mids, COMP_Z_PTS, COMP_W_PTS)
        keep = norm.cdf((m_obs - mlim_o - d50_o) / w_o) > CMIN
        x, sg = m_obs[keep], sig[keep]
        ml_k, d50_k, w_k = mlim_o[keep], d50_o[keep], w_o[keep]
        mg = np.linspace(mlim_sh.min() - 3.0, XHI, 1000)
        gq = np.linspace(-5, 5, 41)

        def nll(t, msf):
            lp_, al_ = t
            pg = phi(mg, msf, lp_, al_)
            lam = 0.0
            for j in range(len(V_sh)):
                Cj = norm.cdf((mg - mlim_sh[j] - d50_s[j]) / w_s[j])
                lam += V_sh[j] * np.trapezoid(pg * np.where(Cj > CMIN, Cj, 0.0), mg)
            mt = x[:, None] + sg[:, None] * gq[None, :]
            Cm = norm.cdf((mt - ml_k[:, None] - d50_k[:, None]) / w_k[:, None])
            integ = (
                phi(mt, msf, lp_, al_)
                * Cm
                * np.exp(-0.5 * gq[None, :] ** 2)
                / (sg[:, None] * np.sqrt(2 * np.pi))
            )
            ll = np.sum(np.log(np.maximum(np.trapezoid(integ, mt, axis=1), 1e-300)))
            return -(-lam + ll)

        for k, msf in enumerate(ms_pins):
            try:
                res = optimize.minimize(
                    nll,
                    [-4.0, -1.6],
                    args=(msf,),
                    method="Nelder-Mead",
                    options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000),
                )
                out[r, k] = res.x[1] - tv_al
            except Exception:
                pass
        done = np.isfinite(out[: r + 1]).all(axis=1).sum()
        print(f"  [{r + 1}/{n_real}] complete realisations: {done}")

    print(f"\n  alpha bias vs pinned M*  ({n_real} realisations)")
    print(f"  {'M* pin':>8} {'mean bias':>10} {'scatter':>9} {'N':>4}")
    for k, msf in enumerate(ms_pins):
        col = out[:, k][np.isfinite(out[:, k])]
        if col.size:
            print(f"  {msf:8.2f} {col.mean():+10.3f} {col.std():9.3f} {col.size:4d}")
    np.savetxt(
        "alpha_bias_calibration.csv",
        out,
        delimiter=",",
        header=",".join(f"ms_{m:.2f}" for m in ms_pins),
        comments="",
    )
    print("  saved alpha_bias_calibration.csv")
    print("  -> put the mean biases in ALPHA_BIAS_MS/ALPHA_BIAS_DA and add the")
    print("     scatter in quadrature to alpha's error bar.")
    return out


def run_selftest(model_kind="marg"):
    """Synthetic recovery, no data files. Confirms the chosen Stan model
    recovers a known MRP.

    simple: masses drawn from the MRP with a SHARP limit, no errors -> pure
            likelihood/normalisation check (expected to PASS).
    marg  : true masses from the MRP, Gaussian errors added, then selected on
            the NOISY mass above the limit -> exercises the boundary + error
            de-biasing. The simple model would be biased on this data; the
            marginalised model should recover it.
    """
    print("=" * 60)
    print(f"  SELF-TEST [{model_kind}]: synthetic recovery (no data files)")
    print("=" * 60)
    rng = np.random.default_rng(7)

    sky_frac = 0.0005  # ~1k groups: enough to catch gross bias, ~10 min numerical
    z_mids, V_sh = shell_volumes(sky_frac)
    mlim_func = lambda z: 12.6 + 1.5 * z
    mlim_sh = mlim_func(z_mids)

    if model_kind == "simple":
        # sharp limit, no errors
        mg = np.arange(10.0, XHI, 0.001)
        phi_g = mrp_phi(mg, **TRUE)
        x_list = []
        for j in range(len(z_mids)):
            mask = mg >= mlim_sh[j]
            n_j = rng.poisson(V_sh[j] * _trapz(phi_g[mask], mg[mask]))
            if n_j == 0:
                continue
            cdf = np.cumsum(phi_g[mask])
            cdf = cdf / cdf[-1]
            x_list.append(np.interp(rng.random(n_j), cdf, mg[mask]))
        x_fit = np.concatenate(x_list)
        sig_fit = sig_sh = None
    else:
        # realistic: draw true masses well below the limit, add errors, select
        # on the noisy mass. A constant sigma keeps the test clean and stresses
        # the boundary.
        sig0 = 0.30
        m_floor = mlim_sh.min() - 5 * sig0  # generate below the limit too
        mg = np.arange(m_floor, XHI, 0.001)
        phi_g = mrp_phi(mg, **TRUE)
        cdf = np.cumsum(phi_g)
        cdf = cdf / cdf[-1]
        xt_list, z_list = [], []
        for j in range(len(z_mids)):
            n_j = rng.poisson(V_sh[j] * _trapz(phi_g, mg))
            if n_j == 0:
                continue
            xt_list.append(np.interp(rng.random(n_j), cdf, mg))
            z_list.append(np.full(n_j, z_mids[j]))
        m_true = np.concatenate(xt_list)
        z_all = np.concatenate(z_list)
        m_obs = m_true + rng.normal(0, sig0, size=m_true.size)
        above = m_obs > mlim_func(z_all)
        x_fit = m_obs[above]
        sig_fit = np.full(x_fit.size, sig0)
        sig_sh = np.full(len(z_mids), sig0)

    print(f"  generated {x_fit.size} groups above the limit")
    data = build_data(model_kind, x_fit, sig_fit, mlim_sh, V_sh, sig_sh=sig_sh)
    map_par, flat = run_stan(
        model_kind, data, warmup=500, sampling=500, show_progress=False
    )
    res = summarise(flat)
    tv = np.array([TRUE[p] for p in PARAMS])
    bias = np.abs((res["median"] - tv) / res["sd"])
    worst = PARAMS[int(np.argmax(bias))]
    if bias.max() < 1.0:
        print(f"\n  max |bias| = {bias.max():.2f} sd ({worst})  -> PASS")
    else:
        print(f"\n  max |bias| = {bias.max():.2f} sd ({worst})  -> FAIL")
        print("  Bias exceeds 1 sd on a clean mock: the model is NOT trustworthy.")
        print("  Do not proceed to coverage/real data until this is understood.")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--model",
        choices=["simple", "marg", "marg_comp", "marg_comp_serr"],
        default="marg",
        help="which likelihood: 'simple' (baseline) or 'marg' "
        "(marginalised + boundary). Default marg.",
    )
    ap.add_argument(
        "--driver-prior",
        action="store_true",
        help="use Driver+22 GSR chains as a multivariate-normal prior",
    )
    ap.add_argument(
        "--driver-prior-inflate",
        type=float,
        default=1.0,
        help="widen the Driver prior covariance by this factor",
    )
    ap.add_argument(
        "--calibrate-alpha",
        action="store_true",
        help="measure the alpha-vs-M* bias over many mock realisations",
    )
    ap.add_argument(
        "--selftest",
        action="store_true",
        help="run synthetic recovery without needing the parquet files",
    )
    ap.add_argument(
        "--coverage",
        action="store_true",
        help="run the N-realisation coverage loop on the GAMA-mock",
    )
    ap.add_argument(
        "--nreal",
        type=int,
        default=20,
        help="number of realisations for --coverage (default 20)",
    )
    ap.add_argument(
        "--realgama",
        action="store_true",
        help="fit the REAL GAMA catalogue (our marg model by default)",
    )
    ap.add_argument(
        "--gama-fits",
        default="../data/G3CFoFGroupv10.fits",
        help="path to the GAMA G3C group FITS file (for --realgama)",
    )
    ap.add_argument(
        "--gama-area",
        type=float,
        default=179.92,
        help="GAMA sky area in deg^2 (for --realgama; default 179.92)",
    )
    ap.add_argument(
        "--gama-model",
        choices=["marg", "marg_comp", "marg_comp_serr", "gama", "simple"],
        default="marg",
        help="model for --realgama: 'marg' (sharp cut), 'marg_comp' "
        "(completeness forward-model), or 'gama' (R port, check only)",
    )
    ap.add_argument(
        "--realsdss",
        action="store_true",
        help="fit a REAL per-object SDSS group catalogue (parquet) with our marg model",
    )
    ap.add_argument(
        "--sdss-parquet",
        default="sdss_groups.parquet",
        help="path to the SDSS group parquet (for --realsdss)",
    )
    ap.add_argument(
        "--sdss-area",
        type=float,
        default=None,
        help="SDSS sky area in deg^2 (for --realsdss; or use --sdss-frac)",
    )
    ap.add_argument(
        "--combined",
        action="store_true",
        help="joint GAMA+SDSS fit with the shared-MRP two-survey model",
    )
    ap.add_argument(
        "--sdss-frac",
        type=float,
        default=0.2126803,
        help="SDSS fractional sky area (for --combined / --realsdss)",
    )
    ap.add_argument(
        "--sdss-zmin",
        type=float,
        default=0.01,
        help="SDSS lower z limit (default 0.01)",
    )
    ap.add_argument(
        "--sdss-zmax",
        type=float,
        default=0.08,
        help="SDSS upper z limit (default 0.08, Driver's SDSS cut)",
    )
    ap.add_argument(
        "--combined-comp",
        action="store_true",
        help="joint GAMA+SDSS with completeness (GAMA ramp measured, SDSS ramp fitted)",
    )
    ap.add_argument(
        "--full-sample",
        action="store_true",
        help="for --combined-comp: also run full MCMC after the MAP",
    )
    ap.add_argument(
        "--reflex",
        action="store_true",
        help="for --combined-comp: add REFLEX II binned points to anchor the cutoff",
    )
    ap.add_argument(
        "--reflex-mmin",
        type=float,
        default=None,
        help="only use REFLEX points above this logM (anchor the cutoff only)",
    )
    ap.add_argument(
        "--fix-offset",
        action="store_true",
        help="pin the REFLEX mass offset at --reflex-dx instead of fitting it",
    )
    ap.add_argument(
        "--reflex-dx",
        type=float,
        default=0.0,
        help="fixed X-ray->dynamical mass offset applied to REFLEX (dex, default 0)",
    )
    args = ap.parse_args()
    USE_DRIVER_PRIOR = args.driver_prior
    DRIVER_PRIOR_INFLATE = args.driver_prior_inflate
    if args.calibrate_alpha:
        calibrate_alpha_bias(n_real=args.nreal)
    elif args.selftest:
        run_selftest(model_kind=args.model)
    elif args.coverage:
        run_coverage(model_kind=args.model, n_real=args.nreal)
    elif args.realgama:
        run_real_gama(
            args.gama_fits,
            sky_area_deg2_val=args.gama_area,
            model_kind=args.gama_model,
            reflex=args.reflex,
            reflex_dX=args.reflex_dx,
            fix_offset=args.fix_offset,
            reflex_mmin=args.reflex_mmin,
        )
    elif args.realsdss:
        run_real_sdss(
            args.sdss_parquet,
            sky_area_deg2_val=args.sdss_area,
            sky_frac=(None if args.sdss_area else args.sdss_frac),
            sdss_zmin=args.sdss_zmin,
            sdss_zmax=args.sdss_zmax,
            model_kind=args.model,
        )
    elif args.combined_comp:
        run_combined_comp(
            args.gama_fits,
            args.sdss_parquet,
            gama_area=args.gama_area,
            sdss_frac=args.sdss_frac,
            sdss_zmin=args.sdss_zmin,
            sdss_zmax=args.sdss_zmax,
            full_sample=args.full_sample,
            reflex=args.reflex,
            reflex_dX=args.reflex_dx,
            fix_offset=args.fix_offset,
            reflex_mmin=args.reflex_mmin,
        )
    elif args.combined:
        run_combined(
            args.gama_fits,
            args.sdss_parquet,
            gama_area=args.gama_area,
            sdss_frac=args.sdss_frac,
            sdss_zmin=args.sdss_zmin,
            sdss_zmax=args.sdss_zmax,
        )
    else:
        run_real_pipeline(model_kind=args.model)
