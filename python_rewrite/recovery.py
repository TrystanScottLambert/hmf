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
# Driver+22 GSR, converted to the h=1 units this pipeline works in:
#   log M*   14.13  - log10(1/h) = 13.958
#   log phi* -3.96  - 3 log10(h) = -3.445
#   alpha, beta unchanged (dimensionless)
TRUE = dict(ms=13.958, lp=-3.445, al=-1.68, be=0.63)
TRUE_DRIVER = dict(ms=14.13, lp=-3.96, al=-1.68, be=0.63)  # as published
PARAMS = ["ms", "lp", "al", "be"]

# Cosmology. h=1, Om=0.25 is the Robotham+11 convention that the GAMA group
# catalogue and Nessie are built in, so the fit is done in those units:
# Rad50 is Mpc/h, masses come out Msun/h, volumes (Mpc/h)^3, phi* h^3 Mpc^-3.
# Note the (100/H0) factor in load_real_gama then becomes exactly 1.
# Driver+22 works at h=0.6737, Om=0.3147 -- use to_driver_cosmology() to
# convert results for comparison rather than changing the fit.
H0, OMEGA_M = 100.0, 0.25
H_DRIVER = 0.6737

# Survey / selection
ZMIN, ZLIMIT = 0.01, 0.25
MULTI = 5  # min members for a detection
# Selection band + limit. GAMA: r_SDSS < 19.8 ; WAVES: Z_VISTA < 21.1 (deeper NIR).
SEL_COL = "total_ap_dust_r_SDSS"
MAG_LIMIT = 19.65
# Dynamical-mass calibration prefactor M = A * sigma^2 R / G. Driver's fiducial
# is 13.9; A=10 is his variant (Fig. A2). Applied to BOTH surveys for a common
# mass scale.
A_SCALE = 10.0

# Completeness ramp C(Delta) = 0.5(1+erf((Delta-D50)/(sqrt2 w))), Delta=m-mlim(z),
# measured from the mock (measure_completeness.py). z-dependent: interpolate
# (D50, w) from these per-z-bin values. CMIN = floor (drop groups below 20%
# completeness; corrections below that are unreliable, cf. Driver+22).
# Measured for r < 19.65 (measure_completeness.py, 1335 detected).
# GLOBAL fit used, not the per-z one: with only ~1300 groups split three ways
# the per-z D50 is noise-dominated -- it reversed sign (-0.148/-0.193/-0.232 ->
# -0.236/-0.163/-0.063) under a 0.15 mag change, while mlim + D50, the mass at
# which completeness is actually 50%, moved by <0.1 dex. The global fit uses
# every group and is far more stable.
COMP_Z_PTS = [0.045, 0.115, 0.20]
COMP_D50_PTS = [-0.100, -0.100, -0.100]
COMP_W_PTS = [0.256, 0.256, 0.256]
# per-z alternative (r<19.65), for the systematic check:
#   COMP_D50_PTS = [-0.236, -0.163, -0.063]
#   COMP_W_PTS   = [0.314, 0.270, 0.223]
CMIN = 0.0

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
# !! measured at h=0.674 -- INVALID at h=1. Re-run --calibrate-alpha.
ALPHA_BIAS_MS = [14.13, 14.35, 14.60]
ALPHA_BIAS_DA = [0.189, 0.013, -0.136]
ALPHA_BIAS_SCATTER = 0.073  # realisation-to-realisation rms (20 mocks)

# The correction's zero-point is the alpha injected into Shark, so a corrected
# value carries that input plus only the real-minus-mock deviation. Off by
# default; enable with --alpha-correction if you want it reported/plotted.
SHOW_ALPHA_CORRECTION = False

# Our own binned points on the publication HMF are raw counts / total volume:
# no 1/Vmax, no completeness correction. They roll over below the limit and sit
# well under every other dataset, which reads as "our data disagrees with our
# fit" when it is really just an uncorrected estimator. Off by default.
SHOW_OWN_POINTS = False

# The dark-green "bias-corrected (halo MF)" curve is the fit with the
# closed-loop NESSIE_BIAS removed, i.e. an estimate of the underlying HALO mass
# function rather than the density of DETECTED groups.  It legitimately does not
# track the plotted points -- those are detected objects, and the two differ by
# the completeness -- but on the figure it just reads as a third unexplained
# line sitting off the data.  The corrected numbers are printed in the
# "closed-loop bias-corrected" table regardless, which is where they belong.
# Off by default; set True (or --show-halo-mf) to draw it.
SHOW_HALO_MF_CURVE = False


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
def to_driver_cosmology(ms, lp, al=None, be=None, h=H_DRIVER):
    """Convert MRP parameters from h=1 units (this fit) to Driver+22's h=0.6737.

        M[Msun]      = M[Msun/h] / h        ->  log10 M*  +0.1715
        phi[Mpc^-3]  = phi[h^3 Mpc^-3] * h^3 ->  log10 phi* -0.5146

    alpha and beta are dimensionless and unchanged. Returns (ms, lp, al, be)."""
    return (ms + np.log10(1.0 / h), lp + 3.0 * np.log10(h), al, be)


def from_driver_cosmology(ms, lp, al=None, be=None, h=H_DRIVER):
    """Inverse of to_driver_cosmology: Driver's published values -> h=1 units.
    Driver's GSR (14.13, -3.96) becomes (13.958, -3.445) at h=1."""
    return (ms - np.log10(1.0 / h), lp - 3.0 * np.log10(h), al, be)


# Driver+22's GAMA-ONLY fit (his "GAMA5" row: N>=5, A=13.9), as published in his
# h=0.6737 units. This is the like-for-like comparison for a GAMA-only analysis;
# the GSR values are a combined GAMA+SDSS+REFLEX fit.
#
# TWO conversions are needed, not one:
#   h    : log M* -log10(1/h),  log phi* -3log10(h)   (as for the GSR values)
#   A    : his A=13.9 vs our A_SCALE -- M ~ A sigma^2 R/G, so log M* shifts by
#          log10(A_SCALE/13.9). phi*, alpha and beta are unaffected by A.
DRIVER_GAMA5 = dict(ms=13.51, lp=-3.19, al=-1.27, be=0.47)  # as published
DRIVER_GAMA5_A = 13.9


def driver_gama5(match_A=True):
    """Driver's GAMA-only fit in this pipeline's units (h=1, and optionally
    rescaled from his A=13.9 to A_SCALE)."""
    ms, lp, al, be = from_driver_cosmology(
        DRIVER_GAMA5["ms"], DRIVER_GAMA5["lp"], DRIVER_GAMA5["al"], DRIVER_GAMA5["be"]
    )
    if match_A:
        ms = ms + np.log10(A_SCALE / DRIVER_GAMA5_A)
    return (ms, lp, al, be)


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
    z_obs,
    m_obs,
    nbin_z=30,
    hist_bw=0.3,
    min_in_bin=20,
    zmin=ZMIN,
    zmax=ZLIMIT,
    form=None,
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

    # `form` forces the functional form. The AIC choice is unstable: a 0.075 dex
    # shift in the masses (rebuilt -> MassA) flipped it linear -> quad and moved
    # mlim by 0.69 dex at low z. C is tabulated against Delta = m - mlim(z) and
    # the mock's mlim is linear, so a flip here silently redefines Delta.
    use_quad = (aic_q < aic_l - 2) if form is None else (form == "quad")
    if use_quad:
        c = beta_q
        func = lambda z: c[0] + c[1] * z + c[2] * z**2
        kind = "quad"
    else:
        c = beta_l
        func = lambda z: c[0] + c[1] * z
        kind = "linear"
    kind += (
        " (forced)"
        if form is not None
        else f" [AIC lin {aic_l:.1f} / quad {aic_q:.1f}]"
    )
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
  ms ~ normal(13.958, 0.42);   // Driver+22 M* (broad -> data-driven)
  lp ~ normal(-3.445, 0.69);   // Driver+22 logphi* (broad -> data-driven)
  al ~ normal(-1.68, 0.22);   // Driver+22 alpha (informative)
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width (+0.25/-0.11)

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
  ms ~ normal(13.958, 0.42);   // Driver+22 M* (broad -> data-driven)
  lp ~ normal(-3.445, 0.69);   // Driver+22 logphi* (broad -> data-driven)
  al ~ normal(-1.68, 0.22);   // Driver+22 alpha (informative)
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width (+0.25/-0.11)

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
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width (+0.25/-0.11)
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
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width (+0.25/-0.11)

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
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width (+0.25/-0.11)
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
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);    // Driver+22 beta, published width
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

# Tabulated-completeness model. C(Delta,z) comes from running Nessie on the
# mock (measure_completeness_nessie.py) and is NOT a parametric ramp -- the
# measured curve saturates below 1 in some z bins and above 1 in others
# (fragmentation), which no erf can represent. Because C does not depend on the
# MRP parameters it is evaluated in Python on the integration grids and passed
# in as matrices, so the sampler does no interpolation and no erf calls.
MARG_TAB_CODE = r"""
data {
  int<lower=1> N;
  int<lower=2> Nint;
  int<lower=1> Nsh;
  int<lower=2> Ng;
  vector[N] x_obs;
  vector<lower=0>[N] sig;
  matrix[N, Nint] mt;          // per-object integration nodes (latent mass)
  matrix[N, Nint] Cobj;        // C at those nodes
  vector[Nsh] V_sh;
  vector[Ng] xg;               // Lambda grid
  matrix[Nsh, Ng] Csh;         // C on the Lambda grid, per shell
  real dx;
}
transformed data {
  real ln10 = log(10.0);
  real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
}
model {
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);

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
    real inv_s = 1.0 / sig[i];
    real dmt = mt[i, 2] - mt[i, 1];
    real sm = 0;
    for (g in 1:Nint) {
      real u = mt[i, g] - ms;
      real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
      real zsc = (x_obs[i] - mt[i, g]) * inv_s;
      real term = phi_g * Cobj[i, g] * exp(-0.5 * zsc * zsc);
      sm += (g == 1 || g == Nint) ? 0.5 * term : term;
    }
    sm *= dmt * inv_s * inv_sqrt2pi;
    target += (sm > 1e-300) ? log(sm) : -300;
  }
}
"""
_STAN["marg_tab"] = MARG_TAB_CODE

# As marg_tab, but with the error scale fitted. The trick that keeps it fast:
# the integration NODES are fixed on a wide grid (+-5*S_HI sigma), so Cobj can
# still be precomputed in Python; s_scale only changes the width of the Gaussian
# kernel evaluated on those nodes. The prior is centred on the mock-measured
# ratio (0.375/0.280 = 1.34) rather than on 1.
MARG_TAB_SERR_CODE = r"""
data {
  int<lower=1> N;
  int<lower=2> Nint;
  int<lower=1> Nsh;
  int<lower=2> Ng;
  vector[N] x_obs;
  vector<lower=0>[N] sig;      // catalogue sigma, UNSCALED
  matrix[N, Nint] mt;          // fixed integration nodes
  matrix[N, Nint] Cobj;
  vector[Nsh] V_sh;
  vector[Ng] xg;
  matrix[Nsh, Ng] Csh;
  real dx;
  real s_mu;                   // prior median for s_scale
  real s_sd;
}
transformed data {
  real ln10 = log(10.0);
  real inv_sqrt2pi = 1.0 / sqrt(2 * pi());
}
parameters {
  real ms;
  real lp;
  real al;
  real<lower=0.1, upper=2.0> be;
  real<lower=0.4, upper=2.6> s_scale;
}
model {
  ms ~ normal(13.958, 0.42);
  lp ~ normal(-3.445, 0.69);
  al ~ normal(-1.68, 0.22);
  be ~ normal(0.63, 0.18);
  s_scale ~ lognormal(log(s_mu), s_sd);

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
    real inv_s = 1.0 / si;
    real dmt = mt[i, 2] - mt[i, 1];
    real sm = 0;
    for (g in 1:Nint) {
      real u = mt[i, g] - ms;
      real phi_g = be * ln10 * pow(10, lp) * pow(10, (al + 1) * u) * exp(-pow(10, be * u));
      real zsc = (x_obs[i] - mt[i, g]) * inv_s;
      real term = phi_g * Cobj[i, g] * exp(-0.5 * zsc * zsc);
      sm += (g == 1 || g == Nint) ? 0.5 * term : term;
    }
    sm *= dmt * inv_s * inv_sqrt2pi;
    target += (sm > 1e-300) ? log(sm) : -300;
  }
}
"""
_STAN["marg_tab_serr"] = MARG_TAB_SERR_CODE

# marg_tab with M* PINNED. M* is exactly degenerate with the mass calibration A
# (M = A sigma^2 R / G), so fixing it asserts the calibration and asks what GAMA
# says about the remaining shape. ms is passed as DATA rather than as a parameter
# with a spike prior -- a near-zero-width prior collapses HMC's step size.
MARG_TAB_PINMS_CODE = (
    MARG_TAB_CODE.replace(
        "  real dx;\n}",
        "  real dx;\n  real ms;          // PINNED, supplied as data\n}",
        1,
    )
    .replace("parameters {\n  real ms;\n", "parameters {\n", 1)
    .replace("  ms ~ normal(13.958, 0.42);\n", "", 1)
)
_STAN["marg_tab_pinms"] = MARG_TAB_PINMS_CODE

# marg_tab with BETA pinned. Unlike M*, beta is not degenerate with the mass
# calibration -- it is genuinely unconstrained because GAMA has almost no groups
# above 10^15 and the completeness table above Delta=+1 rests on 17 halos. Fixing
# it to a literature value is defensible, but note it is correlated with alpha
# (rho ~ -0.8), so pinning shifts alpha and shrinks its error bar by ~40%: the
# resulting precision on alpha is conditional on beta being exactly right.
MARG_TAB_PINBE_CODE = (
    MARG_TAB_CODE.replace(
        "  real dx;\n}",
        "  real dx;\n  real be;          // PINNED, supplied as data\n}",
        1,
    )
    .replace("  real<lower=0.1, upper=2.0> be;\n}", "}", 1)
    .replace("  be ~ normal(0.63, 0.18);\n", "", 1)
)
_STAN["marg_tab_pinbe"] = MARG_TAB_PINBE_CODE

# marg_tab with ALPHA pinned. Note this pins the one parameter the closed-loop
# validation shows is essentially unbiased (+0.011 dex, 0.11 sigma), so it is a
# sensitivity test rather than a preferred configuration.
MARG_TAB_PINAL_CODE = (
    MARG_TAB_CODE.replace(
        "  real dx;\n}",
        "  real dx;\n  real al;          // PINNED, supplied as data\n}",
        1,
    )
    .replace(
        "  real al;\n  real<lower=0.1, upper=2.0> be;",
        "  real<lower=0.1, upper=2.0> be;",
        1,
    )
    .replace("  al ~ normal(-1.68, 0.22);\n", "", 1)
)
_STAN["marg_tab_pinal"] = MARG_TAB_PINAL_CODE


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
    # The chains are in Driver's h=0.6737 units; this pipeline works at h=1,
    # so shift log M* and log phi* (alpha and beta are dimensionless).
    d = d.copy()
    d[:, 0] -= np.log10(1.0 / H_DRIVER)
    d[:, 1] -= 3.0 * np.log10(H_DRIVER)
    mu = d.mean(axis=0)
    cov = np.cov(d.T) * float(inflate) ** 2
    return mu, cov, d


def _with_driver_prior(code):
    """Swap a model's independent-Gaussian prior block for a multivariate normal."""
    old = (
        "  ms ~ normal(13.958, 0.42);\n"
        "  lp ~ normal(-3.445, 0.69);\n"
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
        pin0 = {"_pinms": "ms", "_pinbe": "be", "_pinal": "al"}.get(kind[-6:])
        if pin0:
            init.pop(pin0, None)
        try:
            opt = model.optimize(
                data=data, inits=init, algorithm="lbfgs", iter=20000, seed=seed
            )
            val = float(opt.optimized_params_dict["lp__"])
            if best is None or val > best[0]:
                pin0 = {"_pinms": "ms", "_pinbe": "be", "_pinal": "al"}.get(kind[-6:])
                best = (
                    val,
                    np.array(
                        [
                            float(data[p])
                            if p == pin0
                            else float(opt.optimized_params_dict[p])
                            for p in PARAMS
                        ]
                    ),
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

    pin1 = {"_pinms": "ms", "_pinbe": "be", "_pinal": "al"}.get(kind[-6:])
    if pin1:
        inits = [{k: v for k, v in d.items() if k != pin1} for d in inits]
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

    pinned = {"_pinms": "ms", "_pinbe": "be", "_pinal": "al"}.get(kind[-6:])
    if pinned:
        n_draw = int(chains) * int(sampling)
        flat = np.column_stack(
            [
                np.full(n_draw, float(data[pinned]))
                if p == pinned
                else fit.stan_variable(p)
                for p in PARAMS
            ]
        )
    else:
        flat = np.column_stack([fit.stan_variable(p) for p in PARAMS])

    # diagnostics (column names vary slightly across cmdstanpy versions)
    try:
        summ = fit.summary()
        free = [
            p
            for p in PARAMS
            if p != {"_pinms": "ms", "_pinbe": "be", "_pinal": "al"}.get(kind[-6:])
        ]
        rhat = float(summ.loc[free, "R_hat"].max())
        ess_col = next((c for c in ("ESS_bulk", "N_Eff") if c in summ.columns), None)
        ess = float(summ.loc[free, ess_col].min()) if ess_col else float("nan")
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
    if abs(H0 - 100.0) < 1e-6:
        dms, dlp, _, _ = to_driver_cosmology(med[0], med[1])
        sms = 0.5 * (q84[0] - q16[0])
        slp = 0.5 * (q84[1] - q16[1])
        g5 = driver_gama5(match_A=True)
        print(
            f"  [Driver+22 GAMA-only (GAMA5), h=1 and rescaled to A={A_SCALE:g}: "
            f"{g5[0]:.3f} / {g5[1]:.3f} / {g5[2]:.3f} / {g5[3]:.3f}]"
        )
        print(
            f"  [h=1 units. Converted to Driver's h={H_DRIVER}: "
            f"log M* = {dms:.3f} +/- {sms:.3f}, log phi* = {dlp:.3f} +/- {slp:.3f}"
            f"   (Driver GSR: 14.13, -3.96)]"
        )
    if not SHOW_ALPHA_CORRECTION:
        return dict(median=med, sd=sd, q16=q16, q84=q84)
    ac = corrected_alpha(flat)
    acm = float(np.median(ac))
    tot = float(np.hypot(np.std(ac), ALPHA_BIAS_SCATTER))
    print(
        f"  {'al_corr':9s} {TRUE['al']:8.3f} {acm:9.3f} {tot:7.3f} "
        f"{np.percentile(ac, 16):8.3f} {np.percentile(ac, 84):8.3f} "
        f"{(acm - TRUE['al']) / np.std(ac):+9.2f}   <- mock bias-corrected"
    )
    print(
        f"  [alpha bias at M*={med[0]:.2f} is {float(alpha_bias(med[0])):+.3f} dex; "
        f"sd combines stat {np.std(ac):.3f} + calib {ALPHA_BIAS_SCATTER:.3f} (20 mocks)]"
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
    if mlim_func is not None:
        a.plot(z_plot, mlim_func(z_plot), "r-", lw=2, label="mlim(z)")
    elif COMP_MODE == "mz":
        # No mlim(z) exists in this mode -- the selection is C(m, z).  Draw the
        # 50% and 10% completeness masses instead, which is the same information
        # in the form the model actually uses.
        try:
            _t = np.load(NESSIE_TABLE_MZ)
            _m, _z, _C = _t["m"], _t["z"], _t["C"]
            for _lvl, _ls in ((0.5, "r-"), (0.1, "r--")):
                _mm = [
                    np.interp(_lvl, _C[i], _m) if _C[i].max() >= _lvl else np.nan
                    for i in range(_z.size)
                ]
                a.plot(_z, _mm, _ls, lw=2, label=f"C = {_lvl:g}")
        except Exception as _e:
            print(f"  [C(m,z) overlay skipped: {_e}]")
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
    prior_mu = [13.958, -3.445, -1.68, 0.63]  # must match the Stan model priors
    prior_sd = [0.42, 0.69, 0.22, 0.18]
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
    prior_sd = {"ms": 0.42, "lp": 0.69, "al": 0.22, "be": 0.18}
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


# ---------------------------------------------------------------------------
# Tabulated completeness from measure_completeness_nessie.py
# ---------------------------------------------------------------------------
# The mock says the run.R lookup table understates the mass errors: the measured
# scatter of recovered vs true mass is 0.375 dex against the table's 0.280, a
# ratio of 1.34. The fitted s_scale hyperparameter independently gave
# 1.367 +/- 0.042. Applying the measured factor lets sigma be fixed data, which
# means the completeness can be precomputed and the sampler runs far faster.
MLIM_FORM = None  # None = choose by AIC; 'linear'/'quad' to force
PIN_MS = 13.958  # Driver+22 GSR M* in h=1 units, for marg_tab_pinms
PIN_BE = 0.63  # Driver+22 GSR beta, for marg_tab_pinbe
PIN_AL = -1.68  # Driver+22 GSR alpha, for marg_tab_pinal
SIGMA_SCALE = 1.34

# Closed-loop bias of the marg_tab fit, from run_mock_nessie: the MRP is injected
# by abundance matching, Nessie recovers the groups, the fit is run with the same
# C table, and this is (fitted - injected). Subtracting it turns the fit into an
# estimate of the underlying HALO mass function.
#   ms: -0.229, of which -0.155 is the A=10 dynamical mass under-estimating the
#       true halo mass; the remaining -0.074 dex is 0.34 sigma, i.e. consistent
#       with no residual bias.
#   lp: +0.418, essentially ms riding the rho=-0.97 M*-phi* ridge.
#   al: +0.011 (0.11 sigma) -- alpha is recovered essentially unbiased.
# ONE realisation, so these carry no error bar yet. Jackknifing the mock
# footprint would give one.
NESSIE_BIAS = dict(ms=-0.229, lp=+0.418, al=+0.011, be=-0.037)


def apply_nessie_bias(flat):
    """Subtract the closed-loop bias from every draw, so the covariance is
    carried through rather than the median being shifted on its own."""
    out = np.array(flat, float, copy=True)
    for i, k in enumerate(PARAMS):
        out[:, i] -= NESSIE_BIAS[k]
    return out


NESSIE_TABLE = "nessie_completeness_table.npz"
NESSIE_TABLE_MZ = "nessie_completeness_mz.npz"

# 'delta' = the incumbent C(Delta = m - mlim(z), z), which needs turnover_mlim.
# 'mz'    = C(m, z) keyed on absolute true mass, so mlim(z) drops out of the fit
#           entirely.  Set from --comp-mode in the __main__ block, alongside
#           MLIM_FORM -- run_real_gama reads these as globals at call time.
COMP_MODE = "delta"
COMP_DEF = "entries"  # 'entries' (Poisson intensity) | 'repr' (bounded, check)

# The C(m,z) table is keyed on the TRUE abundance-matched halo mass, while
# x_obs is MassA, a dynamical mass low by ~0.155 dex (A=10 and no h-scaling).
# Under Delta keying that offset partly cancelled, because mlim was fit to
# OBSERVED masses in both mock and data, and the remainder was absorbed into
# NESSIE_BIAS post hoc.  Under absolute keying nothing cancels, so it has to be
# applied explicitly.  N(m_obs | m_t + b, sig) == N(m_obs - b | m_t, sig), so
# shifting x_obs is exact and needs no Stan change.
# measured by measure_completeness_nessie.py: median(log_dyn - log_mass_am)
# over 1132 clean matches = -0.155 dex, i.e. the dynamical mass sits BELOW the
# true abundance-matched mass, so x_obs - MASS_BIAS shifts the data UP onto the
# table's scale.  (The same run gives scatter 0.375/0.280 = 1.34 = SIGMA_SCALE.)
MASS_BIAS = -0.155

# Closed-loop residual for COMP_MODE='mz'.  Still zeros until the mock recovery
# in that mode has actually been run -- apply_nessie_bias warns if it is used.
NESSIE_BIAS_MZ = dict(ms=0.0, lp=0.0, al=0.0, be=0.0)


def load_nessie_table_mz(path=NESSIE_TABLE_MZ, comp_def="entries"):
    """(m, z, C, Nh) from measure_completeness_nessie.tabulate_C_mz.

    ``comp_def='entries'`` returns the expected NUMBER of catalogue entries per
    halo, which is what the Poisson intensity wants and which may exceed 1.
    ``'repr'`` returns the bounded recovered-fraction, for the systematic check.
    """
    t = np.load(path)
    m, z = t["m"], t["z"]
    C = t["C_repr"] if comp_def == "repr" else t["C"]
    Nh = t["Nh"] if "Nh" in t else np.zeros(C.shape, int)
    thin = int((Nh < int(t["min_n"])).sum()) if "min_n" in t else -1
    print(
        f"  [C table] {path} [{comp_def}]: {C.shape[0]} z bins "
        f"({z[0]:.3f}..{z[-1]:.3f}) x {m.size} mass points "
        f"({m[0]:.2f}..{m[-1]:.2f})"
    )
    print(
        f"            plateau {np.round(np.nanmax(C, axis=1), 3)}"
        + (f", {thin}/{Nh.size} cells below min_n" if thin >= 0 else "")
    )
    return m, z, C, Nh


def load_nessie_table(path=NESSIE_TABLE):
    """(d, z, C) from measure_completeness_nessie.py. C is entries-per-halo, so
    it may exceed 1 where the finder fragments a halo into several groups."""
    t = np.load(path)
    d, z, C = t["d"], t["z"], t["C"]
    print(
        f"  [C table] {path}: {C.shape[0]} z bins x {d.size} Delta points, "
        f"C(Delta=0) = {np.round(np.interp(0.0, d, C[C.shape[0] // 2]), 3)}"
    )
    return d, z, C


def eval_C(delta, zval, d, z, C):
    """Bilinear lookup of C(Delta, z). Clamped: 0 below the measured Delta
    range, held flat above it and outside the z range."""
    delta = np.asarray(delta, float)
    zval = np.broadcast_to(np.asarray(zval, float), delta.shape)
    iz = np.clip(np.searchsorted(z, zval) - 1, 0, len(z) - 2)
    t = np.clip((zval - z[iz]) / (z[iz + 1] - z[iz]), 0.0, 1.0)
    out = np.empty(delta.shape)
    for i in range(len(z) - 1):
        m = iz == i
        if not m.any():
            continue
        c0 = np.interp(delta[m], d, C[i], left=0.0, right=C[i][-1])
        c1 = np.interp(delta[m], d, C[i + 1], left=0.0, right=C[i + 1][-1])
        out[m] = (1 - t[m]) * c0 + t[m] * c1
    return np.clip(out, 0.0, None)


def _prep_tab_mz(
    z_obs, m_obs, sigma, z_mids, V_sh, table_path_mz, nint, cmin,
    fit_scale, s_hi, mbias, comp_def,
):
    """prep_tab with the completeness keyed on absolute mass, C(m, z).

    Identical machinery to the Delta path -- same nodes, same Stan data block,
    no new Stan variables -- with two deliberate differences:

    * ``x_obs`` is shifted by ``mbias``. The table is keyed on the TRUE
      abundance-matched mass while the data is MassA, ~0.155 dex lower. Under
      Delta keying that partly cancelled through mlim; here nothing cancels, and
      N(m_obs | m_t + b, sig) == N(m_obs - b | m_t, sig) makes the shift exact.
    * the keep mask is ``Cobj.max(axis=1) > cmin`` -- "the model can explain this
      object somewhere on its node grid" -- rather than C at the observed mass.
      Testing C(x_obs) would reintroduce a hard cut on the observed value, which
      is the exact failure being removed.
    """
    mtab, ztab, C, Nh = load_nessie_table_mz(table_path_mz, comp_def)
    if fit_scale:
        sig = np.asarray(sigma, float)
        nint = max(nint, int(np.ceil(41 * s_hi)))
        span = 5.0 * s_hi
    else:
        sig = np.asarray(sigma, float) * SIGMA_SCALE
        span = 5.0

    xc = np.asarray(m_obs, float) - mbias
    g = np.linspace(-span, span, nint)
    mt_all = xc[:, None] + sig[:, None] * g[None, :]
    C_all = eval_C(
        mt_all, np.repeat(np.asarray(z_obs, float)[:, None], nint, axis=1),
        mtab, ztab, C,
    )
    keep = C_all.max(axis=1) > cmin
    mt, Cobj = mt_all[keep], C_all[keep]

    xg = np.linspace(float(mtab[0]), XHI, NG)
    dx = float(xg[1] - xg[0])
    Csh = np.vstack(
        [eval_C(xg, np.full(NG, z_mids[j]), mtab, ztab, C) for j in range(len(V_sh))]
    )

    dropped = int((~keep).sum())
    beyond = int((xc[keep] > mtab[-1]).sum())
    print(
        f"  tabulated C(m,z): kept {int(keep.sum())}/{xc.size} groups "
        f"(C > {cmin}); mass shifted by {-mbias:+.3f} dex; "
        + ("sigma FITTED" if fit_scale else f"sigma x{SIGMA_SCALE}")
    )
    if dropped:
        print(f"  !! {dropped} groups dropped -- expected ~0 in mz mode; if this "
              f"is large the table does not cover the data")
    if beyond:
        print(
            f"  !! {beyond} groups sit above logM = {mtab[-1]:.2f}, the top of the "
            f"table -- C is held flat there, and that is where M* and beta are set."
        )
    data = dict(
        N=int(keep.sum()), Nint=int(nint), Nsh=int(len(V_sh)), Ng=int(NG),
        x_obs=xc[keep], sig=sig[keep], mt=mt, Cobj=Cobj,
        V_sh=np.asarray(V_sh, float), xg=xg, Csh=Csh, dx=dx,
    )
    if fit_scale:
        data.update(s_mu=float(SIGMA_SCALE), s_sd=0.15)
    return data, keep


def prep_tab(
    z_obs,
    m_obs,
    sigma,
    mlim_func,
    z_mids,
    mlim_sh,
    V_sh,
    table_path=NESSIE_TABLE,
    nint=41,
    cmin=1e-3,
    fit_scale=False,
    s_hi=2.6,
    *,
    mode_key=None,
    mbias=None,
    comp_def=None,
    table_path_mz=None,
):
    """Stan data for marg_tab: integration nodes and the completeness evaluated
    on them, both precomputed. Keeps every group whose own completeness is
    non-negligible (no observed-mass cut: the selection is carried by C).

    ``mode_key='mz'`` switches the completeness from C(Delta = m - mlim(z), z) to
    C(m, z) keyed on absolute true mass, which removes ``turnover_mlim`` -- and
    with it the mode-as-limit problem, the 54% hard cut and the linear/quad
    instability -- from the fit entirely. In that mode ``mlim_func`` and
    ``mlim_sh`` may be None, and ``x_obs`` is shifted by ``mbias`` so the data
    and the table share a mass definition (see MASS_BIAS)."""
    mode_key = COMP_MODE if mode_key is None else mode_key
    mbias = MASS_BIAS if mbias is None else float(mbias)
    comp_def = COMP_DEF if comp_def is None else comp_def
    table_path_mz = NESSIE_TABLE_MZ if table_path_mz is None else table_path_mz
    if mode_key == "mz":
        return _prep_tab_mz(
            z_obs, m_obs, sigma, z_mids, V_sh, table_path_mz, nint, cmin,
            fit_scale, s_hi, mbias, comp_def,
        )
    d, ztab, C = load_nessie_table(table_path)
    if fit_scale:
        # nodes must cover the widest kernel the sampler can reach, so they are
        # laid out at +-5*s_hi sigma and held FIXED; s_scale then only alters the
        # Gaussian width, leaving Cobj precomputable.
        sig = np.asarray(sigma, float)
        nint = max(nint, int(np.ceil(41 * s_hi)))
        span = 5.0 * s_hi
    else:
        sig = np.asarray(sigma, float) * SIGMA_SCALE
        span = 5.0

    mlim_obj = mlim_func(z_obs)
    keep = eval_C(m_obs - mlim_obj, z_obs, d, ztab, C) > cmin
    x, sg, ml, zk = m_obs[keep], sig[keep], mlim_obj[keep], z_obs[keep]

    g = np.linspace(-span, span, nint)
    mt = x[:, None] + sg[:, None] * g[None, :]
    Cobj = eval_C(mt - ml[:, None], np.repeat(zk[:, None], nint, axis=1), d, ztab, C)

    xlo = float(np.min(mlim_sh) + d[0])
    xg = np.linspace(xlo, XHI, NG)
    dx = float(xg[1] - xg[0])
    Csh = np.vstack(
        [
            eval_C(xg - mlim_sh[j], np.full(NG, z_mids[j]), d, ztab, C)
            for j in range(len(V_sh))
        ]
    )

    mode = (
        f"sigma FITTED (nodes +-{span:.0f}sigma, {nint} pts, "
        f"prior median {SIGMA_SCALE})"
        if fit_scale
        else f"sigma x{SIGMA_SCALE} (mock-measured, fixed)"
    )
    beyond = int((m_obs[keep] - mlim_obj[keep] > d.max()).sum())
    print(
        f"  tabulated C: kept {int(keep.sum())}/{m_obs.size} groups "
        f"(C > {cmin}); {mode}"
    )
    if beyond:
        print(
            f"  !! {beyond} groups ({beyond / max(int(keep.sum()), 1):.1%}) sit above "
            f"Delta = {d.max():+.2f}, the top of the measured table -- C is held "
            f"flat there, and that is where M* and beta are set."
        )
    data = dict(
        N=int(keep.sum()),
        Nint=int(nint),
        Nsh=int(len(V_sh)),
        Ng=int(NG),
        x_obs=x,
        sig=sg,
        mt=mt,
        Cobj=Cobj,
        V_sh=np.asarray(V_sh, float),
        xg=xg,
        Csh=Csh,
        dx=dx,
    )
    if fit_scale:
        data.update(s_mu=float(SIGMA_SCALE), s_sd=0.15)
    return data, keep


def check_lambda(data, par, label=""):
    """Lambda at the fitted parameters vs the number of groups actually fitted.
    Lambda is the model's expected count, so Lambda/N ~ 1 means the normalisation
    is self-consistent. It also reports the halo count the same MRP implies with
    no completeness correction -- the gap between the two is how much of phi* is
    being set by the completeness rather than by the data."""
    ms, lp, al, be = [float(v) for v in par[:4]]
    xg = np.asarray(data["xg"], float)
    Csh, V_sh = np.asarray(data["Csh"], float), np.asarray(data["V_sh"], float)
    pg = mrp_phi(xg, ms, lp, al, be)
    lam = float(sum(V_sh[j] * np.trapezoid(pg * Csh[j], xg) for j in range(len(V_sh))))
    lam0 = float(np.sum(V_sh) * np.trapezoid(pg, xg))
    n = int(data["N"])
    print(
        f"  [Lambda{label}] expected detections = {lam:.0f}, fitted N = {n}, "
        f"ratio = {lam / max(n, 1):.3f}"
    )
    print(
        f"            same MRP with C=1: {lam0:.0f} halos "
        f"({lam0 / max(n, 1):.1f}x the catalogue)"
    )
    return lam


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


def load_real_gama(
    fits_path, regions=None, use_veldisp_err=False, dec_cut=None, mass_col=None
):
    """Read the GAMA G3C group catalogue and build observed log-masses exactly
    as run.R does: A=13.9 dynamical mass, multiplicity error model, and the
    empirical masscorr(Nfof) calibration. Returns (log_mass, sigma, z, Nfof).
    Column names/units follow run.R; adjust here if the FITS schema differs."""
    from astropy.io import fits as afits

    parsec, Gnewton, msol = 3.0857e16, 6.67408e-11, 1.988e30
    with afits.open(fits_path) as hdul:
        t = hdul[1].data
    cols = set(t.columns.names)
    Nfof = np.asarray(t["Nfof"], float)
    Zfof = np.asarray(t["Zfof"], float)
    MassAfunc = np.asarray(t["MassAfunc"], float)
    VelDisp = np.asarray(t["VelDisp"], float)
    Rad50 = np.asarray(t["Rad50"], float)
    IterCenDec = np.asarray(t["IterCenDec"], float)
    VelDispErr = np.asarray(t["VelDispErr"], float) if "VelDispErr" in cols else None
    region = np.asarray(t["GAMARegion"]).astype(str) if "GAMARegion" in cols else None

    sel = (Nfof > MULTI - 1) & (Zfof < ZLIMIT) & (Zfof > ZMIN) & (MassAfunc > 1e1)

    # Region selection. The old hardcoded `IterCenDec > -3.5` picks out the
    # equatorial fields (G09/G12/G15) and would DELETE G23 (Dec ~ -35..-30).
    # Prefer an explicit GAMARegion cut when that column exists.
    if regions is not None:
        if region is None:
            raise KeyError("regions= requested but no GAMARegion column in this file")
        sel &= np.isin(region, [str(r) for r in regions])
    elif dec_cut is not None:
        sel &= IterCenDec > float(dec_cut)

    if region is not None:
        import collections

        print(
            f"  regions kept: {dict(sorted(collections.Counter(region[sel]).items()))}"
        )

    Nfof, Zfof, VelDisp, Rad50 = Nfof[sel], Zfof[sel], VelDisp[sel], Rad50[sel]
    if VelDispErr is not None:
        VelDispErr = VelDispErr[sel]

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

    # Optional: propagate the catalogue's own velocity-dispersion error instead
    # of the multiplicity lookup. M ~ A sigma^2 R / G, so ignoring the radius
    # error, sigma_logM = (2/ln10) * (VelDispErr / VelDisp).
    if use_veldisp_err:
        if VelDispErr is None:
            raise KeyError("use_veldisp_err=True but no VelDispErr column here")
        with np.errstate(divide="ignore", invalid="ignore"):
            err_vd = (2.0 / np.log(10.0)) * (VelDispErr / VelDisp)
        bad = ~np.isfinite(err_vd) | (err_vd <= 0)
        err_vd = np.where(bad, err, err_vd)  # fall back to the table
        err_vd = np.maximum(err_vd, 0.05)  # floor
        print(
            f"  [errors] VelDispErr -> med {np.median(err_vd):.3f} dex "
            f"(lookup table med {np.median(err):.3f}); "
            f"{int(bad.sum())}/{bad.size} fell back to the table"
        )
        err = err_vd

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

    # The completeness table is measured against MassA = mass_proxy * A, with no
    # masscorr. Reading that column directly puts the fit on exactly the same
    # mass definition as the mock; the rebuilt path above applies Driver's old
    # masscorr, which the new DMU replaces with its own MassAfunc.
    if mass_col is not None:
        if mass_col not in cols:
            raise KeyError(f"{mass_col!r} not in this catalogue; have {sorted(cols)}")
        v = np.asarray(t[mass_col], float)[sel]
        with np.errstate(divide="ignore", invalid="ignore"):
            lm_cat = np.log10(v)
        f1 = np.isfinite(lm_cat)
        f2 = np.isfinite(log_mass)
        print(
            f"  [mass] catalogue column {mass_col!r}: median "
            f"{np.median(lm_cat[f1]):.3f}   (rebuilt from VelDisp/Rad50 would "
            f"give {np.median(log_mass[f2]):.3f})"
        )
        log_mass = lm_cat

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
    regions=None,
    use_veldisp_err=False,
    dec_cut=None,
    mass_col=None,
):
    """Fit the MRP to the REAL GAMA group catalogue using OUR developed model
    ('marg' by default, with the current tight cutoff/slope priors; 'gama' =
    verbatim R port, kept only for a pure port-check). load_real_gama does the
    run.R data prep (A=13.9 masses, multiplicity error model); the FIT is our
    marginalised model. 'bias vs truth' and the plot 'truth' lines are Driver+22,
    so read them as 'offset vs Driver+22' -- with tight priors, be/al are
    prior-set near Driver by construction; M* and logphi* are the measurements."""
    print(f"Reading GAMA catalogue: {fits_path}")
    log_mass, sigma, z, nfof = load_real_gama(
        fits_path,
        regions=regions,
        use_veldisp_err=use_veldisp_err,
        dec_cut=dec_cut,
        mass_col=mass_col,
    )
    sky_frac = sky_area_deg2_val * (np.pi / 180) ** 2 / (4 * np.pi)
    print(
        f"  N groups: {log_mass.size}   mass {log_mass.min():.2f}..{log_mass.max():.2f} "
        f"(med {np.median(log_mass):.2f})   sigma med {np.median(sigma):.2f}"
    )

    _mz = COMP_MODE == "mz" and model_kind.startswith("marg_tab")
    if _mz:
        # The whole point of mz mode: the selection is carried by C(m, z), so
        # mlim(z) -- the histogram-mode estimator that put 54% of the catalogue
        # below its own "limit" -- is never formed.
        print("mlim(z) NOT used: selection carried by C(m,z) [--comp-mode mz]")
        mlim_func = coefs = tkind = turn_pts = None
        z_mids, V_sh = shell_volumes(sky_frac)
        Vsurvey = float(V_sh.sum())
        mlim_sh = mlim_per = None
    else:
        print("Turnover mlim(z) ...")
        mlim_func, coefs, tkind, turn_pts = turnover_mlim(z, log_mass, form=MLIM_FORM)
        print(
            f"  mlim(z) [{tkind}]: mlim({ZMIN})={mlim_func(ZMIN):.2f} "
            f"mlim({ZLIMIT})={mlim_func(ZLIMIT):.2f}"
        )

        z_mids, V_sh = shell_volumes(sky_frac)
        Vsurvey = float(V_sh.sum())
        mlim_sh = mlim_func(z_mids)

        mlim_per = mlim_func(z)

    if model_kind in (
        "marg_tab",
        "marg_tab_serr",
        "marg_tab_pinms",
        "marg_tab_pinbe",
        "marg_tab_pinal",
    ):
        data, keep = prep_tab(
            z,
            log_mass,
            sigma,
            mlim_func,
            z_mids,
            mlim_sh,
            V_sh,
            fit_scale=(model_kind == "marg_tab_serr"),
        )
        if model_kind == "marg_tab_pinal":
            data["al"] = float(PIN_AL)
            print(
                f"  alpha PINNED at {PIN_AL:.3f}. Note the closed loop shows alpha"
                f" is the LEAST biased\n    parameter (+0.011 dex), so this trades"
                f" away the best-validated result."
            )
        if model_kind == "marg_tab_pinbe":
            data["be"] = float(PIN_BE)
            print(
                f"  beta PINNED at {PIN_BE:.3f}  (Driver GSR 0.63, GAMA-only 0.47,"
                f" LCDM 0.71).\n"
                f"    beta correlates with alpha at rho ~ -0.8, so alpha's error"
                f" bar below is CONDITIONAL on this value."
            )
        if model_kind == "marg_tab_pinms":
            data["ms"] = float(PIN_MS)
            print(
                f"  M* PINNED at {PIN_MS:.3f} (h=1) = "
                f"{to_driver_cosmology(PIN_MS, 0)[0]:.3f} in Driver units;\n"
                f"    equivalent to asserting the mass calibration A"
            )
        # In mz mode data["x_obs"] is log_mass - MASS_BIAS, i.e. the TRUE-mass
        # coordinate the model and xg live in. plot_ppc/emit_publication bin
        # x_fit against a model integral on xg, so using the raw MassA here
        # would shift every figure by MASS_BIAS relative to the fit.
        x_fit = np.asarray(data["x_obs"], float) if _mz else log_mass[keep]
    elif model_kind in ("marg_comp", "marg_comp_serr"):
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
    if model_kind in ("marg_tab", "marg_tab_serr"):
        try:
            check_lambda(data, np.median(flat, axis=0))
        except Exception as e:
            print(f"  [Lambda check failed: {e}]")
        fb = apply_nessie_bias(flat)
        mb, sb = (
            np.median(fb, axis=0),
            0.5 * (np.percentile(fb, 84, axis=0) - np.percentile(fb, 16, axis=0)),
        )
        dms, dlp, _, _ = to_driver_cosmology(mb[0], mb[1])
        print(f"\n  === closed-loop bias-corrected (halo mass function) ===")
        print(f"  {'param':6} {'h=1':>18} {'Driver units':>18} {'Driver+22':>10}")
        print(
            f"  {'ms':6} {mb[0]:9.3f} +/-{sb[0]:5.3f} {dms:9.3f} +/-{sb[0]:5.3f}"
            f" {TRUE_DRIVER['ms']:10.2f}"
        )
        print(
            f"  {'lp':6} {mb[1]:9.3f} +/-{sb[1]:5.3f} {dlp:9.3f} +/-{sb[1]:5.3f}"
            f" {TRUE_DRIVER['lp']:10.2f}"
        )
        print(
            f"  {'al':6} {mb[2]:9.3f} +/-{sb[2]:5.3f} {mb[2]:9.3f} +/-{sb[2]:5.3f}"
            f" {TRUE_DRIVER['al']:10.2f}"
        )
        print(
            f"  {'be':6} {mb[3]:9.3f} +/-{sb[3]:5.3f} {mb[3]:9.3f} +/-{sb[3]:5.3f}"
            f" {TRUE_DRIVER['be']:10.2f}"
        )
        print(f"  [bias from one mock realisation, no error bar on the correction]")
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
    np.savetxt(
        f"gama_{model_kind}_draws.csv",
        flat,
        delimiter=",",
        header="ms,lp,al,be",
        comments="",
    )
    print(f"  saved draws -> gama_{model_kind}_draws.csv")
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
    if model_kind.startswith("marg_tab"):
        try:
            plot_ppc(
                flat,
                data,
                x_fit,
                fname=f"ppc_gama_{model_kind}.pdf",
                title=f"GAMA: observed vs predicted detections [{model_kind}]",
            )
        except Exception as e:
            print(f"  [ppc failed: {e}]")
    emit_publication(
        flat,
        {"GAMA": dict(x_fit=x_fit, Vsurvey=Vsurvey)},
        tag=f"gama_{model_kind}",
        title="GAMA HMF",
        nessie_bias=(
            SHOW_HALO_MF_CURVE and model_kind in ("marg_tab", "marg_tab_serr")
        ),
        mmax_data=float(np.percentile(x_fit, 99.5)),
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
    prior_mu, prior_sd = [13.958, -3.445, -1.68, 0.63], [0.42, 0.69, 0.22, 0.18]
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
# The LCDM curve is normalised so that the integral of M phi(M) dM equals
# Omega_M rho_crit.  That Omega_M must be the one Driver used (0.3147), NOT
# recovery.py's OMEGA_M = 0.25, which is the MOCK's cosmology and is used here
# only for comoving volumes.  Using 0.25 scaled the curve to hold 25% of the
# critical density instead of 31.47%, putting it log10(0.3147/0.25) = 0.100 dex
# low -- which is why the LCDM line sat ~0.12 dex BELOW Driver's GSR fit here
# while the two lie on top of each other in his own figures.
LCDM_OMEGA_M = 0.3147


def lcdm_curve():
    """Murray+21 LCDM MRP, Omega_M-normalised and z=0.1-shifted, exactly as
    Driver+22 allhmf.r. Returns (x, log10 phi) for the curve plus the (M*, logphi*,
    alpha, beta) reference point for the corner plot.

    Normalised with LCDM_OMEGA_M (Driver's 0.3147), not the module OMEGA_M."""
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
        / (LCDM_OMEGA_M * rhocrit)
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
    driver=(13.958, -3.445, -1.68, 0.63),
    data_dir="../data",
    fname="hmf_publication.pdf",
    title="HMF",
    show_corrected=None,
    corrected_band=False,
    nessie_bias=False,
    mmax_data=None,
):
    """Driver-style HMF: our MCMC band (highlighted), the LCDM curve, the Driver+22
    MRP curve, our own survey binned points (grey, 'not fitted'), and external
    comparison data (REFLEX/2PIGG/ELMO) from data_dir."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if show_corrected is None:
        show_corrected = SHOW_ALPHA_CORRECTION
    med = np.median(flat, axis=0)
    mgrid = np.linspace(12.5, 16, 400)
    fig, ax = plt.subplots(figsize=(8, 6))

    # our posterior: a 16-84 credible band, not 400 translucent curves.  The
    # spaghetti version read as noise and hid the comparison lines.
    idx = np.random.default_rng(0).choice(
        flat.shape[0], size=min(400, flat.shape[0]), replace=False
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        _curves = np.array([np.log10(mrp_phi(mgrid, *flat[k])) for k in idx])
    _lo, _hi = np.nanpercentile(_curves, [16, 84], axis=0)
    ax.fill_between(
        mgrid, _lo, _hi, color=(100 / 255, 149 / 255, 237 / 255), alpha=0.30, lw=0
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

    # Closed-loop bias-corrected band: an estimate of the HALO mass function.
    # It will NOT track the comparison points, and should not -- those are
    # densities of DETECTED objects, which differ from the halo density by the
    # completeness (C ~ 0.28 at mlim).
    if nessie_bias:
        fb = apply_nessie_bias(flat)
        for k in idx:
            ax.plot(
                mgrid, np.log10(mrp_phi(mgrid, *fb[k])), color="darkgreen", alpha=0.02
            )
        ax.plot(
            mgrid,
            np.log10(mrp_phi(mgrid, *np.median(fb, axis=0))),
            color="darkgreen",
            lw=2,
            ls="-",
            label="bias-corrected (halo MF)",
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
    g5 = driver_gama5(match_A=True)
    ax.plot(
        mgrid,
        np.log10(mrp_phi(mgrid, *g5)),
        color="darkred",
        lw=1.6,
        ls=(0, (6, 2)),
        label=f"Driver+22 GAMA only (A={A_SCALE:g})",
    )

    # our own survey binned points (grey, not fitted) -- see SHOW_OWN_POINTS
    if my_surveys and SHOW_OWN_POINTS:
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

    # Beyond the most massive group in the sample the curve extrapolates a cutoff
    # that nothing constrains: no groups, the completeness table held flat, and
    # beta only weakly determined. Unmarked, that region reads as a measurement.
    if mmax_data is not None:
        ax.axvspan(float(mmax_data), 16.0, color="0.85", alpha=0.45, zorder=0)
        ax.axvline(float(mmax_data), color="0.55", lw=1, zorder=1)
        ax.annotate(
            "no groups above here",
            xy=(float(mmax_data) + 0.06, -7.7),
            fontsize=7.5,
            color="0.35",
            ha="left",
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
    driver=(13.958, -3.445, -1.68, 0.63),
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
                ax.axvline(driver_gama5(match_A=True)[i], color="darkred", ls=":")
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
                g5 = driver_gama5(match_A=True)
                ax.plot(g5[j], g5[i], marker="P", color="darkred", ms=7, ls="")
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


def plot_ppc(flat, data, x_fit, fname="ppc.pdf", nbin=14, title="Posterior predictive"):
    """Does the model fit the GAMA data?

    The HMF figure shows phi(m), a HALO mass function, while the catalogue holds
    DETECTED groups -- the two differ by the completeness, so the fit is not
    supposed to pass through the observed counts and the figure cannot tell you
    whether it fits.

    What the model actually predicts is the number of catalogue entries per mass
    bin: dN/dm = sum_j V_j phi(m) C(m, z_j). That is compared here with the
    observed histogram. If the points sit in the band, the model fits."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xg = np.asarray(data["xg"], float)
    Csh = np.asarray(data["Csh"], float)
    V_sh = np.asarray(data["V_sh"], float)

    lo = float(np.percentile(x_fit, 0.5))
    hi = float(np.percentile(x_fit, 99.5))
    edges = np.linspace(lo, hi, nbin + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    obs, _ = np.histogram(x_fit, bins=edges)

    # expected detections per bin, for a sample of the posterior
    idx = np.random.default_rng(0).choice(
        flat.shape[0], size=min(300, flat.shape[0]), replace=False
    )
    pred = np.empty((idx.size, nbin))
    for n, k in enumerate(idx):
        pg = mrp_phi(xg, *flat[k])
        dens = np.sum(V_sh[:, None] * Csh * pg[None, :], axis=0)  # dN/dm on the grid
        for b in range(nbin):
            m = (xg >= edges[b]) & (xg < edges[b + 1])
            pred[n, b] = np.trapezoid(dens[m], xg[m]) if m.sum() > 1 else np.nan

    med = np.nanmedian(pred, axis=0)
    q16, q84 = np.nanpercentile(pred, [16, 84], axis=0)

    fig, (a1, a2) = plt.subplots(
        2,
        1,
        figsize=(6.4, 6.0),
        sharex=True,
        gridspec_kw=dict(height_ratios=[2.4, 1], hspace=0.06),
    )
    a1.fill_between(
        cen,
        q16,
        q84,
        color="steelblue",
        alpha=0.35,
        label="model 68% (predicted detections)",
    )
    a1.plot(cen, med, color="steelblue", lw=1.8, label="model median")
    a1.errorbar(
        cen,
        obs,
        yerr=np.sqrt(np.maximum(obs, 1)),
        fmt="o",
        ms=4.5,
        color="k",
        capsize=2.5,
        lw=1,
        label="GAMA groups (observed)",
    )
    a1.set_yscale("log")
    a1.set_ylabel("groups per bin")
    a1.set_title(title)
    a1.legend(fontsize=8)

    with np.errstate(divide="ignore", invalid="ignore"):
        pull = (obs - med) / np.sqrt(np.maximum(med, 1))
    a2.axhline(0, color="0.6", lw=1)
    a2.axhspan(-1, 1, color="0.85", alpha=0.6)
    a2.plot(cen, pull, "o-", color="crimson", ms=4.5, lw=1.2)
    a2.set(
        xlabel=r"$\log_{10}(M_{\rm halo}/M_\odot)$",
        ylabel=r"(obs$-$model)/$\sqrt{\rm model}$",
    )
    lim = max(3.0, float(np.nanmax(np.abs(pull))) * 1.15)
    a2.set_ylim(-lim, lim)

    chi2 = float(np.nansum(pull**2))
    a2.annotate(
        f"$\\chi^2$/bin = {chi2 / nbin:.2f}",
        xy=(0.02, 0.08),
        xycoords="axes fraction",
        fontsize=9,
    )
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    print(
        f"  saved {fname}   (chi2/bin = {chi2 / nbin:.2f}, "
        f"total observed {obs.sum()}, predicted {med.sum():.0f})"
    )
    return fname


def emit_publication(
    flat, my_surveys, tag, title="HMF", nessie_bias=False, mmax_data=None
):
    """Always emit the corner + publication HMF plots for a fit."""
    try:
        plot_corner(flat, fname=f"corner_{tag}.pdf")
    except Exception as e:
        print(f"  [corner_{tag} failed: {e}]")
    try:
        plot_publication(
            flat,
            my_surveys=my_surveys,
            fname=f"hmf_{tag}.pdf",
            title=title,
            nessie_bias=nessie_bias,
            mmax_data=mmax_data,
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


def run_mock_nessie(path="nessie_mock_groups.npz", model_kind="marg_tab"):
    """Closed-loop validation of the tabulated-completeness fit.

    The MRP in TRUE was injected into the Shark halos by abundance matching;
    Nessie then recovered the groups saved by measure_completeness_nessie.py.
    Fitting those groups with the C table measured from the same run must
    return the injected parameters. This is the only test of marg_tab against
    a known answer -- the ordinary mock pipeline uses the membership proxy to
    define detections, which is inconsistent with a Nessie-measured C.

    Expect log M* low by ~0.155 dex: MassA (A=10) under-estimates the true halo
    mass by that much, and the model has no term for it."""
    t = np.load(path)
    log_mass, z = t["log_mass"], t["z"]
    mult = t["multiplicity"]
    area = float(t["area_deg2"])
    sky_frac = area * (np.pi / 180) ** 2 / (4 * np.pi)
    sigma = sigma_from_nfof(mult)

    ok = np.isfinite(log_mass) & np.isfinite(z) & (z > ZMIN) & (z < ZLIMIT)
    log_mass, z, sigma = log_mass[ok], z[ok], sigma[ok]
    print(
        f"  Nessie mock groups: {log_mass.size} over {area:.1f} deg^2 "
        f"({log_mass.size / area:.2f} per deg^2)"
    )
    print(
        f"  mass {log_mass.min():.2f}..{log_mass.max():.2f} "
        f"(med {np.median(log_mass):.2f})"
    )

    _mz = COMP_MODE == "mz" and model_kind.startswith("marg_tab")
    if _mz:
        print("mlim(z) NOT used: selection carried by C(m,z) [--comp-mode mz]")
        mlim_func = mlim_sh = turn_pts = None
        z_mids, V_sh = shell_volumes(sky_frac)
    else:
        mlim_func, _, kind, turn_pts = turnover_mlim(z, log_mass)
        print(f"  mlim(z) [{kind}]: {mlim_func(ZMIN):.2f} -> {mlim_func(ZLIMIT):.2f}")
        z_mids, V_sh = shell_volumes(sky_frac)
        mlim_sh = mlim_func(z_mids)

    data, keep = prep_tab(
        z,
        log_mass,
        sigma,
        mlim_func,
        z_mids,
        mlim_sh,
        V_sh,
        fit_scale=(model_kind == "marg_tab_serr"),
    )
    # In mz mode the fit lives in the MASS_BIAS-shifted true-mass coordinate, so
    # the figures must use the same one or they sit MASS_BIAS dex off the model.
    x_plot = np.asarray(data["x_obs"], float) if _mz else log_mass[keep]
    print(f"\nFitting [{model_kind}] on the Nessie mock catalogue ...")
    map_par, flat = run_stan(model_kind, data)
    res = summarise(flat)
    try:
        check_lambda(data, np.median(flat, axis=0))
    except Exception as e:
        print(f"  [Lambda check failed: {e}]")
    # PPC on the MOCK, so the real-data chi2/bin has a baseline. Without it a
    # value like 15.15 on GAMA is uninterpretable: it could mean the data reject
    # the model, or simply that this statistic is not normalised to ~1.
    try:
        plot_ppc(flat, data, x_plot, fname=f"ppc_nessiemock_{model_kind}.pdf",
                 title="Posterior predictive (Nessie mock)")
    except Exception as e:
        print(f"  [ppc failed: {e}]")
    print("\n  Truth here is the INJECTED MRP, so 'bias(sd)' is genuine recovery.")
    print("  A log M* deficit of ~0.155 dex is expected (the A=10 dynamical mass")
    print("  under-estimates the true halo mass by that much on this mock).")
    emit_publication(
        flat,
        {"Nessie mock": dict(x_fit=x_plot, Vsurvey=float(V_sh.sum()))},
        tag=f"nessiemock_{model_kind}",
        title="Nessie mock HMF",
        mmax_data=float(np.percentile(x_plot, 99.5)),
    )
    return res


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
        "--alpha-correction",
        action="store_true",
        help="report/plot the Shark-calibrated alpha bias correction",
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
        "--mock-nessie",
        action="store_true",
        help="closed-loop validation: fit the Nessie mock catalogue",
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
        "--gama-regions",
        nargs="+",
        default=None,
        help="GAMARegion values to keep, e.g. G09 G12 G15 G23 (default: all)",
    )
    ap.add_argument(
        "--gama-dec-cut",
        type=float,
        default=None,
        help="legacy Dec cut; Driver used -3.5, which EXCLUDES G23",
    )
    ap.add_argument(
        "--pin-al",
        type=float,
        default=None,
        help="value to pin alpha at for marg_tab_pinal; default -1.68 "
        "(Driver GSR). GAMA-only -1.27, LCDM -1.865",
    )
    ap.add_argument(
        "--pin-be",
        type=float,
        default=None,
        help="value to pin beta at for marg_tab_pinbe; default 0.63 "
        "(Driver GSR). GAMA-only 0.47, LCDM 0.71",
    )
    ap.add_argument(
        "--pin-ms",
        type=float,
        default=None,
        help="value to pin log10 M* at (h=1 units) for marg_tab_pinms; "
        "default 13.958 = Driver+22 GSR",
    )
    ap.add_argument(
        "--mlim-form",
        choices=["linear", "quad"],
        default=None,
        help="force the mlim(z) form; the mock used linear and the "
        "AIC choice is unstable. Irrelevant under --comp-mode mz, where "
        "mlim(z) is not formed at all",
    )
    ap.add_argument(
        "--comp-mode",
        choices=["delta", "mz"],
        default="delta",
        help="completeness keying for the marg_tab models. 'delta' is the "
        "incumbent C(m - mlim(z), z), which needs turnover_mlim; 'mz' uses "
        "C(m, z) on absolute true mass and drops mlim(z) entirely",
    )
    ap.add_argument(
        "--comp-def",
        choices=["entries", "repr"],
        default="entries",
        help="which completeness the mz table supplies: 'entries' (expected "
        "number of catalogue entries per halo, the Poisson intensity, may "
        "exceed 1) or 'repr' (bounded recovered fraction, systematic check)",
    )
    ap.add_argument("--comp-table-mz", default=NESSIE_TABLE_MZ)
    ap.add_argument(
        "--mass-bias",
        type=float,
        default=None,
        help="dex offset between the observed mass column and the TRUE mass the "
        "mz table is keyed on; x_obs is shifted by -this. Measured from the "
        "mock by measure_completeness_nessie.py",
    )
    ap.add_argument(
        "--mass-col",
        default=None,
        help="use this catalogue mass column (e.g. MassA) instead of "
        "rebuilding from VelDisp/Rad50; MassA matches the mock",
    )
    ap.add_argument(
        "--veldisp-err",
        action="store_true",
        help="per-group mass errors from the catalogue VelDispErr "
        "instead of the multiplicity lookup table",
    )
    ap.add_argument(
        "--gama-model",
        choices=[
            "marg",
            "marg_comp",
            "marg_comp_serr",
            "marg_tab",
            "marg_tab_serr",
            "marg_tab_pinms",
            "marg_tab_pinbe",
            "marg_tab_pinal",
            "gama",
            "simple",
        ],
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
    SHOW_ALPHA_CORRECTION = args.alpha_correction
    MLIM_FORM = args.mlim_form
    # These are read as globals at call time by run_real_gama / prep_tab, so
    # they have to be assigned in this block -- anywhere else and the default wins.
    COMP_MODE = args.comp_mode
    COMP_DEF = args.comp_def
    NESSIE_TABLE_MZ = args.comp_table_mz
    if args.mass_bias is not None:
        MASS_BIAS = float(args.mass_bias)
    if COMP_MODE == "mz":
        print(
            f"[comp-mode mz] C(m,z) from {NESSIE_TABLE_MZ} [{COMP_DEF}], "
            f"MASS_BIAS = {MASS_BIAS:+.3f} dex, mlim(z) not used"
        )
    if args.pin_ms is not None:
        PIN_MS = args.pin_ms
    if args.pin_be is not None:
        PIN_BE = args.pin_be
    if args.pin_al is not None:
        PIN_AL = args.pin_al
    USE_DRIVER_PRIOR = args.driver_prior
    DRIVER_PRIOR_INFLATE = args.driver_prior_inflate
    if args.mock_nessie:
        run_mock_nessie(model_kind=args.gama_model)
    elif args.calibrate_alpha:
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
            regions=args.gama_regions,
            use_veldisp_err=args.veldisp_err,
            dec_cut=args.gama_dec_cut,
            mass_col=args.mass_col,
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
