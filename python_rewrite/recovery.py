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
TRUE = dict(ms=13.51, lp=-3.19, al=-1.27, be=0.47)
PARAMS = ["ms", "lp", "al", "be"]

# Cosmology
H0, OMEGA_M = 67.37, 0.3147

# Survey / selection
ZMIN, ZLIMIT = 0.01, 0.25
MULTI = 5  # min members for a detection
MAG_LIMIT = 19.8  # r-band limit
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
    """Detection = >= MULTI galaxies brighter than MAG_LIMIT (log Mstar > 8).

    NOTE: log_mstar_total is log10(Mstar/Msun) (~8-12), so the stellar-mass
    floor is `> 8`, i.e. log10(1e8) -- NOT `> 1e8`."""
    gal = galaxies[
        (galaxies["id_fof"] != -1)
        & (galaxies["log_mstar_total"] > 8)
        & (galaxies["total_ap_dust_r_SDSS"] < MAG_LIMIT)
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
        m_obs = det["log_mass_am"].values + rng.normal(0, sigma)
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


def turnover_mlim(z_obs, m_obs, nbin_z=30, hist_bw=0.3, min_in_bin=20):
    z_edges = np.linspace(ZMIN, ZLIMIT, nbin_z + 1)
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
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.3, 0.3);   // tightened: literature slope ~ -1.3, regularises the degeneracy
  be ~ normal(0.47, 0.05);   // TIGHT: cutoff pinned near Driver+22 (0.47); reported as prior-set

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
  ms ~ normal(14.0, 1.5);
  lp ~ normal(-4.0, 2.0);
  al ~ normal(-1.3, 0.3);   // tightened: literature slope ~ -1.3, regularises the degeneracy
  be ~ normal(0.47, 0.05);   // TIGHT: cutoff pinned near Driver+22 (0.47); reported as prior-set

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

_STAN = {"simple": SIMPLE_CODE, "marg": MARG_CODE}
_MODELS = {}

# Compile OUTSIDE any iCloud-synced tree (e.g. ~/Desktop, ~/Documents).
# A compiled Stan binary placed under a synced folder can be evicted mid-run,
# which shows up as "No such file or directory: .../mrp_marg" and retcode -1.
# Override with the HMF_STAN_DIR environment variable if you like.
STAN_BUILD_DIR = os.environ.get(
    "HMF_STAN_DIR", os.path.join(os.path.expanduser("~"), ".cache", "hmf_mrp_stan")
)
os.makedirs(STAN_BUILD_DIR, exist_ok=True)


def get_model(kind):
    """Compile a Stan model once, in a non-synced build dir (cmdstanpy caches)."""
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


def build_data(kind, x_fit, sig_fit, mlim_sh, V_sh, sig_sh=None):
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
    return dict(median=med, sd=sd, q16=q16, q84=q84)


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
    prior_mu = [14.0, -4.0, -1.3, 0.47]  # must match the Stan model priors
    prior_sd = [1.5, 2.0, 0.3, 0.05]
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
    prior_sd = {"ms": 1.5, "lp": 2.0, "al": 0.3, "be": 0.05}
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
def run_real_pipeline(model_kind="marg"):
    print("Reading catalogues ...")
    groups, galaxies = load_catalogues(DATA_DIR)
    sky_frac = (
        sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )

    print("Abundance matching ...")
    gv, Vsurvey = abundance_match(groups, sky_frac)

    print("GAMA selection ...")
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
        choices=["simple", "marg"],
        default="marg",
        help="which likelihood: 'simple' (baseline) or 'marg' "
        "(marginalised + boundary). Default marg.",
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
    args = ap.parse_args()
    if args.selftest:
        run_selftest(model_kind=args.model)
    elif args.coverage:
        run_coverage(model_kind=args.model, n_real=args.nreal)
    else:
        run_real_pipeline(model_kind=args.model)
