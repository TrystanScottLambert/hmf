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
NSH = 20  # redshift shells for the Poisson normalisation


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


def add_mass_errors(gv, rng):
    """Multiplicity-dependent Gaussian errors on log-mass (GAMA-like)."""
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
    det = gv[gv["detected"]].copy()
    sigma = np.interp(det["n_gama"].values, xx, yy)  # np.interp clamps at ends
    sigma = np.maximum(sigma, 0.10)
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
# 4. Stan model -- document-1 `stan_simple`, run via cmdstanpy
# ----------------------------------------------------------
STAN_CODE = r"""
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
"""

STAN_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mrp_simple.stan")
_MODEL = None


def get_model():
    """Compile the Stan model once (cmdstanpy caches the binary)."""
    global _MODEL
    if _MODEL is None:
        from cmdstanpy import CmdStanModel

        with open(STAN_PATH, "w") as f:
            f.write(STAN_CODE)
        _MODEL = CmdStanModel(stan_file=STAN_PATH)
    return _MODEL


def run_stan(
    x_fit,
    mlim_sh,
    V_sh,
    chains=4,
    warmup=2000,
    sampling=2000,
    adapt_delta=0.95,
    max_treedepth=12,
    seed=42,
    show_progress=True,
):
    """MAP (optimize) + posterior (sample) for the simple model.
    Returns (map_par, flat) where flat is an (Ndraws, 4) array in PARAMS order."""
    model = get_model()
    data = dict(
        N=int(np.asarray(x_fit).size),
        x_obs=np.asarray(x_fit, float),
        Nsh=int(np.asarray(V_sh).size),
        V_sh=np.asarray(V_sh, float),
        mlim_sh=np.asarray(mlim_sh, float),
        xhi=float(XHI),
        Ng=int(NG),
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
    print("  MAP:", dict(zip(PARAMS, np.round(map_par, 3))))

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
    print(f"  Rhat={rhat:.3f}  min ESS={ess:.0f}  divergences={ndiv}")
    return map_par, flat


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
    for k, (i, j) in enumerate([(1, 0), (1, 1), (1, 2), (0, 2)]):
        a = ax[i, j]
        a.hist(flat[:, k], bins=40, color="steelblue", density=True)
        a.axvline(tv[k], color="blue", lw=2, ls="--", label="truth")
        a.axvline(med[k], color="red", lw=2, label="median")
        a.set(title=labels[k])
        a.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(fname)
    print(f"  saved {fname}")
    return fname


# ----------------------------------------------------------
# 6. Drivers
# ----------------------------------------------------------
def run_real_pipeline():
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
    print(
        f"  N above mlim: {x_fit.size} / {m_obs.size} ({100 * x_fit.size / m_obs.size:.1f}%)"
    )

    print("\nFitting (cmdstanpy) ...")
    map_par, flat = run_stan(x_fit, mlim_sh, V_sh)
    res = summarise(flat)
    plot_recovery(flat, z_obs, m_obs, x_fit, mlim_func, Vsurvey, turn_pts=turn_pts)
    return res


def run_selftest():
    """Generate a sample DIRECTLY from the MRP under the model's own
    generative assumptions and confirm the Stan model recovers it.
    Isolates the statistical core from the data I/O."""
    print("=" * 60)
    print("  SELF-TEST: synthetic recovery (no data files)")
    print("=" * 60)
    rng = np.random.default_rng(7)

    sky_frac = 0.004
    z_mids, V_sh = shell_volumes(sky_frac)
    mlim_func = lambda z: 12.6 + 1.5 * z
    mlim_sh = mlim_func(z_mids)

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
    print(f"  generated {x_fit.size} detected groups")

    map_par, flat = run_stan(
        x_fit, mlim_sh, V_sh, warmup=1000, sampling=1000, show_progress=False
    )
    res = summarise(flat)
    tv = np.array([TRUE[p] for p in PARAMS])
    bias = np.abs((res["median"] - tv) / res["sd"])
    print(
        f"\n  max |bias| = {bias.max():.2f} sd",
        " -> PASS" if bias.max() < 3 else " -> CHECK",
    )
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--selftest",
        action="store_true",
        help="run synthetic recovery without needing the parquet files",
    )
    args = ap.parse_args()
    if args.selftest:
        run_selftest()
    else:
        run_real_pipeline()
