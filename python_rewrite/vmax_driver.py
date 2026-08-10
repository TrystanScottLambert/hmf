"""
============================================================
vmax_driver.py -- Driver+22's method, applied to the new data
============================================================

Reproduces the Driver+22 analysis (binned 1/Vmax, first-order Eddington
correction, chi^2 fit of the MRP with a Poisson penalty, Monte-Carlo errors
including cosmic variance) on the new GAMA DMU and the Nessie SDSS catalogue.

This is the "same method, better data" cell: everything follows Driver, only the
input catalogues change. Comparing it with the forward-modelled fit isolates the
method; comparing it with Driver's published values isolates the data.

WHAT IS TAKEN VERBATIM FROM allhmf.r
------------------------------------
  * chi^2 in log10 space with sigma_log = f / ln(10) for fractional error f
  * the Poisson penalty  2 * V * sum(phi(m > m_max)) * logbin  per survey,
    which suppresses fits predicting groups above the highest occupied bin
  * logbin = 0.2
  * cosmic variance  cosvar(V, N) = (219.7 - 52.4 logV + 3.21 (logV)^2)/sqrt(N)/100
    evaluated per field (GAMA: V/3 with N=3), 5 per cent for REFLEX
  * the Monte-Carlo loop: perturb each point by cosmic variance in log space,
    then by its fractional error in linear space, refit, repeat
  * REFLEX fractional errors 1/sqrt(20), and 1/sqrt(3) at the two endpoints

WHAT IS RECONSTRUCTED
---------------------
Driver's binning script is not available, so the 1/Vmax estimator and the
Eddington correction are implemented from standard definitions:

  * Vmax: for a group of mass m, the comoving volume out to the redshift where
    mlim(z) = m, capped at the survey limit. phi = sum(1/Vmax) / binwidth.
  * Eddington: iterative first-order correction. Given the current MRP, smear it
    by the per-object mass errors, take the ratio true/smeared per bin, and
    divide the observed phi by it. Two iterations.

These follow the description in Driver+22 but are not his literal code, so
small differences from his published points are expected.

Run:
    python vmax_driver.py --gama-only
    python vmax_driver.py --combined
============================================================
"""

import argparse
import numpy as np
from scipy import optimize
from scipy.stats import norm

import recovery as R

LOGBIN = 0.2
MASS_EDGES = np.arange(10.0, 16.0 + 1e-9, LOGBIN)
MASS_MIDS = 0.5 * (MASS_EDGES[1:] + MASS_EDGES[:-1])
N_MC = 1001  # allhmf.r uses 10001; 1001 is ample and much faster
REFLEX_VOLUME = 1.3e7  # volumeREFLEXII in allhmf.r, Mpc^3


def cosvar(V, N):
    """Cosmic variance, verbatim from allhmf.r."""
    lv = np.log10(V)
    return ((219.7 - 52.4 * lv + 3.21 * lv**2) / np.sqrt(N)) / 100.0


# ----------------------------------------------------------------------
def vmax_hmf(log_mass, z, sigma, mlim_func, sky_frac, zmin, zmax, label="", mmax=None):
    """Binned 1/Vmax mass function.

    For each group the maximum redshift at which it would still lie above
    mlim(z) is found by inverting mlim, and Vmax is the comoving volume out to
    it (capped at the survey limits). Returns (mass, log10 phi, frac_err, N)
    for the occupied bins."""
    # 1/Vmax is defined only for objects that are actually detectable, i.e.
    # above the limit at their own redshift. Including groups below mlim gives
    # them Vmax -> 0 and hence enormous weight, which drags alpha to nonsense.
    if mmax is not None:
        n_hi = int((log_mass > mmax).sum())
        if n_hi:
            print(f"  [{label}] capping at logM < {mmax}: dropping {n_hi} groups")
        keep_hi = log_mass <= mmax
        log_mass, z = log_mass[keep_hi], z[keep_hi]
        if sigma is not None and np.ndim(sigma):
            sigma = np.asarray(sigma)[keep_hi]

    above = log_mass > mlim_func(z)
    n_drop = int((~above).sum())
    log_mass, z = log_mass[above], z[above]
    if sigma is not None and np.ndim(sigma) and len(sigma) == above.size:
        sigma = np.asarray(sigma)[above]

    zg = np.linspace(zmin, zmax, 2001)
    mg = mlim_func(zg)
    # mlim rises with z, so a group drops out once mlim exceeds its mass
    zmax_i = np.interp(log_mass, mg, zg, left=zmin, right=zmax)
    zmax_i = np.clip(zmax_i, zmin, zmax)

    d_hi = R.comoving_distance(zmax_i)
    d_lo = R.comoving_distance(np.full_like(zmax_i, zmin))
    vmax = (4 / 3) * np.pi * (d_hi**3 - d_lo**3) * sky_frac

    # a group barely above the limit has a tiny Vmax and a correspondingly huge
    # weight; drop the worst offenders rather than let one object set a bin
    vfull = (
        (4 / 3)
        * np.pi
        * (R.comoving_distance(np.array([zmax]))[0] ** 3 - d_lo[0] ** 3)
        * sky_frac
    )
    usable = vmax > 1e-4 * vfull
    n_tiny = int((~usable).sum())
    log_mass, vmax = log_mass[usable], vmax[usable]

    idx = np.digitize(log_mass, MASS_EDGES) - 1
    ok = (idx >= 0) & (idx < MASS_MIDS.size)
    phi = np.zeros(MASS_MIDS.size)
    cnt = np.zeros(MASS_MIDS.size, int)
    for i, w in zip(idx[ok], 1.0 / vmax[ok]):
        phi[i] += w
        cnt[i] += 1
    phi /= LOGBIN

    keep = cnt >= 3
    frac = np.where(cnt > 0, 1.0 / np.sqrt(np.maximum(cnt, 1)), 1.0)
    print(
        f"  [{label}] 1/Vmax: {int(cnt.sum())} groups above mlim in "
        f"{int(keep.sum())} bins, logM {MASS_MIDS[keep].min():.2f}-"
        f"{MASS_MIDS[keep].max():.2f}"
    )
    print(f"           dropped {n_drop} below mlim, {n_tiny} with Vmax < 1e-4 Vsurvey")
    with np.errstate(divide="ignore"):
        lp = np.log10(phi)
    return MASS_MIDS[keep], lp[keep], frac[keep], cnt[keep]


def eddington_factor(mass, par, sigma_typ):
    """First-order Eddington correction: ratio of the true MRP to the same MRP
    smeared by the mass errors, evaluated per bin. Multiplying the observed phi
    by this removes the up-scatter excess."""
    fine = np.arange(mass.min() - 6 * sigma_typ, mass.max() + 6 * sigma_typ, 0.02)
    true = R.mrp_phi(fine, *par)
    g = np.arange(-5 * sigma_typ, 5 * sigma_typ + 1e-9, 0.02)
    k = np.exp(-0.5 * (g / sigma_typ) ** 2)
    k /= k.sum()
    smear = np.convolve(true, k, mode="same")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(smear > 0, true / smear, 1.0)
    return np.interp(mass, fine, ratio)


# ----------------------------------------------------------------------
def massfn_chi2(par, sets):
    """chi^2 + Poisson penalty, exactly as allhmf.r's massfn().

    `sets` is a list of (mass, logphi, frac_err, volume) -- one per survey."""
    ms, lp, al, be = par
    if not (0.1 < be < 2.0):
        return 1e12
    tot = 0.0
    for mass, y, f, vol in sets:
        model = np.log10(np.maximum(R.mrp_phi(mass, ms, lp, al, be), 1e-300))
        sig = np.maximum(f, 1e-3) / np.log(10)
        tot += float(np.sum(((y - model) / sig) ** 2))
        # penalty: 2 x expected counts in the ten bins above the highest point
        above = mass.max() + np.arange(1, 11) * LOGBIN
        tot += 2.0 * vol * float(np.sum(R.mrp_phi(above, ms, lp, al, be))) * LOGBIN
    return tot


BOUNDS = [(12.0, 16.0), (-8.0, -1.0), (-2.6, -0.5), (0.15, 1.6)]


def fit_mrp(sets, p0=(13.9, -3.4, -1.4, 0.6)):
    """Bounded so a bad Monte-Carlo realisation cannot run the fit off to
    alpha ~ -5 or phi* > 0, which overflows the MRP."""

    def obj(par):
        for v, (lo, hi) in zip(par, BOUNDS):
            if not (lo <= v <= hi):
                return 1e12
        return massfn_chi2(par, sets)

    r = optimize.minimize(
        obj,
        np.asarray(p0, float),
        method="Nelder-Mead",
        options=dict(xatol=1e-5, fatol=1e-4, maxiter=20000),
    )
    return r.x


def mc_errors(sets, cosvars, n_mc=N_MC, seed=42):
    """Driver's Monte-Carlo: perturb each point by cosmic variance (log space)
    then by its fractional error (linear space), refit, and take the spread."""
    rng = np.random.default_rng(seed)
    out = np.full((n_mc, 4), np.nan)
    for i in range(n_mc):
        pert = []
        for (mass, y, f, vol), cv in zip(sets, cosvars):
            yy = y + np.abs(y) * rng.normal(0.0, cv, size=y.size)
            lin = 10**yy
            lin = lin + lin * rng.normal(0.0, np.clip(f, 0, 0.999), size=y.size)
            good = lin > 0
            pert.append((mass[good], np.log10(lin[good]), f[good], vol))
        try:
            out[i] = fit_mrp(pert)
        except Exception:
            pass
        if (i + 1) % 200 == 0:
            print(f"    MC {i + 1}/{n_mc}")
    return out[np.isfinite(out).all(axis=1)]


# ----------------------------------------------------------------------
def load_gama(args):
    lm, sig, z, nf = R.load_real_gama(args.gama_fits, mass_col=args.mass_col)
    sky_frac = args.gama_area * (np.pi / 180) ** 2 / (4 * np.pi)
    mlim_func, _, kind, _ = R.turnover_mlim(z, lm, form="linear")
    print(
        f"  [GAMA] mlim(z) [{kind}]: {mlim_func(R.ZMIN):.2f} -> {mlim_func(R.ZLIMIT):.2f}"
    )
    vol = float(R.survey_volume(sky_frac))
    m, y, f, n = vmax_hmf(lm, z, sig, mlim_func, sky_frac, R.ZMIN, R.ZLIMIT, "GAMA")
    return dict(
        m=m,
        y=y,
        f=f,
        n=n,
        vol=vol,
        sigma=float(np.median(sig)),
        nfield=4,
        label="GAMA (this work)",
    )


def load_sdss(args):
    lm, sig, z, mult = R.load_sdss_groups(
        args.sdss_parquet, zmin=args.sdss_zmin, zmax=args.sdss_zmax
    )
    sky_frac = args.sdss_frac
    mlim_func, _, kind, _ = R.turnover_mlim(
        z, lm, zmin=args.sdss_zmin, zmax=args.sdss_zmax, form="linear"
    )
    print(
        f"  [SDSS] mlim(z) [{kind}]: {mlim_func(args.sdss_zmin):.2f} -> "
        f"{mlim_func(args.sdss_zmax):.2f}"
    )
    vol = float(R.survey_volume(sky_frac, zmin=args.sdss_zmin, zmax=args.sdss_zmax))
    m, y, f, n = vmax_hmf(
        lm,
        z,
        sig,
        mlim_func,
        sky_frac,
        args.sdss_zmin,
        args.sdss_zmax,
        "SDSS",
        mmax=args.sdss_mmax,
    )
    return dict(
        m=m,
        y=y,
        f=f,
        n=n,
        vol=vol,
        sigma=float(np.median(sig)),
        nfield=1,
        label="SDSS (Nessie, this work)",
    )


def load_reflex(data_dir="../data"):
    comp = R._load_comparison(data_dir)
    rf = comp.get("REFLEX II (Böhringer+17)")
    if rf is None:
        return None
    x, y = np.asarray(rf["x"], float), np.asarray(rf["y"], float)
    f = np.full(x.size, 1 / np.sqrt(20))
    f[0] = f[-1] = 1 / np.sqrt(3)  # allhmf.r
    print(f"  [REFLEX] {x.size} binned points, logM {x.min():.2f}-{x.max():.2f}")
    return dict(
        m=x,
        y=y,
        f=f,
        n=np.full(x.size, 20),
        vol=REFLEX_VOLUME,
        sigma=0.0,
        nfield=1,
        label="REFLEX II (Böhringer+17)",
    )


# ----------------------------------------------------------------------
def run(sets_raw, tag, title, args, eddington=True):
    """Eddington-correct, fit, Monte-Carlo, report and plot."""
    sets = [(d["m"], d["y"].copy(), d["f"], d["vol"]) for d in sets_raw]
    par = fit_mrp(sets)

    if eddington:
        for _ in range(2):
            sets = []
            for d in sets_raw:
                y = d["y"].copy()
                if d["sigma"] > 0.02:
                    y = y + np.log10(eddington_factor(d["m"], par, d["sigma"]))
                sets.append((d["m"], y, d["f"], d["vol"]))
            par = fit_mrp(sets)
        for d, (mm, yy, ff, vv) in zip(sets_raw, sets):
            d["y_corr"] = yy
        print(
            f"  Eddington correction applied (2 iterations, "
            f"sigma = {[round(d['sigma'], 3) for d in sets_raw]})"
        )
        for d in sets_raw:
            if "y_corr" in d:
                sh = float(np.median(d["y_corr"] - d["y"]))
                print(f"    {d['label']}: median shift {sh:+.3f} dex")

    cosvars = [
        cosvar(d["vol"] / max(d["nfield"], 1), d["nfield"])
        if d["label"].startswith(("GAMA", "SDSS"))
        else 0.05
        for d in sets_raw
    ]
    print(f"  cosmic variance per survey: {[round(c, 4) for c in cosvars]}")
    print(f"  running {N_MC} Monte-Carlo refits ...")
    chains = mc_errors(sets, cosvars)

    med = np.median(chains, axis=0)
    lo, hi = np.percentile(chains, [16, 84], axis=0)
    names = ["log M*", "log phi*", "alpha", "beta"]
    print(f"\n  === Driver+22 method on {title} ===")
    print(f"  {'param':>9} {'best fit':>10} {'MC median':>11} {'16-84%':>18}")
    for i, nm in enumerate(names):
        blo, bhi = BOUNDS[i]
        flag = (
            "  <-- AT BOUND"
            if (abs(par[i] - blo) < 1e-3 or abs(par[i] - bhi) < 1e-3)
            else ""
        )
        print(
            f"  {nm:>9} {par[i]:10.3f} {med[i]:11.3f}   "
            f"[{lo[i]:7.3f}, {hi[i]:7.3f}]{flag}"
        )
    at_bound = [
        names[i]
        for i in range(4)
        if abs(par[i] - BOUNDS[i][0]) < 1e-3 or abs(par[i] - BOUNDS[i][1]) < 1e-3
    ]
    if at_bound:
        print(
            f"  !! {', '.join(at_bound)} railed against a bound -- the data are"
            f" pulling outside the physical range, so this fit is not usable."
        )

    dms, dlp, _, _ = R.to_driver_cosmology(par[0], par[1])
    print(
        f"\n  in Driver's h={R.H_DRIVER} units: log M* = {dms:.3f}, "
        f"log phi* = {dlp:.3f}"
    )
    g5 = R.driver_gama5(match_A=True)
    print(
        f"  Driver+22 GAMA-only (h=1, A={R.A_SCALE:g}): "
        f"{g5[0]:.3f} / {g5[1]:.3f} / {g5[2]:.3f} / {g5[3]:.3f}"
    )
    print(
        f"  Driver+22 GSR      (h=1): {R.TRUE['ms']:.3f} / {R.TRUE['lp']:.3f} / "
        f"{R.TRUE['al']:.3f} / {R.TRUE['be']:.3f}"
    )

    np.savetxt(
        f"vmax_{tag}_chains.csv",
        chains,
        delimiter=",",
        header="ms,lp,al,be",
        comments="",
    )
    print(f"  saved vmax_{tag}_chains.csv")
    plot(sets_raw, par, chains, tag, title, args)
    return par, chains


def plot(sets_raw, par, chains, tag, title, args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.size": 10,
            "legend.frameon": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    mg = np.linspace(12.5, 16, 400)

    idx = np.random.default_rng(0).choice(
        chains.shape[0], size=min(400, chains.shape[0]), replace=False
    )
    for k in idx:
        ax.plot(
            mg, np.log10(R.mrp_phi(mg, *chains[k])), color="cornflowerblue", alpha=0.02
        )
    ax.plot(
        mg,
        np.log10(R.mrp_phi(mg, *par)),
        color="royalblue",
        lw=2,
        label="MRP fit (this work, Driver method)",
    )

    lx, ly, _ = R.lcdm_curve()
    ax.plot(lx, ly, "k--", lw=1.8, label=r"$\Lambda$CDM (Murray+21, z=0.1)")
    g5 = R.driver_gama5(match_A=True)
    ax.plot(
        mg,
        np.log10(R.mrp_phi(mg, *g5)),
        color="darkred",
        lw=1.5,
        ls=(0, (6, 2)),
        label="Driver+22 GAMA only",
    )
    ax.plot(
        mg,
        np.log10(R.mrp_phi(mg, R.TRUE["ms"], R.TRUE["lp"], R.TRUE["al"], R.TRUE["be"])),
        color="red",
        lw=1.4,
        ls=":",
        label="Driver+22 GSR",
    )

    styles = {
        "GAMA": ("o", "crimson"),
        "SDSS": ("P", "purple"),
        "REFLEX": ("D", "forestgreen"),
    }
    for d in sets_raw:
        key = next((k for k in styles if d["label"].startswith(k)), "GAMA")
        mk, col = styles[key]
        e_hi = np.log10(1 + np.clip(d["f"], 0, 0.99))
        e_lo = -np.log10(1 - np.clip(d["f"], 0, 0.99))
        yy = d.get("y_corr", d["y"])
        if "y_corr" in d:
            # the raw 1/Vmax points, before the Eddington correction that the
            # fit actually used -- shown faint so the size of the correction is
            # visible and the fit is not compared against uncorrected data
            ax.plot(d["m"], d["y"], mk, ms=4, color=col, alpha=0.25, zorder=3)
        ax.errorbar(
            d["m"],
            yy,
            yerr=[e_lo, e_hi],
            fmt=mk,
            ms=5,
            color=col,
            capsize=2,
            lw=1,
            label=d["label"] + " (Eddington-corr.)",
            zorder=5,
        )

    if args.show_2pigg:
        comp = R._load_comparison(args.data_dir)
        tp = comp.get("2PIGG (Eke+08)")
        if tp is not None:
            ax.errorbar(
                tp["x"],
                tp["y"],
                yerr=[np.abs(tp["elo"]), np.abs(tp["ehi"])],
                fmt="s",
                ms=4,
                color="0.5",
                capsize=2,
                lw=1,
                label="2PIGG (Eke+08, not fitted)",
                zorder=4,
            )

    ax.set(
        xlim=(12.75, 16),
        ylim=(-8, -2),
        xlabel=r"$\log_{10}(M_{\rm halo}/M_\odot)$  [$h=1$]",
        ylabel=r"$\log_{10}$ number density [Mpc$^{-3}$ dex$^{-1}$]",
        title=title,
    )
    ax.legend(fontsize=8, loc="lower left")
    fig.savefig(f"hmf_vmax_{tag}.pdf", bbox_inches="tight")
    print(f"  saved hmf_vmax_{tag}.pdf")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gama-only", action="store_true")
    ap.add_argument("--combined", action="store_true")
    ap.add_argument(
        "--gama-fits",
        default="/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits",
    )
    ap.add_argument("--gama-area", type=float, default=238.11)
    ap.add_argument("--mass-col", default="MassA")
    ap.add_argument(
        "--sdss-parquet",
        default="/Users/00115372/Desktop/my_tools/nessie_tutorials/"
        "python/SDSS/sdss_groups.parquet",
    )
    ap.add_argument("--sdss-frac", type=float, default=0.2126803)
    ap.add_argument("--sdss-zmin", type=float, default=0.01)
    ap.add_argument("--sdss-zmax", type=float, default=0.08)
    ap.add_argument("--data-dir", default="../data")
    ap.add_argument("--show-2pigg", action="store_true", default=True)
    ap.add_argument(
        "--sdss-mmax",
        type=float,
        default=None,
        help="drop SDSS groups above this logM. The Nessie SDSS masses "
        "reach 15.7 and sit ~1.4 dex above REFLEX at the same mass, "
        "which is not credible in 8770 deg^2 out to z=0.08; try 15.0",
    )
    ap.add_argument("--no-eddington", action="store_true")
    a = ap.parse_args()

    if not (a.gama_only or a.combined):
        a.gama_only = a.combined = True

    print("=" * 68)
    print("  Driver+22 method (binned 1/Vmax, chi^2 MRP) on the new data")
    print(
        f"  h = {R.H0 / 100:g}, Om = {R.OMEGA_M}, A = {R.A_SCALE:g}, bin = {LOGBIN} dex"
    )
    print("=" * 68)

    if a.gama_only:
        print("\n--- GAMA only ---")
        g = load_gama(a)
        run([g], "gama", "GAMA only, Driver+22 method", a, eddington=not a.no_eddington)

    if a.combined:
        print("\n--- GAMA + SDSS + REFLEX II (Driver's GSR) ---")
        sets = [load_gama(a), load_sdss(a)]
        rf = load_reflex(a.data_dir)
        if rf is not None:
            sets.append(rf)
        run(
            sets,
            "gsr",
            "GAMA + SDSS + REFLEX II, Driver+22 method",
            a,
            eddington=not a.no_eddington,
        )


if __name__ == "__main__":
    main()
