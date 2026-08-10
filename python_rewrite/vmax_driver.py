"""
============================================================
driver_hmf.py -- a faithful port of Driver+22's gamahmf.r
============================================================

This follows gamahmf.r line by line rather than reconstructing the method.
The pieces that matter, and that earlier guesses got wrong:

 1. Vmax comes from the MEMBER GALAXIES, not from a mass limit.

        zmax_group = the multi-th largest zmax_19p8 among its members

    i.e. the redshift at which the 5th-brightest member would drop below the
    magnitude limit, so the group would fall below N >= 5. Nothing is cut on
    mass, every group counts, and each carries its own volume.

        vmax   = V(zmax) - V(zmin)
        weight = clip(vmax, vlimit/1000, vlimit)

 2. The Eddington correction is a Monte-Carlo ratio of counts, not a
    deconvolution of the model:

        for 1001 trials: mockmass = logM + N(0, log10MassErr)
                         mockcounts = weighted.hist(mockmass, w=1/weight)
        edb   = mean(mockcounts) / counts
        phi   = counts / (logbin * edb)

    Smearing the data increases the counts near the knee, so dividing by edb
    removes that excess. No model is involved.

 3. Errors combine a Monte-Carlo term with Poisson:

        mcerr    = sqrt(quantile((mean - mock)^2, 0.66)) / mean
        rootnerr = 1/sqrt(raw counts)
        f        = sqrt(mcerr^2 + rootnerr^2),  capped at 0.9999

 4. The fit is chi^2 in log space with sigma_log = f/ln(10), plus the Poisson
    penalty 2*vlimit*sum(phi(m > max)) * logbin, over bins with logM > mlimit.

 5. Errors on the parameters come from refitting perturbed realisations, with
    cosmic variance cosvar(vlimit/3, 3) applied per field.

COSMOLOGY. gamahmf.r uses h=0.6737, Om=0.3147 and multiplies the catalogue
masses by 100/ho to put them in Msun/h. recovery.py works natively at h=1, so
the masses are already in those units and no rescaling is applied here. Set
--driver-cosmology to reproduce his numbers exactly instead.

Run:
    python driver_hmf.py --gama-only
    python driver_hmf.py --combined
============================================================
"""

import argparse
import os

import numpy as np
from scipy import optimize

import recovery as R

LOGBIN = 0.2
MASSX = np.arange(10.3, 16.1 + 1e-9, LOGBIN)  # gamahmf.r: seq(10.3,16.1,logbin)
MIDS = 0.5 * (MASSX[1:] + MASSX[:-1])
N_MC = 1001  # as gamahmf.r
MULTI = 5


def cosvar(V, N):
    lv = np.log10(V)
    return ((219.7 - 52.4 * lv + 3.21 * lv**2) / np.sqrt(N)) / 100.0


# ----------------------------------------------------------------------
def group_zmax(
    gig_path, group_ids, nfof, multi=MULTI, zcol="zmax_19p8", idcol="GroupID"
):
    """zmax per group: the multi-th largest member zmax (2nd for pairs).

    Verbatim from gamahmf.r:
        if(Nfof==2) zmax = sort(member zmax, dec)[2]
        else        zmax = sort(member zmax, dec)[multi]
    """
    import pandas as pd

    gig = pd.read_csv(gig_path)
    if zcol not in gig.columns:
        cand = [c for c in gig.columns if "zmax" in c.lower()]
        raise KeyError(
            f"{zcol!r} not in {os.path.basename(gig_path)}; "
            f"zmax-like columns present: {cand}"
        )
    print(f"  member catalogue: {len(gig)} galaxies, using {zcol!r}")

    grp = gig.groupby(idcol)[zcol]
    out = np.full(len(group_ids), np.nan)
    for i, (gid, n) in enumerate(zip(group_ids, nfof)):
        try:
            v = np.sort(grp.get_group(gid).values)[::-1]
        except KeyError:
            continue
        k = 2 if int(n) == 2 else int(multi)
        if v.size >= k:
            out[i] = v[k - 1]
    miss = int(np.isnan(out).sum())
    if miss:
        print(f"  !! {miss}/{out.size} groups have no usable member zmax")
    return out


def vmax_weights(zmax, zfof, sky_frac, zmin, zlimit):
    """vmax = V(zmax) - V(zmin), clipped to [vlimit/1000, vlimit] as gamahmf.r."""
    zmax = np.where(np.isnan(zmax), zfof, zmax)
    zmax = np.clip(zmax, zfof, zlimit)  # ifelse(...) lines
    d_max = R.comoving_distance(zmax)
    d_min = R.comoving_distance(np.full_like(zmax, zmin))
    vmax = (4 / 3) * np.pi * (d_max**3 - d_min**3) * sky_frac
    d_lim = R.comoving_distance(np.array([zlimit]))[0]
    vlimit = (4 / 3) * np.pi * (d_lim**3 - d_min[0] ** 3) * sky_frac
    w = np.clip(vmax, vlimit / 1000.0, vlimit)
    return vmax, w, float(vlimit)


def hmf_with_edb(log_mass, weights, mass_err, n_mc=N_MC, seed=0):
    """Weighted 1/Vmax histogram with the Monte-Carlo Eddington correction and
    the combined MC + Poisson fractional error, exactly as gamahmf.r."""
    rng = np.random.default_rng(seed)
    counts_raw, _ = np.histogram(log_mass, bins=MASSX)  # gamahmf
    counts_w, _ = np.histogram(log_mass, bins=MASSX, weights=1.0 / weights)  # gamahmf2

    mock = np.empty((n_mc, MIDS.size))
    for i in range(n_mc):
        mm = log_mass + rng.normal(0.0, mass_err)
        mock[i], _ = np.histogram(mm, bins=MASSX, weights=1.0 / weights)
    meanc = mock.mean(axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        edb = meanc / counts_w
    edb[~np.isfinite(edb)] = 1.0

    with np.errstate(divide="ignore", invalid="ignore"):
        mcerr = np.sqrt(np.quantile((meanc[None, :] - mock) ** 2, 0.66, axis=0)) / meanc
        rootn = np.where(
            counts_raw > 0, 1.0 / np.sqrt(np.maximum(counts_raw, 1)), np.inf
        )
        y = counts_w / (LOGBIN * edb)
    mcerr = np.nan_to_num(mcerr, nan=0.0, posinf=1.0, neginf=0.0)
    rootn = np.nan_to_num(rootn, nan=1.0, posinf=1.0, neginf=1.0)
    f = np.sqrt(mcerr**2 + rootn**2)
    f[~np.isfinite(f)] = 0.9999
    f = np.where(f >= 1, 0.9999, f)

    ok = np.isfinite(y) & (y > 0) & (counts_raw > 0)
    print(
        f"  Eddington (MC): median edb = {np.nanmedian(edb[ok]):.3f} "
        f"(>1 means smearing inflates the counts)"
    )
    return MIDS, y, f, counts_raw, edb, ok


# ----------------------------------------------------------------------
def massfn(par, allx, ally, allf, vlimit, use_penalty=True):
    """gamahmf.r's massfn(): chi^2 in log space plus the Poisson penalty."""
    mstar, phi, alpha, beta = par
    if phi <= 0 or beta <= 0:
        return 1e12
    model = np.log10(
        beta
        * np.log(10)
        * np.exp(-(10 ** (beta * (allx - mstar))))
        * (phi * (10**allx / 10**mstar) ** (alpha + 1))
    )
    chi2 = float(np.sum(((ally - model) / (allf / np.log(10))) ** 2))
    if use_penalty:
        xx = allx.max() + np.arange(1, 11) * LOGBIN
        pen = (
            2
            * vlimit
            * float(
                np.sum(
                    np.log(10)
                    * beta
                    * np.exp(-(10 ** (beta * (xx - mstar))))
                    * (phi * (10**xx / 10**mstar) ** (alpha + 1))
                )
            )
            * LOGBIN
        )
        chi2 += pen
    return chi2 if np.isfinite(chi2) else 1e12


# sane ranges; the chi^2 surface has flat directions and Nelder-Mead will
# happily wander to beta ~ 0 or phi* ~ 1 without them
BOUNDS = [(11.5, 16.5), (-8.0, -0.5), (-2.5, -0.3), (0.15, 1.8)]


def fit(allx, ally, allf, vlimit, p0, use_penalty=True):
    """gamahmf.r calls optim(..., parscale=c(1,1,1,0.5)), which rescales the
    parameters internally. Without an equivalent, scipy would be optimising
    mstar ~ 14 alongside phi ~ 4e-4 -- four orders of magnitude apart -- and the
    simplex collapses. Fitting log10(phi) instead puts every parameter at O(1),
    which is the same trick and numerically identical.

    p0 is given with phi LINEAR, as in gamahmf.r; it is converted here."""
    q0 = np.array([p0[0], np.log10(max(p0[1], 1e-12)), p0[2], p0[3]], float)
    q0 = np.array([np.clip(v, lo, hi) for v, (lo, hi) in zip(q0, BOUNDS)])

    def obj(q):
        for v, (lo, hi) in zip(q, BOUNDS):
            if not (lo <= v <= hi):
                return 1e12
        return massfn(
            [q[0], 10 ** q[1], q[2], q[3]], allx, ally, allf, vlimit, use_penalty
        )

    best, bq = np.inf, q0
    # a few restarts, since the surface is not convex
    for jitter in (0.0, 0.15, -0.15):
        start = q0 + np.array([jitter, jitter, 0.5 * jitter, 0.1 * jitter])
        start = np.array([np.clip(v, lo, hi) for v, (lo, hi) in zip(start, BOUNDS)])
        r = optimize.minimize(
            obj,
            start,
            method="Nelder-Mead",
            options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-8),
        )
        if r.fun < best:
            best, bq = r.fun, r.x
    return np.array([bq[0], 10 ** bq[1], bq[2], bq[3]])


def mc_params(
    allx, ally, allf, vlimit, p0, cosvariance, n=200, seed=1, use_penalty=True
):
    """Parameter errors: perturb the points by their fractional error and by
    cosmic variance, refit each time."""
    rng = np.random.default_rng(seed)
    out = np.full((n, 4), np.nan)
    for i in range(n):
        lin = 10**ally
        lin = lin + lin * rng.normal(0.0, np.clip(allf, 0, 0.999), size=lin.size)
        lin = lin * 10 ** rng.normal(0.0, cosvariance)
        g = lin > 0
        try:
            out[i] = fit(allx[g], np.log10(lin[g]), allf[g], vlimit, p0, use_penalty)
        except Exception:
            pass
    good = np.isfinite(out).all(axis=1)
    good &= (out[:, 1] > 0) & (out[:, 3] > 0)
    return out[good]


# ----------------------------------------------------------------------
def build_gama(a):
    from astropy.io import fits as afits

    with afits.open(a.gama_fits) as h:
        t = h[1].data
    cols = set(t.columns.names)
    nfof = np.asarray(t["Nfof"], float)
    zfof = np.asarray(t["Zfof"], float)
    mafunc = np.asarray(t["MassAfunc"], float)
    gid = np.asarray(t["GroupID"])

    sel = (nfof > MULTI - 1) & (zfof < R.ZLIMIT) & (zfof > R.ZMIN) & (mafunc > 1e1)
    if a.regions and "GAMARegion" in cols:
        reg = np.asarray(t["GAMARegion"]).astype(str)
        sel &= np.isin(reg, a.regions)
    print(f"  groups after selection: {int(sel.sum())}")

    mass = np.asarray(t[a.mass_col], float)[sel]
    nfof, zfof, gid = nfof[sel], zfof[sel], gid[sel]
    log_mass = np.log10(mass)

    # gamahmf.r's multiplicity -> sigma table
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
    err = np.interp(nfof, xx, yy, left=np.nan, right=np.nan)
    err = np.where(np.isfinite(err), err, 0.03)
    err = np.where(err < 0.1, 0.1, err)

    sky_frac = a.gama_area * (np.pi / 180) ** 2 / (4 * np.pi)
    zmax = group_zmax(a.gig, gid, nfof, multi=MULTI, zcol=a.zmax_col)
    have = np.isfinite(zmax)
    if (~have).any():
        print(
            f"  dropping {int((~have).sum())} groups with no member zmax "
            f"(they would otherwise be given zmax = Zfof and a huge weight)"
        )
        log_mass, nfof, zfof, err, zmax = (
            log_mass[have],
            nfof[have],
            zfof[have],
            err[have],
            zmax[have],
        )
    vmax, w, vlimit = vmax_weights(zmax, zfof, sky_frac, R.ZMIN, R.ZLIMIT)
    print(
        f"  vlimit = {vlimit:.4e} Mpc^3 (h=1); median vmax/vlimit = "
        f"{np.median(vmax) / vlimit:.3f}"
    )

    x, y, f, cnt, edb, ok = hmf_with_edb(log_mass, w, err)
    keep = ok & (x > a.mlimit)
    print(f"  bins kept: {int(keep.sum())} with logM > {a.mlimit}")
    return dict(
        x=x[keep],
        y=np.log10(y[keep]),
        f=f[keep],
        cnt=cnt[keep],
        vlimit=vlimit,
        nfield=len(a.regions) if a.regions else 4,
        label="GAMA (this work)",
    )


def build_reflex(data_dir):
    comp = R._load_comparison(data_dir)
    rf = comp.get("REFLEX II (Böhringer+17)")
    if rf is None:
        return None
    x, y = np.asarray(rf["x"], float), np.asarray(rf["y"], float)
    f = np.full(x.size, 1 / np.sqrt(20))
    f[0] = f[-1] = 1 / np.sqrt(3)
    return dict(
        x=x,
        y=y,
        f=f,
        cnt=np.full(x.size, 20),
        vlimit=1.3e7,
        nfield=1,
        label="REFLEX II (Böhringer+17)",
    )


# ----------------------------------------------------------------------
def report(name, par, chains, sets):
    print(f"\n  === {name} ===")
    lp = np.log10(chains[:, 1])
    med = [
        np.median(chains[:, 0]),
        np.median(lp),
        np.median(chains[:, 2]),
        np.median(chains[:, 3]),
    ]
    q16 = [
        np.percentile(chains[:, 0], 16),
        np.percentile(lp, 16),
        np.percentile(chains[:, 2], 16),
        np.percentile(chains[:, 3], 16),
    ]
    q84 = [
        np.percentile(chains[:, 0], 84),
        np.percentile(lp, 84),
        np.percentile(chains[:, 2], 84),
        np.percentile(chains[:, 3], 84),
    ]
    best = [par[0], np.log10(par[1]), par[2], par[3]]
    for i, nm in enumerate(["log M*", "log phi*", "alpha", "beta"]):
        print(
            f"  {nm:>9} {best[i]:9.3f}   MC {med[i]:8.3f} "
            f"[{q16[i]:7.3f}, {q84[i]:7.3f}]"
        )
    dms, dlp, _, _ = R.to_driver_cosmology(best[0], best[1])
    print(f"  in Driver's h={R.H_DRIVER}: log M* = {dms:.3f}, log phi* = {dlp:.3f}")
    g5 = R.driver_gama5(match_A=True)
    print(
        f"  Driver+22 GAMA-only (h=1, A={R.A_SCALE:g}): "
        f"{g5[0]:.3f} / {g5[1]:.3f} / {g5[2]:.3f} / {g5[3]:.3f}"
    )
    return np.column_stack([chains[:, 0], lp, chains[:, 2], chains[:, 3]])


def plot(sets, par, chains, tag, title, data_dir):
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
        chains.shape[0], size=min(300, chains.shape[0]), replace=False
    )
    for k in idx:
        p = chains[k]
        ax.plot(
            mg,
            np.log10(R.mrp_phi(mg, p[0], np.log10(max(p[1], 1e-30)), p[2], p[3])),
            color="cornflowerblue",
            alpha=0.02,
        )
    ax.plot(
        mg,
        np.log10(R.mrp_phi(mg, par[0], np.log10(par[1]), par[2], par[3])),
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

    style = {
        "GAMA": ("o", "crimson"),
        "SDSS": ("P", "purple"),
        "REFLEX": ("D", "forestgreen"),
    }
    for d in sets:
        mk, col = next(
            (v for k, v in style.items() if d["label"].startswith(k)), ("o", "crimson")
        )
        ehi = np.log10(1 + np.clip(d["f"], 0, 0.99))
        elo = -np.log10(1 - np.clip(d["f"], 0, 0.99))
        ax.errorbar(
            d["x"],
            d["y"],
            yerr=[elo, ehi],
            fmt=mk,
            ms=5,
            color=col,
            capsize=2,
            lw=1,
            label=d["label"],
            zorder=5,
        )

    comp = R._load_comparison(data_dir)
    tp = comp.get("2PIGG (Eke+08)")
    if tp is not None:
        ax.errorbar(
            tp["x"],
            tp["y"],
            yerr=[np.abs(tp["elo"]), np.abs(tp["ehi"])],
            fmt="s",
            ms=4,
            color="0.55",
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
    fig.savefig(f"hmf_driver_{tag}.pdf", bbox_inches="tight")
    print(f"  saved hmf_driver_{tag}.pdf")


def run(sets, tag, title, a):
    allx = np.concatenate([d["x"] for d in sets])
    ally = np.concatenate([d["y"] for d in sets])
    allf = np.concatenate([d["f"] for d in sets])
    vlimit = sets[0]["vlimit"]

    p0 = (13.9, 10**-3.4, -1.4, 0.6)  # phi LINEAR, as gamahmf.r
    par = fit(allx, ally, allf, vlimit, p0, use_penalty=not a.no_penalty)
    c2 = massfn(par, allx, ally, allf, vlimit, use_penalty=not a.no_penalty)
    print(
        f"  best-fit chi2 = {c2:.1f} over {allx.size} bins "
        f"({c2 / max(allx.size - 4, 1):.2f} per dof)"
    )
    cv = cosvar(vlimit / max(sets[0]["nfield"], 1), sets[0]["nfield"])
    print(f"  cosmic variance = {cv:.4f}")
    print(f"  {a.n_mc} Monte-Carlo refits ...")
    chains = mc_params(
        allx, ally, allf, vlimit, par, cv, n=a.n_mc, use_penalty=not a.no_penalty
    )
    print(f"  {chains.shape[0]}/{a.n_mc} converged")
    out = report(title, par, chains, sets)
    np.savetxt(
        f"driver_{tag}_chains.csv",
        out,
        delimiter=",",
        header="ms,lp,al,be",
        comments="",
    )
    plot(sets, par, chains, tag, title, a.data_dir)
    return par, chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gama-only", action="store_true")
    ap.add_argument("--combined", action="store_true")
    ap.add_argument(
        "--gama-fits",
        default="/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits",
    )
    ap.add_argument(
        "--gig",
        default="../data/GAMAGalsInGroups.csv",
        help="member-galaxy catalogue holding the per-galaxy zmax",
    )
    ap.add_argument(
        "--zmax-col",
        default="zmax_19p8",
        help="per-galaxy zmax column; the new DMU may name it differently",
    )
    ap.add_argument("--gama-area", type=float, default=238.11)
    ap.add_argument("--mass-col", default="MassA")
    ap.add_argument(
        "--mlimit",
        type=float,
        default=12.7,
        help="lowest mass bin used in the fit (gamahmf.r's mlimit)",
    )
    ap.add_argument("--regions", nargs="+", default=None)
    ap.add_argument("--data-dir", default="../data")
    ap.add_argument("--n-mc", type=int, default=200)
    ap.add_argument("--no-penalty", action="store_true")
    a = ap.parse_args()
    if not (a.gama_only or a.combined):
        a.gama_only = a.combined = True

    print("=" * 68)
    print("  Driver+22 gamahmf.r, ported")
    print(
        f"  h = {R.H0 / 100:g}, Om = {R.OMEGA_M}, mass = {a.mass_col}, "
        f"bin = {LOGBIN} dex, mlimit = {a.mlimit}"
    )
    print("=" * 68)

    g = build_gama(a)
    if a.gama_only:
        run([g], "gama", "GAMA only (Driver method)", a)
    if a.combined:
        sets = [g]
        rf = build_reflex(a.data_dir)
        if rf is not None:
            sets.append(rf)
        run(sets, "gsr", "GAMA + REFLEX II (Driver method)", a)


if __name__ == "__main__":
    main()
