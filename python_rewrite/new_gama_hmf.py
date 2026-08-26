#!/usr/bin/env python3
"""Driver et al. (2022) figure 4, applied to the new GAMA group catalogue.

Step 2 of the exercise: run the *identical* 1/Vmax + Eddington + MRP method of
``gamahmf.r`` on the new GAMA DMU (Nessie group finder, ProFound photometry,
four regions including G23, r < 19.65), and plot it against Driver's original
result on the old data.  Everything numerical is imported from
``driver_recovery.py``, which is verified bit-for-bit against Driver's R script,
so any difference in the output is a difference in the *data*, not the method.

The one thing the new DMU does not provide is a per-galaxy ``zmax``, which the
whole Vmax calculation rests on.  It is computed here by inverting the DMU's own
selection:

    apparent(z) = AbsoluteMagR + distmod(z) + ke(z)  ==  19.65

using exactly the k+e polynomial and the h = 1, OmegaM = 0.25 cosmology that
``make_gama_dmu`` used to build ``AbsoluteMagR`` in the first place (verified to
3e-6 mag).  That makes zmax self-consistent with the catalogue.

Usage
-----
    python new_gama_hmf.py                    # comparison figure
    python new_gama_hmf.py --zmax-control     # also run old data with recomputed
                                              # zmax, to separate the zmax method
                                              # from the catalogue change
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

import driver_recovery as dr
from driver_recovery import (
    G, HO, LOGBIN, MAGICA, MASSCORR, MLIMIT, MSOL, MULTI, NFOF_XX, NFOF_YY,
    PARSEC, ZLIMIT, ZMIN, bin_hmf, co_vol, fit_mrp, fit_selection, lcdm_curve,
    monte_carlo_fits, mrp, r_approx, survey_volume,
)

# ---------------------------------------------------------------------------
# New-DMU constants.  From make_gama_dmu/config.py and the DMU notes.
# ---------------------------------------------------------------------------

NEW_AREA = 238.11                  # deg^2, four regions (fractional 0.005771988)
NEW_MAGLIM = 19.65                 # config.APPARENT_MAG_LIM  (was 19.8)
NEW_REGIONS = ("g09", "g12", "g15", "g23")

# The group finder's cosmology -- config.py: FlatCosmology(1.0, 0.25).
# This is the frame AbsoluteMagR lives in, so zmax must be inverted in it.
DMU_COSMO = FlatLambdaCDM(H0=100, Om0=0.25)

# make_gama_dmu.main.calc_ke_correction
KCORR = np.array([0.20848, 1.0226, 0.52366, 3.5902, 2.3843])

GROUPS_NEW = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"
GALS_NEW = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CGal.fits"


def ke_correction(z):
    """make_gama_dmu's k+e correction (main.py line 112)."""
    z = np.asarray(z, dtype=float)
    return ((z[..., None] - 0.2) ** np.arange(len(KCORR))) @ KCORR - 1.75 * z


def _appmag_offset_grid(ngrid=6000, zhi=2.0):
    """Tabulate distmod(z) + ke(z), which is monotonic, for inversion."""
    zg = np.logspace(np.log10(1e-4), np.log10(zhi), ngrid)
    off = DMU_COSMO.distmod(zg).value + ke_correction(zg)
    if not np.all(np.diff(off) > 0):
        raise RuntimeError("apparent-magnitude offset is not monotonic in z")
    return zg, off


def compute_zmax(absmag, maglim=NEW_MAGLIM):
    """The redshift at which a galaxy of this absolute magnitude hits the limit.

    Inverts ``apparent(z) = absmag + distmod(z) + ke(z) = maglim``.  This is the
    replacement for the old member file's ``zmax_19p8`` column, which the new
    DMU does not emit and which would in any case be wrong here: it was computed
    at r < 19.8, so it would hand every group too large a volume.
    """
    zg, off = _appmag_offset_grid()
    return np.interp(np.asarray(maglim, float) - np.asarray(absmag, float), off, zg)


# ---------------------------------------------------------------------------
# Group building
# ---------------------------------------------------------------------------

def _masses_and_errors(g3c, masserr="driver"):
    """gamahmf.r lines 263-289, myoption="GAMA".

    The mass is rebuilt from the velocity dispersion and Rad50 rather than read
    from the catalogue's MassAfunc, exactly as Driver does.  The (100/ho) is
    already inside mymass -- no further h conversion.  Both catalogues store
    Rad50 in Mpc/h and VelDisp in km/s at h = 1, so the formula transfers
    unchanged.
    """
    g3c["mymass"] = (MAGICA * (g3c.VelDisp * 1000) ** 2 * g3c.Rad50 * PARSEC * 1e6
                     / (G * MSOL) * (100 / HO))

    err, mc = dr.mass_error_and_bias(g3c.Nfof.values, masserr)
    g3c["log10MassErr"] = err
    g3c["masscorr"] = mc
    g3c["MassAfunc"] = g3c.mymass / 10 ** g3c.masscorr
    return g3c


def _vmax_from_members(g3c, zmax_by_group, area, verbose=True, vmax_floor_frac=1e-3):
    """gamahmf.r lines 295-308, with survey area and Vmax floor as parameters."""
    zmax = np.full(len(g3c), np.nan)
    for i, (gid, nfof) in enumerate(zip(g3c.GroupID.values, g3c.Nfof.values)):
        arr = zmax_by_group.get(gid)
        if arr is None:
            continue
        k = 2 if nfof == 2 else MULTI
        if len(arr) >= k:
            zmax[i] = arr[k - 1]
    n_na = int(np.isnan(zmax).sum())

    zmax = np.where(zmax < g3c.Zfof.values, g3c.Zfof.values, zmax)     # line 302
    zmax = np.where(zmax > ZLIMIT, ZLIMIT, zmax)                       # line 303
    g3c["zmax"] = zmax

    vmax = survey_volume(zmax, area) - survey_volume(np.array([ZMIN]), area)[0]
    g3c["vmax"] = vmax
    vlimit = survey_volume(np.array([ZLIMIT]), area)[0]
    vlimitmin = vlimit * vmax_floor_frac

    # Lines 307-308 as written: the second assignment overwrites the first and
    # discards the upper clip.  Reproduced deliberately.
    _ = np.where(vmax > vlimit, vlimit, vmax)
    g3c["weightszlimit"] = np.where(vmax < vlimitmin, vlimitmin, vmax)

    if verbose:
        print(f"  groups with NA zmax    : {n_na}")
        print(f"  vlimit                 : {vlimit:.6e} Mpc^3")
        print(f"  median vmax/vlimit     : {np.median(vmax) / vlimit:.4f}")
        print(f"  fraction at zlimit     : {np.mean(zmax >= ZLIMIT - 1e-12):.4f}")
    return g3c, vlimit


def build_groups_new(group_file=GROUPS_NEW, gal_file=GALS_NEW, verbose=True,
                     vmax_floor_frac=1e-3, masserr="driver"):
    """The new DMU, put through Driver's selection and Vmax construction."""
    g = Table.read(group_file).to_pandas()
    g["GAMARegion"] = g.GAMARegion.str.decode("utf-8") if g.GAMARegion.dtype == object \
        else g.GAMARegion

    # gamahmf.r line 258, minus the IterCenDec > -3.5 cut.  That cut selects the
    # three equatorial fields and would delete G23 entirely; the region column
    # is the correct selector here.
    sel = ((g.Nfof > MULTI - 1) & (g.Zfof < ZLIMIT) & (g.Zfof > ZMIN)
           & (g.MassAfunc > 1e1) & (g.GAMARegion.isin(NEW_REGIONS)))
    g3c = g[sel].copy().reset_index(drop=True)

    g3c = _masses_and_errors(g3c, masserr)

    # NOTE: gamahmf.r line 310 forces GroupID 100622 to 1e9 -- a known-bad object
    # in the *old* catalogue.  GroupID 100622 also exists here (g09 offset is
    # 1e5) but is an entirely different group, so that line is NOT applied.

    gal = Table.read(gal_file).to_pandas()
    gal = gal[gal.GroupID != 0]                     # 0 is the ungrouped sentinel
    gal = gal[gal.GroupID.isin(set(g3c.GroupID.values.tolist()))]
    gal = gal[np.isfinite(gal.AbsoluteMagR) & (gal.AbsoluteMagR > -40)
              & (gal.AbsoluteMagR < 0)]
    gal["zmax"] = compute_zmax(gal.AbsoluteMagR.values)

    zmax_by_group = {gid: np.sort(sub.zmax.values)[::-1]
                     for gid, sub in gal.groupby("GroupID")}

    if verbose:
        print(f"  groups after selection : {len(g3c)}")
        by = g3c.GAMARegion.value_counts().reindex(NEW_REGIONS, fill_value=0)
        print("  by region              : "
              + ", ".join(f"{k} {int(v)}" for k, v in by.items()))
    return _vmax_from_members(g3c, zmax_by_group, NEW_AREA, verbose,
                              vmax_floor_frac)


def build_groups_old_recomputed_zmax(verbose=True):
    """Control: the OLD catalogue with zmax recomputed the same way as the new.

    The old member file ships an official ``zmax_19p8`` built with per-galaxy
    k-corrections; the new data has to use a mean k+e polynomial instead, which
    sits ~+0.01 in z high.  Running the old data through the polynomial isolates
    how much of any old/new difference is that methodological change rather than
    the catalogue itself.
    """
    g3c, _ = dr.build_groups("../data/G3CFoFGroupv10.fits",
                             "../data/GAMAGalsInGroups.csv", verbose=False)
    gig = pd.read_csv("../data/GAMAGalsInGroups.csv")
    gig = gig[(gig.GroupID != 0) & (gig.Z > 0)]
    gig = gig[gig.GroupID.isin(set(g3c.GroupID.values.tolist()))]
    mabs = (gig.Rpetro.values - DMU_COSMO.distmod(gig.Z.values).value
            - ke_correction(gig.Z.values))
    gig = gig.assign(zmax=compute_zmax(mabs, maglim=19.8))     # old limit
    zmax_by_group = {gid: np.sort(sub.zmax.values)[::-1]
                     for gid, sub in gig.groupby("GroupID")}
    if verbose:
        print(f"  groups after selection : {len(g3c)}")
    return _vmax_from_members(g3c, zmax_by_group, dr.AREA, verbose)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

OLD_C = (100 / 255, 149 / 255, 237 / 255)      # Driver's cornflower
PUBLISHED_C = "#c8781e"                        # Driver+22's published fit
# his abstract / table 2 GSR row -- one fixed reference line on every figure
DRIVER_ABSTRACT_FIT = (14.13, -3.96, -1.68, 0.63)
NEW_C = "#c2185b"


MSUN_R = 4.65      # absolute r-band magnitude of the Sun, for M/L only


def ml_excess(g3c, nmin=15, edges=None):
    """dex by which each group's mass-to-light ratio exceeds the running median
    at the same mass.

    A group's M/L is the one quantity that says whether its *mass* is credible
    independently of its Vmax.  Real clusters sit near log10(M/L_r) ~ 2.5-3.5
    and the ratio rises smoothly with mass, so the residual against a running
    median is a clean, nearly mass-independent outlier statistic.

    This exists because of GAMA 205509: N_fof = 5 at z = 0.019 with
    sigma = 482 km/s, giving logM = 14.26 from one M_r = -19.5 galaxy and four
    dwarfs spread over 900 km/s.  At that redshift GAMA reaches M_r ~ -15, so a
    genuine 10^14.3 halo would have hundreds of members, not five.  Its
    log10(M/L) is 4.38 against a median of 3.22 for its mass -- a factor of 14.
    Sitting at vmax/vlimit = 2.4e-3 it then carried **38% of the 14.2 bin**.

    Only the *relative* excess is used, so the absolute M/L zero point (and
    hence ``MSUN_R``) does not matter.
    """
    lm = np.log10(g3c.MassAfunc.values)
    lum = 10 ** (-0.4 * (g3c.TotRmag.values - MSUN_R))
    with np.errstate(divide="ignore", invalid="ignore"):
        lml = np.log10(g3c.MassAfunc.values / lum)
    if edges is None:
        edges = np.arange(12.4, 15.9, 0.2)
    cen, med = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (lm >= a) & (lm < b) & np.isfinite(lml)
        if m.sum() > nmin:
            cen.append(0.5 * (a + b))
            med.append(np.median(lml[m]))
    if len(cen) < 2:
        return np.zeros(len(g3c))
    return lml - np.interp(lm, cen, med)


def apply_ml_cut(g3c, thresh, verbose=True):
    """Drop groups whose M/L exceeds the running median by more than ``thresh``.

    ``thresh = 1.0`` is the recommended value: an order of magnitude more mass
    per unit light than comparable systems, well beyond the 97.5th percentile
    of the excess distribution (+0.60), so it removes only extreme outliers --
    5 of 1833 groups.  It is a statement about the mass being wrong, which is
    why it is preferred to a Vmax floor: the floor caps the *weight* of a group
    whose mass is still believed, and is a biased estimator, while this removes
    objects whose mass is not credible in the first place.
    """
    if thresh is None:
        return g3c, None
    ex = ml_excess(g3c)
    bad = ex > thresh
    if verbose:
        print(f"  M/L cut > {thresh} dex           : {int(bad.sum())} of "
              f"{len(g3c)} groups dropped")
    flagged = g3c[bad].copy()
    flagged["ml_excess"] = ex[bad]
    return g3c[~bad].copy(), flagged


def plot_comparison(results, mrpx, mrpy, factor, outfile):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import plotstyle
    plotstyle.apply()
    outfile = plotstyle.as_png(outfile)
    from matplotlib.collections import LineCollection

    fig = plt.figure(figsize=(7.4, 5.4))
    ax_top = fig.add_axes([0.115, 0.795, 0.875, 0.185])
    ax = fig.add_axes([0.115, 0.088, 0.875, 0.690])
    xfit = np.arange(0.0, 20 + 1e-9, dr.FITBINWID)

    for res in results:
        b = res["binned"]
        ax_top.step(b.gamax - 0.5 * LOGBIN, b.raw_counts, where="post",
                    color=res["colour"], lw=1.1, label=res["short"])
    ax_top.set_xlim(12, 16)
    ax_top.set_ylabel("Number", fontsize=9)
    ax_top.set_xticklabels([])
    ax_top.tick_params(direction="in", top=True, right=True, labelsize=8)
    ax_top.legend(fontsize=7.5, frameon=False, loc="upper right", ncol=1,
                  handlelength=1.4)

    # Monte-Carlo bands first, so the points sit on top
    for res in results:
        if not res.get("mc"):
            continue
        segs = []
        for par in res["mc"]["curves"][:600]:
            with np.errstate(divide="ignore", invalid="ignore"):
                yy = np.log10(mrp(xfit, *par))
            m = np.isfinite(yy) & (yy > -9) & (xfit > 11.5) & (xfit < 16.3)
            if m.sum() > 2:
                segs.append(np.column_stack([xfit[m], yy[m]]))
        ax.add_collection(LineCollection(segs, colors=[res["colour"]],
                                         linewidths=0.6, alpha=0.012))

    # Driver+22's published GAMA5 fit (table 2), behind the data for reference.
    # This is his printed answer, not our refit of his catalogue.
    # the same fixed reference on every figure: his abstract / table 2 GSR row
    _p = DRIVER_ABSTRACT_FIT
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(mrp(xfit, _p[0], 10 ** _p[1], _p[2], _p[3])),
                color=PUBLISHED_C, lw=3.2, alpha=0.9, zorder=2,
                label="Driver+22 published GSR fit")

    with np.errstate(divide="ignore"):
        ax.plot(mrpx - 0.08, np.log10(mrpy) - np.log10(factor) + 0.08,
                ls="--", color="black", lw=1.8, zorder=3, label="LCDM prediction")

    for res in results:
        b, c = res["binned"], res["colour"]
        with np.errstate(divide="ignore", invalid="ignore"):
            ly = np.log10(b.gamay)
            lo = np.abs(np.log10(1 - b.gamaf))
            hi = np.abs(np.log10(1 + b.gamaf))
        good = np.isfinite(ly)
        ax.errorbar(b.gamax[good], ly[good], yerr=[lo[good], hi[good]], fmt="none",
                    ecolor=c, elinewidth=0.9, capsize=1.5, alpha=0.9, zorder=4)
        himask = b.gamax > MLIMIT
        ax.plot(b.gamax[himask], ly[himask], ls="none", marker=res["marker"],
                color=c, ms=5, zorder=6)
        ax.plot(b.gamax[~himask], ly[~himask], ls="none", marker=res["marker"],
                mfc="none", mec=c, ms=5, zorder=6)
        with np.errstate(divide="ignore", invalid="ignore"):
            ax.plot(xfit, np.log10(mrp(xfit, *res["par"])), color=c, lw=2,
                    zorder=7, label=res["label"])

    ax.set_xlim(12, 16)
    ax.set_ylim(-8, -2)
    ax.set_xlabel(r"log$_{10}$(Halo Mass / M$_\odot$)", fontsize=10)
    ax.set_ylabel(r"log$_{10}$(number density) [Mpc$^{-3}$ dex$^{-1}$]", fontsize=10)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=9)
    ax.minorticks_on()
    ax.legend(fontsize=8.5, frameon=False, loc="lower left", handlelength=1.8)

    fig.savefig(outfile, dpi=240)
    plt.close(fig)
    print(f"\n  wrote {outfile}")


def run_one(g3c, vlimit, phimrp, seed, nmc_edb, nmc_fit, label, short, colour,
            marker, verbose=True, nboot=0, mmax=None):
    b, rng = bin_hmf(g3c, vlimit, seed=seed, nmc=nmc_edb, verbose=verbose,
                     nboot=nboot)
    allx, ally, allf = fit_selection(b, mmax=mmax)
    par, value, fncount, conv = fit_mrp(allx, ally, allf, vlimit, phimrp)
    mc = None
    if nmc_fit > 0:
        print(f"  {short}: Monte-Carlo error band, {nmc_fit} refits ...")
        mc = monte_carlo_fits(b, rng, phimrp, nmc=nmc_fit, mmax=mmax)
    return dict(binned=b, par=par, value=value, conv=conv, mc=mc, nfit=len(allx),
                label=label, short=short, colour=colour, marker=marker,
                ngroups=len(g3c), vlimit=vlimit, mfit=(allx.min(), allx.max()))


def report(results):
    print("\n" + "=" * 78)
    print("  MRP fits -- identical method, different data")
    print("=" * 78)
    hdr = f"  {'':<26}" + "".join(f"{r['short']:>17}" for r in results)
    print(hdr)
    rows = [("N groups", lambda r: f"{r['ngroups']}"),
            ("vlimit [Mpc^3]", lambda r: f"{r['vlimit']:.3e}"),
            ("bins fitted", lambda r: f"{r['nfit']}"),
            ("fitted range", lambda r: f"{r['mfit'][0]:.1f}-{r['mfit'][1]:.1f}"),
            ("log10(M*)", lambda r: f"{r['par'][0]:.3f}"),
            ("log10(phi*)", lambda r: f"{np.log10(r['par'][1]):.3f}"),
            ("alpha", lambda r: f"{r['par'][2]:.3f}"),
            ("beta", lambda r: f"{r['par'][3]:.3f}"),
            ("chi2", lambda r: f"{r['value']:.2f}")]
    for name, fn in rows:
        print(f"  {name:<26}" + "".join(f"{fn(r):>17}" for r in results))
    for r in results:
        if r.get("mc"):
            m = r["mc"]
            q = lambda a: (np.quantile(a[np.isfinite(a)], 0.16),
                           np.quantile(a[np.isfinite(a)], 0.84))
            lo, hi = q(m["mstar"])
            print(f"\n  {r['short']}: MC 16-84 log10(M*) {lo:.2f} to {hi:.2f}, "
                  f"alpha {q(m['alphastar'])[0]:.2f} to {q(m['alphastar'])[1]:.2f}")
    print("\n  Reminder: optim(maxit=500) is a truncated, non-converged walk, so "
          "each\n  fit is one draw from a wide spread (see driver_recovery.py "
          "--seed-scan).\n  Compare the binned points, which are stable, before "
          "the fitted parameters.")


def print_binned_comparison(results):
    print("\n  Binned HMF, log10(phi_corr)")
    print(f"  {'logM':>6}" + "".join(f"{r['short']:>20}" for r in results)
          + f"{'diff':>9}")
    b0 = results[0]["binned"]
    for i in range(len(b0.gamax))[::-1]:
        if b0.gamax[i] < 12.0 or b0.gamax[i] > 15.6:
            continue
        cells, vals = "", []
        for r in results:
            b = r["binned"]
            with np.errstate(divide="ignore", invalid="ignore"):
                v = np.log10(b.gamay[i])
            vals.append(v)
            n = int(b.raw_counts[i])
            cells += f"{v:>12.3f} ({n:4d})" if np.isfinite(v) else f"{'--':>12} ({n:4d})"
        d = vals[-1] - vals[0]
        print(f"  {b0.gamax[i]:6.1f}" + cells
              + (f"{d:>+9.3f}" if np.isfinite(d) else f"{'':>9}"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ml-cut", type=float, default=None,
                   help="drop groups whose mass-to-light ratio exceeds the "
                        "running median at their mass by more than this many "
                        "dex.  1.0 is the recommended value; it removes 5 of "
                        "1833 groups and fixes the 14.2 bin.  Applied to BOTH "
                        "catalogues so the comparison stays fair.")
    p.add_argument("--out", default="hmf_old_vs_new.png")
    p.add_argument("--mass-err", default="driver",
                   help="mass-error curve: 'driver' (his hardcoded arrays, the "
                        "default) or a vuvuzela CSV such as gama_masserr.csv, "
                        "which also replaces MASSCORR with the measured "
                        "multiplicity debiasing.")
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc-edb", type=int, default=1001)
    p.add_argument("--nmc-fit", type=int, default=2001)
    p.add_argument("--nboot", type=int, default=2000,
                   help="group-bootstrap errors (0 = Driver's Poisson term)")
    p.add_argument("--fit-max", type=float, default=15.5,
                   help="upper edge of the fitted mass range. 15.5 matches the "
                        "range Driver's catalogue reaches, so the two are "
                        "compared like with like; use 99 to disable.")
    p.add_argument("--zmax-control", action="store_true",
                   help="also run the old data with zmax recomputed the new way")
    args = p.parse_args()

    mrpx, mrpy, factor, phimrp = lcdm_curve()
    results = []

    print("OLD GAMA (G3C v10, 3 equatorial fields, 179.92 deg^2, r<19.8)")
    g_old, v_old = dr.build_groups("../data/G3CFoFGroupv10.fits",
                                   "../data/GAMAGalsInGroups.csv")
    g_old, _ = apply_ml_cut(g_old, args.ml_cut)
    results.append(run_one(g_old, v_old, phimrp, args.seed, args.nmc_edb,
                           args.nmc_fit, "Driver+22, old GAMA (v10)",
                           "old (Driver+22)", OLD_C, "o",
                           nboot=args.nboot, mmax=args.fit_max))

    if args.zmax_control:
        print("\nCONTROL: old GAMA, zmax recomputed with the k+e polynomial")
        g_c, v_c = build_groups_old_recomputed_zmax()
        results.append(run_one(g_c, v_c, phimrp, args.seed, args.nmc_edb, 0,
                               "old GAMA, recomputed zmax", "old (poly zmax)",
                               "#888888", "^", nboot=args.nboot,
                               mmax=args.fit_max))

    print(f"\nNEW GAMA (Nessie DMU, 4 regions incl. G23, {NEW_AREA} deg^2, "
          f"r<{NEW_MAGLIM})")
    g_new, v_new = build_groups_new(masserr=args.mass_err)
    g_new, flagged_ml = apply_ml_cut(g_new, args.ml_cut)
    if flagged_ml is not None and len(flagged_ml):
        flagged_ml.to_csv("flagged_ml_outliers.csv", index=False)
        print(f"  wrote flagged_ml_outliers.csv ({len(flagged_ml)} groups)")
    results.append(run_one(g_new, v_new, phimrp, args.seed, args.nmc_edb,
                           args.nmc_fit, "This work, new GAMA DMU",
                           "new (Nessie DMU)", NEW_C, "s",
                           nboot=args.nboot, mmax=args.fit_max))

    print_binned_comparison(results)
    report(results)
    plot_comparison(results, mrpx, mrpy, factor, args.out)


if __name__ == "__main__":
    main()
