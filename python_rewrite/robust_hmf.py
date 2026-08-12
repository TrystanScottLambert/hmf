#!/usr/bin/env python3
"""Robust GAMA HMF: bootstrap errors, a Vmax floor, and the offending groups.

Driver's 1/Vmax estimator is dominated by a handful of nearby groups sitting
exactly at the multiplicity threshold.  Their Vmax is set by their *faintest*
member -- by construction, since a group leaves the sample when its 5th
brightest galaxy does -- and at low redshift that member is an intrinsically
feeble dwarf that disappears almost immediately.  The result is a group with a
Vmax of ~0.1 per cent of the survey volume carrying several hundred times the
typical weight.

Three things are done about it.

1. **Bootstrap errors** (``--nboot``).  Driver's Poisson term ``1/sqrt(N)``
   counts groups as if they contributed equally.  Resampling groups with
   replacement measures the true sampling variance of the weighted sum.  It
   reduces to 1/sqrt(N) for equal weights, so it is a strict generalisation.
   In practice it inflates the broken bins by 3-6x and leaves the healthy ones
   at 1.0-1.1x -- it finds the bad bins by itself.

2. **A Vmax floor** (``--floor``).  Driver already has one, ``vlimitmin =
   vlimit/1000`` (gamahmf.r line 229), set so permissively it never binds.
   Raising it caps the weight any one group can carry.  ``--scan`` shows the
   result as a function of the floor.

   *The floor is a biased estimator and the bias is mass dependent.*  It bites
   hardest in the low-mass bins, where nearby low-Vmax groups are legitimate --
   a 10^12.8 group with 5 members genuinely can only be found nearby.  Flooring
   flattens the faint-end slope: alpha moves from about -1.0 to -0.4 going from
   floor 0.001 to 0.02.  The bootstrap achieves the same protection with no
   bias, so the floor is best used as a robustness check rather than a default.
   Both are reported side by side.

3. **A flagged-group table** (``--flagged``), written to CSV for the paper.

Separately, and more seriously: the four-parameter MRP is **not constrained by
GAMA alone**.  Run to convergence the fit runs away to logM* ~ 10.5 (old) and
~6.3 (new) with alpha > 0; bounded, it rails at whatever lower limit is imposed.
Driver's published 13.51 is an artefact of stopping Nelder-Mead at
``maxit = 500`` (see ``driver_recovery.py``).  Both fits are reported here so the
difference is visible, but neither should be quoted as a measurement of M*
without external high-mass data (SDSS + REFLEX, as in ``allhmf.r``) or a prior.

Usage
-----
    python robust_hmf.py --scan       # floor scan
    python robust_hmf.py              # figure + flagged CSV + both error models
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

import driver_recovery as dr
import new_gama_hmf as ng
from driver_recovery import (LOGBIN, MASSX, MLIMIT, PUBLISHED_FIT, bin_hmf,
                             bin_index, fit_mrp, fit_selection, lcdm_curve,
                             make_massfn)

DEFAULT_FLOOR = 0.02          # from --scan; see module docstring
DRIVER_FLOOR = 0.001          # gamahmf.r line 229
NBOOT = 2000

# Physical bounds for the bounded fit.  M* is bounded below at the lowest fitted
# bin: below that the MRP is unconstrained by construction, because the whole
# fitted range then sits in the exponential tail.
BOUNDS = [(12.5, 15.5), (-6.0, -1.0), (-2.5, 0.0), (0.10, 1.50)]


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

def fit_bounded(allx, ally, allf, vlimit, seed=0, x0=None):
    """Global bounded fit of the MRP, in log10(phi).

    Driver's ``optim(maxit=500)`` is a truncated, non-converged walk whose answer
    is set by where the simplex happened to be when the budget ran out.  This is
    an actual minimum of the same objective within physical bounds, so it is
    reproducible and seed-stable -- but see the docstring: the minimum lies at
    the M* boundary, which is the real result.
    """
    fn = make_massfn(allx, ally, allf, vlimit)
    obj = lambda t: fn(np.array([t[0], 10 ** t[1], t[2], t[3]]))
    if x0 is None:
        r = differential_evolution(obj, BOUNDS, seed=seed, tol=1e-10,
                                   maxiter=4000, polish=True)
        t, v = r.x, r.fun
    else:
        r = minimize(obj, x0, method="L-BFGS-B", bounds=BOUNDS)
        t, v = r.x, r.fun
    return np.array([t[0], 10 ** t[1], t[2], t[3]]), v, t


def bounded_param_errors(b, vlimit, x0, nmc=2000, seed=0):
    """Parameter uncertainties by perturbing the binned points and refitting.

    Same perturbation as gamahmf.r lines 408-413 (fractional error plus cosmic
    variance), but refitting with the bounded minimiser started from the global
    solution, so the spread reflects the data rather than the optimiser.
    """
    rng = np.random.default_rng(seed)
    n = len(b.gamaf)
    out = []
    for _ in range(nmc):
        cv = b.gamay * rng.normal(0.0, b.cosvariance, n)
        mock = b.gamay + b.gamay * rng.normal(0.0, np.where(b.gamaf > 0, b.gamaf, 0), n)
        mock = mock + cv
        ax, ay, af = fit_selection(b, gamay=mock)
        if len(ax) < 5:
            continue
        par, _, _ = fit_bounded(ax, ay, af, vlimit, x0=x0)
        out.append([par[0], np.log10(par[1]), par[2], par[3]])
    return np.array(out)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def effective_n(logm, w, breaks):
    """N_eff = (sum w)^2 / sum w^2 per bin -- how many groups a bin is worth."""
    idx = bin_index(logm, breaks)
    nbin = len(breaks) - 1
    neff = np.zeros(nbin)
    raw = np.zeros(nbin, dtype=int)
    for b in range(nbin):
        ws = w[idx == b]
        raw[b] = len(ws)
        neff[b] = ws.sum() ** 2 / (ws ** 2).sum() if len(ws) and ws.sum() > 0 else 0.0
    return neff, raw


def flag_problem_groups(g3c, vlimit, floor_frac, bin_share_thresh=0.10,
                        min_bin_n=10):
    """The groups that break the estimator, by two independent criteria.

    * ``floored``  -- raw Vmax below the floor, so the weight had to be clipped.
    * ``dominant`` -- carries more than ``bin_share_thresh`` of the total 1/Vmax
      weight of its own mass bin, *and* that bin holds at least ``min_bin_n``
      groups.  The population requirement matters: without it the single most
      massive group in the survey is flagged merely for being alone in its bin,
      which is not a pathology but the measurement working as intended.

    ``Nfof`` is carried through because every one of these is a multiplicity-5
    group: its velocity dispersion, and hence its mass, rests on five redshifts.
    For the high-mass offenders that is the deeper problem -- a genuine 10^14.3
    halo at z = 0.02 would have far more than five GAMA members, so the mass is
    likely inflated by an interloper as well as the weight being extreme.
    """
    logm = np.log10(g3c.MassAfunc.values)
    vmax = g3c.vmax.values
    w = 1.0 / g3c.weightszlimit.values
    idx = bin_index(logm, MASSX)

    bin_total = np.array([w[idx == b].sum() for b in range(len(MASSX) - 1)])
    bin_n = np.array([int((idx == b).sum()) for b in range(len(MASSX) - 1)])
    with np.errstate(invalid="ignore", divide="ignore"):
        share = np.where(idx >= 0, w / np.where(bin_total[idx] > 0, bin_total[idx],
                                                np.nan), np.nan)
    n_in_bin = np.where(idx >= 0, bin_n[idx], 0)

    floored = vmax < floor_frac * vlimit
    dominant = (share > bin_share_thresh) & (n_in_bin >= min_bin_n)
    sel = floored | dominant

    out = pd.DataFrame({
        "GroupID": g3c.GroupID.values[sel].astype(int),
        "Nfof": g3c.Nfof.values[sel].astype(int),
        "Zfof": g3c.Zfof.values[sel],
        "VelDisp": g3c.VelDisp.values[sel],
        "Rad50": g3c.Rad50.values[sel],
        "log10M": logm[sel],
        "zmax": g3c.zmax.values[sel],
        "vmax_over_vlimit": vmax[sel] / vlimit,
        "weight_vs_median": (w / np.median(w))[sel],
        "mass_bin": np.where(idx[sel] >= 0,
                             MASSX[np.clip(idx[sel], 0, None)] + 0.5 * LOGBIN, np.nan),
        "bin_N": n_in_bin[sel],
        "bin_weight_share": share[sel],
        "floored": floored[sel],
        "dominant": dominant[sel],
    })
    if "GAMARegion" in g3c.columns:
        out.insert(1, "region", g3c.GAMARegion.values[sel])
    return out.sort_values("bin_weight_share", ascending=False).reset_index(drop=True)


def scan_floor(builder, phimrp, floors, seed, nboot, nmc_edb, label):
    print(f"\n  --- Vmax floor scan: {label} ---")
    print("  (fit columns are the bounded global fit; logM* rails at the bound "
          "in every case)")
    print(f"  {'floor':>7}{'Nfloor':>8}{'medNeff/N':>11}{'phi(14.2)':>11}"
          f"{'phi(12.8)':>11}{'alpha':>8}{'chi2':>9}")
    for f in floors:
        g3c, vlimit = builder(f)
        logm = np.log10(g3c.MassAfunc.values)
        w = 1.0 / g3c.weightszlimit.values
        neff, raw = effective_n(logm, w, MASSX)
        mids = MASSX[:-1] + 0.5 * LOGBIN
        good = (mids > MLIMIT) & (raw >= 20)
        b, _ = bin_hmf(g3c, vlimit, seed=seed, nmc=nmc_edb, verbose=False, nboot=nboot)
        par, val, _ = fit_bounded(*fit_selection(b), vlimit)
        with np.errstate(divide="ignore", invalid="ignore"):
            ly = np.log10(b.gamay)
        p142 = ly[int(np.argmin(np.abs(b.gamax - 14.2)))]
        p128 = ly[int(np.argmin(np.abs(b.gamax - 12.8)))]
        print(f"  {f:7.4f}{int((g3c.vmax.values < f * vlimit).sum()):8d}"
              f"{np.median((neff / np.maximum(raw, 1))[good]):11.2f}"
              f"{p142:11.3f}{p128:11.3f}{par[2]:8.3f}{val:9.2f}")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_final(results, mrpx, mrpy, factor, outfile, floor):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.4, 5.6))
    ax_top = fig.add_axes([0.115, 0.800, 0.875, 0.180])
    ax = fig.add_axes([0.115, 0.085, 0.875, 0.695])
    xfit = np.arange(11.0, 16.5, dr.FITBINWID)

    for r in results:
        ax_top.step(r["binned"].gamax - 0.5 * LOGBIN, r["binned"].raw_counts,
                    where="post", color=r["colour"], lw=1.1, label=r["short"])
    ax_top.set_xlim(12, 16)
    ax_top.set_ylabel("Number", fontsize=9)
    ax_top.set_xticklabels([])
    ax_top.tick_params(direction="in", top=True, right=True, labelsize=8)
    ax_top.legend(fontsize=7.5, frameon=False, loc="upper right", handlelength=1.4)

    with np.errstate(divide="ignore"):
        ax.plot(mrpx - 0.08, np.log10(mrpy) - np.log10(factor) + 0.08, ls="--",
                color="black", lw=1.8, zorder=3, label="LCDM prediction")
        ax.plot(xfit, np.log10(dr.mrp(xfit, PUBLISHED_FIT["logmstar"],
                10 ** PUBLISHED_FIT["logphistar"], PUBLISHED_FIT["alpha"],
                PUBLISHED_FIT["beta"])), color="0.45", lw=1.6, ls=(0, (5, 2)),
                zorder=3, label="Driver+22 published MRP")

    for r in results:
        b, c = r["binned"], r["colour"]
        with np.errstate(divide="ignore", invalid="ignore"):
            ly = np.log10(b.gamay)
            lo = np.abs(np.log10(np.clip(1 - b.gamaf, 1e-6, None)))
            hi = np.abs(np.log10(1 + b.gamaf))
        good = np.isfinite(ly)
        ax.errorbar(b.gamax[good], ly[good], yerr=[lo[good], hi[good]], fmt="none",
                    ecolor=c, elinewidth=0.9, capsize=1.5, alpha=0.9, zorder=4)
        hm = b.gamax > MLIMIT
        ax.plot(b.gamax[hm], ly[hm], ls="none", marker=r["marker"], color=c, ms=5,
                zorder=6)
        ax.plot(b.gamax[~hm], ly[~hm], ls="none", marker=r["marker"], mfc="none",
                mec=c, ms=5, zorder=6)
        with np.errstate(divide="ignore", invalid="ignore"):
            ax.plot(xfit, np.log10(dr.mrp(xfit, *r["par"])), color=c, lw=2,
                    zorder=7, label=r["label"])

    ax.set_xlim(12, 16)
    ax.set_ylim(-8, -2)
    ax.set_xlabel(r"log$_{10}$(Halo Mass / M$_\odot$)", fontsize=10)
    ax.set_ylabel(r"log$_{10}$(number density) [Mpc$^{-3}$ dex$^{-1}$]", fontsize=10)
    ax.tick_params(direction="in", top=True, right=True, which="both", labelsize=9)
    ax.minorticks_on()
    ax.legend(fontsize=8, frameon=False, loc="lower left", handlelength=2.0)
    ax.text(0.985, 0.955, f"bootstrap errors, Vmax floor = {floor}"
            r"$\,V_{\rm lim}$", transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color="0.35")
    fig.savefig(outfile, dpi=240)
    plt.close(fig)
    print(f"\n  wrote {outfile}")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--floor", type=float, default=DEFAULT_FLOOR)
    p.add_argument("--nboot", type=int, default=NBOOT)
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc-edb", type=int, default=1001)
    p.add_argument("--nmc-fit", type=int, default=2000)
    p.add_argument("--scan", action="store_true")
    p.add_argument("--out", default="hmf_robust.pdf")
    p.add_argument("--flagged", default="flagged_groups.csv")
    args = p.parse_args()

    mrpx, mrpy, factor, phimrp = lcdm_curve()
    old_b = lambda f: dr.build_groups("../data/G3CFoFGroupv10.fits",
                                      "../data/GAMAGalsInGroups.csv",
                                      verbose=False, vmax_floor_frac=f)
    new_b = lambda f: ng.build_groups_new(verbose=False, vmax_floor_frac=f)

    if args.scan:
        floors = [0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.12]
        scan_floor(new_b, phimrp, floors, args.seed, args.nboot, args.nmc_edb,
                   "new GAMA DMU")
        scan_floor(old_b, phimrp, floors, args.seed, args.nboot, args.nmc_edb,
                   "old GAMA (Driver+22)")
        return

    cats = [("old (Driver+22)", old_b, ng.OLD_C, "o"),
            ("new (Nessie DMU)", new_b, ng.NEW_C, "s")]

    # --- error-model and floor diagnostics -----------------------------------
    print("Bootstrap vs Poisson fractional error, and the effect of the floor")
    for short, build, _, _ in cats:
        g0, v0 = build(DRIVER_FLOOR)
        bp, _ = bin_hmf(g0, v0, seed=args.seed, nmc=args.nmc_edb, verbose=False)
        bb, _ = bin_hmf(g0, v0, seed=args.seed, nmc=args.nmc_edb, verbose=False,
                        nboot=args.nboot)
        neff, raw = effective_n(np.log10(g0.MassAfunc.values),
                                1.0 / g0.weightszlimit.values, MASSX)
        print(f"\n  --- {short}")
        print(f"  {'logM':>6}{'N':>6}{'Neff/N':>9}{'1/sqrtN':>10}{'boot':>8}"
              f"{'inflation':>11}")
        for i in range(len(bp.gamax))[::-1]:
            if bp.gamax[i] < 12.7 or bp.gamax[i] > 15.5 or bp.raw_counts[i] == 0:
                continue
            print(f"  {bp.gamax[i]:6.1f}{bp.raw_counts[i]:6.0f}"
                  f"{neff[i] / max(raw[i], 1):9.2f}{bp.rootnerr[i]:10.3f}"
                  f"{bb.rootnerr[i]:8.3f}{bb.rootnerr[i] / bp.rootnerr[i]:10.1f}x")

    # --- flagged groups ------------------------------------------------------
    frames = []
    for short, build, _, _ in cats:
        g0, v0 = build(DRIVER_FLOOR)
        f = flag_problem_groups(g0, v0, args.floor)
        f.insert(0, "catalogue", "new_DMU" if "new" in short else "old_v10")
        if "region" not in f.columns:
            f.insert(2, "region", "equatorial")
        frames.append(f)
    flagged = pd.concat(frames, ignore_index=True)
    flagged.to_csv(args.flagged, index=False)
    print(f"\n  Flagged groups -> {args.flagged}  "
          f"({len(frames[1])} new, {len(frames[0])} old)")
    cols = ["GroupID", "region", "Nfof", "Zfof", "VelDisp", "log10M", "zmax",
            "vmax_over_vlimit", "weight_vs_median", "mass_bin", "bin_N",
            "bin_weight_share", "floored", "dominant"]
    fmt = {"Zfof": "{:.4f}".format, "VelDisp": "{:.0f}".format,
           "log10M": "{:.2f}".format, "zmax": "{:.4f}".format,
           "vmax_over_vlimit": "{:.2e}".format, "weight_vs_median": "{:.0f}".format,
           "mass_bin": "{:.1f}".format, "bin_weight_share": "{:.3f}".format}
    print("\n  New catalogue -- groups to flag in the paper (top 12 by bin share):")
    print(frames[1][cols].head(12).to_string(index=False, formatters=fmt))
    print("\n  Old catalogue (Driver+22) -- the same pathology is already there:")
    print(frames[0][cols].drop(columns=["region"]).head(8)
          .to_string(index=False, formatters=fmt))

    # --- fits, both error models and both floors -----------------------------
    print("\n" + "=" * 96)
    print("  MRP fits.  'maxit=500' is Driver's truncated walk; 'bounded' is a "
          "true minimum")
    print("  of the same objective within physical bounds.")
    print("=" * 96)
    print(f"  {'catalogue':<18}{'floor':>7}{'errors':>11}{'logM*':>9}"
          f"{'logphi*':>10}{'alpha':>8}{'beta':>7}{'chi2':>9}   note")
    results = []
    floor_list = sorted({DRIVER_FLOOR, args.floor})
    for short, build, colour, marker in cats:
        for floor in floor_list:
            g3c, vlimit = build(floor)
            b, _ = bin_hmf(g3c, vlimit, seed=args.seed, nmc=args.nmc_edb,
                           verbose=False, nboot=args.nboot)
            ax_, ay_, af_ = fit_selection(b)
            p5, v5, _, _ = fit_mrp(ax_, ay_, af_, vlimit, phimrp)
            pb, vb, t0 = fit_bounded(ax_, ay_, af_, vlimit)
            for nm, par, val, note in (("maxit=500", p5, v5, "not converged"),
                                       ("bounded", pb, vb, "M* at bound")):
                print(f"  {short:<18}{floor:>7.3f}{'boot':>11}{par[0]:9.3f}"
                      f"{np.log10(abs(par[1])):10.3f}{par[2]:8.3f}{par[3]:7.3f}"
                      f"{val:9.2f}   {note}")
            if floor == args.floor:
                errs = bounded_param_errors(b, vlimit, t0, nmc=args.nmc_fit,
                                            seed=args.seed)
                results.append(dict(binned=b, par=pb, value=vb, errs=errs,
                                    short=short, colour=colour, marker=marker,
                                    label=f"{short}, bounded MRP fit"))
        print()

    print("  Bounded-fit parameter spread (16-84 per cent of the perturbed refits):")
    for r in results:
        e = r["errs"]
        q = lambda k: (np.quantile(e[:, k], .16), np.quantile(e[:, k], .84))
        print(f"    {r['short']:<18} logM* {q(0)[0]:6.2f}-{q(0)[1]:5.2f}   "
              f"alpha {q(2)[0]:6.2f}-{q(2)[1]:5.2f}   beta {q(3)[0]:5.2f}-{q(3)[1]:4.2f}")

    print("\n  WARNING: logM* sits at the lower bound in every configuration. "
          "GAMA alone\n  does not constrain M*; do not quote it without SDSS + "
          "REFLEX or a prior.")

    plot_final(results, mrpx, mrpy, factor, args.out, args.floor)


if __name__ == "__main__":
    main()
