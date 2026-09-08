#!/usr/bin/env python3
"""Shape-only HMF: rescale each leg's y-axis to a common integral, then fit.

Each leg's density is multiplied by a constant so that

    Int_{10^14}^{10^15} phi dlog10M

matches the anchor leg's.  The absolute normalisation of every leg is then
deliberately discarded and the MRP is fitted to the rescaled points, so the fit
carries only the SHAPE (M*, alpha, beta) and a single common amplitude.

Rationale: the legs sit on different mass calibrations (GAMA's A = 13.9, the
SDSS NFW scale, REFLEX's L_X-M), and their absolute amplitudes are therefore not
mutually comparable.  Throwing the normalisations away and asking what shape the
data collectively prefer is a well-posed question that does not require the
calibrations to be reconciled first.

Caveat, stated because it bounds the interpretation: a mass-calibration offset
does NOT produce a constant vertical shift.  It produces `slope x offset`, which
over 14-15 alone varies by ~0.8 dex.  So this procedure removes the *mean*
disagreement between legs but not its mass dependence, and the residual scatter
after rescaling is a measure of how much that matters.  See `leg_offsets.py`
and `mass_definition_audit.py` for the mass-axis treatment.

The fractional errors are unchanged by a constant rescale.  The penalty term's
volumes ARE rescaled: data scaled by s correspond to an effective volume V/s,
since phi = N/(V dlogM).

Usage
-----
    uv run python normalised_shape_fit.py
    uv run python normalised_shape_fit.py --lo 14 --hi 15 --anchor sdss_nessie
"""

from __future__ import annotations

import argparse

import numpy as np

import combined_hmf as ch
import driver_recovery as dr

LEG_STYLE = {
    "gama_nessie": ("o", "red", "GAMA (Nessie)"),
    "sdss_nessie": ("P", "purple", r"SDSS (Nessie, NFW + rms $\sigma$)"),
    "sdss_tempel": ("s", "darkorange", "SDSS (Tempel+14)"),
    "reflex":      ("D", "forestgreen", "REFLEX II (Bohringer+17)"),
}


def integral(x, y, lo, hi, n=400):
    """Int phi dlog10M over [lo, hi] from binned log10 phi, linear interp."""
    o = np.argsort(x)
    xs, ys = np.asarray(x)[o], np.asarray(y)[o]
    grid = np.linspace(lo, hi, n)
    lg = np.interp(grid, xs, ys, left=np.nan, right=np.nan)
    if np.isnan(lg).any():
        return np.nan
    return np.trapezoid(10.0 ** lg, grid)


def _rescale_curve(par, legs, args, n=400):
    """Rescale an MRP's amplitude so its 10^lo-10^hi integral matches the anchor.

    Lets a MODEL be plotted on the same rescaled axis as the data: the shape is
    untouched and only phi* moves, which is exactly the operation applied to
    each leg.
    """
    grid = np.linspace(args.lo, args.hi, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        I = np.trapezoid(dr.mrp(grid, *par), grid)
    ia = integral(legs[args.anchor]["x"], legs[args.anchor]["y"],
                  args.lo, args.hi)
    out = np.array(par, dtype=float)
    out[1] = out[1] * (ia / I)
    return out


def build_legs(seed, nmc, ml_cut, masserr, sdss_masserr, fit_min, verbose):
    """label -> dict(x, y, f, vol) for every leg, on the fitted mass range."""
    import sdss_hmf, nessie_sdss_hmf as nsd, new_gama_hmf as ng
    out = {}

    g, v = ng.build_groups_new(verbose=verbose, masserr=masserr)
    g, _ = ng.apply_ml_cut(g, ml_cut, verbose=False)
    gset, _ = ch.gama_set(g, v, ng.NEW_AREA, seed, nmc, 0, ch.MLIMIT_GAMA,
                          fit_min=fit_min)
    out["gama_nessie"] = dict(x=np.asarray(gset[0]), y=np.asarray(gset[1]),
                              f=np.asarray(gset[2]), vol=gset[3])

    lo_s = max(ch.MLIMIT_SDSS, fit_min or ch.MLIMIT_SDSS)
    n, _, _ = nsd.build(seed=seed, nmc=nmc, verbose=verbose,
                        mass_mode="tempel_rms", masserr=sdss_masserr)
    n = n[(n.V2 > 0) & np.isfinite(n.V4) & (n.V1 > lo_s)]
    out["sdss_nessie"] = dict(x=n.V1.values, y=n.V4.values, f=n.V8.values,
                              vol=ch.survey_volume(ch.ZLIMITSDSS, ch.AREA_SDSS))

    s, _, _ = sdss_hmf.build(seed=seed, nmc=nmc, verbose=verbose)
    s = s[(s.V2 > 0) & np.isfinite(s.V4) & (s.V1 > lo_s)]
    out["sdss_tempel"] = dict(x=s.V1.values, y=s.V4.values, f=s.V8.values,
                              vol=ch.survey_volume(ch.ZLIMITSDSS, ch.AREA_SDSS))

    rx, ry, rf = ch.load_reflex()
    out["reflex"] = dict(x=rx, y=ry, f=rf, vol=ch.VOLUME_REFLEXII)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lo", type=float, default=14.0)
    p.add_argument("--hi", type=float, default=15.0)
    p.add_argument("--anchor", default="sdss_nessie")
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc-edb", type=int, default=1001)
    p.add_argument("--ml-cut", type=float, default=1.0)
    p.add_argument("--mass-err", default="gama_masserr.csv")
    p.add_argument("--sdss-mass-err", default="sdss_masserr.csv")
    p.add_argument("--fit-min", type=float, default=12.9)
    p.add_argument("--maxit", type=int, default=5000)
    p.add_argument("--fit-lo", type=float, default=None,
                   help="fit only points above this logM, ALL legs including "
                        "REFLEX.  The rescaling is only valid where the legs "
                        "were matched, so --fit-lo 14 restricts the fit to "
                        "where it holds -- at the cost of alpha, which is "
                        "measured from the low-mass end.")
    p.add_argument("--mcmc", action="store_true",
                   help="sample the posterior on the same rescaled objective "
                        "and write a corner plot")
    p.add_argument("--mcmc-steps", type=int, default=25000)
    p.add_argument("--from-chain", action="store_true",
                   help="redraw from the saved chain_<name>.npz instead of "
                        "re-sampling.  Restyling a figure is a redraw, not a "
                        "fresh 25k-step run.")
    p.add_argument("--mcmc-walkers", type=int, default=64)
    p.add_argument("--name", default="shape_normalised",
                   help="basename for hmf_<name>.png, corner_<name>.png and "
                        "chain_<name>.npz")
    p.add_argument("--out", default=None)
    args = p.parse_args()
    if args.out is None:
        args.out = f"hmf_{args.name}.png"

    print("Shape-only HMF: y-rescaled to a common integral over "
          f"10^{args.lo:g}-10^{args.hi:g}")
    legs = build_legs(args.seed, args.nmc_edb, args.ml_cut, args.mass_err,
                      args.sdss_mass_err, args.fit_min, verbose=False)

    ia = integral(legs[args.anchor]["x"], legs[args.anchor]["y"], args.lo, args.hi)
    print(f"\n  anchor = {args.anchor},  integral = {ia:.4e} Mpc^-3\n")
    print(f"  {'leg':<14s} {'integral':>12s} {'scale':>9s} {'dex':>8s}")
    for lab, d in legs.items():
        I = integral(d["x"], d["y"], args.lo, args.hi)
        d["scale"] = ia / I
        d["ynorm"] = d["y"] + np.log10(d["scale"])
        # phi = N/(V dlogM), so scaling phi by s is an effective volume V/s
        d["volnorm"] = d["vol"] / d["scale"]
        print(f"  {lab:<14s} {I:12.4e} {d['scale']:9.3f} "
              f"{np.log10(d['scale']):+8.3f}")

    # --- fit the MRP to the rescaled points -------------------------------
    order = ["gama_nessie", args.anchor if args.anchor != "sdss_nessie"
             else "sdss_nessie", "reflex"]
    order = [l for i, l in enumerate(order) if l not in order[:i]]
    def _cut(d):
        if args.fit_lo is None:
            return np.ones(len(d["x"]), bool)
        return d["x"] > args.fit_lo
    allx = np.concatenate([legs[l]["x"][_cut(legs[l])] for l in order])
    ally = np.concatenate([legs[l]["ynorm"][_cut(legs[l])] for l in order])
    allf = np.concatenate([legs[l]["f"][_cut(legs[l])] for l in order])
    allf = np.where(allf >= 1.0, 0.9999, allf)
    vols = np.array([legs[l]["volnorm"] for l in order])

    mrpx, mrpy, factor, phimrp = dr.lcdm_curve()
    fn = ch.make_massfn(allx, ally, allf, vols, False)
    par, val, nfe, conv = dr.r_optim_nm(
        np.array([dr.MSTARMRP, phimrp, dr.ALPHAMRP, dr.BETAMRP]), fn,
        maxit=args.maxit, reltol=1e-8, parscale=np.array([1.0, 1.0, 1.0, 0.1]))
    print(f"\n  fit to {len(allx)} rescaled points from {', '.join(order)}")
    for k, nm in enumerate(["log10(M*)", "log10(phi*)", "alpha", "beta"]):
        v = par[k] if k != 1 else np.log10(abs(par[1]))
        print(f"    {nm:<12s} {v:8.3f}")
    print(f"    chi2 = {val:.2f}   fevals = {nfe}   convergence = {conv}")

    par_plot, chain = par, None
    if args.mcmc or args.from_chain:
        import mcmc_hmf as mh
        nmv = np.array([par[0], np.log10(abs(par[1])), par[2], par[3]])
        # The MRP we EXPECT: Murray+21 LCDM, rescaled by the same integral
        # matching so it is comparable with the rescaled data.
        mur = np.array([dr.MSTARMRP, phimrp, dr.ALPHAMRP, dr.BETAMRP])
        murn = _rescale_curve(mur, legs, args)
        murv = np.array([murn[0], np.log10(murn[1]), murn[2], murn[3]])
        if args.from_chain:
            d = np.load(f"chain_{args.name}.npz", allow_pickle=True)
            chain, best = d["chain"], d["best"]
            print(f"\n  --- reusing chain_{args.name}.npz "
                  f"({len(chain)} samples), no re-sampling ---")
            mh.summarise(chain, "posterior, rescaled shape fit", nm=nmv)
        else:
            print(f"\n  --- MCMC on the same rescaled objective ---")
            chain, best, info = mh.run_emcee(fn, par, nwalkers=args.mcmc_walkers,
                                             nsteps=args.mcmc_steps,
                                             burn=args.mcmc_steps // 4,
                                             seed=args.seed)
            mh.summarise(chain, "posterior, rescaled shape fit", nm=nmv)
            if not info["converged"]:
                print("  WARNING: n_eff < 50 for at least one parameter")
            np.savez(f"chain_{args.name}.npz", chain=chain,
                     median=np.median(chain, axis=0), best=best, nm=nmv,
                     murray_rescaled=murv,
                     scales=np.array([legs[l]["scale"] for l in order]),
                     legs=np.array(order, dtype=object), allow_pickle=True)
        par_plot = np.array([best[0], 10.0 ** best[1], best[2], best[3]])
        med = np.median(chain, axis=0)
        markers = [
            (murv, "Murray+21 LCDM (rescaled)", "white"),
            (med, "This work (posterior median)", "#d4af37"),
        ]
        mh.corner_plot([(chain, "Rescaled shape fit", "#4878a8")],
                       f"corner_{args.name}.png", markers=markers,
                       title="MRP posterior, y-rescaled to a common integral")

    # residual scatter after rescaling, the honest diagnostic
    print("\n  residual scatter of each leg about the fit (rescaled):")
    xfit = np.arange(12.0, 16.01, 0.01)
    with np.errstate(divide="ignore", invalid="ignore"):
        yfit = np.log10(dr.mrp(xfit, *par))
    for lab, d in legs.items():
        r = d["ynorm"] - np.interp(d["x"], xfit, yfit)
        print(f"    {lab:<14s} median {np.median(r):+.3f}   rms {r.std():.3f} dex")

    plot(legs, par_plot, order, args, chain)


def plot(legs, par, order, args, chain=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotting

    fig = plotting.start_plot(r"log$_{10}$(Halo Mass / M$_\odot$)",
                              r"log$_{10}$(number density) [Mpc$^{-3}$ dex$^{-1}$]",
                              figsize=plotting.FULL_PAGE)
    ax = plt.gca()

    xfit = np.arange(12.0, 16.01, 0.01)

    # the MRP we EXPECT: Murray+21 LCDM, rescaled the same way
    murn = _rescale_curve(np.array([dr.MSTARMRP, dr.lcdm_curve()[3],
                                    dr.ALPHAMRP, dr.BETAMRP]), legs, args)
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(dr.mrp(xfit, *murn)), color="k", lw=1.8,
                ls=(0, (6, 3)), alpha=0.9, zorder=4,
                label="Murray+21 LCDM MRP, rescaled")

    if chain is not None:
        from matplotlib.collections import LineCollection
        rng = np.random.default_rng(0)
        idx = rng.choice(len(chain), size=min(600, len(chain)), replace=False)
        segs = []
        for c in chain[idx]:
            with np.errstate(divide="ignore", invalid="ignore"):
                yy = np.log10(dr.mrp(xfit, c[0], 10.0 ** c[1], c[2], c[3]))
            m = np.isfinite(yy) & (yy > -9)
            if m.sum() > 2:
                segs.append(np.column_stack([xfit[m], yy[m]]))
        ax.add_collection(LineCollection(segs, colors=["#d62728"],
                                         linewidths=0.6, alpha=0.02,
                                         rasterized=True, zorder=6))
    with np.errstate(divide="ignore", invalid="ignore"):
        ax.plot(xfit, np.log10(dr.mrp(xfit, *par)), color="#d62728", lw=2,
                ls=(0, (1, 2.2)), zorder=8,
                label="MRP fit to the rescaled data"
                      + (" (MCMC)" if chain is not None else ""))

    for lab, d in legs.items():
        mk, c, name = LEG_STYLE[lab]
        f = np.where(d["f"] >= 1.0, 0.9999, d["f"])
        ax.errorbar(d["x"], d["ynorm"],
                    yerr=[np.abs(np.log10(1 - f)), np.abs(np.log10(1 + f))],
                    fmt=mk, ms=4.5, color=c, mfc=c if lab in order else "none",
                    ecolor=c, elinewidth=0.8, capsize=0, ls="none", zorder=5,
                    label=name + ("" if lab in order else "  (not fitted)"))

    ax.axvspan(args.lo, args.hi, color="0.5", alpha=0.2, zorder=0, lw=0)
    ax.set_xlim(12.75, 16)
    ax.set_ylim(-8, -2)
    ax.legend(loc="upper right", fontsize=7.5, frameon=False,
              handletextpad=0.5, labelspacing=0.35)
    plotting.end_plot(args.out)
    plt.close(fig)
    print(f"\n  wrote {args.out}")


if __name__ == "__main__":
    main()
