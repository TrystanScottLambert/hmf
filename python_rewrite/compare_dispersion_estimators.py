#!/usr/bin/env python3
"""NFW mass from the gapper dispersion vs from Tempel's rms dispersion.

A referee will reasonably ask how much of the Nessie/Tempel mass difference is
the choice of velocity-dispersion estimator rather than anything physical.  This
isolates exactly that: **same groups, same members, same sigma_sky, same NFW
kappa(M)** -- the only thing that changes is how sigma_v is measured.

* Tempel+2014 eq. 3 is the plain rms,
  ``sigma_v^2 = 1/[(1+z_m)^2 (n-1)] sum (v_i - v_mean)^2``.  Recomputing it from
  the members reproduces his published ``col12`` to +0.0003 dex, so the
  identification is certain.
* Nessie uses the **gapper** estimator, which is more robust for small N but is
  not the same statistic.

The velocity-error term is deliberately shown to be irrelevant: subtracting
Nessie's constant ``sigma_err_squared = 50`` changes sigma by 0.0001 dex,
because 50 (km/s)^2 against a typical sigma^2 of 40 000 is 0.06%.

Only groups whose membership is **byte-identical** in both catalogues are used,
so the group finder cannot contribute.

Usage
-----
    python compare_dispersion_estimators.py
    python compare_dispersion_estimators.py --ngroups 6000
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

import nessie_sdss_hmf as nsd
import vuvuzela as V

C_KMS = 299792.458


def rms_dispersion(zs):
    """Tempel+2014 eq. 3 -- the plain rms, in km/s."""
    n = len(zs)
    zm = np.mean(zs)
    v = np.asarray(zs, float) * C_KMS
    return np.sqrt(np.sum((v - v.mean()) ** 2) / (n - 1)) / (1.0 + zm)


def nfw_mass(sigma_v1d, sigma_sky, iters=40):
    """Tempel eq. 8 with the NFW kappa, iterated (kappa depends on the mass).

    ``sigma_v = sqrt(3) sigma_1D`` per section 4, hence the factor 3.
    """
    m = np.full(np.shape(sigma_v1d), 1e13, dtype=float)
    for _ in range(iters):
        with np.errstate(divide="ignore", invalid="ignore"):
            k = np.interp(np.log10(np.where(m > 0, m, np.nan)),
                          nsd.KAPPA_NFW_LOGM, nsd.KAPPA_NFW_VAL)
        m = 2.325e12 * k * sigma_sky * 3.0 * (sigma_v1d / 100.0) ** 2
    return m


def matched_groups(ngroups, verbose=True):
    """Groups whose member set is identical in both catalogues."""
    gal = pd.read_parquet(V.NESSIE_GALS)
    gal = gal[(gal.group_id != -1) & (gal.GroupID != -1)].copy()
    gal = V._absmag_by_join(gal, verbose=False) if False else gal
    gal["gk"] = (np.round(gal.RAJ2000, 6).astype(str) + "_"
                 + np.round(gal.DEJ2000, 6).astype(str))
    nes = gal.groupby("group_id")["gk"].apply(frozenset)
    tem = gal.groupby("GroupID")["gk"].apply(frozenset)
    tmap = {v: k for k, v in tem.items()}
    match = {g: tmap[s] for g, s in nes.items() if s in tmap}

    ng = pd.read_parquet(V.NESSIE_GROUPS).set_index("group_id")
    t2 = Table.read(V.TEMPEL_GROUPS)
    tg = pd.DataFrame({"idcl": np.asarray(t2["col1"], float),
                       "sig": np.asarray(t2["col12"], float)}).set_index("idcl")

    rows = []
    for g, t in match.items():
        if len(rows) >= ngroups:
            break
        if g not in ng.index or float(t) not in tg.index:
            continue
        m = gal[gal.group_id == g]
        if len(m) < 5:
            continue
        a = ng.loc[g]
        zs = m.zobs.values
        sd = V.sky_distribution(m.RAJ2000.values, m.DEJ2000.values, zs,
                                a.iter_ra, a.iter_dec)
        if not np.isfinite(sd) or sd <= 0:
            continue
        s_gap = V.velocity_dispersion_gapper(zs, 50.0)
        s_gap0 = V.velocity_dispersion_gapper(zs, 0.0)
        s_rms = rms_dispersion(zs)
        if min(s_gap, s_rms) <= 0:
            continue
        rows.append((len(m), s_gap, s_gap0, s_rms, sd, tg.loc[float(t)].sig))
    d = pd.DataFrame(rows, columns=["N", "s_gap", "s_gap0", "s_rms",
                                    "sky", "s_tempel"])
    if verbose:
        print(f"  identical-membership groups, N >= 5 : {len(d)}")
    return d


def plot(d, outfile):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    lm_r, lm_g = np.log10(d.m_rms), np.log10(d.m_gap)
    dif = lm_g - lm_r

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = GridSpec(2, 1, height_ratios=[2.6, 1.0], hspace=0.06,
                  left=0.115, right=0.845, top=0.965, bottom=0.082)
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)

    # multiplicity is a magnitude, so a sequential (perceptually uniform,
    # CVD-safe) ramp rather than categorical hues
    norm = matplotlib.colors.LogNorm(5, max(60, np.percentile(d.N, 99)))
    cmap = plt.get_cmap("viridis")

    pad = 0.12
    lo = float(min(lm_r.min(), lm_g.min())) - pad
    hi = float(max(lm_r.max(), lm_g.max())) + pad
    ax.plot([lo, hi], [lo, hi], color="0.25", lw=1.4, ls="--", zorder=4)
    sc = ax.scatter(lm_r, lm_g, c=d.N, cmap=cmap, norm=norm, s=9,
                    linewidths=0, alpha=0.65, zorder=3)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_ylabel(r"log$_{10}$ M$_{\rm NFW}$ from the **gapper** [M$_\odot$]"
                  .replace("**", ""), fontsize=10)
    ax.tick_params(direction="in", top=True, right=True, which="both",
                   labelsize=9, labelbottom=False)
    ax.minorticks_on()
    ax.grid(alpha=0.12, lw=0.6)

    med = np.median(dif)
    ax.text(0.032, 0.965,
            "Same groups, same members, same $\\sigma_{\\rm sky}$,\n"
            "same NFW $\\kappa$(M) — only $\\sigma_v$ differs",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            color="0.25")
    ax.text(0.975, 0.04,
            f"N = {len(d)} identical-membership groups\n"
            f"median offset = {med:+.3f} dex",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=9,
            color="0.25")

    cax = fig.add_axes([0.865, 0.082, 0.028, 0.883])
    cb = fig.colorbar(sc, cax=cax)
    # a galaxy count wants plain integers, not 6x10^0 / 2x10^1
    ticks = [t for t in (5, 7, 10, 15, 20, 30, 50, 100, 200)
             if norm.vmin <= t <= norm.vmax]
    cb.set_ticks(ticks)
    cb.ax.set_yticklabels([str(t) for t in ticks])
    cb.ax.minorticks_off()
    cb.set_label("multiplicity $N_{\\rm fof}$", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    axr.axhline(0, color="0.25", lw=1.4, ls="--", zorder=4)
    axr.scatter(lm_r, dif, c=d.N, cmap=cmap, norm=norm, s=7, linewidths=0,
                alpha=0.5, zorder=3)
    # running median, the only summary line that earns its place here
    edges = np.arange(12.5, 15.81, 0.25)
    cen, run = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (lm_r >= a) & (lm_r < b)
        if m.sum() > 20:
            cen.append(0.5 * (a + b))
            run.append(np.median(dif[m]))
    axr.plot(cen, run, color="crimson", lw=2.0, zorder=5,
             label="running median")
    lim = float(np.quantile(np.abs(dif - np.median(dif)), 0.995)) + 0.06
    axr.set_ylim(np.median(dif) - lim, np.median(dif) + lim)
    axr.set_xlabel(r"log$_{10}$ M$_{\rm NFW}$ from Tempel's rms (eq. 3) "
                   r"[M$_\odot$]", fontsize=10)
    axr.set_ylabel(r"$\Delta\,$log$_{10}$M", fontsize=10)
    axr.tick_params(direction="in", top=True, right=True, which="both",
                    labelsize=9)
    axr.minorticks_on()
    axr.grid(alpha=0.12, lw=0.6)
    axr.legend(fontsize=8.5, frameon=False, loc="upper right")

    fig.savefig(outfile, dpi=240)
    plt.close(fig)
    print(f"\n  wrote {outfile}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ngroups", type=int, default=4000)
    p.add_argument("--out", default="dispersion_estimator_comparison.pdf")
    p.add_argument("--csv", default="dispersion_estimator_comparison.csv")
    args = p.parse_args()

    print("NFW mass: gapper dispersion vs Tempel's rms (eq. 3)")
    d = matched_groups(args.ngroups)

    # validate the rms port against his published column before using it
    r = np.log10(d.s_rms / d.s_tempel)
    print(f"  eq. 3 rms vs his col12      : median {np.median(r):+.4f} dex "
          f"(scatter {r.std():.4f}) -- the identification check")
    rg = np.log10(d.s_gap / d.s_tempel)
    print(f"  gapper vs his col12         : median {np.median(rg):+.4f} dex")
    r0 = np.log10(d.s_gap0 / d.s_gap.replace(0, np.nan))
    print(f"  effect of the velocity-error term : "
          f"{np.nanmedian(r0):+.5f} dex -- negligible")

    d["m_gap"] = nfw_mass(d.s_gap.values, d.sky.values)
    d["m_rms"] = nfw_mass(d.s_rms.values, d.sky.values)
    d = d[(d.m_gap > 0) & (d.m_rms > 0)].copy()

    dif = np.log10(d.m_gap / d.m_rms)
    print(f"\n  log10(M_gapper / M_rms)     : median {np.median(dif):+.4f}   "
          f"16/84 {np.quantile(dif, .16):+.4f}/{np.quantile(dif, .84):+.4f}")
    print("\n  by multiplicity:")
    for lo, hi in [(5, 6), (7, 9), (10, 14), (15, 24), (25, 400)]:
        m = (d.N >= lo) & (d.N <= hi)
        if m.sum() < 10:
            continue
        print(f"    N {lo:3d}-{hi:<3d} n={m.sum():5d}   "
              f"{np.median(dif[m]):+.4f} dex")

    d.to_csv(args.csv, index=False)
    print(f"\n  wrote {args.csv}")
    plot(d, args.out)


if __name__ == "__main__":
    main()
