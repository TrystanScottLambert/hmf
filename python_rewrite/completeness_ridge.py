#!/usr/bin/env python3
"""Is `recovery.turnover_mlim`'s completeness line the detection boundary?

Plots halo mass against redshift for the Nessie GAMA groups and overlays

  * the per-z-bin MODE of the mass histogram -- the quantity `turnover_mlim`
    actually fits;
  * low percentiles of the mass distribution -- the empirical lower envelope,
    i.e. the ridge a detection limit ought to trace;
  * the fitted `mlim(z)`, in both the linear and quadratic forms.

The point of the figure is the gap between the second and third.  `mlim(z)` is
built from a histogram mode, so it lands where groups are most COMMON rather
than where they stop being detectable, and `recovery.py`'s `marg` model then
integrates every object from `mlim` upward -- discarding roughly half the
catalogue.  See CLAUDE.md, "recovery.py -- audit of the hierarchical Bayesian
route".

The 1/Vmax pipeline has no mass limit at all: each group carries its own zmax
from its 5th-brightest member, which is why it is immune to this.

Usage
-----
    uv run python completeness_ridge.py
"""

from __future__ import annotations

import argparse

import numpy as np


def per_bin(z, m, edges, hist_bw=0.3, min_in_bin=20):
    """(z centre, N, mode, p2, p5, p16) for each z bin that is populated."""
    me = np.arange(10.0, 16.0 + 1e-9, hist_bw)
    mid = 0.5 * (me[1:] + me[:-1])
    out = []
    for i in range(len(edges) - 1):
        s = (z >= edges[i]) & (z < edges[i + 1])
        if s.sum() <= min_in_bin:
            continue
        mm = m[s]
        mode = mid[np.argmax(np.histogram(mm, bins=me)[0])]
        p2, p5, p16 = np.percentile(mm, [2, 5, 16])
        out.append((0.5 * (edges[i] + edges[i + 1]), s.sum(), mode, p2, p5, p16))
    return np.array(out)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ml-cut", type=float, default=1.0)
    p.add_argument("--mass-err", default="gama_masserr.csv")
    p.add_argument("--nbin-z", type=int, default=24)
    p.add_argument("--out", default="completeness_ridge.png")
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import new_gama_hmf as ng
    import recovery as rc
    import plotting

    g, _ = ng.build_groups_new(verbose=False, masserr=args.mass_err)
    g, _ = ng.apply_ml_cut(g, args.ml_cut, verbose=False)
    z = g.Zfof.values
    m = np.log10(g.MassAfunc.values)
    print(f"  {len(z)} groups")

    fits = {}
    for form in ("linear", "quad"):
        f, c, kind, _ = rc.turnover_mlim(z, m, form=form)
        fits[form] = f
        print(f"  mlim {kind:20s}  below the line: {(m < f(z)).mean():.1%}")

    edges = np.linspace(z.min(), z.max(), args.nbin_z + 1)
    b = per_bin(z, m, edges)
    zc, nb, mode, p2, p5, p16 = b.T

    fig, (ax, axf) = plt.subplots(
        2, 1, figsize=(plotting.FULL_PAGE[0], 5.4), dpi=600, sharex=True,
        gridspec_kw=dict(height_ratios=[3.2, 1.0], hspace=0.06))

    zg = np.linspace(z.min(), z.max(), 200)
    ax.scatter(z, m, s=5, c="0.55", lw=0, alpha=0.55, zorder=1,
               label=f"Nessie GAMA groups ({len(z)})")
    ax.plot(zc, mode, "o", ms=5, mfc="white", mec="k", mew=1.2, zorder=5,
            label="mode of the mass histogram (what mlim fits)")
    ax.plot(zc, p16, "-", color="#2c7fb8", lw=1.8, zorder=4,
            label="16th percentile")
    ax.plot(zc, p5, "-", color="#41ab5d", lw=1.8, zorder=4,
            label="5th percentile  (the detection ridge)")
    ax.plot(zg, fits["linear"](zg), "-", color="#d62728", lw=2.4, zorder=6,
            label=r"$m_{\rm lim}(z)$, linear")
    ax.plot(zg, fits["quad"](zg), "--", color="#d62728", lw=2.0, zorder=6,
            label=r"$m_{\rm lim}(z)$, quadratic")
    ax.fill_between(zg, 10.4, fits["linear"](zg), color="#d62728", alpha=0.08,
                    lw=0, zorder=0)
    ax.text(0.985, 0.055, "groups here are discarded by the sharp-cut model",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="#d62728")
    ax.set_ylim(10.4, 16.0)
    ax.legend(loc="upper left", fontsize=7.5, frameon=False,
              handletextpad=0.6, labelspacing=0.32)
    plotting.prettify_axis(ax, "", r"log$_{10}$(Halo Mass / M$_\odot$)")

    frac = np.array([(m[(z >= edges[i]) & (z < edges[i + 1])]
                      < fits["linear"](0.5 * (edges[i] + edges[i + 1]))).mean()
                     for i in range(len(edges) - 1)
                     if ((z >= edges[i]) & (z < edges[i + 1])).sum() > 20])
    axf.plot(zc, 100 * frac, "-o", color="#d62728", ms=4, lw=1.6)
    axf.axhline(50, color="0.4", lw=1.0, ls=":")
    axf.set_ylim(0, 100)
    plotting.prettify_axis(axf, "Redshift",
                           "below\n" + r"$m_{\rm lim}$ [%]")
    plt.savefig(args.out, bbox_inches="tight", transparent=False)
    plt.close(fig)
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
