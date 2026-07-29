"""
============================================================
Old vs new GAMA catalogue: what actually changed?
============================================================

alpha moved from -1.806 +/- 0.079 (old catalogue, 179.92 deg^2, z-dependent
ramp) to -2.185 +/- 0.071 (new catalogue, 238.11 deg^2, global ramp) -- 0.38 dex,
about 5 sigma. Three things changed at once, so this pulls them apart.

Compares, per unit area:
  - number of groups, and their mass / multiplicity / redshift distributions
  - the turnover mlim(z) each catalogue produces
  - the completeness C(m,z) that mlim implies, which is what the fit sees

No Stan, no MCMC. Seconds.

Run:  python compare_catalogues.py OLD.fits NEW.fits [--old-area 179.92] [--new-area 238.11]
============================================================
"""

import argparse
import numpy as np
from scipy.stats import norm
import recovery as R


def summarise_cat(path, area, label, **load_kw):
    lm, sig, z, nfof = R.load_real_gama(path, **load_kw)
    mlim_func, _, kind, _ = R.turnover_mlim(z, lm)
    print(f"\n  === {label} ===")
    print(f"  file        {path}")
    print(f"  area        {area:.2f} deg^2")
    print(f"  N groups    {lm.size}   ({lm.size / area:.2f} per deg^2)")
    print(
        f"  mass        {lm.min():.2f} .. {lm.max():.2f}   median {np.median(lm):.2f}"
    )
    print(
        f"  Nfof        median {np.median(nfof):.0f}   "
        f"{(nfof >= 10).sum()} with >=10 members"
    )
    print(f"  sigma       median {np.median(sig):.3f}")
    print(f"  z           median {np.median(z):.3f}")
    print(
        f"  mlim(z)     [{kind}] {mlim_func(R.ZMIN):.2f} -> {mlim_func(R.ZLIMIT):.2f}"
    )
    # counts in fixed mass bins, per unit area -- where did groups appear/vanish?
    edges = np.arange(12.5, 15.5001, 0.25)
    cnt, _ = np.histogram(lm, bins=edges)
    return dict(
        lm=lm,
        sig=sig,
        z=z,
        nfof=nfof,
        mlim=mlim_func,
        kind=kind,
        area=area,
        edges=edges,
        cnt=cnt,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--old-area", type=float, default=179.92)
    ap.add_argument("--new-area", type=float, default=238.11)
    ap.add_argument(
        "--old-dec-cut",
        type=float,
        default=-3.5,
        help="the legacy Dec cut the old catalogue needs",
    )
    a = ap.parse_args()

    O = summarise_cat(a.old, a.old_area, "OLD", dec_cut=a.old_dec_cut)
    N = summarise_cat(a.new, a.new_area, "NEW")

    print("\n  === surface density by mass (groups per deg^2 per 0.25 dex) ===")
    print(f"  {'logM':>12} {'old':>8} {'new':>8} {'new/old':>8}")
    cen = 0.5 * (O["edges"][:-1] + O["edges"][1:])
    for i, c in enumerate(cen):
        o = O["cnt"][i] / O["area"]
        n = N["cnt"][i] / N["area"]
        r = f"{n / o:8.2f}" if o > 0 else "       -"
        print(f"  {c - 0.125:5.2f}-{c + 0.125:5.2f} {o:8.3f} {n:8.3f} {r}")
    print("  ratio ~1 everywhere -> same survey, more area.")
    print("  ratio <1 at low mass only -> the new catalogue is losing small groups")
    print("    (expected for a shallower limit; the question is how much).")

    print("\n  === what the fit actually sees: 50% completeness mass ===")
    print(f"  {'z':>6} {'old m50':>9} {'new m50':>9} {'shift':>7}")
    # old used the z-dependent r<19.8 ramp; new uses the global r<19.65 one
    old_d50 = np.interp(
        [0.05, 0.10, 0.15, 0.20], [0.045, 0.115, 0.20], [-0.148, -0.193, -0.232]
    )
    for k, zz in enumerate([0.05, 0.10, 0.15, 0.20]):
        m50_o = float(O["mlim"](zz)) + old_d50[k]
        m50_n = float(N["mlim"](zz)) + R.COMP_D50_PTS[0]
        print(f"  {zz:6.2f} {m50_o:9.2f} {m50_n:9.2f} {m50_n - m50_o:+7.2f}")
    print("  A positive shift means the model now thinks MORE low-mass groups are")
    print("  missing, so it inflates the faint end -> steeper alpha.")

    print("\n  === how much of that is mlim vs the ramp choice? ===")
    for zz in (0.05, 0.20):
        dm = float(N["mlim"](zz)) - float(O["mlim"](zz))
        dd = R.COMP_D50_PTS[0] - float(
            np.interp(zz, [0.045, 0.115, 0.20], [-0.148, -0.193, -0.232])
        )
        print(
            f"  z={zz:.2f}:  mlim contributes {dm:+.2f} dex, "
            f"D50 choice contributes {dd:+.2f} dex"
        )


if __name__ == "__main__":
    main()
