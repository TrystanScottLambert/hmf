#!/usr/bin/env python3
"""Per-leg mass-scale offsets, measured by matching integrals over 10^14-10^15.

The three legs of the combined fit are on different mass calibrations:

  GAMA (Nessie)   Robotham `A sigma^2 R50 / G`, A = 13.9 -- an empirical
                  calibration with no declared halo definition
  SDSS (Nessie)   Tempel eq. 8, NFW, M200c   <- the anchor
  SDSS (Tempel)   his col15, M200c
  REFLEX          M200c via the L_X-M relation (Bohringer sec 2)

`mass_definition_audit.py` measures the GAMA-vs-M200c offset directly, per
group, on identical memberships: **+0.543 dex**.  This module measures the same
thing a completely different way -- from the *binned* mass functions, by finding
the shift in log10 M that makes each leg's integral over 10^14-10^15 match the
anchor's.  Agreement between the two routes is the check that either means
anything.

**The shift is in MASS, not density.**  On a mass function this steep a mass
shift produces a density offset of `slope x shift`, which within 14.0-15.0 alone
varies by 0.84 dex -- so no single vertical rescale can represent it.

Usage
-----
    uv run python leg_offsets.py
    uv run python leg_offsets.py --lo 14.0 --hi 15.0 --anchor sdss_nessie
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.optimize import brentq

import combined_hmf as ch


def _integral(x, y, lo, hi, delta=0.0, n=400):
    """Int phi dlogM over [lo, hi] with the leg shifted by ``delta`` in log M.

    Shifting the DATA by +delta is the same as evaluating the unshifted data at
    logM - delta, so we interpolate the leg's log10 phi at (grid - delta).
    Linear interpolation in log-log, which is what the binned points support.
    """
    o = np.argsort(x)
    xs, ys = np.asarray(x)[o], np.asarray(y)[o]
    grid = np.linspace(lo, hi, n)
    lg = np.interp(grid - delta, xs, ys, left=np.nan, right=np.nan)
    if np.isnan(lg).any():
        return np.nan
    return np.trapezoid(10.0 ** lg, grid)


def offset_to_match(x, y, xa, ya, lo, hi, span=1.5):
    """Shift in log10 M that makes this leg's integral match the anchor's."""
    target = _integral(xa, ya, lo, hi)
    f = lambda d: _integral(x, y, lo, hi, d) - target
    a, b = -span, span
    fa, fb = f(a), f(b)
    if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
        # fall back to a scan when the leg does not span the shifted range
        ds = np.linspace(-span, span, 601)
        v = np.array([f(d) for d in ds])
        k = np.isfinite(v)
        if k.sum() < 2 or (v[k] > 0).all() or (v[k] < 0).all():
            return np.nan
        return float(ds[k][np.argmin(np.abs(v[k]))])
    return float(brentq(f, a, b, xtol=1e-4))


def legs(seed=10, nmc=1001, ml_cut=1.0, masserr="gama_masserr.csv",
         sdss_masserr="sdss_masserr.csv", verbose=False):
    """The four binned legs as (label -> (x, y))."""
    import sdss_hmf, nessie_sdss_hmf as nsd, new_gama_hmf as ng
    out = {}

    g, v = ng.build_groups_new(verbose=verbose, masserr=masserr)
    g, _ = ng.apply_ml_cut(g, ml_cut, verbose=False)
    b, _ = ch.bin_hmf(g, v, seed=seed, nmc=nmc, verbose=verbose)
    k = (b.gamay > 0) & np.isfinite(b.gamay)
    out["gama_nessie"] = (b.gamax[k], np.log10(b.gamay[k]))

    n, _, _ = nsd.build(seed=seed, nmc=nmc, verbose=verbose,
                        mass_mode="tempel_rms", masserr=sdss_masserr)
    n = n[(n.V2 > 0) & np.isfinite(n.V4)]
    out["sdss_nessie"] = (n.V1.values, n.V4.values)

    s, _, _ = sdss_hmf.build(seed=seed, nmc=nmc, verbose=verbose)
    s = s[(s.V2 > 0) & np.isfinite(s.V4)]
    out["sdss_tempel"] = (s.V1.values, s.V4.values)

    rx, ry, _ = ch.load_reflex()
    out["reflex"] = (rx, ry)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lo", type=float, default=14.0)
    p.add_argument("--hi", type=float, default=15.0)
    p.add_argument("--anchor", default="sdss_nessie")
    args = p.parse_args()

    d = legs()
    xa, ya = d[args.anchor]
    print(f"\n  Mass-scale offsets, integrals matched over "
          f"{args.lo} <= logM <= {args.hi}")
    print(f"  anchor: {args.anchor} (fixed at 0 by definition)\n")
    print(f"  {'leg':<16s} {'delta logM':>11s}   {'meaning':<s}")
    res = {}
    for lab, (x, y) in d.items():
        dd = 0.0 if lab == args.anchor else offset_to_match(x, y, xa, ya,
                                                            args.lo, args.hi)
        res[lab] = dd
        note = ("anchor" if lab == args.anchor else
                ("masses are HIGHER than the anchor" if dd > 0 else
                 "masses are LOWER than the anchor"))
        print(f"  {lab:<16s} {dd:>+11.3f}   {note}")

    print("\n  cross-check against mass_definition_audit.py")
    print("    that measures Robotham A=13.9 vs Tempel M200c per group,")
    print("    on identical memberships:            +0.543 dex")
    print(f"    this measures GAMA vs SDSS(Nessie), binned:  "
          f"{res.get('gama_nessie', float('nan')):+.3f} dex")
    print("    -> two independent routes; they should agree.")


if __name__ == "__main__":
    main()
