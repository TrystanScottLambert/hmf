#!/usr/bin/env python3
"""What halo mass definition does the GAMA leg's A = 13.9 actually correspond to?

The three legs of the combined fit sit at very different offsets from
Tinker+2008 M200c at Planck (REFLEX -0.18 dex, SDSS +0.26, GAMA +0.55), and the
halo mass DEFINITION alone moves Tinker by up to 0.34 dex (M200m vs M200c).  So
before any of those offsets can be interpreted, the GAMA leg needs a declared
definition -- and it does not have one:

  REFLEX  M200c  -- Bohringer+2017 sec 2, r200 wrt the CRITICAL density
  SDSS    M200c  -- Tempel's NFW col15; our kappa table is calibrated to it
  GAMA    ?      -- Robotham+2011 `A * sigma^2 * R50 / G` with A = 13.9, an
                    empirical calibration, not a geometric overdensity

The test uses the SDSS groups whose membership is **byte-identical** in the
Nessie and Tempel catalogues (the same set `vuvuzela.py --mass-audit` uses), so
the group finder is excluded by construction.  For each we compute the
GAMA-style Robotham mass from the same members and compare with Tempel's
`col15`, which is known M200c.

  flat offset  -> A = 13.9 is equivalent to some fixed overdensity, and the
                  offset says which
  sloped       -> the definition drifts with mass and the GAMA leg cannot be
                  placed on a single definition at all

Usage
-----
    uv run python mass_definition_audit.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

import driver_recovery as dr
import vuvuzela as vz


def identical_membership():
    """Groups with byte-identical membership in both SDSS catalogues."""
    gal = pd.read_parquet(vz.NESSIE_GALS)
    gal = gal[(gal.group_id != -1) & (gal.GroupID != -1)].copy()
    gal["gk"] = (np.round(gal.RAJ2000, 6).astype(str) + "_"
                 + np.round(gal.DEJ2000, 6).astype(str))
    nes = gal.groupby("group_id")["gk"].apply(frozenset)
    tem = gal.groupby("GroupID")["gk"].apply(frozenset)
    tmap = {v: k for k, v in tem.items()}
    return gal, {g: tmap[s] for g, s in nes.items() if s in tmap}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nmin", type=int, default=5)
    args = p.parse_args()

    gal, match = identical_membership()
    ng = pd.read_parquet(vz.NESSIE_GROUPS).set_index("group_id")
    t2 = Table.read(vz.TEMPEL_GROUPS)
    tg = pd.DataFrame({"idcl": np.asarray(t2["col1"], float),
                       "sig": np.asarray(t2["col12"], float),
                       "mass": np.asarray(t2["col15"], float)}).set_index("idcl")

    rows = []
    for g, t in match.items():
        if g not in ng.index or float(t) not in tg.index:
            continue
        a, b = ng.loc[g], tg.loc[float(t)]
        if a.multiplicity < args.nmin or b.sig <= 0 or b.mass <= 0 or a.r50 <= 0:
            continue
        # The GAMA leg's mass, built exactly as new_gama_hmf._masses_and_errors
        # does it: A * sigma^2 * R50 / G, with Nessie's own gapper sigma and r50.
        m_rob = (dr.MAGICA * (a.velocity_dispersion_gap * 1000.0) ** 2
                 * a.r50 * dr.PARSEC * 1e6 / (dr.G * dr.MSOL))
        rows.append((a.multiplicity, m_rob, b.mass * 1e12))
    d = pd.DataFrame(rows, columns=["N", "m_rob", "m_t"])
    d = d[(d.m_rob > 0) & (d.m_t > 0)]
    d["lg_t"] = np.log10(d.m_t)
    d["diff"] = np.log10(d.m_rob) - d.lg_t

    print(f"\n  Robotham A=13.9 mass vs Tempel col15 (M200c), identical membership")
    print(f"  {len(d)} groups with N >= {args.nmin}\n")
    print(f"  median log10(M_robotham / M200c) = {np.median(d['diff']):+.3f} dex"
          f"   16/84 {np.percentile(d['diff'],16):+.3f}/"
          f"{np.percentile(d['diff'],84):+.3f}")

    print("\n  is it flat?  (a fixed overdensity would give a constant offset)")
    print(f"  {'logM200c':>10s} {'n':>6s} {'median diff':>13s}")
    edges = [12.5, 13.0, 13.5, 14.0, 14.5, 15.5]
    cen, med = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (d.lg_t >= lo) & (d.lg_t < hi)
        if k.sum() < 20:
            continue
        cen.append(0.5 * (lo + hi)); med.append(np.median(d['diff'][k]))
        print(f"  {lo:4.1f}-{hi:4.1f} {k.sum():6d} {med[-1]:+13.3f}")
    if len(cen) > 2:
        slope = np.polyfit(cen, med, 1)[0]
        print(f"\n  slope = {slope:+.3f} dex per dex of mass"
              f"   (range {max(med)-min(med):.3f} dex across the bins)")
        print("  -> " + ("FLAT: consistent with a fixed overdensity"
                         if abs(slope) < 0.05 else
                         "SLOPED: the definition drifts with mass"))

    print("\n  also by multiplicity, since sigma's bias is N-dependent:")
    for lo, hi in ((5,6),(7,9),(10,14),(15,24),(25,10**6)):
        k = (d.N >= lo) & (d.N <= hi)
        if k.sum() < 20:
            continue
        print(f"    N {lo:3d}-{hi if hi < 10**5 else 999:<4d} {k.sum():6d}"
              f"   {np.median(d['diff'][k]):+.3f}")


if __name__ == "__main__":
    main()
