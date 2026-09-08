#!/usr/bin/env python3
"""A real LCDM halo mass function: Tinker+2008 at Planck 2018.

Why this exists
---------------
``driver_recovery.lcdm_curve`` is the Murray+21 MRP with its amplitude set by
``factor`` so that integrating it to zero mass returns Omega_M = 0.3147 -- i.e.
so that *all* matter sits in haloes.  No halo mass function attains that.
Integrating Tinker+2008 at Planck 2018 down to 10^6 Msun gives Omega_M = 0.188,
60 per cent of Planck's value; the rest is in unbound material and haloes below
any sensible resolution.

The imposed normalisation therefore lifts the curve by
``log10(0.3147/0.188) = +0.223 dex``, against the +0.25 dex measured directly
against this module.  That is why the dashed reference line in the figures sat
well above where a real LCDM mass function lies, and why the REFLEX II points
appeared to fall far below "LCDM" once their own 0.638 dex normalisation error
was corrected.

``lcdm_curve`` is NOT changed: its fourth return value, ``phimrp``, is the
starting value the optimiser uses for phi, and moving it would perturb every
Nelder-Mead result including the bit-exact reproductions.  This module supplies
a curve for plotting only.

Definitions match what we compare against: **M200 with respect to the critical
density**, as in Bohringer+2017 sec 2 and in our ``tempel_nfw`` / ``tempel_rms``
SDSS masses.

Usage
-----
    uv run python lcdm_tinker.py --check     # regression + Omega_M integral
    from lcdm_tinker import tinker_curve
    lgm, lgphi = tinker_curve(z=0.1)
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

# Planck 2018 (TT,TE,EE+lowE+lensing), the cosmology used throughout the project
H0 = 67.37
OM0 = 0.3147
SIGMA8 = 0.811
NS = 0.9665
OB0H2 = 0.02237
TCMB = 2.7255

CACHE = "lcdm_tinker.csv"


MDEFS = {"200c": ("SOCritical", {"overdensity": 200}),
         "200m": ("SOMean", {"overdensity": 200}),
         "vir":  ("SOVirial", {})}


def _compute(H0, Om0, sigma_8, ns, z, lgm_min, lgm_max, dlog10m, mdef="200c"):
    """Tinker+2008 M200c via hmf + CAMB.  Slow (~1-2 min), hence the cache."""
    import astropy.units as u
    from astropy.cosmology import FlatLambdaCDM
    from hmf import MassFunction

    h = H0 / 100.0
    # Tcmb0 must be set explicitly or hmf's CAMB transfer model refuses to run.
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=OB0H2 / h ** 2, Tcmb0=TCMB * u.K)
    mf = MassFunction(cosmo_model=cosmo, sigma_8=sigma_8, n=ns,
                      hmf_model="Tinker08",
                      mdef_model=MDEFS[mdef][0], mdef_params=MDEFS[mdef][1],
                      Mmin=lgm_min + np.log10(h), Mmax=lgm_max + np.log10(h),
                      dlog10m=dlog10m, z=z)
    # hmf works in h-units throughout: M is Msun/h, dndlog10m is (Mpc/h)^-3.
    return np.log10(mf.m / h), np.log10(mf.dndlog10m * h ** 3)


def tinker_curve(z=0.1, H0=H0, Om0=OM0, sigma_8=SIGMA8, ns=NS,
                 lgm_min=6.0, lgm_max=16.5, dlog10m=0.02, cache=CACHE,
                 rebuild=False, verbose=False, mdef="200c"):
    """log10 M [Msun] and log10 dn/dlog10M [Mpc^-3 dex^-1] for Tinker+2008.

    Cached per (z, cosmology) in ``cache``; delete the file or pass
    ``rebuild=True`` to recompute.
    """
    key = f"{z:.4f}_{H0:.4f}_{Om0:.5f}_{sigma_8:.4f}_{ns:.5f}_{mdef}"
    if cache and os.path.exists(cache) and not rebuild:
        d = pd.read_csv(cache)
        sub = d[d.key == key]
        if len(sub):
            if verbose:
                print(f"  tinker_curve: {len(sub)} points from {cache}")
            return sub.lgm.values, sub.lgphi.values
    if verbose:
        print(f"  tinker_curve: computing (CAMB, ~1-2 min) for z={z} ...")
    lgm, lgphi = _compute(H0, Om0, sigma_8, ns, z, lgm_min, lgm_max, dlog10m,
                          mdef)
    if cache:
        new = pd.DataFrame({"key": key, "lgm": lgm, "lgphi": lgphi})
        if os.path.exists(cache):
            old = pd.read_csv(cache)
            new = pd.concat([old[old.key != key], new], ignore_index=True)
        new.to_csv(cache, index=False)
        if verbose:
            print(f"  tinker_curve: cached to {cache}")
    return lgm, lgphi


def omega_in_haloes(lgm, lgphi, above=None, H0=H0):
    """Integrate M phi dlogM / rho_crit -- the matter fraction in haloes."""
    import driver_recovery as dr
    k = np.ones(len(lgm), bool) if above is None else lgm > above
    rho_c = dr.RHOCRIT * (1e6 * dr.PARSEC) ** 3 / dr.MSOL
    return np.trapezoid(10.0 ** lgphi[k] * 10.0 ** lgm[k], lgm[k]) / rho_c


# Regression targets, measured 2026-09-04 with hmf+CAMB at z=0.
CHECK_Z0 = {13.5: -3.858, 14.0: -4.428, 14.5: -5.153, 15.0: -6.220}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true")
    p.add_argument("--rebuild", action="store_true")
    p.add_argument("--z", type=float, default=0.1)
    args = p.parse_args()

    if args.check:
        lgm, lgphi = tinker_curve(z=0.0, rebuild=args.rebuild, verbose=True)
        print("\n  z = 0 regression against the values that diagnosed the problem")
        print(f"  {'logM':>6s} {'this run':>10s} {'expected':>10s} {'diff':>8s}")
        ok = True
        for L, want in CHECK_Z0.items():
            got = float(np.interp(L, lgm, lgphi))
            ok &= abs(got - want) < 0.01
            print(f"  {L:6.1f} {got:10.3f} {want:10.3f} {got - want:+8.3f}")
        print(f"  -> {'PASS' if ok else 'FAIL'}")

        print("\n  matter fraction in Tinker haloes (Planck Omega_M = 0.3147):")
        for lo in (13, 12, 11, 10, 8, 6, None):
            om = omega_in_haloes(lgm, lgphi, above=lo)
            lab = f"logM > {lo}" if lo is not None else "all masses"
            print(f"    {lab:<12s} Omega_M = {om:.4f}  ({100*om/0.3147:5.1f}%)")
        tot = omega_in_haloes(lgm, lgphi)
        print(f"\n  driver_recovery.lcdm_curve forces this to 0.3147, i.e. "
              f"{np.log10(0.3147/tot):+.3f} dex too high.")
        return

    lgm, lgphi = tinker_curve(z=args.z, rebuild=args.rebuild, verbose=True)
    print(f"  z = {args.z}: {len(lgm)} points, logM {lgm.min():.2f} .. {lgm.max():.2f}")
    for L in (13.0, 13.5, 14.0, 14.5, 15.0, 15.5):
        print(f"    logM {L:5.2f}   log10 dn/dlogM = {np.interp(L, lgm, lgphi):+.3f}")


if __name__ == "__main__":
    main()
