#!/usr/bin/env python3
"""A mock-free selection function for the hierarchical GAMA fit.

WHY
---
The hierarchical fit's completeness ``C(m, z)`` currently comes from running
Nessie on the Shark lightcone.  That mock does not reproduce GAMA's group
population: it under-produces N = 5-6 groups by **1.87x** (the table in
``measure_completeness_nessie.py``'s own header).  Because the discrepancy runs
monotonically with richness it runs with mass, and a mass-dependent completeness
error maps straight onto alpha -- measured at **0.35 dex, 3.5x the statistical
error**.  That is not fixable inside ``recovery.py``; it is a property of the
mock.

So derive the selection from the data instead.

THE IDEA
--------
GAMA's selection is exact and knowable: a group is in the catalogue iff it has
>= MULTI members brighter than the flux limit.  Equivalently, writing ``M5`` for
the **5th brightest member's absolute magnitude**, the group is detected at
redshift z iff

    M5 <= Mlim(z) = maglim - distmod(z) - ke(z)

There is no modelling in that -- it is the definition of the catalogue, and it
is the same construction ``new_gama_hmf.compute_zmax`` already uses per group.

So the completeness is

    C(m, z) = P(M5 <= Mlim(z) | m)

and all we need is the intrinsic distribution of M5 at fixed halo mass.

THE TRUNCATION, AND WHY A NAIVE FIT FAILS
-----------------------------------------
``P(M5 | m)`` cannot be read off the observed groups directly: a group only
appears in the catalogue **if** M5 <= Mlim(z), so the observed M5 distribution
is truncated, and increasingly so with redshift.  It shows up plainly in the
data -- median M5 runs -17.08 at z < 0.05 to -19.95 at z > 0.15.

The obvious fix, calibrating on a low-z volume-limited sample, does not work
either: at z < 0.05 the limit is only M_r ~ -16.3, and there are just 59 groups
in that volume.

But the truncation is **exactly known** -- Mlim(z) is analytic -- so instead of
avoiding it, model it.  Fit

    M5 | m  ~  Normal(a + b (logM - 13.5),  s)

to ALL groups by maximum likelihood, with each group's contribution divided by
its own truncation probability ``Phi((Mlim(z_i) - mu_i) / s)``.  That uses all
1833 groups rather than 59, and the truncation being modelled is the very
selection we are trying to measure, so the calibration is self-consistent.

WHAT THIS BUYS AND WHAT IT COSTS
--------------------------------
Buys: no Shark, no HOD, no mock fidelity assumption.  The selection comes from
the same flux limit that defines the catalogue.

Costs: it assumes the M5-mass relation does not evolve over z < 0.25, and that a
Gaussian describes its scatter.  Both are checkable and neither involves galaxy
formation physics.  It also cannot supply ``NESSIE_BIAS`` -- but the closed loop
that produced that number can be replaced by a synthetic one built from this
same relation (``--closed-loop``), which needs no mock either.

Usage
-----
    python selfcal.py                 # fit the relation, report, write the table
    python selfcal.py --closed-loop   # synthetic injection-recovery of the relation
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.optimize import minimize
from scipy.special import erf, log_ndtr, ndtr

import new_gama_hmf as ng
import recovery as R

MULTI = 5
MPIV = 13.5  # pivot for the M5-mass relation, near the median group mass
OUT = "selfcal_selection.npz"


def mlim_of_z(z, maglim=None):
    """The faintest absolute magnitude visible at redshift z.

    The exact inverse of ``new_gama_hmf.compute_zmax``: both are built on the
    same tabulated ``distmod(z) + ke(z)`` and the same DMU cosmology, so the
    selection here is identical to the one the 1/Vmax pipeline applies.
    """
    maglim = ng.NEW_MAGLIM if maglim is None else maglim
    zg, off = ng._appmag_offset_grid()
    return maglim - np.interp(np.asarray(z, float), zg, off)


def load_groups_with_m5(verbose=True):
    """Group table plus M5, the 5th-brightest member's absolute magnitude."""
    g, vlimit = ng.build_groups_new(verbose=False)
    gal = Table.read(ng.GALS_NEW).to_pandas()
    gal = gal[gal.GroupID.isin(set(g.GroupID.values))]
    gal = gal[np.isfinite(gal.AbsoluteMagR) & (gal.AbsoluteMagR > -40)
              & (gal.AbsoluteMagR < 0)]
    gal = gal.sort_values(["GroupID", "AbsoluteMagR"], kind="mergesort")
    gal["rk"] = gal.groupby("GroupID").cumcount()
    m5 = gal[gal.rk == MULTI - 1].set_index("GroupID").AbsoluteMagR

    g = g.set_index("GroupID").copy()
    g["M5"] = m5
    g = g[np.isfinite(g.M5)].copy()
    # Mass MUST be on recovery.py's definition, not new_gama_hmf's: the former
    # uses A_SCALE = 10 and H0 = 100 (so the (100/H0) factor is 1), the latter
    # MAGICA = 13.9 and HO = 67.37.  That is a constant 0.315 dex, and since the
    # fit evaluates C at ITS masses, a table keyed on the other definition would
    # be read 0.315 dex off.
    _ratio = np.log10((ng.MAGICA / R.A_SCALE) * ((100.0 / ng.HO) / (100.0 / R.H0)))
    g["logM"] = np.log10(g.MassAfunc.values) - _ratio
    g["Mlim"] = mlim_of_z(g.Zfof.values)
    if verbose:
        viol = int((g.M5 > g.Mlim + 1e-6).sum())
        print(f"  mass rebased to recovery.py convention: -{_ratio:.3f} dex "
              f"(A {ng.MAGICA:g}->{R.A_SCALE:g}, h {ng.HO:g}->{R.H0:g})")
        print(f"  groups with M5 defined     : {len(g)}")
        print(f"  M5 fainter than its own Mlim: {viol} "
              f"(should be ~0 -- it is the detection condition)")
    return g, vlimit


def fit_m5_relation(logM, M5, Mlim, verbose=True):
    """Truncated ML fit of  M5 | m ~ N(a + b (logM - MPIV), s).

    Each group is observed only because M5 <= Mlim(z), so its likelihood is the
    Gaussian density divided by the truncation mass Phi((Mlim - mu)/s).  Without
    that division the fit chases the truncation and the relation comes out far
    too bright at high z.
    """
    x = np.asarray(logM, float) - MPIV
    y = np.asarray(M5, float)
    L = np.asarray(Mlim, float)

    def nll(p):
        a, b, ls = p
        s = np.exp(ls)
        mu = a + b * x
        zsc = (y - mu) / s
        # log density minus log truncation probability
        logpdf = -0.5 * zsc**2 - ls - 0.5 * np.log(2 * np.pi)
        logcdf = log_ndtr((L - mu) / s)
        return float(-np.sum(logpdf - logcdf))

    p0 = np.array([np.median(y), 0.0, np.log(1.0)])
    res = minimize(nll, p0, method="Nelder-Mead",
                   options=dict(maxiter=20000, fatol=1e-8, xatol=1e-8))
    a, b, s = res.x[0], res.x[1], float(np.exp(res.x[2]))
    if verbose:
        print(f"  M5 | logM  = {a:.3f} {b:+.3f} (logM - {MPIV})   scatter {s:.3f} mag"
              f"   [ok={res.success}, nll={res.fun:.1f}, n={len(y)}]")
        # the naive (untruncated) fit, to show the size of the bias being removed
        A = np.vstack([np.ones_like(x), x]).T
        cn = np.linalg.lstsq(A, y, rcond=None)[0]
        print(f"  naive (no truncation term): {cn[0]:.3f} {cn[1]:+.3f} -- biased "
              f"bright by {cn[0] - a:+.2f} mag at the pivot")
    return a, b, s


def completeness(logM, z, par):
    """C(m, z) = P(M5 <= Mlim(z) | m). Pure data-derived selection."""
    a, b, s = par
    mu = a + b * (np.asarray(logM, float) - MPIV)
    return ndtr((mlim_of_z(z) - mu) / s)


def build_table(par, m_grid=None, z_grid=None):
    m = np.arange(11.0, 15.81, 0.05) if m_grid is None else m_grid
    zz = np.linspace(R.ZMIN, R.ZLIMIT, 40) if z_grid is None else z_grid
    C = np.array([completeness(m, np.full(m.size, zi), par) for zi in zz])
    return m, zz, C


def closed_loop(par, vlimit, n_real=8, seed=0, truth=None):
    """Synthetic injection-recovery using ONLY the fitted relation.

    Draws masses from a known MRP over the survey volume, gives each object an
    M5 from the fitted relation, applies the exact flux selection, and checks
    that the resulting sample's mass distribution is the MRP times C.  This is
    the mock-free replacement for the Shark closed loop: it validates the
    selection model without assuming anything about galaxy formation.
    """
    rng = np.random.default_rng(seed)
    a, b, s = par
    truth = truth or (R.TRUE["ms"], R.TRUE["lp"], R.TRUE["al"], R.TRUE["be"])
    ms, lp, al, be = truth
    mg = np.arange(11.0, 16.0, 0.01)
    dens = R.mrp_phi(mg, ms, lp, al, be) * vlimit * 0.01
    print(f"  injected MRP: ms={ms:.3f} lp={lp:.3f} al={al:.3f} be={be:.3f}")
    print(f"  expected halos in the volume above logM 11: {dens.sum():.0f}")

    out = []
    for it in range(n_real):
        n = rng.poisson(dens)
        mass = np.repeat(mg, n)
        if mass.size == 0:
            continue
        # place uniformly in comoving volume
        u = rng.random(mass.size)
        zed = np.interp(u, np.linspace(0, 1, 400),
                        np.linspace(R.ZMIN, R.ZLIMIT, 400) ** 1.0)
        d = R.comoving_distance(np.linspace(R.ZMIN, R.ZLIMIT, 400))
        cdf = d**3 - d[0] ** 3
        cdf = cdf / cdf[-1]
        zed = np.interp(u, cdf, np.linspace(R.ZMIN, R.ZLIMIT, 400))
        m5 = a + b * (mass - MPIV) + s * rng.standard_normal(mass.size)
        det = m5 <= mlim_of_z(zed)
        out.append((mass, zed, det))
        if it == 0:
            print(f"  realisation 1: {mass.size} halos -> {int(det.sum())} detected "
                  f"({100 * det.mean():.1f}%)")

    # does the detected mass distribution equal MRP x C, bin by bin?
    edges = np.arange(12.0, 15.5, 0.25)
    cen = 0.5 * (edges[1:] + edges[:-1])
    obs = np.zeros(cen.size)
    for mass, zed, det in out:
        h, _ = np.histogram(mass[det], bins=edges)
        obs += h
    obs /= len(out)

    # volume-weighted mean completeness per mass bin
    zs = np.linspace(R.ZMIN, R.ZLIMIT, 60)
    w = np.gradient(R.comoving_distance(zs) ** 3)
    w = w / w.sum()
    pred = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mm = np.arange(lo, hi, 0.01)
        Cbar = np.array(
            [np.sum(w * completeness(np.full(zs.size, x), zs, par)) for x in mm]
        )
        pred.append(float(np.sum(R.mrp_phi(mm, ms, lp, al, be) * Cbar * vlimit * 0.01)))
    pred = np.array(pred)
    print()
    print("  logM     detected   MRPxC pred   ratio")
    for c, o, pv in zip(cen, obs, pred):
        if o > 0 or pv > 1:
            rr = f"{o / pv:8.3f}" if pv > 0.5 else "       -"
            print(f"  {c:5.2f} {o:10.1f} {pv:12.1f} {rr}")
    ok = (pred > 5)
    if ok.sum():
        r = obs[ok] / pred[ok]
        print(f"\n  median ratio over well-populated bins: {np.median(r):.4f} "
              f"(1.000 = the selection model is self-consistent)")
    return obs, pred


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--closed-loop", action="store_true")
    p.add_argument("--out", default=OUT)
    args = p.parse_args()

    print("Mock-free selection function for GAMA (selfcal)")
    g, vlimit = load_groups_with_m5()
    par = fit_m5_relation(g.logM.values, g.M5.values, g.Mlim.values)

    print()
    print("  C(m, z) implied by the fitted relation:")
    print("     logM " + "".join(f"{z:8.3f}" for z in (0.02, 0.05, 0.10, 0.15, 0.25)))
    for lm in (12.5, 13.0, 13.5, 14.0, 14.5):
        row = "".join(
            f"{completeness(lm, np.array([z]), par)[0]:8.3f}"
            for z in (0.02, 0.05, 0.10, 0.15, 0.25))
        print(f"  {lm:7.1f}" + row)

    m, zz, C = build_table(par)
    np.savez(args.out, m=m, z=zz, C=C, a=par[0], b=par[1], s=par[2],
             maglim=ng.NEW_MAGLIM, multi=MULTI, mpiv=MPIV, vlimit=vlimit)
    print(f"\n  wrote {args.out}  (m, z, C + the relation)")

    if args.closed_loop:
        print("\n  === synthetic closed loop (no mock) ===")
        closed_loop(par, vlimit)


if __name__ == "__main__":
    main()
