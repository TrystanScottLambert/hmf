#!/usr/bin/env python3
"""The Nessie SDSS catalogue through Driver's ``sdsshmf.r`` method.

This is the SDSS analogue of ``new_gama_hmf.py``: same 1/Vmax + Eddington
recipe, same binning, same optimiser, new group catalogue.  Only the group
construction differs from ``sdss_hmf.py``; the binning, table and fit are
imported from it unchanged so the two legs cannot drift apart.

Three things had to be resolved before this could be built.  All three turned
out to be answerable from the data rather than by choosing a convention:

* **No k-correction needs inventing.**  The Nessie galaxy file is Tempel+14's
  table 1 -- all 584,447 rows match ``sdssdr10table1.fits`` on (RA, Dec) with
  *identical* ``rmag`` -- so the absolute magnitudes come across by join, and
  ``sdsshmf.r`` line 269's analytic ``dmax`` is used verbatim.  Nessie's
  ``zobs`` differs from Tempel's ``col9`` by up to 0.0012 (frame convention);
  Tempel's redshift is used inside ``dmax`` to keep the formula self-consistent.
* **The two ID columns are two different catalogues.**  ``GroupID`` is Tempel's
  (max 88662, matching ``sdssdr10table2.fits``) and ``group_id`` is Nessie's
  (max 84890, matching ``sdss_groups.parquet``).  Both use ``-1`` for ungrouped.
  Only ``group_id`` joins the parquet pair; joining on ``GroupID`` would silently
  mix the two group finders.
* **The mass range is not the problem it looked like.**  Over the whole
  catalogue ``estimated_mass`` reaches 10^16.1, but that is driven by pairs and
  triples at high z.  Within the selection actually used (multiplicity >= 5,
  z < 0.08) it runs to 10^15.64 against Tempel's 10^15.08, with 19 groups above
  10^15 against Tempel's 4.  ``mass_proxy`` is *not* the analogue -- it sits
  ~0.9 dex low and is uncalibrated.  See ``--fit-max`` for the penalty-range
  check that group 300223 taught us to run.

Nessie's own cosmology is H0 = 70, Omega_m = 0.30 (recovered from its
``co_dist`` column to 4e-8), so masses carry the ``(70/ho)`` conversion that
mirrors Driver's ``(67.8/ho)`` for Tempel.

Note that ``sdsshmf.r``'s hand-fix of group ``idcl == 81455`` is *not* applied
here: that is a Tempel object, and the same lesson as GAMA's ``GroupID 100622``
applies -- the ID exists in the new catalogue but is a different group.

Usage
-----
    python nessie_sdss_hmf.py                  # binned table + fit
    python nessie_sdss_hmf.py --fit-max 15.0   # cap the penalty range
    python nessie_sdss_hmf.py --compare        # old vs new, side by side
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

import driver_recovery as dr
import sdss_hmf
from driver_recovery import co_dist, lcdm_curve, r_approx, r_optim_nm
from sdss_hmf import (LOGBIN, MAGLIM, MASSX, MLIMIT, MULTI, ZLIMIT, ZMIN,
                      _native, bin_hmf, survey_volume, table)

NESSIE_DIR = "/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS"
NESSIE_GROUPS = f"{NESSIE_DIR}/sdss_groups.parquet"
NESSIE_GALS = f"{NESSIE_DIR}/sdss_galaxies.parquet"

NESSIE_H = 70.0      # recovered from the co_dist column, with Omega_m = 0.30
UNGROUPED = -1       # sentinel in *both* GroupID and group_id


def _absmag_by_join(gal, gal_file=sdss_hmf.GALS, verbose=True):
    """Attach Tempel's ``col31`` absolute magnitude to the Nessie galaxies.

    Matched on (RA, Dec) rounded to 1e-6 deg.  This is exact, not approximate:
    the Nessie file is a row subset of Tempel's, so every key is present and the
    apparent magnitudes agree to 0.0.  The join is verified rather than assumed
    -- an unmatched row would mean the files had drifted apart, and is fatal.
    """
    t = Table.read(gal_file)
    tp = pd.DataFrame({
        "ra_t": _native(t["col13"]), "dec_t": _native(t["col14"]),
        "z_t": _native(t["col9"]), "rmag_t": _native(t["col26"]),
        "absmag_r": _native(t["col31"])})

    def key(ra, dec):
        return np.round(ra, 6).astype(str) + "_" + np.round(dec, 6).astype(str)

    tp["_k"] = key(tp.ra_t.values, tp.dec_t.values)
    tp = tp.drop_duplicates("_k")
    gal = gal.copy()
    gal["_k"] = key(gal.RAJ2000.values, gal.DEJ2000.values)
    out = gal.merge(tp, on="_k", how="left")

    n_bad = int(out.absmag_r.isna().sum())
    if n_bad:
        raise RuntimeError(
            f"{n_bad} Nessie galaxies have no Tempel counterpart; the two "
            "catalogues are no longer the same galaxy sample and the absolute "
            "magnitudes cannot be recovered by join.")
    dmag = float(np.abs(out.rmag - out.rmag_t).max())
    if verbose:
        print(f"  Tempel join                 : {len(out)}/{len(gal)} matched, "
              f"max |d rmag| = {dmag:.3g}")
    return out.drop(columns="_k")


def group_masses(grp, mass_mode="tempel_eq8", mass_shift=0.0):
    """The mass column, with the estimator made explicit.

    Nessie carries more than one estimator and they are **not** interchangeable:

    * ``tempel_eq8`` -- ``estimated_mass``, i.e. Tempel+2014 eq. 8 with
      ``grav_rad = 4.582 * sky_disp``.  This is the catalogue's own default and
      the analogue of Tempel's published masses, so it mirrors what Driver does
      for his SDSS leg (take the catalogue's native mass).
    * ``robotham`` -- ``mass_proxy * MAGICA``, i.e. the Robotham+2011 form
      ``A * R50 * sigma^2 / G`` with A = 13.9, then divided by ``10^masscorr``.
      This is **the estimator the GAMA leg uses** (``new_gama_hmf`` line 104),
      so it is the internally consistent choice when comparing the two Nessie
      legs to each other rather than to their parent catalogues.
    * ``shift`` -- ``tempel_eq8`` displaced by ``mass_shift`` dex, for testing
      the +0.285 dex offset against Tempel directly.

    Driver himself mixes estimators (Robotham for GAMA, Tempel's own for SDSS),
    so ``tempel_eq8`` is the faithful default; ``robotham`` is the consistent
    one.  Which is right depends on the question being asked.
    """
    if mass_mode == "tempel_eq8":
        m = grp.estimated_mass.astype(float).values
    elif mass_mode == "shift":
        m = grp.estimated_mass.astype(float).values * 10.0 ** mass_shift
    elif mass_mode == "robotham":
        m = grp.mass_proxy.astype(float).values * dr.MAGICA
        nf = grp.nrich.values.astype(int)
        mc = np.where(nf <= len(dr.MASSCORR),
                      dr.MASSCORR[np.clip(nf, 1, len(dr.MASSCORR)) - 1], np.nan)
        m = m / 10.0 ** np.where(np.isnan(mc), 0.0, mc)
    else:
        raise ValueError(f"unknown mass_mode {mass_mode!r}")
    return m * (NESSIE_H / sdss_hmf.HO)


def build_groups(group_file=NESSIE_GROUPS, gal_file=NESSIE_GALS, verbose=True,
                 mass_mode="tempel_eq8", mass_shift=0.0):
    """The Nessie analogue of ``sdss_hmf.build_groups`` (sdsshmf.r 263-307)."""
    gal = pd.read_parquet(gal_file)
    gal = gal[gal.group_id != UNGROUPED].copy()
    gal = _absmag_by_join(gal, verbose=verbose)

    # sdsshmf.r line 268 selects Tempel's precomputed rank == multi.  Tempel's
    # rank is by *absolute* magnitude ascending (verified exactly on his own
    # groups), so it is reproduced here inside each Nessie group.
    gal = gal.sort_values(["group_id", "absmag_r"], kind="mergesort")
    rank = gal.groupby("group_id").cumcount() + 1
    fifth = gal[rank == MULTI].copy()

    # line 269, verbatim.  Tempel's own redshift keeps the relation
    # self-consistent with the absolute magnitude it is paired with.
    fifth["dmax"] = 10 ** (0.2 * (MAGLIM - fifth.absmag_r - 25
                                  - 1.5 * np.log10(1 + fifth.z_t)))

    grp = pd.read_parquet(group_file)
    grp = grp[(grp.median_redshift < ZLIMIT) & (grp.median_redshift > ZMIN)
              & (grp.multiplicity > MULTI - 1)].copy().reset_index(drop=True)
    grp = grp.rename(columns={"group_id": "idcl", "multiplicity": "nrich",
                              "median_redshift": "zcl"})
    grp["mass"] = group_masses(grp, mass_mode, mass_shift)

    err = r_approx(grp.nrich.values.astype(float), dr.NFOF_XX, dr.NFOF_YY)
    err = np.where(np.isnan(err), 0.03, err)
    grp["log10MassErr"] = np.where(err < 0.1, 0.1, err)

    # lines 294-297: the same deliberately coarse dz = 0.001 grid
    tryz = np.arange(1, 1001) / 1000.0
    tryd = co_dist(tryz) * (1 + tryz)
    fifth["zmax"] = tryz[np.argmin(
        np.abs(fifth.dmax.values[:, None] - tryd[None, :]), axis=1)]

    grp["zmax"] = grp.idcl.map(dict(zip(fifth.group_id.values,
                                        fifth.zmax.values)))
    n_na = int(grp.zmax.isna().sum())

    grp["zmax"] = np.where(grp.zmax < grp.zcl, 1.1 * grp.zcl, grp.zmax)  # line 300
    grp["zmax"] = np.where(grp.zmax > ZLIMIT, ZLIMIT, grp.zmax)

    vmax = survey_volume(grp.zmax.values) - survey_volume(np.array([ZMIN]))[0]
    grp["vmax"] = vmax
    volumesdss = (survey_volume(np.array([ZLIMIT]))[0]
                  - survey_volume(np.array([ZMIN]))[0])
    volumesdssmin = volumesdss / 1000.0

    # lines 304-305, overwrite bug included; line 306's hand-fix deliberately not
    _ = np.where(vmax > volumesdss, volumesdss, vmax)
    grp["weightszlimit"] = np.where(vmax < volumesdssmin, volumesdssmin, vmax)

    if verbose:
        lm = np.log10(grp["mass"].values)
        print(f"  Nessie SDSS groups          : {len(grp)}")
        print(f"  rank-{MULTI} galaxies matched   : {len(grp) - n_na} "
              f"({n_na} unmatched)")
        print(f"  logM range                  : {lm.min():.2f} - {lm.max():.2f}"
              f"  (N > 15: {(lm > 15).sum()})")
        print(f"  volumesdss                  : {volumesdss:.6e} Mpc^3")
        print(f"  median vmax/volumesdss      : "
              f"{np.median(vmax) / volumesdss:.4f}")
    return grp, volumesdss


def fit(b, volumesdss, phimrp, maxit=500, fit_max=None):
    """``sdss_hmf.fit`` with an optional cap on the fitted mass range.

    ``fit_max`` exists because ``make_massfn`` integrates the penalty over
    ``max(allx) + 1..10 bins``: a single high-mass group therefore sets the
    range the penalty sees.  On GAMA that was worth 0.7-1.9 dex.
    """
    sel = (b["sdssy"] > 0) & (b["sdssx"] > MLIMIT)
    if fit_max is not None:
        sel &= b["sdssx"] <= fit_max
    with np.errstate(divide="ignore"):
        allx, ally, allf = (b["sdssx"][sel], np.log10(b["sdssy"][sel]),
                            b["sdssf"][sel])
    fn = sdss_hmf.make_massfn(allx, ally, allf, volumesdss)
    par, val, nfe, conv = r_optim_nm(
        np.array([dr.MSTARMRP, phimrp, dr.ALPHAMRP, dr.BETAMRP]), fn,
        maxit=maxit, reltol=1e-8, parscale=np.array([1.0, 1.0, 1.0, 0.1]))
    return par, val, nfe, conv, len(allx)


def build(seed=10, nmc=1001, verbose=True, nboot=0, mass_mode="tempel_eq8",
          mass_shift=0.0):
    grp, volumesdss = build_groups(verbose=verbose, mass_mode=mass_mode,
                                   mass_shift=mass_shift)
    b = bin_hmf(grp, volumesdss, seed=seed, nmc=nmc, verbose=verbose,
                nboot=nboot)
    return table(b), volumesdss, b


def _report(b, par, val, nfe, conv, n, label, fit_max):
    print(f"\n  --- {label} MRP fit ({n} bins with logM > {MLIMIT}"
          + (f", <= {fit_max}" if fit_max is not None else "") + ") ---")
    got = [par[0], np.log10(abs(par[1])), par[2], par[3]]
    for nm, v in zip(["log10(M*)", "log10(phi*)", "alpha", "beta"], got):
        print(f"  {nm:<13}{v:9.3f}")
    print(f"  chi2 = {val:.3f}   fevals = {nfe}   convergence = {conv}")
    return got


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc", type=int, default=1001)
    p.add_argument("--nboot", type=int, default=0,
                   help="bootstrap errors instead of Poisson")
    p.add_argument("--fit-max", type=float, default=None,
                   help="cap the fitted/penalty mass range, e.g. 15.0")
    p.add_argument("--mass-mode", default="tempel_eq8",
                   choices=["tempel_eq8", "robotham", "shift"],
                   help="which mass estimator; 'robotham' matches the GAMA leg")
    p.add_argument("--mass-shift", type=float, default=0.0,
                   help="dex shift, with --mass-mode shift")
    p.add_argument("--compare", action="store_true",
                   help="also run Driver's Tempel leg for a side-by-side")
    p.add_argument("--out", default="sdsshmfNessie5.csv")
    args = p.parse_args()

    _, _, _, phimrp = lcdm_curve()
    print("Nessie SDSS HMF (sdsshmf.r method, new group catalogue)")
    t, volumesdss, b = build(seed=args.seed, nmc=args.nmc, nboot=args.nboot,
                             mass_mode=args.mass_mode,
                             mass_shift=args.mass_shift)

    print("\n  logM      N   log10(wc)  log10(phi)     edb  rootnerr   mcerr"
          "   sdssf")
    with np.errstate(divide="ignore", invalid="ignore"):
        lw, ly = np.log10(b["wcounts"]), np.log10(b["sdssy"])
    for i in range(len(b["sdssx"]))[::-1]:
        if b["sdssx"][i] < MLIMIT - 0.5 or b["raw"][i] == 0:
            continue
        print(f"  {b['sdssx'][i]:5.2f} {b['raw'][i]:6.0f} {lw[i]:10.3f} "
              f"{ly[i]:11.3f} {b['edb'][i]:7.3f} {b['rootnerr'][i]:9.3f} "
              f"{b['mcerr'][i]:7.3f} {b['sdssf'][i]:7.3f}")

    par, val, nfe, conv, n = fit(b, volumesdss, phimrp, fit_max=args.fit_max)
    new = _report(b, par, val, nfe, conv, n, "Nessie SDSS", args.fit_max)

    if args.compare:
        print("\n" + "=" * 62)
        print("Driver's Tempel+14 leg, identical settings")
        og, ov = sdss_hmf.build_groups(verbose=False)
        ob = bin_hmf(og, ov, seed=args.seed, nmc=args.nmc, verbose=False,
                     nboot=args.nboot)
        opar, oval, onfe, oconv, on = sdss_hmf.fit(ob, ov, phimrp)
        old = _report(ob, opar, oval, onfe, oconv, on, "Tempel", None)
        print("\n  parameter        Tempel     Nessie      diff")
        for nm, o, nw in zip(["log10(M*)", "log10(phi*)", "alpha", "beta"],
                             old, new):
            print(f"  {nm:<13}{o:9.3f}{nw:11.3f}{nw - o:10.3f}")

    t.to_csv(args.out, index=False)
    print(f"\n  wrote {args.out}  ({len(t)} rows, Driver's V1..V8 order)")


if __name__ == "__main__":
    main()
