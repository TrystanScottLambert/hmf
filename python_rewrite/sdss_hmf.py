#!/usr/bin/env python3
"""Python port of ``sdsshmf.r``: Driver's binned SDSS halo mass function.

This is the pipeline that produces ``sdsshmf5.csv``, the missing SDSS input to
``allhmf.r``.  It runs on the Tempel et al. (2014) SDSS DR10 group catalogue
(``sdssdr10table1.fits`` galaxies, ``sdssdr10table2.fits`` groups) and follows
the same 1/Vmax + Eddington-bias recipe as ``gamahmf.r``, with several
deliberate differences that are easy to miss:

* ``logbin = 0.1``, not 0.2 -- the SDSS HMF is binned twice as finely.
* ``mlimit = 12.9``, not 12.7.
* Masses are Tempel's own, ``mass * 1e12 * (67.8/ho)``, not rebuilt from a
  velocity dispersion.  There is no A factor and no multiplicity debiasing.
* ``rootnerr = sqrt(N)/N``.  Algebraically 1/sqrt(N), but at N = 0 it evaluates
  to 0/0 = NaN rather than Inf, so empty bins end up with ``sdssf = 0.9999``
  instead of 0.0.  The distinction only shows up in the written table.
* ``zmax`` comes from the *rank = 5* member's absolute magnitude via a
  luminosity-distance grid, and when it falls below the group redshift it is set
  to ``1.1 * zcl`` (``gamahmf.r`` uses ``Zfof`` itself).
* ``volumesdss`` here subtracts the ``zmin`` volume (line 248); ``allhmf.r``
  line 264 does not.  Both are reproduced where they belong.

Usage
-----
    python sdss_hmf.py                 # binned table + fit, writes sdsshmf5.csv
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

import driver_recovery as dr
from driver_recovery import (FULLSKY, MSOL, PARSEC, RRandom, co_dist, co_vol,
                             cosvar, lcdm_curve, r_approx, r_maghist,
                             r_optim_nm, r_weighted_hist)

HO = 67.37
OMEGAM = 0.3147
LOGBIN = 0.1                 # line 236 -- NOT 0.2
MULTI = 5                    # line 245
MLIMIT = 12.9                # line 246
ZMIN = 0.015                 # line 247
ZLIMIT = 0.08                # line 239
AREA = 7221.0                # line 248
MAGLIM = 17.77               # line 269, SDSS spectroscopic limit
TEMPEL_H = 67.8              # line 278
BAD_IDCL = 81455             # line 306, hand-fixed group
MASSX = np.round(np.arange(10.3, 16.1 + 1e-9, LOGBIN), 10)

GALS = "../data/sdssdr10table1.fits"
GROUPS = "../data/sdssdr10table2.fits"

# Column indices from lines 264-267 and 273-276 (R is 1-based, and astropy
# names the columns col1..colN, so the numbers carry over unchanged).
GAL_COLS = {"idcl": "col4", "rank": "col6", "redshift": "col9",
            "absmag_r": "col31"}
GRP_COLS = {"idcl": "col1", "nrich": "col2", "zcl": "col10", "mass": "col15"}


def _native(a):
    """FITS columns come back big-endian; pandas needs native byte order."""
    a = np.asarray(a)
    return a.astype(a.dtype.newbyteorder("=")) if a.dtype.byteorder == ">" else a


def survey_volume(z):
    return AREA / FULLSKY * 1e9 * co_vol(np.atleast_1d(z))


def build_groups(gal_file=GALS, group_file=GROUPS, verbose=True,
                 vmax_floor_frac=1e-3):
    """sdsshmf.r lines 263-307."""
    g = Table.read(gal_file)
    gal = pd.DataFrame({k: _native(g[c]) for k, c in GAL_COLS.items()})
    # line 268: only the rank-5 member matters -- it is the one whose
    # disappearance drops the group below the multiplicity cut
    gal = gal[(gal["rank"] == MULTI) & (gal.redshift < ZLIMIT + 0.05)].copy()
    # line 269: luminosity distance at which this galaxy reaches the limit
    gal["dmax"] = 10 ** (0.2 * (MAGLIM - gal.absmag_r - 25
                                - 1.5 * np.log10(1 + gal.redshift)))

    t = Table.read(group_file)
    grp = pd.DataFrame({k: _native(t[c]) for k, c in GRP_COLS.items()})
    grp = grp[(grp.zcl < ZLIMIT) & (grp.zcl > ZMIN)
              & (grp.nrich > MULTI - 1)].copy().reset_index(drop=True)
    grp["mass"] = grp["mass"].astype(float) * 1e12 * (TEMPEL_H / HO)  # line 278

    err = r_approx(grp.nrich.values.astype(float), dr.NFOF_XX, dr.NFOF_YY)
    err = np.where(np.isnan(err), 0.03, err)
    grp["log10MassErr"] = np.where(err < 0.1, 0.1, err)

    # lines 294-297: zmax by nearest point on a 1000-step luminosity-distance
    # grid.  The grid is coarse (dz = 0.001); it is reproduced rather than
    # inverted exactly, because the coarseness is part of the original result.
    tryz = np.arange(1, 1001) / 1000.0
    tryd = co_dist(tryz) * (1 + tryz)                    # flat: D_L = D_C (1+z)
    gal["zmax"] = tryz[np.argmin(np.abs(gal.dmax.values[:, None] - tryd[None, :]),
                                 axis=1)]

    zmax_by = dict(zip(gal.idcl.values, gal.zmax.values))
    grp["zmax"] = grp.idcl.map(zmax_by)
    n_na = int(grp.zmax.isna().sum())

    # line 300: note 1.1*zcl, not zcl as in gamahmf.r line 302
    grp["zmax"] = np.where(grp.zmax < grp.zcl, 1.1 * grp.zcl, grp.zmax)
    grp["zmax"] = np.where(grp.zmax > ZLIMIT, ZLIMIT, grp.zmax)

    vmax = survey_volume(grp.zmax.values) - survey_volume(np.array([ZMIN]))[0]
    grp["vmax"] = vmax
    volumesdss = (survey_volume(np.array([ZLIMIT]))[0]
                  - survey_volume(np.array([ZMIN]))[0])          # line 248
    volumesdssmin = volumesdss * vmax_floor_frac

    # lines 304-305: the same overwrite bug as gamahmf.r 307-308
    _ = np.where(vmax > volumesdss, volumesdss, vmax)
    w = np.where(vmax < volumesdssmin, volumesdssmin, vmax)
    w = np.where(grp.idcl.values == BAD_IDCL, volumesdss, w)      # line 306
    grp["weightszlimit"] = w

    if verbose:
        print(f"  SDSS groups after selection : {len(grp)}")
        print(f"  rank-{MULTI} galaxies matched   : {len(grp) - n_na} "
              f"({n_na} unmatched)")
        print(f"  volumesdss                  : {volumesdss:.6e} Mpc^3")
        print(f"  median vmax/volumesdss      : "
              f"{np.median(vmax) / volumesdss:.4f}")
    return grp, volumesdss


def bin_hmf(grp, volumesdss, seed, nmc=1001, verbose=True, nboot=0):
    """sdsshmf.r lines 312-370.

    ``nboot > 0`` swaps Driver's Poisson term for a group bootstrap.  SDSS needs
    this as much as GAMA does: its low-mass bins run at N_eff/N ~ 0.07-0.11, so
    the nominal 1/sqrt(N) understates the sampling error by a factor of ~3.
    """
    logm = np.log10(grp["mass"].values)
    w = 1.0 / grp.weightszlimit.values
    masserr = grp.log10MassErr.values

    raw = r_maghist(logm, MASSX)
    wtd = r_weighted_hist(logm, w, MASSX)
    cosvariance = cosvar(volumesdss, 1)                  # line 317, N=1 not 3

    rng = RRandom(seed)
    nbin = len(raw["mids"])
    mock = np.zeros((nmc, nbin))
    for i in range(nmc):
        mock[i] = r_weighted_hist(logm + rng.norm(len(masserr), masserr), w,
                                  MASSX)["counts"]

    meancounts = mock.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        edb = meancounts / wtd["counts"]
    edb[~np.isfinite(edb)] = 1.0

    mcerr = np.empty(nbin)
    for i in range(nbin):
        q = np.quantile((meancounts[i] - mock[:, i]) ** 2, 0.66, method="linear")
        with np.errstate(divide="ignore", invalid="ignore"):
            mcerr[i] = np.sqrt(q) / meancounts[i]

    # line 363: sqrt(N)/N, which is NaN (not Inf) at N = 0
    with np.errstate(divide="ignore", invalid="ignore"):
        rootnerr = np.sqrt(raw["counts"]) / raw["counts"]
    if nboot > 0:
        rootnerr = dr.bootstrap_error(logm, w, MASSX, nboot, seed, wtd["counts"])
    with np.errstate(divide="ignore", invalid="ignore"):
        sdssy = wtd["counts"] / (LOGBIN * edb)
        sdssf = np.sqrt(mcerr ** 2 + rootnerr ** 2)
    sdssf = np.where(np.isnan(sdssf), 0.9999, sdssf)
    sdssf = np.where(np.isinf(sdssf), 0.0, sdssf)
    sdssf = np.where(sdssf >= 1.0, 0.9999, sdssf)

    if verbose:
        print(f"  cosvariance                 : {cosvariance:.6f}")
    return dict(sdssx=wtd["mids"], raw=raw["counts"], wcounts=wtd["counts"],
                edb=edb, sdssy=sdssy, rootnerr=rootnerr, mcerr=mcerr,
                sdssf=sdssf, cosvariance=cosvariance)


def table(b):
    """The eight columns of sdsshmf.r line 443, in its reversed order."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return pd.DataFrame({
            "V1": b["sdssx"][::-1], "V2": b["raw"][::-1],
            "V3": np.log10(b["wcounts"])[::-1], "V4": np.log10(b["sdssy"])[::-1],
            "V5": b["rootnerr"][::-1], "V6": b["mcerr"][::-1],
            "V7": np.full(len(b["sdssx"]), b["cosvariance"]),
            "V8": b["sdssf"][::-1]})


def make_massfn(allx, ally, allf, volumesdss):
    """sdsshmf.r lines 200-209 -- one volume in the penalty, unlike allhmf.r."""
    allxxx = np.max(allx) + np.arange(1, 11) * LOGBIN
    ln10 = np.log(10.0)
    sig = allf / ln10

    def fn(p):
        mstar, phi, alpha, beta = p
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            penalty = (2 * volumesdss * np.sum(
                ln10 * beta * np.exp(-10 ** (beta * (allxxx - mstar)))
                * (phi * (10 ** allxxx / 10 ** mstar) ** (alpha + 1))) * LOGBIN)
            model = np.log10(beta * ln10 * np.exp(-10 ** (beta * (allx - mstar)))
                             * (phi * (10 ** allx / 10 ** mstar) ** (alpha + 1)))
            return np.sum(((ally - model) / sig) ** 2) + penalty

    return fn


def fit(b, volumesdss, phimrp, maxit=500):
    """sdsshmf.r line 386.  parscale=c(1,1,1,0.1), as in allhmf.r."""
    sel = (b["sdssy"] > 0) & (b["sdssx"] > MLIMIT)
    with np.errstate(divide="ignore"):
        allx, ally, allf = (b["sdssx"][sel], np.log10(b["sdssy"][sel]),
                            b["sdssf"][sel])
    fn = make_massfn(allx, ally, allf, volumesdss)
    par, val, nfe, conv = r_optim_nm(
        np.array([dr.MSTARMRP, phimrp, dr.ALPHAMRP, dr.BETAMRP]), fn, maxit=maxit,
        reltol=1e-8, parscale=np.array([1.0, 1.0, 1.0, 0.1]))
    return par, val, nfe, conv, len(allx)


def build(seed=10, nmc=1001, verbose=True, nboot=0, vmax_floor_frac=1e-3):
    """Convenience: everything, returning the V1..V8 table and the volume."""
    grp, volumesdss = build_groups(verbose=verbose,
                                   vmax_floor_frac=vmax_floor_frac)
    b = bin_hmf(grp, volumesdss, seed=seed, nmc=nmc, verbose=verbose, nboot=nboot)
    return table(b), volumesdss, b


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seed", type=int, default=10)
    p.add_argument("--nmc", type=int, default=1001)
    p.add_argument("--out", default="sdsshmf5.csv")
    args = p.parse_args()

    _, _, _, phimrp = lcdm_curve()
    print("SDSS HMF (sdsshmf.r), Tempel+14 DR10 groups")
    t, volumesdss, b = build(seed=args.seed, nmc=args.nmc)

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

    par, val, nfe, conv, n = fit(b, volumesdss, phimrp)
    print(f"\n  --- SDSS-only MRP fit ({n} bins with logM > {MLIMIT}) ---")
    for nm, v in zip(["log10(M*)", "log10(phi*)", "alpha", "beta"],
                     [par[0], np.log10(abs(par[1])), par[2], par[3]]):
        print(f"  {nm:<13}{v:9.3f}")
    print(f"  chi2 = {val:.3f}   fevals = {nfe}   convergence = {conv}")
    pub = (13.38, -3.00, -1.57, 0.47)          # Driver+22 table 2, "SDSS5"
    got = [par[0], np.log10(abs(par[1])), par[2], par[3]]
    print("  vs published SDSS5: " + "  ".join(
        f"{n}={p:.2f}({g - p:+.2f})" for n, p, g in
        zip(["M*", "phi*", "alpha", "beta"], pub, got)))

    t.to_csv(args.out, index=False)
    print(f"\n  wrote {args.out}  ({len(t)} rows, Driver's V1..V8 order)")


if __name__ == "__main__":
    main()
