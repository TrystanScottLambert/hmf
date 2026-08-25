#!/usr/bin/env python3
"""Driver+22's figure 3 -- the "vuvuzela" -- measured on the Nessie GAMA groups.

Driver's figure 3 calibrates the group mass error against multiplicity, and the
resulting arrays are hardcoded into ``gamahmf.r`` (line 270, ``xx``/``yy``) and
reused verbatim for SDSS in ``sdsshmf.r`` (line 281).  The derivation code was
never published.  ``vuvuzela.py`` measures the curve for the Nessie SDSS
catalogue; this script does the same for the Nessie GAMA catalogue, so that each
leg of the HMF can carry an error function derived from its own groups and its
own mass estimator rather than from Driver's v10 GAMA sample.

**The estimator differs from the SDSS one, and that is the point.**  The GAMA
leg of the HMF does not use Tempel+2014 eq. 8; ``new_gama_hmf._masses_and_errors``
rebuilds Driver's Robotham+2011 mass with the A factor,

    M = A sigma^2 R_50 / G ,     A = 13.9,

so the tracks here recompute ``VelDisp`` and ``Rad50`` from the surviving
members at each step.  Both are ported from ``fof/src/group_properties.rs`` --
``velocity_dispersion_gapper`` and ``calculate_radius`` -- and validated against
the DMU's own columns at full multiplicity by ``--validate``:

    VelDisp   max fractional difference   3e-16
    Rad50     max fractional difference   1.4e-07   (astropy vs Nessie D_C)
    Zfof      == median(z members)        exactly
    iterative centre == IterCenUberID     60/60

Members are removed **faintest-first in apparent magnitude**, following the
figure 3 caption's recipe, because that is what happens to a group receding
through a flux-limited survey.

Usage
-----
    uv run python vuvuzela_gama.py               # writes vuvuzela_gama.png
    uv run python vuvuzela_gama.py --validate    # the port checks above
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

import driver_recovery as dr
import vuvuzela as vz

GROUPS = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"
GALS = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CGal.fits"

# The DMU's own cosmology (make_gama_dmu/config.py lines 9-10).  Rad50 was built
# with it, so it has to be used here or the radii will not reproduce.
DMU_COSMO = FlatLambdaCDM(H0=100, Om0=0.25)

# make_gama_dmu/main.py passes vel_errors = 50 for every galaxy, and the Rust
# gapper subtracts their mean; hence velocity_dispersion_gap_err = sqrt(50).
SIGMA_ERR_SQ = 50.0


# ---------------------------------------------------------------------------
# The GAMA mass estimator, as the HMF actually computes it
# ---------------------------------------------------------------------------

def _cartesian_scaled(ra, dec, distance):
    """``convert_equitorial_to_cartesian_scaled``."""
    r, d = np.radians(ra), np.radians(dec)
    return np.column_stack([distance * np.cos(d) * np.cos(r),
                            distance * np.cos(d) * np.sin(r),
                            distance * np.sin(d)])


def _quantile_interpolated(sorted_a, q):
    """``stats.rs::quantile_interpolated`` -- R's type-7 quantile."""
    n = len(sorted_a)
    h = (n - 1) * q
    i = int(np.floor(h))
    frac = h - i
    return sorted_a[i] * (1.0 - frac) + sorted_a[i + 1] * frac if i + 1 < n \
        else sorted_a[i]


def group_radii(ra, dec, zs, ra0, dec0):
    """``Group::calculate_radius`` -- R50, Rsigma, R100 in Mpc/h.

    Every member is projected on to the sphere at the *group's* comoving
    distance, so this is a projected radius, and it inherits the DMU's h = 1.
    """
    zmed = float(np.median(zs))
    dist = float(DMU_COSMO.comoving_distance(zmed).value)
    pos = _cartesian_scaled(ra, dec, dist)
    cen = _cartesian_scaled(np.array([ra0]), np.array([dec0]), dist)[0]
    d = np.sort(np.linalg.norm(pos - cen, axis=1))
    return (_quantile_interpolated(d, 0.50), _quantile_interpolated(d, 0.68),
            d[-1])


def gama_mass(ra, dec, zs, absmag):
    """Driver's ``mymass`` (gamahmf.r 263) rebuilt from the surviving members.

    Rad50 is in Mpc/h and VelDisp in km/s at h = 1, so the ``(100/ho)`` inside
    the formula is the only h conversion needed -- exactly as in
    ``new_gama_hmf._masses_and_errors``.
    """
    if len(ra) < 3:
        return np.nan
    k = vz.iterative_centre_idx(ra, dec, absmag)
    r50, _, _ = group_radii(ra, dec, zs, ra[k], dec[k])
    sigma = vz.velocity_dispersion_gapper(zs, sigma_err_sq=SIGMA_ERR_SQ)
    if sigma <= 0 or r50 <= 0:
        return np.nan
    return (dr.MAGICA * (sigma * 1000.0) ** 2 * r50 * dr.PARSEC * 1e6
            / (dr.G * dr.MSOL) * (100.0 / dr.HO))


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def _to_frame(path):
    t = Table.read(path)
    out = {}
    for c in t.colnames:
        a = np.asarray(t[c])
        if a.dtype.byteorder == ">":
            a = a.astype(a.dtype.newbyteorder("="))
        out[c] = a
    return pd.DataFrame(out)


def load_members(verbose=True):
    """Nessie GAMA members, in the standard column names ``build_tracks`` wants.

    NOTE the ungrouped sentinel is region dependent: ``0`` in g09 but
    ``offset - 1`` elsewhere (199999, 299999, 499999), so ``GroupID != 0``
    leaves 87 602 ungrouped galaxies masquerading as three enormous groups.
    Selecting on membership of the *group table* is the safe filter.
    """
    grp = _to_frame(GROUPS)
    gal = _to_frame(GALS)
    real = set(grp.GroupID.values.tolist())
    gal = gal[gal.GroupID.isin(real)]
    gal = gal[np.isfinite(gal.AbsoluteMagR) & (gal.AbsoluteMagR > -40)
              & (gal.AbsoluteMagR < 0)]
    if verbose:
        print(f"  members loaded              : {len(gal)}")
    return grp, pd.DataFrame({"gid": gal.GroupID.values, "ra": gal.RAcen.values,
                              "dec": gal.Deccen.values, "z": gal.Z.values,
                              "appmag": gal.ApparentMagR.values,
                              "absmag": gal.AbsoluteMagR.values})


def validate(grp, gal, ntest=60):
    """Reproduce the DMU's own VelDisp and Rad50 at full multiplicity."""
    g = grp.set_index("GroupID")
    sizes = gal.groupby("gid").size()
    rows = []
    for gid in sizes[sizes > 20].index[:ntest]:
        s = gal[gal.gid == gid]
        k = vz.iterative_centre_idx(s.ra.values, s.dec.values, s.absmag.values)
        r50, _, r100 = group_radii(s.ra.values, s.dec.values, s.z.values,
                                   s.ra.values[k], s.dec.values[k])
        sig = vz.velocity_dispersion_gapper(s.z.values, sigma_err_sq=SIGMA_ERR_SQ)
        row = g.loc[gid]
        rows.append((row.VelDisp, sig, row.Rad50, r50, row.Rad100, r100,
                     row.Zfof, float(np.median(s.z.values))))
    d = pd.DataFrame(rows, columns=["vd_cat", "vd_new", "r50_cat", "r50_new",
                                    "r100_cat", "r100_new", "zfof", "zmed"])
    print(f"\n  Port validation on {len(d)} groups with N > 20:")
    print(f"    VelDisp  max |frac diff| : "
          f"{np.nanmax(np.abs(d.vd_new / d.vd_cat - 1)):.2e}")
    print(f"    Rad50    max |frac diff| : "
          f"{np.nanmax(np.abs(d.r50_new / d.r50_cat - 1)):.2e}")
    print(f"    Rad100   max |frac diff| : "
          f"{np.nanmax(np.abs(d.r100_new / d.r100_cat - 1)):.2e}")
    print(f"    Zfof - median(z)         : {np.nanmax(np.abs(d.zfof - d.zmed)):.2e}")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-mult", type=int, default=20)
    p.add_argument("--floor", type=int, default=3)
    p.add_argument("--out", default="vuvuzela_gama.png")
    p.add_argument("--csv", default="gama_masserr.csv")
    p.add_argument("--xmax", type=float, default=100,
                   help="multiplicity axis cap (default 100)")
    p.add_argument("--validate", action="store_true",
                   help="check the ported VelDisp/Rad50 against the DMU columns")
    args = p.parse_args()

    print("GAMA mass error vs multiplicity (Driver+22 figure 3, Nessie GAMA)")
    grp, gal = load_members()

    if args.validate:
        validate(grp, gal)
        return

    tracks = vz.build_tracks(gal, min_mult=args.min_mult, floor=args.floor,
                             mass_fn=gama_mass)
    ns, q16, q50, q84, cnt = vz.quantile_curve(tracks)
    err = vz.error_function(ns, q16, q84)

    # q50 is the analogue of gamahmf.r's `masscorr` (the multiplicity debiasing)
    # and `err` of its `yy` (the scatter).  Both of Driver's are hardcoded from
    # his v10 GAMA sample; these are measured from ours.
    print("\n     N   tracks     q16     q50     q84   sigma |  Driver yy  masscorr")
    for i, n in enumerate(ns):
        if not np.isfinite(q50[i]) or n > 40:
            continue
        g = np.interp(n, vz.GAMA_XX, vz.GAMA_YY, left=np.nan, right=np.nan)
        m = dr.MASSCORR[n - 1] if n <= len(dr.MASSCORR) else np.nan
        print(f"  {n:4d} {cnt[i]:8d} {q16[i]:7.3f} {q50[i]:7.3f} {q84[i]:7.3f}"
              f" {err[i]:7.3f} | {g:10.3f} {m:9.3f}")

    sel = (ns >= 3) & (ns <= 22) & np.isfinite(q50)
    gg = np.interp(ns[sel], vz.GAMA_XX, vz.GAMA_YY)
    mm = np.array([dr.MASSCORR[int(n) - 1] for n in ns[sel]])
    skew = np.median(np.abs(q16[sel]) / np.abs(q84[sel]))
    print(f"\n  sigma  ours/Driver ratio : {np.median(err[sel] / gg):.3f}"
          "   (1.0 => Driver's v10 curve transfers to Nessie GAMA)")
    print(f"  sigma  ours-Driver, dex  : {np.median(err[sel] - gg):+.4f}")
    print(f"  bias   ours-masscorr, dex: {np.median(q50[sel] - mm):+.4f}"
          "   (q50 is the analogue of masscorr; both -> 0 at large N,")
    print("                             so a ratio is meaningless there)")
    print(f"  skew   |q16|/|q84|       : {skew:.2f}"
          "   (1.0 would be symmetric; Driver's MC assumes a Gaussian)")

    out = pd.DataFrame({"N": ns, "ntracks": cnt, "q16": q16, "q50": q50,
                        "q84": q84, "sigma_log10M": err})
    out.to_csv(args.csv, index=False)
    print(f"\n  wrote {args.csv}")
    vz.plot_vuvuzela(tracks, ns, q16, q50, q84, args.out,
                     panel_label=f"Nessie GAMA ($N>{args.min_mult}$)",
                     ref_label="Driver et al. 2022 (v10 GAMA)", xmax=args.xmax)


if __name__ == "__main__":
    main()
