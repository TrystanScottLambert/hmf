#!/usr/bin/env python3
"""The Nessie GAMA and SDSS group catalogues on one schema -- G3C v10's.

Writes, into ``--outdir``:

    G3CFoFGroup_nessie.fits   / .par     GAMA (G09, G12, G15, G23)
    SDSSFoFGroup_nessie.fits  / .par     SDSS DR10 (Tempel+14 galaxies)
    catalogue_columns.tex                every column, which catalogue has it
    catalogue_dropped_v10.tex            the v10 columns not provided, and why

Both group files have identical columns in identical order.  Every column that
exists in ``G3CFoFGroupv10.fits`` keeps its v10 name, so pipelines written for
v10 run unchanged.  ``-999`` marks a value a catalogue cannot provide.

**One deliberate departure.**  v10's ``IterCenCATAID`` and ``BCGCATAID`` become
``IterCenID`` and ``BCGID``: neither new catalogue has CATAIDs (GAMA carries
UberID, SDSS Tempel's galid), so keeping the v10 name would let an old pipeline
join them against CATAID tables and match the wrong galaxies without error.

**Everything is re-derived through Nessie on the stored memberships.**  Setting
``RedshiftCatalog.group_ids`` directly and calling ``calculate_group_table``
reproduces both inputs bit-for-bit when given their original inputs (asserted
below for GAMA; verified for SDSS when this was written).  That is how the two
things the inputs do not carry are obtained:

* GAMA's ``estimated_mass``, which the DMU computes and then drops.
* SDSS in v10's h = 1 units.  ``sdss_groups.parquet`` was built at h = 0.7 with
  ``absmag = rmag - distmod`` and no k-correction; here it is rebuilt at h = 1
  with Tempel's k+e corrected magnitudes.  The groups themselves are identical.

**Tempel's DR10 extract is in physical h = 0.678 units, not the paper's h^-1.**
Its comoving distance (``col11``) fits Planck 2013 (h = 0.678, Om = 0.308) with
zero scatter, and ``col31`` then implies a k+e of 0.009 mag.  So ``col31`` is
translated with Tempel's distance modulus swapped for ours, not by 5 log h.

The masses:

* ``MassHern`` -- Tempel+14 eq. 8, Hernquist.  Nessie's ``estimated_mass``, times
  3^(1/3) only if it came from pre-fix Nessie (cbrt(3) where eq. 8 needs sqrt(3);
  fixed in nessie-py 1.1.1) -- detected from the masses, see hernquist_mass ("Hernquist
  vs NFW").
* ``MassNFW`` -- Tempel+14 eq. 8, NFW: kappa(M) calibrated on his ``mass_nfw``,
  his eq. 3 rms dispersion, sigma_sky from the members.  The HMF's
  ``--mass-mode tempel_rms`` recipe, in h = 1 units.

The galaxy catalogues wait on the v10 G3CGal schema.

Usage
-----
    uv run python build_catalogues.py
    uv run python build_catalogues.py --only sdss
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from dataclasses import dataclass
from datetime import date

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from nessie import FlatCosmology, RedshiftCatalog
from scipy.integrate import quad
from scipy.interpolate import interp1d

import nessie_sdss_hmf as nsd
from sdss_hmf import MAGLIM as SDSS_MAGLIM
from vuvuzela_gama import _to_frame

DMU_DIR = "/Users/00115372/Desktop/my_tools/make_gama_dmu"
GAMA_GROUPS = f"{DMU_DIR}/G3CFoFGroup.fits"
GAMA_GALS = f"{DMU_DIR}/G3CGal.fits"
V10_GROUPS = "../data/G3CFoFGroupv10.fits"
TEMPEL_GALS = "../data/sdssdr10table1.fits"

GAMA_H, GAMA_OM = 1.0, 0.25          # make_gama_dmu/config.py
SDSS_H, SDSS_OM = 1.0, 0.30          # h = 1 for v10 units; Nessie's own Om
TEMPEL_H, TEMPEL_OM = 0.678, 0.308   # fitted to his col11, zero scatter
VEL_ERR = 50.0                       # km/s, as both inputs were built
NULL = -999.0


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

@dataclass
class Col:
    name: str
    unit: str
    ucd: str
    v10: str        # "Y", or what it replaces
    gama: bool
    sdss: bool
    description: str


SCHEMA = [
    Col("GroupID", "-", "meta.id;meta.main", "Y", True, True,
        "Unique group ID, as listed in the associated galaxy catalogue. GAMA: "
        "Nessie ID plus a region offset (1e5 G09, 2e5 G12, 3e5 G15, 5e5 G23)"),
    Col("Nfof", "-", "-", "Y", True, True, "Group multiplicity"),
    Col("IterCenID", "-", "meta.id", "IterCenCATAID", True, True,
        "ID of the iterative central galaxy (preferred position of the group "
        "centre): GAMA III UberID for GAMA, Tempel+14 galid for SDSS. Not a CATAID"),
    Col("IterCenRA", "deg", "pos.eq.ra", "Y", True, True,
        "RA of the iterative central galaxy (J2000)"),
    Col("IterCenDec", "deg", "pos.eq.dec", "Y", True, True,
        "Dec of the iterative central galaxy (J2000)"),
    Col("IterCenZ", "-", "src.redshift", "Y", True, True,
        "Redshift of the iterative central galaxy"),
    Col("Zfof", "-", "src.redshift", "Y", True, True,
        "Median redshift of the group"),
    Col("Rad50", "Mpc/h", "phys.size.radius", "Y", True, True,
        "Group radius defined by the 50th percentile group member, based on the "
        "projected comoving distance from IterCen"),
    Col("Rad1Sig", "Mpc/h", "phys.size.radius", "Y", True, True,
        "Group radius defined by the 68th percentile group member, based on the "
        "projected comoving distance from IterCen"),
    Col("Rad100", "Mpc/h", "phys.size.radius", "Y", True, True,
        "Group radius defined by the most distant group member, based on the "
        "projected comoving distance from IterCen"),
    Col("VelDisp", "km/s", "phys.veloc", "Y", True, True,
        "Gapper velocity dispersion with VelErr subtracted in quadrature; set to "
        "zero if VelErr exceeds the raw dispersion"),
    Col("VelErr", "km/s", "phys.veloc", "Y", True, True,
        "Total group velocity error subtracted from VelDisp (input galaxy "
        "velocity errors of 50 km/s)"),
    Col("TotRmag", "mag-5log(h)", "phys.magAbs", "Y", True, True,
        "k+e corrected r-band absolute magnitude of the group, from the summed "
        "member luminosities (GAMA: DMU k+e polynomial; SDSS: Tempel+14 KCORRECT)"),
    Col("TotFluxProxy", "-", "-", "Y", True, True,
        "Proxy for group total r-band luminosity down to M_r - 5 log h = -14 in "
        "solar luminosities, from the observed flux and the fraction of the GAMA "
        "LF observable at the group redshift (B = 1 in Eq. 22 of Robotham et al. "
        "2011). SDSS uses the r < 17.77 limit"),
    Col("MassProxy", "Msun/h", "-", "Y", True, True,
        "Dynamical mass proxy Rad50 times VelDisp squared over G, i.e. A = 1 in "
        "Eq. 18 of Robotham et al. (2011)"),
    Col("BCGID", "-", "meta.id", "BCGCATAID", True, True,
        "ID of the brightest group galaxy: GAMA III UberID for GAMA, Tempel+14 "
        "galid for SDSS. Not a CATAID"),
    Col("BCGRA", "deg", "pos.eq.ra", "Y", True, True, "RA of BCG (J2000)"),
    Col("BCGDec", "deg", "pos.eq.dec", "Y", True, True, "Dec of BCG (J2000)"),
    Col("CenRA", "deg", "pos.eq.ra", "Y", True, True,
        "RA of the r-band luminosity weighted centre of the group (J2000)"),
    Col("CenDec", "deg", "pos.eq.dec", "Y", True, True,
        "Dec of the r-band luminosity weighted centre of the group (J2000)"),
    Col("MassA", "Msun/h", "phys.mass", "Y", True, False,
        "MassProxy times the global A factor, A = 10.0 (Robotham et al. 2011 "
        "section 4.3). Calibrated on GAMA mocks: GAMA only"),
    Col("MassAfunc", "Msun/h", "phys.mass", "Y", True, False,
        "MassProxy times the functional A factor, a function of Nfof and Zfof "
        "(Robotham et al. 2011 section 4.3). Calibrated on GAMA mocks: GAMA only"),
    Col("LumB", "Lsun/h^2", "phys.luminosity", "Y", True, False,
        "TotFluxProxy times the global B factor, B = 1.04 (Robotham et al. 2011 "
        "section 4.4). Calibrated on GAMA mocks: GAMA only"),
    Col("LumBfunc", "Lsun/h^2", "phys.luminosity", "Y", True, False,
        "TotFluxProxy times the functional B factor, a function of Nfof and Zfof "
        "(Robotham et al. 2011 section 4.4). Calibrated on GAMA mocks: GAMA only"),
    Col("BCGZ", "-", "src.redshift", "-", True, True, "Redshift of BCG"),
    Col("MassHern", "Msun/h", "phys.mass", "-", True, True,
        "Dynamical mass from Eq. 8 of Tempel et al. (2014) for a Hernquist "
        "profile: Rg = 4.582 sigma_sky, sigma_3D = sqrt(3) VelDisp, with "
        "sigma_sky the projected rms member distance from IterCen. Zero where "
        "VelDisp is zero"),
    Col("MassNFW", "Msun/h", "phys.mass", "-", True, True,
        "Dynamical mass from Eq. 8 of Tempel et al. (2014) for an NFW profile: "
        "Rg = kappa(M) sigma_sky with kappa calibrated on the Tempel et al. "
        "(2014) NFW masses, sigma_sky the projected rms member distance from "
        "IterCen, sigma_3D = sqrt(3) times the rms dispersion of their Eq. 3. "
        "Zero where all members share one redshift"),
    Col("Region", "-", "-", "-", True, True,
        "Survey region: g09, g12, g15 or g23 for GAMA; SDSS for SDSS"),
]

# v10 columns this schema does not carry, grouped by why.
DROPPED_V10 = [
    ("Need the randoms and redshift-completeness maps, which were not "
     "propagated to Nessie", ["RelDen", "Zcomp", "GroupEdge"]),
    ("FoF link bookkeeping, not exposed by Nessie",
     ["LinkTot", "LinkStrenProj", "LinkStrenRad", "LinkStrenComb"]),
    ("Shape and kinematic statistics, not implemented in Nessie",
     ["RadKurt", "AxAng", "AxRat", "SymKS", "VelSkew", "VelKurt", "Modality",
      "d3vol", "d3area", "d2radec", "d2raz", "d2decz", "distortcirc",
      "distortellip", "distortsph"]),
    ("Derivable: VelDispRaw is the quadrature sum of VelDisp and VelErr where "
     "VelDisp > 0; Rgap from the galaxy catalogue", ["VelDispRaw", "Rgap"]),
    ("Renamed: the new catalogues carry no CATAIDs (see IterCenID, BCGID)",
     ["IterCenCATAID", "BCGCATAID"]),
    ("No image service", ["URL_IMG"]),
]


def check_schema_covers_v10():
    """Every v10 column is either kept under its name or listed as dropped."""
    v10 = set(fits.getdata(V10_GROUPS, 1).columns.names)   # no unit parsing
    kept = {c.name for c in SCHEMA if c.v10 == "Y"}
    dropped = {n for _, names in DROPPED_V10 for n in names}
    missing = v10 - kept - dropped
    extra = (kept | dropped) - v10
    if missing or extra or kept & dropped:
        raise RuntimeError(f"schema vs v10: missing {missing}, not in v10 {extra}, "
                           f"both kept and dropped {kept & dropped}")


# ---------------------------------------------------------------------------
# make_gama_dmu helpers
# ---------------------------------------------------------------------------

def _dmu():
    """make_gama_dmu/main.py, loaded by path: python_rewrite has a main.py too."""
    if DMU_DIR not in sys.path:
        sys.path.append(DMU_DIR)
    spec = importlib.util.spec_from_file_location("dmu_main", f"{DMU_DIR}/main.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def integrated_lf(path=f"{DMU_DIR}/gama_lf/lf.dat", cut=-14):
    """``luminosity_function.build_integrated_lf``, ported.

    The original uses ``read_csv(delim_whitespace=True)``, which this pandas no
    longer accepts; everything else is unchanged.
    """
    df = pd.read_csv(path, sep=r"\s+")
    df = df[df["mag_bin_centre"] < cut]
    mags = df["mag_bin_centre"].values
    phi_lum = df["luminosity_function"].values * 10 ** (-0.4 * mags)
    func = interp1d(mags, phi_lum, bounds_error=False,
                    fill_value=(phi_lum[0], phi_lum[-1]))
    integrals = np.array([quad(func, -30, m, limit=1000)[0] for m in mags])
    integrals[integrals == 0] = integrals[integrals > 0].min()
    return interp1d(mags, integrals, bounds_error=False,
                    fill_value=(integrals[0], integrals[-1]))


# ---------------------------------------------------------------------------
# Group properties
# ---------------------------------------------------------------------------

def nessie_group_table(ra, dec, z, absmag, group_ids, cosmo):
    """Nessie's group table for a fixed membership (no FoF run)."""
    rc = RedshiftCatalog(np.asarray(ra, float), np.asarray(dec, float),
                         np.asarray(z, float), lambda zz: np.ones_like(zz), cosmo)
    rc.group_ids = np.asarray(group_ids).astype(int)
    tab = rc.calculate_group_table(np.asarray(absmag, float),
                                   np.full(len(rc.group_ids), VEL_ERR))
    return pd.DataFrame(tab).set_index("group_id")


def sky_dispersion(gid, ra, dec, tab, cosmo):
    """Tempel+14 eq. 4 about the iterative centre, comoving over (1+z).

    ``Group::calculate_sky_distribution``: agrees with the sigma_sky implied by
    Nessie's own ``estimated_mass`` to a median 3e-9 dex.  Defined for every
    group, including those with VelDisp = 0 where the inversion is not.
    """
    cen = tab.reindex(np.asarray(gid))
    dra = ((np.asarray(ra) - cen.iter_ra.values + 180.0) % 360.0 - 180.0) \
        * np.cos(np.radians(cen.iter_dec.values))
    ddec = np.asarray(dec) - cen.iter_dec.values
    dc = cosmo.comoving_distance(cen.median_redshift.values)
    r2 = (np.radians(np.hypot(dra, ddec)) * dc) ** 2
    s = pd.DataFrame({"g": np.asarray(gid), "r2": r2}).groupby("g")["r2"] \
        .agg(["sum", "size"])
    zm = tab.median_redshift.reindex(s.index)
    return np.sqrt(s["sum"] / (2.0 * s["size"] * (1.0 + zm) ** 2)).reindex(tab.index)


def nfw_mass(sigma_sky, sigma_v):
    """Tempel+14 eq. 8 with the NFW kappa, iterated because kappa depends on M.

    In h = 1 units (sigma_sky in Mpc/h) the mass is in Msun/h.  kappa was
    calibrated on Tempel's physical h = 0.678 masses, so it is looked up there.
    """
    m0 = 2.325e12 * np.asarray(sigma_sky) * 3.0 * (np.asarray(sigma_v) / 100.0) ** 2
    m = m0 * nsd.KAPPA_NFW
    with np.errstate(divide="ignore", invalid="ignore"):
        for _ in range(30):
            k = np.interp(np.log10(m / TEMPEL_H), nsd.KAPPA_NFW_LOGM,
                          nsd.KAPPA_NFW_VAL)
            m = m0 * k
    return np.where(m0 > 0, m, 0.0)


def hernquist_mass(estimated_mass, sigma_gap, sky_disp):
    """Nessie's ``estimated_mass`` as a correct sqrt(3) Hernquist mass.

    Pre-fix Nessie (nessie-py <= 1.0.x) coded cbrt(3) where eq. 8 needs sqrt(3);
    nessie-py 1.1.1 has the fix.  Which one ran is read off the masses
    themselves, by comparing with the cbrt(3) formula evaluated on this group's
    own gapper sigma and sky dispersion, so the correction is never applied
    twice.
    """
    em = np.asarray(estimated_mass, float)
    old = (2.325e12 * nsd.KAPPA_HERNQUIST * np.asarray(sky_disp, float)
           * (nsd.CBRT3 * np.asarray(sigma_gap, float) / 100.0) ** 2)
    ok = (em > 0) & (old > 0)
    return em * nsd.sqrt3_factor_from_ratio(np.log10(em[ok] / old[ok]),
                                            "re-derived estimated_mass")


def _decode(a):
    a = np.asarray(a)
    return np.char.decode(a) if a.dtype.kind == "S" else a.astype(str)


# ---------------------------------------------------------------------------
# GAMA
# ---------------------------------------------------------------------------

# DMU column -> Nessie key, for the bit-for-bit reproduction check.
_DMU_REPRO = {"Nfof": "multiplicity", "IterCenRA": "iter_ra",
              "Zfof": "median_redshift", "Rad50": "r50", "VelDisp":
              "velocity_dispersion_gap", "MassProxy": "mass_proxy",
              "TotRmag": "total_absolute_mag", "BCGRA": "bcg_ras"}


def build_gama(verbose=True):
    grp = _to_frame(GAMA_GROUPS)
    gal = _to_frame(GAMA_GALS)
    # The DMU's ungrouped sentinel is 0 in G09 but offset - 1 elsewhere
    # (199999, 299999, 499999), so membership of the group table is the filter.
    grouped = gal.GroupID.isin(grp.GroupID).values
    cosmo = FlatCosmology(GAMA_H, GAMA_OM)
    tab = nessie_group_table(gal.RAcen, gal.Deccen, gal.Z, gal.AbsoluteMagR,
                             np.where(grouped, gal.GroupID, -1), cosmo)
    tab = tab.loc[grp.GroupID.values]
    for dmu_col, key in _DMU_REPRO.items():
        if not np.array_equal(grp[dmu_col].values, tab[key].values):
            raise RuntimeError(f"Nessie no longer reproduces the DMU's {dmu_col}")

    mem = gal[grouped]
    sky = sky_dispersion(mem.GroupID.values, mem.RAcen.values, mem.Deccen.values,
                         tab, cosmo)
    sig = nsd.rms_dispersion_by_group(mem, zcol="Z", gcol="GroupID") \
        .reindex(tab.index)

    out = pd.DataFrame({
        "GroupID": grp.GroupID, "Nfof": grp.Nfof, "IterCenID": grp.IterCenUberID,
        "IterCenRA": grp.IterCenRA, "IterCenDec": grp.IterCenDec,
        "IterCenZ": grp.IterCenZ, "Zfof": grp.Zfof, "Rad50": grp.Rad50,
        "Rad1Sig": grp.Rad1Sig, "Rad100": grp.Rad100, "VelDisp": grp.VelDisp,
        "VelErr": grp.VelDispErr, "TotRmag": grp.TotRmag,
        "TotFluxProxy": grp.TotFluxProxy, "MassProxy": grp.MassProxy,
        "BCGID": grp.BCGUberID, "BCGRA": grp.BCGRA, "BCGDec": grp.BCGDec,
        "CenRA": grp.CenRA, "CenDec": grp.CenDec, "MassA": grp.MassA,
        "MassAfunc": grp.MassAfunc, "LumB": grp.LumB, "LumBfunc": grp.LumBfunc,
        "BCGZ": grp.BCGZ,
        "MassHern": hernquist_mass(tab.estimated_mass.values,
                                   tab.velocity_dispersion_gap.values, sky.values),
        "MassNFW": nfw_mass(sky.values, sig.values),
        "Region": _decode(grp.GAMARegion),
    })
    if verbose:
        print(f"  GAMA groups                 : {len(out)} (DMU reproduced)")
    return out


# ---------------------------------------------------------------------------
# SDSS
# ---------------------------------------------------------------------------

def _tempel_join(gal):
    """Tempel's galid, redshift and absolute magnitude, on (RA, Dec)."""
    t = Table.read(TEMPEL_GALS)
    nat = lambda c: np.asarray(t[c]).astype(np.asarray(t[c]).dtype.newbyteorder("="))
    tp = pd.DataFrame({"galid": nat("col1"), "ra_t": nat("col13"),
                       "dec_t": nat("col14"), "z_t": nat("col9"),
                       "absmag_t": nat("col31")})
    key = lambda r, d: np.round(r, 6).astype(str) + "_" + np.round(d, 6).astype(str)
    tp["_k"] = key(tp.ra_t.values, tp.dec_t.values)
    gal = gal.assign(_k=key(gal.RAJ2000.values, gal.DEJ2000.values))
    out = gal.merge(tp.drop(columns=["ra_t", "dec_t"]), on="_k", how="left")
    if len(out) != len(gal) or out.galid.isna().any():
        raise RuntimeError("Nessie and Tempel galaxies no longer match one-to-one")
    return out.drop(columns="_k")


def build_sdss(verbose=True):
    ref = pd.read_parquet(nsd.NESSIE_GROUPS).set_index("group_id")
    gal = _tempel_join(pd.read_parquet(nsd.NESSIE_GALS))
    cosmo = FlatCosmology(SDSS_H, SDSS_OM)
    # Tempel's k+e corrected magnitude with his distance modulus swapped for ours.
    absmag = (gal.absmag_t.values
              + FlatCosmology(TEMPEL_H, TEMPEL_OM).dist_mod(gal.z_t.values)
              - cosmo.dist_mod(gal.z_t.values))
    tab = nessie_group_table(gal.RAJ2000, gal.DEJ2000, gal.zobs, absmag,
                             gal.group_id, cosmo)
    tab = tab.loc[ref.index]
    for key in ("multiplicity", "median_redshift", "velocity_dispersion_gap"):
        if not np.array_equal(ref[key].values, tab[key].values):
            raise RuntimeError(f"re-derived {key} differs from sdss_groups.parquet")

    dmu = _dmu()
    ob_limit = (SDSS_MAGLIM - cosmo.dist_mod(tab.median_redshift.values)
                - dmu.calc_ke_correction(tab.median_redshift.values))
    lf = integrated_lf()
    lum_factor = lf(dmu.AB_CUT) / lf(ob_limit)

    mem = gal[gal.group_id != nsd.UNGROUPED]
    sky = sky_dispersion(mem.group_id.values, mem.RAJ2000.values,
                         mem.DEJ2000.values, tab, cosmo)
    sig = nsd.rms_dispersion_by_group(mem, zcol="z_t", gcol="group_id") \
        .reindex(tab.index)
    galid = gal.galid.values

    out = pd.DataFrame({
        "GroupID": tab.index.values, "Nfof": tab.multiplicity.values,
        "IterCenID": galid[tab.iter_idx.values],
        "IterCenRA": tab.iter_ra.values, "IterCenDec": tab.iter_dec.values,
        "IterCenZ": tab.iter_redshift.values, "Zfof": tab.median_redshift.values,
        "Rad50": tab.r50.values, "Rad1Sig": tab.rsigma.values,
        "Rad100": tab.r100.values, "VelDisp": tab.velocity_dispersion_gap.values,
        "VelErr": tab.velocity_dispersion_gap_err.values,
        "TotRmag": tab.total_absolute_mag.values,
        "TotFluxProxy": tab.flux_proxies.values * lum_factor,
        "MassProxy": tab.mass_proxy.values,
        "BCGID": galid[tab.bcg_idxs.values], "BCGRA": tab.bcg_ras.values,
        "BCGDec": tab.bcg_decs.values,
        "CenRA": tab.center_of_light_ras.values % 360.0,
        "CenDec": tab.center_of_light_decs.values,
        "MassA": NULL, "MassAfunc": NULL, "LumB": NULL, "LumBfunc": NULL,
        "BCGZ": tab.bcg_redshifts.values,
        "MassHern": hernquist_mass(tab.estimated_mass.values,
                                   tab.velocity_dispersion_gap.values, sky.values),
        "MassNFW": nfw_mass(sky.values, sig.values),
        "Region": "SDSS",
    })
    if verbose:
        same = np.mean(tab.iter_idx.values == ref.iter_idx.values)
        print(f"  SDSS groups                 : {len(out)} (membership reproduced)")
        print(f"  iterative centre unchanged  : {same:.1%} "
              "(the rest move with the k+e corrected magnitudes)")
    return out


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_fits(df, path):
    df = df[[c.name for c in SCHEMA]]
    Table.from_pandas(df).write(path, overwrite=True)


def write_par(title, abstract, path):
    """The make_gama_dmu ``.par`` layout, with this catalogue's own header."""
    rows = [(c.name, str(i + 1), c.unit, c.ucd, c.description)
            for i, c in enumerate(SCHEMA)]
    w = [max(len(r[j]) for r in rows) for j in range(4)]
    with open(path, "w", encoding="utf8") as fh:
        fh.write(f"# GroupFinding\n# {title}\n# {date.today()}\n#\n# {abstract}\n#\n")
        for r in rows:
            fh.write("  ".join(f"{r[j]:<{w[j]}}" for j in range(4)) + f"  {r[4]}\n")


TEX_UNITS = {"-": "--", "deg": "deg", "km/s": r"km\,s$^{-1}$",
             "Mpc/h": r"$h^{-1}$\,Mpc", "mag-5log(h)": r"mag\,$-\,5\log h$",
             "Msun/h": r"$h^{-1}\,\mathrm{M}_\odot$",
             "Lsun/h^2": r"$h^{-2}\,\mathrm{L}_\odot$"}


def _tex(s):
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("_", r"\_"), ("#", r"\#"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def write_tex_columns(path):
    yes = r"$\checkmark$"
    lines = [r"\begin{tabular}{llcccp{7.5cm}}", r"\hline",
             r"Column & Units & v10 & GAMA & SDSS & Description \\", r"\hline"]
    for c in SCHEMA:
        v10 = yes if c.v10 == "Y" else ("--" if c.v10 == "-" else _tex(c.v10))
        lines.append(" & ".join([_tex(c.name), TEX_UNITS[c.unit], v10,
                                 yes if c.gama else "--", yes if c.sdss else "--",
                                 _tex(c.description)]) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_tex_dropped(path):
    lines = [r"\begin{tabular}{p{6cm}p{8cm}}", r"\hline",
             r"v10 columns & Reason not provided \\", r"\hline"]
    for reason, names in DROPPED_V10:
        lines.append(_tex(", ".join(names)) + " & " + _tex(reason) + r" \\")
    lines += [r"\hline", r"\end{tabular}"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="catalogues")
    p.add_argument("--only", choices=["gama", "sdss"])
    args = p.parse_args()

    check_schema_covers_v10()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)

    if args.only in (None, "gama"):
        write_fits(build_gama(), out("G3CFoFGroup_nessie.fits"))
        write_par("G3CFoFGroup (Nessie, G3C v10 schema)",
                  "GAMA III groups for G09, G12, G15 and G23, built with Nessie.",
                  out("G3CFoFGroup_nessie.par"))
    if args.only in (None, "sdss"):
        write_fits(build_sdss(), out("SDSSFoFGroup_nessie.fits"))
        write_par("SDSSFoFGroup (Nessie, G3C v10 schema)",
                  "SDSS DR10 groups built with Nessie on the Tempel et al. (2014) "
                  "galaxies.", out("SDSSFoFGroup_nessie.par"))
    write_tex_columns(out("catalogue_columns.tex"))
    write_tex_dropped(out("catalogue_dropped_v10.tex"))
    print(f"  wrote {args.outdir}/")


if __name__ == "__main__":
    main()
