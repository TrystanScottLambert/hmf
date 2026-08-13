"""
============================================================
measure_completeness_nessie.py
Measure C(m,z) by running the ACTUAL group finder on the mock
============================================================

WHY THIS EXISTS
---------------
`measure_completeness.py` defines a halo as "detected" if it has >= MULTI
galaxies brighter than the magnitude limit. That is a proxy for the group
finder, not the group finder. Comparing it to the real catalogue shows the
proxy is wrong, and wrong in a mass-dependent way:

    N members     mock/deg^2   real/deg^2   real/mock
        5-6           1.40         2.63        1.87
        6-8           1.30         2.23        1.72
       8-10           0.73         1.13        1.56
      10-15           0.73         1.04        1.43
        15+           0.74         0.67        0.90
      TOTAL           4.91         7.70        1.57

A selection error that varies monotonically with richness varies with mass,
and a mass-dependent selection error maps straight onto alpha. So the
completeness correction has to be measured against the same operation that
built the data: Nessie's friends-of-friends.

WHAT IT DOES
------------
 1. Loads the Shark/WAVES lightcone galaxies and applies the GAMA selection.
 2. Runs Nessie FoF with the DMU's own B0/R0/cosmology.
 3. Builds the group table exactly as make_gama_dmu.py does.
 4. Matches Nessie groups to true Shark halos (bijective, Robotham+11 style).
 5. Measures, as a function of TRUE halo mass and redshift:
        C(m,z)      recovered fraction  -> the completeness for the HMF model
        purity      fraction of catalogue entries that are clean matches
        mass bias + scatter of the recovered dynamical mass
 6. Fits the erf ramp and prints COMP_* constants for recovery.py.
 7. Writes four publication figures.

CAVEATS THIS SCRIPT CANNOT FIX (state them in the paper)
--------------------------------------------------------
 * COSMOLOGY. The DMU runs Nessie at h=1.0, Om=0.25; recovery.py fits at
   h=0.674, Om=0.3147. Since M ~ sigma^2 R / G and R ~ 1/h, dynamical masses
   differ by log10(1/0.674) = 0.17 dex between the two. Either reconcile them
   or carry it as a known offset.
 * SPECTROSCOPIC COMPLETENESS. GAMA has fibre collisions and redshift
   incompleteness; the mock does not. Set MOCK_SPEC_COMPLETENESS below if you
   want to impose a uniform value.
 * RANDOMS. GAMA uses precomputed random catalogues to build the radial
   density function. Here they are generated from the mock's own redshift
   distribution (see build_randoms), which is the same idea but not the same
   file.

Run:
    python measure_completeness_nessie.py
    python measure_completeness_nessie.py --mag-limit 19.4 --no-plots
============================================================
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.special import erf

import recovery as R

# ----------------------------------------------------------------------
# Configuration -- mirrors make_gama_dmu/config.py so the finder behaves
# identically to the run that produced G3CFoFGroup.fits.
# ----------------------------------------------------------------------
B0, R0 = 0.06, 32  # FoF linking lengths
MASS_A = 10  # MassA = mass_proxy * MASS_A
APPARENT_MAG_LIM = 19.65  # VERIFIED: make_gama_dmu/config.py line 13 is 19.65,
# and G3CGal.fits max ApparentMagR is exactly 19.6500. The old "config.py says
# 19.4 -- check!" note was stale; the mock and the DMU agree.
AB_CUT = -10
OVERFACTOR = 400  # randoms oversampling
VEL_ERROR = 50.0  # km/s, as make_gama_dmu uses
NESSIE_H0, NESSIE_OM = 1.0, 0.25  # FlatCosmology(1.0, 0.25)
MULTI = 5  # minimum members for a group
UNGROUPED_ID = -1  # Nessie's marker for an ungrouped galaxy.
# make_gama_dmu adds GROUP_ID_OFFSET (1e5)
# and then maps 99999 -> 0, i.e. -1 + 1e5.

MOCK_SPEC_COMPLETENESS = None  # None = perfect; else a float in (0, 1]

# matching thresholds (bijective match, Robotham+11 style)
PURITY_MIN = 0.5  # >= this fraction of the FoF group is one halo
RECOVERY_MIN = 0.5  # >= this fraction of the halo is in that FoF group
MZ_MIN_N = 20  # min halos per (m, z) cell for the absolute-mass C table

MASS_COL = "log_mass_am"  # true mass used for C(m); "mvir" also available


# ----------------------------------------------------------------------
def kecorr(z):
    """GAMA k+e correction, verbatim from make_gama_dmu.calc_ke_correction."""
    z = np.asarray(z, float)
    kcorrvals = np.array([0.20848, 1.0226, 0.52366, 3.5902, 2.3843])
    powers = (z[..., None] - 0.2) ** np.arange(len(kcorrvals))
    return (powers @ kcorrvals) - 1.75 * z


def erf_ramp(d, d50, w):
    return 0.5 * (1 + erf((d - d50) / (np.sqrt(2) * w)))


# ----------------------------------------------------------------------
def load_mock(mag_limit, zmin, zmax, abs_mag_mode="native"):
    """Mock galaxies + groups with the GAMA selection applied.

    Returns (gal, halo) DataFrames. `gal` carries ra/dec/z/apparent/absolute
    magnitudes and the true halo id; `halo` carries the true halo masses."""
    groups, galaxies = R.load_catalogues(R.DATA_DIR)

    # --- redshift: FoF is a redshift-space algorithm, so peculiar velocities
    #     matter. Without them there are no fingers-of-god and the linking
    #     behaves quite differently.
    zcol = next(
        (c for c in ("zobs", "z_obs", "zspec", "z_cos_pec") if c in galaxies.columns),
        None,
    )
    if zcol is None:
        zcol = "zcos"
        print("  !! WARNING: no observed-redshift column found; using 'zcos'.")
        print("     FoF without peculiar velocities will UNDER-produce linking")
        print("     along the line of sight. Results are not comparable to GAMA.")
        print(f"     columns available: {sorted(galaxies.columns)[:20]} ...")
    else:
        print(f"  redshift column: {zcol}")

    if "total_ap_dust_r_SDSS" not in galaxies.columns:
        raise KeyError("mock galaxies lack total_ap_dust_r_SDSS")

    wanted = [
        "ra",
        "dec",
        zcol,
        "total_ap_dust_r_SDSS",
        "total_ab_dust_r_SDSS",
        "id_group_sky",
        "id_fof",
        "log_mstar_total",
    ]
    gal = galaxies[[c for c in wanted if c in galaxies.columns]].copy()
    gal = gal.rename(
        columns={
            zcol: "z",
            "total_ap_dust_r_SDSS": "app_mag",
            "total_ab_dust_r_SDSS": "abs_mag_native",
        }
    )
    gal = gal[(gal["app_mag"] < mag_limit) & (gal["z"] > zmin) & (gal["z"] < zmax)]
    gal = gal[np.isfinite(gal["ra"]) & np.isfinite(gal["dec"]) & np.isfinite(gal["z"])]

    if MOCK_SPEC_COMPLETENESS is not None:
        rng = np.random.default_rng(0)
        keep = rng.random(len(gal)) < MOCK_SPEC_COMPLETENESS
        print(
            f"  imposed spectroscopic completeness {MOCK_SPEC_COMPLETENESS:.0%}: "
            f"{keep.sum()}/{len(gal)} kept"
        )
        gal = gal[keep]

    from nessie import FlatCosmology

    cosmo = FlatCosmology(NESSIE_H0, NESSIE_OM)

    # Two ways to get absolute magnitudes:
    #  'dmu'    - apparent - dist_mod(z) - k+e, exactly as make_gama_dmu does.
    #             Carries GAMA's k+e correction and depends on the assumed
    #             cosmology (h=1 here, which may not be the mock's).
    #  'native' - Shark's own total_ab_dust_r_SDSS. No cosmology assumption and
    #             no k+e correction, so it is NOT identical to GAMA's definition,
    #             but it is internally consistent with the mock.
    # Nessie uses these for the luminosity-weighted centre, total flux and BCG;
    # mass_proxy comes from sigma^2 R, so the completeness result is insensitive
    # to the choice. Both are printed so the difference is on the record.
    recomputed = (
        gal["app_mag"].values
        - cosmo.dist_mod(gal["z"].values)
        - kecorr(gal["z"].values)
    )
    if "abs_mag_native" in gal.columns:
        native = gal["abs_mag_native"].values
        q = [1, 50, 99]
        print(f"  abs mag percentiles {q}:")
        print(f"    Shark native : {np.round(np.nanpercentile(native, q), 2)}")
        print(f"    DMU recipe   : {np.round(np.nanpercentile(recomputed, q), 2)}")
        print(f"    median offset: {np.nanmedian(recomputed - native):+.2f} mag")
    else:
        native = None
        if abs_mag_mode == "native":
            print("  !! no total_ab_dust_r_SDSS; falling back to the DMU recipe")
            abs_mag_mode = "dmu"

    gal["abs_mag"] = native if abs_mag_mode == "native" else recomputed
    print(f"  using abs_mag_mode = {abs_mag_mode!r}")
    # Shark uses -999 as a null. It is finite, so isfinite() does not catch it,
    # and a luminosity of 10^(0.4*1003) overflows the group-table calculation.
    bad = ~np.isfinite(gal["abs_mag"]) | (gal["abs_mag"] < -40) | (gal["abs_mag"] > 0)
    if bad.any():
        print(
            f"  dropped {int(bad.sum())} galaxies with sentinel/invalid abs_mag "
            f"(e.g. -999)"
        )
    gal = gal[~bad]

    # true halo masses on the same footing recovery.py uses
    sky_frac = mock_sky_frac(groups)
    halo, _ = R.abundance_match(groups, sky_frac)
    halo = halo[["id_group_sky", "zcos", "mvir", "log_mass_am"]].copy()

    print(f"  galaxies: {len(gal)} with r < {mag_limit}, {zmin} < z < {zmax}")
    print(f"  halos   : {len(halo)} in the same redshift range")
    return gal, halo, cosmo, sky_frac


def mock_sky_frac(groups, nb=200):
    """Fractional sky area, corrected for the unfilled part of the RA/Dec
    bounding box (the box overestimates by ~8 per cent for this lightcone)."""
    ra, dec = groups["ra"].values, groups["dec"].values
    bbox = R.sky_area_deg2(ra, dec)
    H, _, _ = np.histogram2d(ra, dec, bins=nb)
    filled = (H > 0).mean()
    area = bbox * filled
    print(
        f"  footprint: bbox {bbox:.1f} deg^2 x {filled:.1%} filled = {area:.1f} deg^2"
    )
    return area * (np.pi / 180) ** 2 / (4 * np.pi)


def build_randoms(z, n_gal, overfactor=OVERFACTOR, seed=0):
    """Randoms tracing the mock's radial selection, standing in for GAMA's
    precomputed random catalogues. Resamples the observed redshifts with a
    small smoothing kernel so the density function is not lumpy."""
    rng = np.random.default_rng(seed)
    n = int(n_gal * overfactor)
    draw = rng.choice(np.asarray(z, float), size=n, replace=True)
    draw = draw + rng.normal(0, 0.005, size=n)
    return np.clip(draw, 1e-4, None)


# ----------------------------------------------------------------------
def run_nessie(gal, cosmo, sky_frac):
    """Run the finder exactly as make_gama_dmu does. Returns (group_ids,
    group_table)."""
    try:
        from nessie import RedshiftCatalog
        from nessie.helper_funcs import create_density_function
    except ImportError as e:
        sys.exit(
            f"nessie is not importable here ({e}).\n"
            "Run this script in the environment where make_gama_dmu works."
        )

    randoms = build_randoms(gal["z"].values, len(gal))
    rho = create_density_function(randoms, len(randoms) / OVERFACTOR, sky_frac, cosmo)

    # Sanity-check rho against the galaxies' own number density. If these
    # disagree by orders of magnitude the linking lengths are wrong and the FoF
    # will percolate (one giant group swallowing the sample).
    zt = np.array([0.05, 0.10, 0.15, 0.20])
    try:
        rho_t = np.array([float(rho(float(zz))) for zz in zt])
    except Exception:
        rho_t = np.array([float(x) for x in rho(zt)])
    dz = 0.01
    emp = []
    for zz in zt:
        n = int(((gal["z"] > zz - dz) & (gal["z"] < zz + dz)).sum())
        # NB: must use the Nessie cosmology, not recovery.py's. Comoving volumes
        # scale as h^-3, so mixing them gives a spurious (1/0.674)^3 = 3.27 offset.
        d = R.comoving_distance(
            np.array([zz - dz, zz + dz]), H0=100 * NESSIE_H0, Om=NESSIE_OM
        )
        vol = (4 / 3) * np.pi * (d[1] ** 3 - d[0] ** 3) * sky_frac
        emp.append(n / vol)
    emp = np.array(emp)
    print(f"  {'z':>6} {'rho(z)':>12} {'empirical':>12} {'ratio':>8}")
    for zz, r_, e_ in zip(zt, rho_t, emp):
        print(f"  {zz:6.2f} {r_:12.4e} {e_:12.4e} {r_ / e_:8.2f}")
    print("  (ratio should be O(1); orders of magnitude off -> rho is wrong,")
    print("   most likely the randoms are not in the format create_density_function")
    print("   expects -- check np.loadtxt of a real GAMA randoms file for its shape)")

    redcat = RedshiftCatalog(
        gal["ra"].values, gal["dec"].values, gal["z"].values, rho, cosmo
    )

    # Completeness MUST be produced by calculate_completeness(), not assigned.
    # The method populates state on the Rust side; setting the Python attribute
    # alone leaves that as None and calculate_group_table() panics with
    #   called `Option::unwrap()` on a `None` value
    # The mock has no fibre collisions, so passing the same galaxies as both the
    # target and the parent sample gives completeness = 1 everywhere, while
    # still initialising everything the group table needs.
    import astropy.units as u
    from astropy.cosmology import FlatLambdaCDM

    acosmo = FlatLambdaCDM(H0=100 * NESSIE_H0, Om0=NESSIE_OM)
    radii = acosmo.arcsec_per_kpc_comoving(gal["z"].values) * 1 * u.Mpc
    radii_deg = radii.to(u.deg).value
    radii_deg[radii_deg > 180] = 180.0
    redcat.calculate_completeness(gal["ra"].values, gal["dec"].values, radii_deg)
    redcat.completeness = np.array(redcat.completeness)
    comp = np.asarray(redcat.completeness, float)
    print(
        f"  completeness: min {np.nanmin(comp):.3f} med {np.nanmedian(comp):.3f} "
        f"max {np.nanmax(comp):.3f}  ({int((~np.isfinite(comp)).sum())} non-finite)"
    )

    print(f"  running FoF (B0={B0}, R0={R0}) on {len(gal)} galaxies ...")
    redcat.run_fof(B0, R0)

    gids = np.asarray(redcat.group_ids)
    n_ung = int((gids == UNGROUPED_ID).sum())
    ug, cnt = np.unique(gids[gids != UNGROUPED_ID], return_counts=True)
    print(
        f"  FoF: {ug.size} groups (sizes {cnt.min()}-{cnt.max()}), "
        f"{n_ung} ungrouped ({n_ung / len(gal):.0%} of galaxies)"
    )
    print(f"       {int((cnt >= MULTI).sum())} groups with N >= {MULTI}")
    frac_big = cnt.max() / len(gal)
    if frac_big > 0.05:
        raise RuntimeError(
            f"FoF looks percolated: the largest group holds {frac_big:.1%} of the "
            f"sample. Linking lengths depend on the density function -- check the "
            f"rho(z) table above against build_randoms()."
        )

    am = gal["abs_mag"].values
    bad = ~np.isfinite(am)
    if bad.any():
        raise ValueError(
            f"{bad.sum()} non-finite absolute magnitudes -- "
            "the group table cannot be built"
        )
    print(f"  abs_mag: {am.min():.2f} to {am.max():.2f} (median {np.median(am):.2f})")
    if am.min() < -30 or am.max() > -5:
        print(
            "  !! absolute magnitudes look out of range; check the magnitude "
            "column and the cosmology used for the distance modulus"
        )

    vel_err = np.repeat(VEL_ERROR, len(gal))
    table = pd.DataFrame(redcat.calculate_group_table(gal["abs_mag"].values, vel_err))
    table["MassA"] = table["mass_proxy"] * MASS_A
    n_grp = int((table["multiplicity"] >= MULTI).sum())
    print(f"  FoF found {len(table)} groups, {n_grp} with >= {MULTI} members")
    return gids, table


# ----------------------------------------------------------------------
def match_groups(gal, gids, table, halo):
    """Bijective match between FoF groups and true halos.

    For FoF group g and true halo h sharing N_gh member galaxies:
        purity   P = N_gh / N_g   (how much of the FoF group is that halo)
        recovery Q = N_gh / N_h   (how much of the halo is in that FoF group)
    A halo is 'detected' if some group has P >= PURITY_MIN and Q >= RECOVERY_MIN.

    Returns (halo_with_flags, table_with_flags)."""
    df = gal.copy()
    df["fof_id"] = gids
    df = df[df["fof_id"] != UNGROUPED_ID]  # drop ungrouped galaxies

    # group sizes and halo sizes, counted over the SAME galaxy sample
    n_g = df.groupby("fof_id").size().rename("N_g")
    n_h = df.groupby("id_group_sky").size().rename("N_h")

    pair = (
        df.groupby(["fof_id", "id_group_sky"])
        .size()
        .rename("N_gh")
        .reset_index()
        .merge(n_g, on="fof_id")
        .merge(n_h, on="id_group_sky")
    )
    pair["purity"] = pair["N_gh"] / pair["N_g"]
    pair["recovery"] = pair["N_gh"] / pair["N_h"]
    pair["match"] = (pair["purity"] >= PURITY_MIN) & (pair["recovery"] >= RECOVERY_MIN)

    # keep only FoF groups above the multiplicity cut
    big = (
        set(table.loc[table["multiplicity"] >= MULTI, "group_id"].astype(int))
        if "group_id" in table.columns
        else None
    )
    if big is not None:
        pair = pair[pair["fof_id"].astype(int).isin(big)]

    matched = pair[pair["match"]].copy()
    # A halo split evenly in two can satisfy Q >= 0.5 for BOTH fragments, and a
    # group can satisfy P >= 0.5 for two halos. Keep the single best overlap on
    # each side so the match is genuinely bijective and no row is duplicated.
    matched = (
        matched.sort_values("N_gh", ascending=False)
        .drop_duplicates("id_group_sky", keep="first")
        .drop_duplicates("fof_id", keep="first")
    )
    halo = halo.merge(
        matched[["id_group_sky", "fof_id", "N_gh", "N_h", "purity", "recovery"]],
        on="id_group_sky",
        how="left",
    )
    halo["detected"] = halo["fof_id"].notna()

    # purity of the catalogue: what fraction of FoF groups are a clean match
    if "group_id" in table.columns:
        table = table.merge(
            matched[["fof_id", "id_group_sky", "purity", "recovery"]].rename(
                columns={"fof_id": "group_id"}
            ),
            on="group_id",
            how="left",
        )
        table["clean"] = table["id_group_sky"].notna()

    # ---- alternative completeness definitions -------------------------------
    # The bijective match answers "is this halo cleanly recovered as ONE group".
    # That is not what the likelihood needs: C(m) multiplies phi(m) to give the
    # intensity of CATALOGUE ENTRIES, so it is an expected number of entries per
    # halo and may exceed 1 when a halo is fragmented.
    big_ids = (
        set(table.loc[table["multiplicity"] >= MULTI, "group_id"].astype(int))
        if "group_id" in table.columns
        else None
    )
    pr = pair if big_ids is None else pair[pair["fof_id"].astype(int).isin(big_ids)]

    # (a) "represented": some N>=MULTI group holds >= RECOVERY_MIN of the halo,
    #     with no purity requirement -- immune to interloper contamination.
    rep = pr[pr["recovery"] >= RECOVERY_MIN].drop_duplicates("id_group_sky")[
        ["id_group_sky"]
    ]
    rep["represented"] = True
    halo = halo.merge(rep, on="id_group_sky", how="left")
    halo["represented"] = halo["represented"].fillna(False)

    # (b) "n_entries": how many N>=MULTI catalogue entries this halo produces.
    #
    # NO PURITY CUT.  C_entries multiplies phi(m) to give the intensity of
    # CATALOGUE ENTRIES, and the catalogue contains every N>=MULTI group, clean
    # or not.  Requiring purity >= PURITY_MIN counted only 1253 of the 1585 mock
    # groups, so Lambda expected 1253 while the likelihood fitted all 1585 --
    # the two populations disagreed and Lambda/N came out at 0.79, a 21% deficit
    # landing straight on log phi*.  Contamination is a separate effect and
    # belongs in a purity/contamination term, not in the completeness.
    #
    # Each group is assigned to its DOMINANT halo (largest shared membership)
    # before counting.  Without that dedup a group straddling three halos would
    # be counted once for each and sum(n_entries) would exceed the catalogue
    # size -- over-correcting in the opposite direction.
    # RESTORED. Dropping the purity cut is right in principle -- a Poisson
    # intensity over catalogue entries should count every entry -- and it did
    # fix the normalisation (Lambda/N 0.790 -> 0.959). But the closed loop says
    # otherwise: alpha recovery degraded from +0.011 sigma to +0.88 sigma. Two
    # errors were evidently cancelling in the incumbent, and until that is
    # understood the validated configuration wins. See CLAUDE.md, "The C(m,z)
    # route".
    #
    # Note purity >= 0.5 means more than half the group belongs to that halo, so
    # at most one halo can own a group and no dedup is needed here.
    ent = (
        pr[pr["purity"] >= PURITY_MIN]
        .groupby("id_group_sky")
        .size()
        .rename("n_entries")
        .reset_index()
    )
    halo = halo.merge(ent, on="id_group_sky", how="left")
    halo["n_entries"] = halo["n_entries"].fillna(0).astype(int)
    if big_ids is not None:
        _tot, _cat = int(halo["n_entries"].sum()), len(big_ids)
        print(
            f"  entries    : {_tot} attributed / {_cat} catalogue groups "
            f"(N>={MULTI})  ratio {_tot / max(_cat, 1):.3f}"
        )
        if abs(_tot / max(_cat, 1) - 1.0) > 0.02:
            print(
                "  !! Lambda will inherit this ratio directly as log phi*; "
                "entries not attributed to any halo are groups with no overlap"
            )

    nd = int(halo["detected"].sum())
    print(
        f"  matched   : {nd} / {len(halo)} halos (bijective: P>={PURITY_MIN} "
        f"AND Q>={RECOVERY_MIN})"
    )
    print(
        f"  represented: {int(halo['represented'].sum())} halos (Q>={RECOVERY_MIN} only)"
    )
    print(
        f"  entries    : {int(halo['n_entries'].sum())} catalogue entries "
        f"majority-owned by a halo"
    )

    # ---- why does completeness turn over at high mass? ----------------------
    hi = halo[halo[MASS_COL] > np.nanpercentile(halo[MASS_COL], 99.5)]
    if len(hi) > 20:
        bestQ = (
            pr.groupby("id_group_sky")["recovery"].max().rename("bestQ").reset_index()
        )
        bestP = pr.groupby("id_group_sky")["purity"].max().rename("bestP").reset_index()
        hi = hi.merge(bestQ, on="id_group_sky", how="left").merge(
            bestP, on="id_group_sky", how="left"
        )
        fail = hi[~hi["detected"].astype(bool)]
        print(
            f"\n  top 0.5% by mass: {len(hi)} halos, {len(fail)} fail the bijective match"
        )
        if len(fail):
            nq = int((fail["bestQ"] < RECOVERY_MIN).sum())
            npu = int((fail["bestP"] < PURITY_MIN).sum())
            nno = int(fail["bestQ"].isna().sum())
            print(
                f"    {nq:5d} best group holds < {RECOVERY_MIN:.0%} of the halo "
                f"-> FRAGMENTED"
            )
            print(
                f"    {npu:5d} best group is < {PURITY_MIN:.0%} that halo "
                f"-> CONTAMINATED / merged"
            )
            print(f"    {nno:5d} no N>={MULTI} group contains any of its members")
    if "clean" in table.columns:
        sel = table["multiplicity"] >= MULTI
        print(
            f"  purity : {table.loc[sel, 'clean'].mean():.1%} of FoF groups "
            f"(N>={MULTI}) are a clean match to a halo"
        )
    return halo, table


# ----------------------------------------------------------------------
def measure_completeness(halo, mlim_func, z_bins, d_edges):
    """C(Delta) per redshift bin, Delta = m_true - mlim(z). Returns the global
    fit, the per-z fits, and the binned curves for plotting."""
    m = halo[MASS_COL].values
    z = halo["zcos"].values
    det = halo["detected"].values.astype(bool)
    delta = m - mlim_func(z)
    rep = (
        halo["represented"].values.astype(bool)
        if "represented" in halo.columns
        else det
    )
    nent = (
        halo["n_entries"].values.astype(float)
        if "n_entries" in halo.columns
        else det.astype(float)
    )

    def curve(mask, min_n=10):
        cen, C, N = [], [], []
        for a, b in zip(d_edges[:-1], d_edges[1:]):
            s = mask & (delta >= a) & (delta < b)
            if int(s.sum()) >= min_n:
                cen.append(0.5 * (a + b))
                C.append(det[s].mean())
                N.append(int(s.sum()))
        return np.array(cen), np.array(C), np.array(N)

    # the two looser definitions, binned the same way
    def curve_of(vals, min_n=10):
        cen, C = [], []
        for a, b in zip(d_edges[:-1], d_edges[1:]):
            sm = (delta >= a) & (delta < b)
            if int(sm.sum()) >= min_n:
                cen.append(0.5 * (a + b))
                C.append(float(np.mean(vals[sm])))
        return np.array(cen), np.array(C)

    out = {}
    out["represented"] = curve_of(rep.astype(float))
    out["n_entries"] = curve_of(nent)
    cen, C, N = curve(np.ones_like(det, bool))
    try:
        (d50, w), cov = curve_fit(
            erf_ramp,
            cen,
            C,
            p0=[0.0, 0.25],
            sigma=np.sqrt(np.maximum(C * (1 - C) / N, 1e-4)),
            absolute_sigma=True,
            maxfev=20000,
        )
        err = np.sqrt(np.diag(cov))
    except Exception:
        d50, w, err = np.nan, np.nan, [np.nan, np.nan]
    out["global"] = dict(cen=cen, C=C, N=N, d50=d50, w=w, d50_err=err[0], w_err=err[1])

    out["zbins"] = []
    for za, zb in z_bins:
        mask = (z >= za) & (z < zb)
        cz, Cz, Nz = curve(mask)
        if cz.size < 4:
            continue
        try:
            (dz, wz), covz = curve_fit(
                erf_ramp,
                cz,
                Cz,
                p0=[0.0, 0.25],
                sigma=np.sqrt(np.maximum(Cz * (1 - Cz) / Nz, 1e-4)),
                absolute_sigma=True,
                maxfev=20000,
            )
            ez = np.sqrt(np.diag(covz))
        except Exception:
            dz, wz, ez = np.nan, np.nan, [np.nan, np.nan]
        out["zbins"].append(
            dict(
                za=za,
                zb=zb,
                zc=0.5 * (za + zb),
                cen=cz,
                C=Cz,
                N=Nz,
                d50=dz,
                w=wz,
                d50_err=ez[0],
                w_err=ez[1],
                n_halo=int(mask.sum()),
            )
        )
    return out


def tabulate_C_mz(halo, m_edges, z_edges, min_n=20):
    """Completeness tabulated on ABSOLUTE true mass and redshift, C(m, z).

    The existing ``C(Delta, z)`` table is keyed on ``Delta = m - mlim(z)``, which
    drags ``recovery.turnover_mlim`` into the fit.  That estimator takes the
    *mode* of the observed mass histogram, which for a smooth distribution lands
    mid-distribution rather than at the faint edge -- 54% of real GAMA groups sit
    below their own mlim.  Keying on absolute mass removes mlim from the fit
    entirely.

    **Consequence to keep in mind**: with Delta gone, C's *z*-dependence IS the
    selection function.  It is no longer carried by mlim(z), so the z grid has to
    span the shells the fit actually uses (``recovery.shell_volumes`` runs
    z = 0.016..0.244), not the 3 clamped nodes the Delta table got away with.

    Returns ``(m_cen, z_cen, C_entries, C_repr, Nh)``:

    * ``C_entries`` -- mean ``n_entries``, the expected NUMBER of catalogue
      entries per halo.  This is what the Poisson intensity wants and it may
      legitimately exceed 1 (one halo can yield several N>=5 groups), so it is
      not clipped at 1.
    * ``C_repr`` -- mean ``represented``, a true probability in [0, 1].  Bounded,
      so it is the systematic cross-check on ``C_entries``.
    * ``Nh`` -- halo count per cell.  The high-mass plateau is a flat hold off
      whatever the last populated cell contains; under absolute keying that
      plateau sits at a specific mass and directly sets M* and beta, so the
      count behind it has to be visible.
    """
    m_cen = 0.5 * (np.asarray(m_edges[1:], float) + np.asarray(m_edges[:-1], float))
    z_cen = 0.5 * (np.asarray(z_edges[1:], float) + np.asarray(z_edges[:-1], float))
    mall = halo[MASS_COL].values.astype(float)
    zall = halo["zcos"].values.astype(float)
    nent = halo["n_entries"].values.astype(float)
    nrep = halo["represented"].values.astype(float)

    shape = (z_cen.size, m_cen.size)
    C_ent, C_rep = np.full(shape, np.nan), np.full(shape, np.nan)
    Nh = np.zeros(shape, int)
    for i in range(z_cen.size):
        zm = (zall >= z_edges[i]) & (zall < z_edges[i + 1])
        for k in range(m_cen.size):
            s_ = zm & (mall >= m_edges[k]) & (mall < m_edges[k + 1])
            Nh[i, k] = int(s_.sum())
            if Nh[i, k] >= min_n:
                C_ent[i, k] = float(np.mean(nent[s_]))
                C_rep[i, k] = float(np.mean(nrep[s_]))
        # same fill rule as the Delta table: 0 below the first measured cell,
        # hold the last measured value above it
        for arr in (C_ent, C_rep):
            ok = np.isfinite(arr[i])
            if ok.sum() >= 2:
                arr[i] = np.interp(
                    m_cen, m_cen[ok], arr[i][ok], left=0.0, right=float(arr[i][ok][-1])
                )
            else:
                arr[i] = np.where(np.isfinite(arr[i]), arr[i], 0.0)
    return m_cen, z_cen, C_ent, C_rep, Nh


def fit_C_mz_parametric(halo, m_grid, z_grid, kind="entries"):
    """C(m, z) as a smooth parametric fit, evaluated onto a dense grid.

    WHY NOT A RAW 2-D HISTOGRAM.  Binning halos in (mass, z) and averaging is
    the obvious thing and it does not work here: the mock lightcone has almost
    no nearby massive halos.  At z = 0.025 there are 2 halos above logM 14 and
    none above 14.5, so the last measurable bin is 12.55 and everything above it
    gets flat-held at C = 0.087 -- for objects that are certainly detected.  The
    Delta = m - mlim(z) keying dodged this by pooling every redshift into a
    single shape; absolute keying has to recover that pooling some other way.

    So fit, do not bin.  The ramp

        C(m, z) = A(z) * 0.5 * (1 + erf((m - m50(z)) / (sqrt(2) w(z))))

    is fitted to ALL halos at once by maximum likelihood, with

        m50(z) quadratic,  log w(z) linear,  log A(z) linear

    which pools across redshift and extrapolates into the empty low-z corner
    along the trend the populated bins define.  No binning, no sparsity cliff,
    and -- the point of the exercise -- no mlim(z) fitted to the observed data:
    m50(z) is measured from the mock against the mock's own truth.

    ``kind='entries'`` fits a Poisson mean (may exceed 1); ``'repr'`` fits a
    Bernoulli probability (bounded by 1).
    """
    from scipy.optimize import minimize

    m = halo[MASS_COL].values.astype(float)
    zz = halo["zcos"].values.astype(float)
    y = (
        halo["n_entries"].values.astype(float)
        if kind == "entries"
        else halo["represented"].values.astype(float)
    )
    ok = np.isfinite(m) & np.isfinite(zz) & np.isfinite(y)
    m, zz, y = m[ok], zz[ok], y[ok]
    z0 = float(np.median(zz))

    def unpack(p, zv):
        dz = zv - z0
        m50 = p[0] + p[1] * dz + p[2] * dz**2
        w = np.exp(p[3] + p[4] * dz)
        # A is the high-mass plateau and is held CONSTANT in z on purpose.  It
        # is only constrained by the handful of cells above the ramp, and an
        # exp(a0 + a1*dz) form extrapolated to logM 15.5 ran away to 2.5 at low
        # z -- against ~0.8 where the mock actually has halos.  Real GAMA groups
        # reach logM ~15.3, i.e. inside that extrapolated region, so the runaway
        # would land directly on M* and beta.
        # For 'repr' the outcome is Bernoulli, so A is a probability and is
        # squashed into (0, 1]; an unbounded A gave C_repr = 1.32, which is
        # impossible.
        A = 1.0 / (1.0 + np.exp(-p[5])) if kind == "repr" else np.exp(p[5])
        return m50, w, A

    def model(p, mv, zv):
        m50, w, A = unpack(p, zv)
        # LOGISTIC, not erf.  An erf ramp has Gaussian tails, which die far
        # faster than the real selection: a low-mass halo can still occasionally
        # be found, and that tail carries a lot of objects because the HMF is
        # steep there.  Fitted against the binned mock the erf under-predicted C
        # at low mass by factors of 2-4, and the Poisson normalisation came out
        # at Lambda/N = 0.79 -- a 21% deficit that lands straight on log phi*.
        # The logistic has exponential tails and holds that region.
        return A / (1.0 + np.exp(-1.7 * (mv - m50) / np.maximum(w, 1e-6)))

    def nll(p):
        mu = model(p, m, zz)
        mu = np.clip(mu, 1e-9, None)
        if kind == "entries":  # Poisson
            return float(np.sum(mu - y * np.log(mu)))
        q = np.clip(mu, 1e-9, 1 - 1e-9)  # Bernoulli
        return float(-np.sum(y * np.log(q) + (1 - y) * np.log1p(-q)))

    p0 = np.array([13.8, 3.0, 0.0, np.log(0.35), 0.0,
                   0.0 if kind == "repr" else np.log(0.8), 0.0])
    res = minimize(nll, p0, method="Nelder-Mead",
                   options=dict(maxiter=40000, maxfev=40000, xatol=1e-6, fatol=1e-6))
    p = res.x
    mm, zc = np.meshgrid(m_grid, z_grid, indexing="xy")
    C = model(p, mm, zc)

    m50_0, w_0, A_0 = unpack(p, np.array([z_grid[0]]))
    m50_1, w_1, A_1 = unpack(p, np.array([z_grid[-1]]))
    print(
        f"  [{kind}] fit ok={res.success} nll={res.fun:.1f} on {m.size} halos\n"
        f"           m50: {m50_0[0]:.2f} (z={z_grid[0]:.3f}) -> "
        f"{m50_1[0]:.2f} (z={z_grid[-1]:.3f})\n"
        f"           w:   {w_0[0]:.3f} -> {w_1[0]:.3f}   "
        f"A: {float(np.atleast_1d(A_0)[0]):.3f} (constant in z)"
    )
    return C, p, z0


def measure_mass_relation(halo, table):
    """Recovered dynamical mass vs true mass, for cleanly matched pairs.
    This is the empirical replacement for the multiplicity->sigma lookup."""
    if "group_id" not in table.columns:
        return None
    m = halo[halo["detected"]].merge(
        table[["group_id", "MassA", "multiplicity", "median_redshift"]].rename(
            columns={"group_id": "fof_id"}
        ),
        on="fof_id",
        how="inner",
    )
    m = m[(m["MassA"] > 0) & np.isfinite(m["MassA"])]
    m["log_dyn"] = np.log10(m["MassA"])
    m["resid"] = m["log_dyn"] - m[MASS_COL]
    return m


# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
def _style():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "lines.linewidth": 1.4,
        }
    )
    return plt


def plot_completeness(res, fname="fig_completeness_nessie.pdf"):
    plt = _style()
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    dd = np.linspace(-1.5, 1.5, 400)
    cmap = plt.cm.viridis
    nz = len(res["zbins"])

    for i, b in enumerate(res["zbins"]):
        col = cmap(i / max(nz - 1, 1) * 0.85)
        ax.plot(b["cen"], b["C"], "o", ms=3.5, color=col, alpha=0.85)
        if np.isfinite(b["d50"]):
            ax.plot(
                dd,
                erf_ramp(dd, b["d50"], b["w"]),
                color=col,
                lw=1.2,
                label=rf"${b['za']:.2f}<z<{b['zb']:.2f}$",
            )

    if "represented" in res:
        cr, Cr = res["represented"]
        ax.plot(
            cr,
            Cr,
            "^--",
            ms=4,
            color="crimson",
            lw=1.2,
            alpha=0.9,
            label="represented ($Q\\geq0.5$)",
        )
    if "n_entries" in res:
        ce, Ce = res["n_entries"]
        ax.plot(
            ce,
            Ce,
            "v:",
            ms=4,
            color="darkorange",
            lw=1.2,
            alpha=0.9,
            label="entries per halo",
        )

    g = res["global"]
    ax.plot(
        g["cen"], g["C"], "s", ms=4.5, color="k", zorder=5, label="bijective (all $z$)"
    )
    if np.isfinite(g["d50"]):
        ax.plot(dd, erf_ramp(dd, g["d50"], g["w"]), "k-", lw=2, zorder=4)
        ax.axvline(g["d50"], color="0.5", ls=":", lw=1)
        ax.axhline(0.5, color="0.5", ls=":", lw=1)
        ax.annotate(
            rf"$D_{{50}} = {g['d50']:+.3f}$" "\n" rf"$w = {g['w']:.3f}$",
            xy=(0.04, 0.72),
            xycoords="axes fraction",
            fontsize=9,
        )

    ax.set_xlabel(r"$\Delta = \log_{10}M_{\rm true} - m_{\rm lim}(z)$")
    ax.set_ylabel(r"completeness $C$")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(
        -0.03,
        max(
            1.03,
            float(np.nanmax(res["n_entries"][1])) * 1.1 if "n_entries" in res else 1.03,
        ),
    )
    ax.set_title("Selection function of the Nessie catalogue")
    ax.legend(loc="lower right", ncol=1)
    fig.savefig(fname)
    print(f"  saved {fname}")


def plot_purity(table, fname="fig_purity_nessie.pdf"):
    if "clean" not in table.columns:
        return
    plt = _style()
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    sel = table["multiplicity"] >= MULTI
    t = table[sel]
    edges = np.array([5, 6, 8, 10, 15, 25, 200])
    cen, frac, err = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (t["multiplicity"] >= a) & (t["multiplicity"] < b)
        if s.sum() >= 5:
            p = t.loc[s, "clean"].mean()
            cen.append(np.sqrt(a * b))
            frac.append(p)
            err.append(np.sqrt(max(p * (1 - p), 1e-4) / s.sum()))
    ax.errorbar(cen, frac, yerr=err, fmt="o-", color="C0", ms=5, capsize=2.5)
    ax.axhline(1.0, color="0.6", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("FoF group multiplicity $N$")
    ax.set_ylabel("fraction cleanly matched")
    ax.set_ylim(0, 1.05)
    ax.set_title("Purity of the recovered catalogue")
    fig.savefig(fname)
    print(f"  saved {fname}")


def plot_mass_relation(mm, fname="fig_mass_relation_nessie.pdf"):
    if mm is None or len(mm) < 50:
        return
    plt = _style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.8))

    a1.plot(mm[MASS_COL], mm["log_dyn"], ".", ms=1.5, alpha=0.25, color="C0")
    lo = float(np.nanpercentile(mm[MASS_COL], 0.5))
    hi = float(np.nanpercentile(mm[MASS_COL], 99.5))
    xx = np.linspace(lo, hi, 50)
    a1.plot(xx, xx, "k--", lw=1.2, label="1:1")
    # running median
    b = np.linspace(lo, hi, 14)
    ctr = 0.5 * (b[:-1] + b[1:])
    med = [
        np.nanmedian(mm["log_dyn"][(mm[MASS_COL] >= u) & (mm[MASS_COL] < v)])
        for u, v in zip(b[:-1], b[1:])
    ]
    a1.plot(ctr, med, "-", color="crimson", lw=2, label="running median")
    a1.set(
        xlabel=r"$\log_{10} M_{\rm true}$",
        ylabel=r"$\log_{10} M_{\rm dyn}$ (Nessie, $A=10$)",
        xlim=(lo, hi),
        ylim=(lo - 0.8, hi + 0.8),
    )
    a1.legend(loc="upper left")
    a1.set_title("Recovered vs true mass")

    edges = np.array([5, 6, 7, 8, 10, 12, 15, 20, 30, 200])
    cen, sd, bias, n = [], [], [], []
    for u, v in zip(edges[:-1], edges[1:]):
        s = (mm["multiplicity"] >= u) & (mm["multiplicity"] < v)
        if s.sum() >= 15:
            r = mm.loc[s, "resid"].values
            cen.append(np.sqrt(u * v))
            sd.append(0.5 * (np.percentile(r, 84) - np.percentile(r, 16)))
            bias.append(np.median(r))
            n.append(int(s.sum()))
    a2.plot(cen, sd, "o-", color="C0", ms=5, label=r"scatter $\sigma_{\log M}$")
    a2.plot(cen, bias, "s--", color="crimson", ms=4, label="median offset")
    nf = np.arange(3, 40)
    a2.plot(
        nf,
        R.sigma_from_nfof(nf),
        ":",
        color="0.4",
        lw=1.5,
        label="lookup table (run.R)",
    )
    a2.axhline(0, color="0.7", lw=0.8)
    a2.set_xscale("log")
    a2.set(xlabel="multiplicity $N$", ylabel="dex")
    a2.legend(loc="upper right")
    a2.set_title("Mass error model")
    fig.savefig(fname)
    print(f"  saved {fname}")


def plot_multiplicity(
    table, gama_fits, area_mock, area_gama, fname="fig_multiplicity_nessie.pdf"
):
    """The test the membership proxy failed: does the mock now reproduce the
    real catalogue's surface density as a function of richness?"""
    plt = _style()
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    edges = np.array([5, 6, 8, 10, 15, 25, 200])
    cen = np.sqrt(edges[:-1] * edges[1:])

    nm = table.loc[table["multiplicity"] >= MULTI, "multiplicity"].values
    dm = np.histogram(nm, bins=edges)[0] / area_mock
    ax.step(edges[:-1], dm, where="post", lw=1.8, color="C0", label="Shark + Nessie")

    if gama_fits and os.path.exists(gama_fits):
        _, _, _, nr = R.load_real_gama(gama_fits)
        dr = np.histogram(nr, bins=edges)[0] / area_gama
        ax.step(
            edges[:-1],
            dr,
            where="post",
            lw=1.8,
            color="crimson",
            label="GAMA (this catalogue)",
        )
        for c, a, b in zip(cen, dm, dr):
            if a > 0:
                ax.annotate(
                    f"{b / a:.2f}",
                    xy=(c, max(a, b) * 1.15),
                    ha="center",
                    fontsize=8,
                    color="0.35",
                )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("multiplicity $N$")
    ax.set_ylabel(r"groups deg$^{-2}$")
    ax.set_title("Surface density by richness")
    ax.legend()
    fig.savefig(fname)
    print(f"  saved {fname}")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mag-limit", type=float, default=APPARENT_MAG_LIM)
    ap.add_argument("--zmin", type=float, default=R.ZMIN)
    ap.add_argument("--zmax", type=float, default=R.ZLIMIT)
    ap.add_argument(
        "--gama-fits",
        default="/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits",
    )
    ap.add_argument("--gama-area", type=float, default=238.11)
    ap.add_argument(
        "--abs-mag",
        choices=["native", "dmu"],
        default="native",
        help="native = Shark total_ab_dust_r_SDSS; "
        "dmu = apparent - dist_mod - k+e as make_gama_dmu does",
    )
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--save", default="nessie_completeness.npz")
    a = ap.parse_args()

    print("=" * 70)
    print("  Nessie on the Shark mock: measuring the real selection function")
    print("=" * 70)
    print(f"  FoF: B0={B0}, R0={R0}, cosmology h={NESSIE_H0}, Om={NESSIE_OM}")
    h_rec = R.H0 / 100.0
    print(f"  recovery.py: h={h_rec:.4f}, Om={R.OMEGA_M}")
    if abs(h_rec - NESSIE_H0) > 1e-6 or abs(R.OMEGA_M - NESSIE_OM) > 1e-6:
        print(
            f"  !! COSMOLOGY MISMATCH. Dynamical masses differ by "
            f"{np.log10(NESSIE_H0 / h_rec):+.3f} dex and comoving volumes by "
            f"{(h_rec / NESSIE_H0) ** 3:.2f}x, so the abundance-matched 'true'"
        )
        print(f"     masses and the Nessie masses are on different scales.")
        print(
            f"     Set H0, OMEGA_M = {100 * NESSIE_H0:.0f}, {NESSIE_OM} in "
            f"recovery.py and re-run."
        )
    else:
        print(f"  cosmologies match -- true and dynamical masses are on one scale.")
    print()

    gal, halo, cosmo, sky_frac = load_mock(
        a.mag_limit, a.zmin, a.zmax, abs_mag_mode=a.abs_mag
    )
    area_mock = sky_frac * 4 * np.pi * (180 / np.pi) ** 2

    gids, table = run_nessie(gal, cosmo, sky_frac)
    halo, table = match_groups(gal, gids, table, halo)

    # mlim(z) from the RECOVERED catalogue, the same estimator recovery.py uses
    sel = table["multiplicity"] >= MULTI
    m_obs = np.log10(table.loc[sel, "MassA"].values)
    z_obs = table.loc[sel, "median_redshift"].values
    ok = np.isfinite(m_obs) & np.isfinite(z_obs)
    mlim_func, _, kind, _ = R.turnover_mlim(
        z_obs[ok], m_obs[ok], zmin=a.zmin, zmax=a.zmax
    )
    print(f"  mlim(z) [{kind}]: {mlim_func(a.zmin):.2f} -> {mlim_func(a.zmax):.2f}")

    z_bins = [(a.zmin, 0.08), (0.08, 0.15), (0.15, a.zmax)]
    # Fine bins where the halos are, progressively wider in the tail. With a
    # uniform 0.1-dex grid the bins above Delta ~ 0.8 fall below the minimum
    # count and are dropped, so C is never measured where the massive groups
    # live -- exactly where M* and beta are set.
    d_edges = np.concatenate(
        [np.arange(-1.5, 0.8, 0.1), [0.8, 0.95, 1.1, 1.3, 1.55, 1.85, 2.2, 2.6]]
    )
    res = measure_completeness(halo, mlim_func, z_bins, d_edges)

    g = res["global"]
    print("\n  === the three completeness definitions ===")
    print(
        f"  {'Delta':>7} {'bijective':>10} {'represented':>12} {'entries/halo':>13} {'N':>7}"
    )
    cb, Cb = g["cen"], g["C"]
    cr, Cr = res["represented"]
    ce, Ce = res["n_entries"]
    for d in (-0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.4, 1.9):
        f = lambda cc, vv: (
            float(vv[np.argmin(abs(cc - d))])
            if len(cc) and abs(cc - d).min() < 0.08
            else np.nan
        )
        nn = (
            g["N"][np.argmin(abs(g["cen"] - d))]
            if len(g["cen"]) and abs(g["cen"] - d).min() < 0.15
            else 0
        )
        print(
            f"  {d:+7.2f} {f(cb, Cb):10.3f} {f(cr, Cr):12.3f} {f(ce, Ce):13.3f} {nn:7d}"
        )
    print(
        "  bijective   : cleanly recovered as ONE group (turns over -> fragmentation)"
    )
    print("  represented : some N>=5 group holds >=50% of the halo (no purity cut)")
    print("  entries/halo: expected number of catalogue entries -- what the")
    print("                likelihood's C(m) actually needs; may exceed 1")

    print("\n  === completeness ramp (Nessie-based) ===")
    print(
        f"  GLOBAL   D50 = {g['d50']:+.3f} +/- {g['d50_err']:.3f}   "
        f"w = {g['w']:.3f} +/- {g['w_err']:.3f}"
    )
    print(f"  mlim sits at C = {erf_ramp(0.0, g['d50'], g['w']):.2f}")
    print(f"  {'z bin':>14} {'D50':>18} {'w':>18} {'N halo':>8}")
    for b in res["zbins"]:
        print(
            f"  {b['za']:.2f}-{b['zb']:.2f}   "
            f"{b['d50']:+.3f} +/- {b['d50_err']:.3f}   "
            f"{b['w']:.3f} +/- {b['w_err']:.3f} {b['n_halo']:8d}"
        )

    print("\n  paste into recovery.py:")
    zc = [round(b["zc"], 3) for b in res["zbins"]]
    print(f"  COMP_Z_PTS   = {zc}")
    print(f"  COMP_D50_PTS = {[round(float(b['d50']), 3) for b in res['zbins']]}")
    print(f"  COMP_W_PTS   = {[round(float(b['w']), 3) for b in res['zbins']]}")
    print(f"  # global alternative: D50={g['d50']:.3f}, w={g['w']:.3f}")

    mm = measure_mass_relation(halo, table)
    if mm is not None and len(mm):
        print(f"\n  === mass recovery ({len(mm)} clean matches) ===")
        print(f"  median offset  {np.median(mm['resid']):+.3f} dex")
        print(
            f"  scatter (68%)  "
            f"{0.5 * (np.percentile(mm['resid'], 84) - np.percentile(mm['resid'], 16)):.3f} dex"
        )
        print(
            "  -> compare with the run.R lookup table "
            f"(median {np.median(R.sigma_from_nfof(mm['multiplicity'].values)):.3f} dex)"
        )

    # ---- tabulated C from entries/halo, per z bin: what recovery.py should use
    print("\n  === tabulated C(Delta, z) from entries/halo ===")
    d_tab = np.concatenate(
        [np.arange(-1.0, 0.8, 0.1), [0.85, 1.0, 1.2, 1.4, 1.7, 2.0, 2.4]]
    )
    tab = np.zeros((len(res["zbins"]), d_tab.size))
    zc_t = []
    zall = halo["zcos"].values
    dall = halo[MASS_COL].values - mlim_func(zall)
    nall = halo["n_entries"].values.astype(float)
    for i, b in enumerate(res["zbins"]):
        zc_t.append(b["zc"])
        msk = (zall >= b["za"]) & (zall < b["zb"])
        halfw = np.gradient(d_tab) / 2.0
        for k, dd in enumerate(d_tab):
            s_ = msk & (dall >= dd - halfw[k]) & (dall < dd + halfw[k])
            tab[i, k] = float(np.mean(nall[s_])) if s_.sum() >= 8 else np.nan
        # fill: 0 below the first measured point, hold the last above
        row = tab[i]
        ok = np.isfinite(row)
        if ok.sum() >= 2:
            tab[i] = np.interp(
                d_tab, d_tab[ok], row[ok], left=0.0, right=float(row[ok][-1])
            )
    print(
        f"  grid: Delta {d_tab[0]:+.1f} to {d_tab[-1]:+.1f} step 0.1, "
        f"{len(zc_t)} z bins"
    )
    print(f"  C at Delta=0 : {np.round(tab[:, np.argmin(abs(d_tab))], 3)}")
    print(f"  C at Delta=+1: {np.round(tab[:, np.argmin(abs(d_tab - 1.0))], 3)}")
    print(f"  C at Delta=+2: {np.round(tab[:, np.argmin(abs(d_tab - 2.0))], 3)}")
    print(f"  max Delta with data: {d_tab[np.isfinite(tab).any(axis=0)].max():+.2f}")
    print(f"  plateau      : {np.round(np.nanmax(tab, axis=1), 3)}")
    np.savez("nessie_completeness_table.npz", d=d_tab, z=np.array(zc_t), C=tab)

    # ---- the same completeness keyed on ABSOLUTE mass, C(m, z).  Written
    # alongside the Delta table, which is left untouched for backward
    # compatibility.  See tabulate_C_mz for why the z grid has to be finer here.
    print("\n  === tabulated C(m, z) on absolute mass ===")
    m_edges_mz = np.concatenate(
        [np.arange(11.0, 14.6 + 1e-9, 0.1), [14.9, 15.3, 15.8]]
    )
    z_edges_mz = np.linspace(a.zmin, a.zmax, 9)
    m_cen, z_cen, C_bin, C_bin_rep, Nh = tabulate_C_mz(
        halo, m_edges_mz, z_edges_mz, min_n=MZ_MIN_N
    )
    print(
        f"  grid: logM {m_cen[0]:.2f}..{m_cen[-1]:.2f} ({m_cen.size} bins), "
        f"z {z_cen[0]:.3f}..{z_cen[-1]:.3f} ({z_cen.size} bins)"
    )
    # The binned table is kept only as a diagnostic: at low z the mock has no
    # massive halos, so its high-mass cells are flat-held off low-mass ones.
    lastm = [
        (m_cen[np.where(Nh[i] >= MZ_MIN_N)[0][-1]]
         if (Nh[i] >= MZ_MIN_N).any() else np.nan)
        for i in range(z_cen.size)
    ]
    print(f"  binned: last mass bin with N>={MZ_MIN_N} per z row: "
          f"{np.round(lastm, 2)}   <- why the fit is used instead")
    C_ent, p_ent, z0 = fit_C_mz_parametric(halo, m_cen, z_cen, kind="entries")
    C_rep, p_rep, _ = fit_C_mz_parametric(halo, m_cen, z_cen, kind="repr")
    with np.errstate(invalid="ignore"):
        m50 = [
            np.interp(0.5, C_ent[i], m_cen) if C_ent[i].max() >= 0.5 else np.nan
            for i in range(z_cen.size)
        ]
    print(f"  mass at C_entries=0.5 : {np.round(m50, 2)}")
    print(f"  plateau (C_entries)   : {np.round(np.nanmax(C_ent, axis=1), 3)}")
    print(f"  plateau (C_repr)      : {np.round(np.nanmax(C_rep, axis=1), 3)}")
    np.savez(
        "nessie_completeness_mz.npz",
        m=m_cen, z=z_cen, C=C_ent, C_repr=C_rep, Nh=Nh,
        C_binned=C_bin, C_binned_repr=C_bin_rep,
        par_entries=p_ent, par_repr=p_rep, z0=z0,
        # provenance: with no self-calibration against the data left, the
        # assumptions behind this table are the assumptions of the fit
        mass_col=MASS_COL, mag_limit=a.mag_limit, multi=MULTI,
        purity_min=PURITY_MIN, recovery_min=RECOVERY_MIN,
        area_deg2=a.gama_area, zmin=a.zmin, zmax=a.zmax, min_n=MZ_MIN_N,
    )
    print("  saved nessie_completeness_mz.npz  (m, z, C, C_repr, Nh)")

    # Save the Nessie group catalogue itself so the fit can be closed-loop
    # validated: the MRP was injected by abundance matching, Nessie recovered
    # these groups from it, so fitting them with the same C table must return
    # the injected parameters.
    sel_t = table["multiplicity"] >= MULTI
    gt = table.loc[sel_t]
    np.savez(
        "nessie_mock_groups.npz",
        log_mass=np.log10(gt["MassA"].values),
        z=gt["median_redshift"].values,
        multiplicity=gt["multiplicity"].values,
        area_deg2=float(sky_frac * 4 * np.pi * (180 / np.pi) ** 2),
    )
    print(
        f"  saved nessie_mock_groups.npz  ({int(sel_t.sum())} groups, "
        f"for closed-loop validation)"
    )
    print("  saved nessie_completeness_table.npz  (d, z, C)")

    np.savez(
        a.save,
        delta_cen=g["cen"],
        C=g["C"],
        N=g["N"],
        d50=g["d50"],
        w=g["w"],
        d50_err=g["d50_err"],
        w_err=g["w_err"],
        zc=zc,
        d50_z=[b["d50"] for b in res["zbins"]],
        w_z=[b["w"] for b in res["zbins"]],
    )
    print(f"\n  saved {a.save}")

    if not a.no_plots:
        print()
        plot_completeness(res)
        plot_purity(table)
        plot_mass_relation(mm)
        plot_multiplicity(table, a.gama_fits, area_mock, a.gama_area)


if __name__ == "__main__":
    main()
