#!/usr/bin/env python3
"""An SDSS-specific version of Driver+22's figure 3 -- the "vuvuzela" plot.

Figure 3 determines the group mass error as a function of multiplicity.  Its
caption gives the recipe: for every group with more than 20 members, remove
members one at a time, recompute the mass at each multiplicity, and take the
16/50/84 per cent quantiles of the resulting tracks.  The trails are coloured by
the original halo mass.

**Driver applies the GAMA-derived curve to SDSS unchanged.**  ``sdsshmf.r``
line 281 uses exactly the same hardcoded ``xx``/``yy`` as ``gamahmf.r`` line
270.  That is a reasonable default when you have no SDSS mocks, but we have the
SDSS group members, so the curve can be measured directly for this catalogue and
for the mass estimator actually in use.

That last point matters more than it looks.  Nessie's ``estimated_mass`` and
Tempel's mass are **not** the same estimator: on 7258 groups with byte-identical
membership, Nessie's masses sit +0.285 dex higher (see ``--mass-audit``).  An
error function derived for one does not automatically describe the other.

The derivation code was never published.  ``gamahmf.r`` and ``sdsshmf.r`` only
carry the resulting arrays; grepping the machine finds nothing but copies of
them.  Nessie's Rust has the right idea in ``fof/src/group_properties.rs`` --
``calculate_error_function`` documents this exact procedure -- but
``create_mass_error_track`` is still ``todo!()``.  So this is a fresh
implementation, with the mass estimator ported from that same Rust so the tracks
use precisely the masses the HMF uses.

Members are removed **faintest-first in apparent magnitude**, because that is
what physically happens to a group as it recedes through a flux-limited survey.

Usage
-----
    python vuvuzela.py                    # Nessie SDSS, writes vuvuzela_sdss.pdf
    python vuvuzela.py --mass-audit       # why Nessie and Tempel masses differ
    python vuvuzela.py --min-mult 20 --out vuvuzela_sdss.pdf
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

import driver_recovery as dr

NESSIE_DIR = "/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS"
NESSIE_GROUPS = f"{NESSIE_DIR}/sdss_groups.parquet"
NESSIE_GALS = f"{NESSIE_DIR}/sdss_galaxies.parquet"
TEMPEL_GALS = "../data/sdssdr10table1.fits"
TEMPEL_GROUPS = "../data/sdssdr10table2.fits"

C_KMS = 299792.458
NESSIE_H0, NESSIE_OM = 70.0, 0.30      # recovered from Nessie's own co_dist
SIGMA_ERR_SQ = 50.0                    # velocity_dispersion_gap_err == sqrt(50)
GRAV_RAD_FACTOR = 4.582                # group_properties.rs: grav_rad = 4.582 * sky_disp
TEMPEL_RG_FACTOR = 3.9322              # recovered by inverting Tempel's own Eq 8

# Driver+22's GAMA curve (gamahmf.r line 270) and the Robotham+2011 estimate it
# replaced (the "OrigErr" branch, line 273 -- the dotted line in figure 3).
GAMA_XX = np.arange(3, 23, dtype=float)
GAMA_YY = np.array([
    0.68389355, 0.38719116, 0.40325591, 0.32696735, 0.27680685, 0.24018684,
    0.20226682, 0.18645475, 0.17437005, 0.14271506, 0.13922450, 0.13482418,
    0.13741619, 0.11715141, 0.12134983, 0.10078830, 0.09944761, 0.09913166,
    0.08590223, 0.07588408])
ROBOTHAM_XX = np.array([3.0, 4.3, 6.25, 8.1, 12.1, 19.6, 31.0, 47.5, 68.9,
                        81.9, 100.0])
ROBOTHAM_YY = np.array([0.866, 0.763, 0.716, 0.679, 0.572, 0.419, 0.310, 0.246,
                        0.189, 0.171, 0.126])


# ---------------------------------------------------------------------------
# Nessie's mass estimator, ported from fof/src/group_properties.rs
# ---------------------------------------------------------------------------

def _comoving_distance(z, h0=NESSIE_H0, om=NESSIE_OM, n=2048):
    """Flat LCDM comoving distance, vectorised over z."""
    z = np.atleast_1d(np.asarray(z, float))
    zz = np.linspace(0.0, z.max() if z.max() > 0 else 1e-6, n)
    e = 1.0 / np.sqrt(om * (1 + zz) ** 3 + (1 - om))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(zz))])
    return np.interp(z, zz, cum) * C_KMS / h0


def velocity_dispersion_gapper(zs, sigma_err_sq=SIGMA_ERR_SQ):
    """``Group::velocity_dispersion_gapper``.

    Note the peculiar-velocity convention: the Rust uses ``z * c / (1 + z_med)``
    on the raw redshift, not ``(z - z_med) * c / (1 + z_med)``.  Since only
    *gaps* between sorted velocities are used, the constant offset cancels and
    the two agree exactly; it is written the same way here so the port can be
    diffed against the original.
    """
    n = len(zs)
    if n < 2:
        return 0.0
    zmed = np.median(zs)
    v = np.sort(np.asarray(zs, float) * C_KMS / (1.0 + zmed))
    gaps = np.diff(v)
    i = np.arange(1, n)
    w = i * (n - i)
    sigma_gap = (np.sqrt(np.pi) / (n * (n - 1.0))) * np.sum(w * gaps)
    raw_sq = (n * sigma_gap ** 2) / (n - 1.0)
    return np.sqrt(raw_sq - sigma_err_sq) if raw_sq > sigma_err_sq else 0.0


def _to_cartesian(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def iterative_centre_idx(ra, dec, absmag):
    """``Group::calculate_iterative_center_idx``.

    Flux-weighted centre; repeatedly drop the object furthest from it until two
    remain, then return the brighter survivor's original index.
    """
    xyz = _to_cartesian(ra, dec)
    flux = 10.0 ** (-0.4 * np.asarray(absmag, float))
    idx = np.arange(len(ra))
    while len(idx) > 2:
        f = flux[idx]
        centre = (xyz[idx] * f[:, None]).sum(axis=0) / f.sum()
        d = np.linalg.norm(xyz[idx] - centre, axis=1)
        idx = np.delete(idx, np.argmax(d))
    return idx[np.argmax(flux[idx])]


def _angsep_small(ra, dec, ra0, dec0):
    """``angular_separation_small_angle``, in degrees."""
    dra = (np.asarray(ra, float) - ra0) * np.cos(np.radians(dec0))
    ddec = np.asarray(dec, float) - dec0
    return np.hypot(dra, ddec)


def sky_distribution(ra, dec, zs, ra0, dec0):
    """``Group::calculate_sky_distribution``, in Mpc.

    kpc/arcsec comoving = D_C * 1000 * (pi/180/3600), so the arcsecond
    round-trip in the Rust cancels to a plain small-angle projected radius in
    comoving Mpc, then divided by (1+z) via the (1+z)^2 under the square root.
    """
    zmed = float(np.median(zs))
    dc = float(_comoving_distance(zmed)[0])
    proj = np.radians(_angsep_small(ra, dec, ra0, dec0)) * dc      # Mpc comoving
    n = len(ra)
    return np.sqrt(np.sum(proj ** 2) / (n * 2.0 * (1.0 + zmed) ** 2))


def total_mass(grav_rad, sigma):
    """``calculate_total_mass`` -- Tempel+2014 equation 8."""
    return 2.325e12 * grav_rad * ((3.0 ** (1.0 / 3.0)) * sigma / 100.0) ** 2


def nessie_mass(ra, dec, zs, absmag):
    """The full ``estimated_mass`` chain for one set of members."""
    if len(ra) < 2:
        return np.nan
    k = iterative_centre_idx(ra, dec, absmag)
    sd = sky_distribution(ra, dec, zs, ra[k], dec[k])
    sig = velocity_dispersion_gapper(zs)
    return total_mass(GRAV_RAD_FACTOR * sd, sig)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_members(verbose=True):
    """Nessie SDSS members with Tempel's absolute magnitudes joined on."""
    gal = pd.read_parquet(NESSIE_GALS)
    gal = gal[gal.group_id != -1].copy()
    t = Table.read(TEMPEL_GALS)

    def native(a):
        a = np.asarray(a)
        return a.astype(a.dtype.newbyteorder("=")) if a.dtype.byteorder == ">" else a

    tp = pd.DataFrame({"ra_t": native(t["col13"]), "dec_t": native(t["col14"]),
                       "absmag_r": native(t["col31"])})
    key = lambda r, d: np.round(r, 6).astype(str) + "_" + np.round(d, 6).astype(str)
    tp["_k"] = key(tp.ra_t.values, tp.dec_t.values)
    tp = tp.drop_duplicates("_k")
    gal["_k"] = key(gal.RAJ2000.values, gal.DEJ2000.values)
    gal = gal.merge(tp[["_k", "absmag_r"]], on="_k", how="left").drop(columns="_k")
    if gal.absmag_r.isna().any():
        raise RuntimeError("unmatched galaxies; the two catalogues have drifted")
    if verbose:
        print(f"  members loaded              : {len(gal)}")
    return gal


# ---------------------------------------------------------------------------
# The tracks
# ---------------------------------------------------------------------------

def build_tracks(gal, min_mult=20, floor=3, verbose=True):
    """Remove the apparently faintest member repeatedly, recomputing the mass.

    Returns a list of (multiplicity array, delta-log10-mass array, log10 M_full)
    -- one entry per group, which is one "trail" of the vuvuzela.
    """
    sizes = gal.groupby("group_id").size()
    big = sizes[sizes > min_mult].index
    if verbose:
        print(f"  groups with N > {min_mult}          : {len(big)}")
    tracks = []
    for j, gid in enumerate(big):
        m = gal[gal.group_id == gid]
        order = np.argsort(m.rmag.values)          # brightest first
        ra = m.RAJ2000.values[order]
        dec = m.DEJ2000.values[order]
        zs = m.zobs.values[order]
        am = m.absmag_r.values[order]
        n_full = len(ra)
        ns, ms = [], []
        for n in range(n_full, floor - 1, -1):     # keep the n brightest
            mm = nessie_mass(ra[:n], dec[:n], zs[:n], am[:n])
            if np.isfinite(mm) and mm > 0:
                ns.append(n)
                ms.append(np.log10(mm))
        if len(ns) < 5:
            continue
        ns, ms = np.array(ns), np.array(ms)
        full = ms[0]                                # mass at full multiplicity
        tracks.append((ns, ms - full, full))
        if verbose and (j + 1) % 25 == 0:
            print(f"    {j + 1}/{len(big)}", end="\r", flush=True)
    if verbose:
        print(" " * 30, end="\r")
        print(f"  usable tracks               : {len(tracks)}")
    return tracks


def quantile_curve(tracks, nmax=None):
    """16/50/84 per cent quantiles of delta-log10-mass at each multiplicity."""
    alln = np.concatenate([t[0] for t in tracks])
    alld = np.concatenate([t[1] for t in tracks])
    ns = np.arange(int(alln.min()), int(nmax or alln.max()) + 1)
    q16, q50, q84, cnt = [], [], [], []
    for n in ns:
        d = alld[alln == n]
        if len(d) < 5:
            q16.append(np.nan); q50.append(np.nan); q84.append(np.nan)
            cnt.append(len(d))
            continue
        q16.append(np.quantile(d, 0.16))
        q50.append(np.quantile(d, 0.50))
        q84.append(np.quantile(d, 0.84))
        cnt.append(len(d))
    return (ns, np.array(q16), np.array(q50), np.array(q84), np.array(cnt))


def error_function(ns, q16, q84):
    """Half the 16-84 spread -- the 1-sigma log10 mass error at each N."""
    return 0.5 * (np.array(q84) - np.array(q16))


# ---------------------------------------------------------------------------
# Mass audit
# ---------------------------------------------------------------------------

def mass_audit():
    """Why Nessie's SDSS masses sit above Tempel's for the very same groups."""
    gal = pd.read_parquet(NESSIE_GALS)
    gal = gal[(gal.group_id != -1) & (gal.GroupID != -1)].copy()
    gal["gk"] = (np.round(gal.RAJ2000, 6).astype(str) + "_"
                 + np.round(gal.DEJ2000, 6).astype(str))
    nes = gal.groupby("group_id")["gk"].apply(frozenset)
    tem = gal.groupby("GroupID")["gk"].apply(frozenset)
    tmap = {v: k for k, v in tem.items()}
    match = {g: tmap[s] for g, s in nes.items() if s in tmap}

    ng = pd.read_parquet(NESSIE_GROUPS).set_index("group_id")
    t2 = Table.read(TEMPEL_GROUPS)
    tg = pd.DataFrame({
        "idcl": np.asarray(t2["col1"], float), "sig": np.asarray(t2["col12"], float),
        "skydisp": np.asarray(t2["col13"], float),
        "mass": np.asarray(t2["col15"], float)}).set_index("idcl")

    rows = []
    for g, t in match.items():
        if g in ng.index and float(t) in tg.index:
            a, b = ng.loc[g], tg.loc[float(t)]
            if a.multiplicity >= 5 and b.sig > 0 and b.skydisp > 0:
                rows.append((a.multiplicity, a.velocity_dispersion_gap, b.sig,
                             a.estimated_mass, b.mass * 1e12))
    d = pd.DataFrame(rows, columns=["N", "sig_n", "sig_t", "m_n", "m_t"])
    ls = np.log10(d.sig_n / d.sig_t)
    lm = np.log10(d.m_n / d.m_t)
    rad = np.log10(GRAV_RAD_FACTOR / TEMPEL_RG_FACTOR)

    print(f"\n  Mass audit -- groups with byte-identical membership: {len(d)}")
    print("  Both catalogues use Tempel+2014 eq. 8, so the difference is in its")
    print("  inputs, not its form.\n")
    print(f"  log10(M_nessie / M_tempel)      median {np.median(lm):+.4f}"
          f"   scatter {lm.std():.4f}")
    print(f"    from sigma (gapper vs Tempel) {2 * np.median(ls):+.4f}"
          f"   (sigma ratio {np.median(ls):+.4f}, mass goes as sigma^2)")
    print(f"    from the R_g constant         {rad:+.4f}"
          f"   ({GRAV_RAD_FACTOR} vs {TEMPEL_RG_FACTOR} x sky dispersion)")
    print(f"    unexplained residual          {np.median(lm) - 2 * np.median(ls) - rad:+.4f}"
          f"   (the sky-dispersion measure itself)")
    print(f"\n  For reference log10(1/h) at h=0.7 = {np.log10(1 / 0.7):+.4f}, so this"
          "\n  is not an h-convention slip.")
    return d


# ---------------------------------------------------------------------------

def plot(tracks, ns, q16, q50, q84, outfile, min_mult):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    fig, ax = plt.subplots(figsize=(7.1, 5.4))
    masses = np.array([t[2] for t in tracks])
    lo, hi = np.quantile(masses, 0.02), np.quantile(masses, 0.98)
    cmap = plt.get_cmap("winter")
    norm = plt.Normalize(lo, hi)

    segs = [np.column_stack([t[0], t[1]]) for t in tracks]
    ax.add_collection(LineCollection(segs, colors=[cmap(norm(m)) for m in masses],
                                     linewidths=0.7, alpha=0.35))

    good = np.isfinite(q50)
    ax.plot(ns[good], q50[good], color="red", lw=2, zorder=5)
    ax.plot(ns[good], q16[good], color="red", lw=2, ls="-", zorder=5)
    ax.plot(ns[good], q84[good], color="red", lw=2, ls="-", zorder=5)

    ax.plot(GAMA_XX, GAMA_YY, color="black", ls=(0, (6, 2)), lw=1.8, zorder=6)
    ax.plot(GAMA_XX, -GAMA_YY, color="black", ls=(0, (6, 2)), lw=1.8, zorder=6)
    ax.plot(ROBOTHAM_XX, ROBOTHAM_YY, color="black", ls=":", lw=1.8, zorder=6)
    ax.plot(ROBOTHAM_XX, -ROBOTHAM_YY, color="black", ls=":", lw=1.8, zorder=6)

    ax.axhline(0, color="grey", lw=0.6, zorder=1)
    ax.set_xlim(2, max(60, int(np.nanmax(ns[good])) + 2))
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("Multiplicity $N_{fof}$", fontsize=11)
    ax.set_ylabel(r"$\Delta\,\log_{10}$(Mass)", fontsize=11)
    ax.tick_params(direction="in", top=True, right=True, which="both")
    ax.minorticks_on()

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, pad=0.015)
    cb.set_label(r"$\log_{10}$(original halo mass / M$_\odot$)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    ax.plot([], [], color="red", lw=2, label="16, 50, 84 per cent (this work, SDSS)")
    ax.plot([], [], color="black", ls=(0, (6, 2)), lw=1.8,
            label="Driver+22 GAMA curve (used for SDSS in sdsshmf.r)")
    ax.plot([], [], color="black", ls=":", lw=1.8, label="Robotham+2011")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.set_title(f"Nessie SDSS, groups with $N>{min_mult}$, faintest-first removal",
                 fontsize=9)

    fig.tight_layout()
    fig.savefig(outfile, dpi=240)
    plt.close(fig)
    print(f"\n  wrote {outfile}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-mult", type=int, default=20)
    p.add_argument("--floor", type=int, default=3)
    p.add_argument("--out", default="vuvuzela_sdss.pdf")
    p.add_argument("--csv", default="sdss_masserr.csv")
    p.add_argument("--mass-audit", action="store_true")
    args = p.parse_args()

    print("SDSS mass error vs multiplicity (Driver+22 figure 3, SDSS-specific)")
    if args.mass_audit:
        mass_audit()
        return

    gal = load_members()
    tracks = build_tracks(gal, min_mult=args.min_mult, floor=args.floor)
    ns, q16, q50, q84, cnt = quantile_curve(tracks)
    err = error_function(ns, q16, q84)

    # q50 is the direct analogue of gamahmf.r's `masscorr` (the multiplicity
    # debiasing), and `err` of its `yy` (the scatter).  Print both comparisons.
    print("\n     N   tracks     q16     q50     q84   sigma |  GAMA yy  masscorr")
    for i, n in enumerate(ns):
        if not np.isfinite(q50[i]) or n > 30:
            continue
        g = np.interp(n, GAMA_XX, GAMA_YY, left=np.nan, right=np.nan)
        m = dr.MASSCORR[n - 1] if n <= len(dr.MASSCORR) else np.nan
        print(f"  {n:4d} {cnt[i]:8d} {q16[i]:7.3f} {q50[i]:7.3f} {q84[i]:7.3f}"
              f" {err[i]:7.3f} | {g:8.3f} {m:9.3f}")

    sel = (ns >= 3) & (ns <= 22) & np.isfinite(q50)
    gg = np.interp(ns[sel], GAMA_XX, GAMA_YY)
    mm = np.array([dr.MASSCORR[int(n) - 1] for n in ns[sel]])
    skew = np.median(np.abs(q16[sel]) / np.abs(q84[sel]))
    print(f"\n  sigma  SDSS/GAMA ratio : {np.median(err[sel] / gg):.3f}"
          "   (1.0 => Driver's reuse of the GAMA curve is justified)")
    print(f"  bias   SDSS/GAMA ratio : {np.median(q50[sel] / mm):.3f}"
          "   (q50 is the analogue of masscorr)")
    print(f"  skew   |q16|/|q84|     : {skew:.2f}"
          "   (1.0 would be symmetric; Driver's MC assumes a Gaussian)")

    out = pd.DataFrame({"N": ns, "ntracks": cnt, "q16": q16, "q50": q50,
                        "q84": q84, "sigma_log10M": err})
    out.to_csv(args.csv, index=False)
    print(f"\n  wrote {args.csv}")
    plot(tracks, ns, q16, q50, q84, args.out, args.min_mult)


if __name__ == "__main__":
    main()
