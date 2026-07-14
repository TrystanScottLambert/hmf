"""
============================================================
Measure the SDSS completeness ramp C(Delta, z) from the Shark lightcone.
============================================================

The Shark lightcone is deeper/wider than SDSS, so imposing the SDSS selection
on it gives a valid SDSS-like mock. Completeness is a property of the SELECTION,
so as long as the cuts match, the measured ramp applies to the real SDSS sample.

SDSS selection imposed here:
   r_SDSS < 17.77   (Main Galaxy Sample limit; GAMA uses 19.8)
   z < 0.08         (Driver's SDSS cut)
   >= 5 members     (same multiplicity as the Nessie SDSS catalogue)

CAVEAT (state it in the paper): 'detected' here means '>=5 members brighter
than 17.77', i.e. magnitude-limit-driven incompleteness. It does NOT model
Nessie's group-finding losses (fragmentation, mergers, missed groups). That
would require running Nessie on this cut lightcone -- the higher-fidelity
version. This captures the dominant effect.

Outputs SDSS_D50_PTS / SDSS_W_PTS to paste into recovery.py, and a plot.

Run:  python measure_completeness_sdss.py
============================================================
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
import recovery as R

# --- SDSS selection ---
SDSS_MAG_LIMIT = 17.77
SDSS_SEL_COL = "total_ap_dust_r_SDSS"
SDSS_ZMIN, SDSS_ZMAX = 0.01, 0.08
SDSS_MULTI = 5


def erf_ramp(d, d50, w):
    return 0.5 * (1 + erf((d - d50) / (np.sqrt(2) * w)))


def completeness_curve(delta, detected, edges, min_n=15):
    cen, C = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (delta >= a) & (delta < b)
        if int(sel.sum()) >= min_n:
            cen.append(0.5 * (a + b))
            C.append(detected[sel].mean())
    return np.array(cen), np.array(C)


def sdss_select(gv, galaxies):
    """Detection = >= SDSS_MULTI galaxies brighter than SDSS_MAG_LIMIT."""
    gal = galaxies[
        (galaxies["id_fof"] != -1)
        & (galaxies["log_mstar_total"] > 8)
        & (galaxies[SDSS_SEL_COL] < SDSS_MAG_LIMIT)
    ]
    counts = gal.groupby("id_group_sky").size().rename("n_sdss")
    gv = gv.merge(counts, left_on="id_group_sky", right_index=True, how="left")
    gv["n_sdss"] = gv["n_sdss"].fillna(0).astype(int)
    gv["detected"] = gv["n_sdss"] >= SDSS_MULTI
    return gv


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, _ = R.abundance_match(groups, sky_frac)
    gv = sdss_select(gv, galaxies)

    # restrict to the SDSS redshift range
    inz = (gv["zcos"].values > SDSS_ZMIN) & (gv["zcos"].values < SDSS_ZMAX)
    gv = gv[inz]
    m_true = gv["log_mass_am"].values
    det = gv["detected"].values.astype(bool)
    z = gv["zcos"].values

    print(
        f"  SDSS selection: {SDSS_SEL_COL} < {SDSS_MAG_LIMIT}, "
        f"z {SDSS_ZMIN}-{SDSS_ZMAX}, >= {SDSS_MULTI} members"
    )
    print(f"  {det.sum()} detected / {det.size} groups in z-range")
    if det.sum() < 200:
        print("  WARNING: few detections -- the ramp will be noisy.")

    # mlim(z) from the OBSERVED detected masses (as the pipeline does)
    rng = np.random.default_rng(42)
    sig = R.sigma_from_nfof(gv["n_sdss"].values[det])
    m_obs_det = m_true[det] + rng.normal(0, sig)
    z_det = z[det]
    mlim_func, _, tkind, _ = R.turnover_mlim(
        z_det, m_obs_det, zmin=SDSS_ZMIN, zmax=SDSS_ZMAX
    )
    print(
        f"  mlim(z) [{tkind}]: {mlim_func(SDSS_ZMIN):.2f} -> {mlim_func(SDSS_ZMAX):.2f}"
    )

    delta = m_true - mlim_func(z)
    edges = np.arange(-1.5, 1.5 + 1e-9, 0.1)

    cen, C = completeness_curve(delta, det, edges)
    (d50, w), _ = curve_fit(erf_ramp, cen, C, p0=[0.0, 0.25], maxfev=10000)
    print(
        f"\n  GLOBAL erf: D50 = {d50:+.3f},  w = {w:.3f}   "
        f"(mlim at C={erf_ramp(0.0, d50, w):.2f})"
    )

    print("\n  Per-redshift (universality check):")
    z_bins = [(0.01, 0.03), (0.03, 0.05), (0.05, 0.08)]
    d50s, ws, zcs, curves = [], [], [], []
    for za, zb in z_bins:
        m = (z >= za) & (z < zb)
        if m.sum() < 100:
            print(f"  {za:.2f}-{zb:.2f}   too few groups")
            continue
        cz, Cz = completeness_curve(delta[m], det[m], edges)
        try:
            (dz, wz), _ = curve_fit(erf_ramp, cz, Cz, p0=[0.0, 0.25], maxfev=10000)
            print(f"  {za:.2f}-{zb:.2f}   D50={dz:+.3f}  w={wz:.3f}  N={int(m.sum())}")
            d50s.append(dz)
            ws.append(wz)
            zcs.append(0.5 * (za + zb))
            curves.append((f"{za:.2f}-{zb:.2f}", cz, Cz))
        except Exception:
            print(f"  {za:.2f}-{zb:.2f}   fit failed")

    if d50s:
        print("\n  Paste into recovery.py:")
        print(f"  SDSS_COMP_Z_PTS = {[round(v, 3) for v in zcs]}")
        print(f"  SDSS_COMP_D50_PTS = {[round(v, 3) for v in d50s]}")
        print(f"  SDSS_COMP_W_PTS = {[round(v, 3) for v in ws]}")

    print("\n  Compare to GAMA ramp:")
    print(f"    GAMA  D50 {R.COMP_D50_PTS}   w {R.COMP_W_PTS}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    dd = np.linspace(-1.5, 1.5, 300)
    ax.plot(
        dd,
        erf_ramp(dd, d50, w),
        "k-",
        lw=2,
        label=f"SDSS global erf (D50={d50:+.2f}, w={w:.2f})",
    )
    ax.scatter(cen, C, c="k", s=35, zorder=5, label="SDSS measured")
    for name, cz, Cz in curves:
        ax.plot(cz, Cz, "o-", ms=4, alpha=0.7, label=f"z {name}")
    ax.plot(
        dd,
        erf_ramp(dd, R.COMP_D50_PTS[1], R.COMP_W_PTS[1]),
        "r--",
        lw=1.5,
        label="GAMA ramp (for comparison)",
    )
    ax.axvline(0, color="grey", ls=":")
    ax.set(
        xlabel=r"$\Delta = m_{\rm true} - m_{\rm lim}(z)$",
        ylabel="completeness C",
        xlim=(-1.5, 1.5),
        ylim=(-0.02, 1.02),
        title="SDSS-selected mock completeness ramp",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("completeness_ramp_sdss.pdf")
    print("\n  saved completeness_ramp_sdss.pdf")


if __name__ == "__main__":
    main()
