"""
============================================================
STEP 1: measure the completeness ramp C(m,z) from the mock.
============================================================

C(m,z) = P(detected | TRUE mass m, redshift z), measured straight from the
mock's detection flags (we have true AM mass + detected flag for ALL groups).
Measured against TRUE mass because that is what Lambda integrates over.

Checks the key simplification: does C collapse to a near-universal 1D ramp in
   Delta = m_true - mlim(z)
(i.e. same shape at all z, just shifted by the turnover)? If yes, C(m,z) is one
1D curve C(Delta) and we fit it with an erf:
   C(Delta) = 0.5 * (1 + erf((Delta - D50) / (sqrt(2) * w)))
giving the 50%-completeness offset D50 (relative to mlim) and the ramp width w.

Outputs the erf parameters (-> straight into the Lambda model) and a plot
overlaying the per-redshift curves against the global erf fit.

Run:  python measure_completeness.py
============================================================
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
import recovery as R


def erf_ramp(d, d50, w):
    return 0.5 * (1 + erf((d - d50) / (np.sqrt(2) * w)))


def completeness_curve(delta, detected, edges):
    cen, C, Ntot = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (delta >= a) & (delta < b)
        n = int(sel.sum())
        if n >= 20:
            cen.append(0.5 * (a + b))
            C.append(detected[sel].mean())
            Ntot.append(n)
    return np.array(cen), np.array(C), np.array(Ntot)


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)

    # mlim(z) from the observed detected masses, exactly as the pipeline builds it
    z_obs, m_obs, sigma_obs, _ = R.add_mass_errors(gv, np.random.default_rng(42))
    mlim_func, _, tkind, _ = R.turnover_mlim(z_obs, m_obs, zmin=R.ZMIN, zmax=R.ZLIMIT)

    m_true = gv["log_mass_am"].values
    det = gv["detected"].values.astype(bool)
    z = gv["zcos"].values
    delta = m_true - mlim_func(z)

    print(f"  mlim(z) [{tkind}]: {mlim_func(R.ZMIN):.2f} -> {mlim_func(R.ZLIMIT):.2f}")
    print(f"  {det.sum()} detected / {det.size} total")

    edges = np.arange(-1.5, 1.5 + 1e-9, 0.1)

    # global curve + erf fit
    cen, C, Nt = completeness_curve(delta, det, edges)
    (d50, w), _ = curve_fit(erf_ramp, cen, C, p0=[0.0, 0.2], maxfev=10000)
    print(f"\n  GLOBAL erf fit:  D50 = {d50:+.3f} dex (rel. to mlim),  w = {w:.3f} dex")
    print(f"  -> mlim sits at completeness C = {erf_ramp(0.0, d50, w):.2f}")
    print(
        f"  -> C reaches 90% at Delta = {d50 + 1.2816 * w:+.2f}, "
        f"50% at {d50:+.2f}, 10% at {d50 - 1.2816 * w:+.2f}"
    )

    # per-redshift universality check
    print("\n  Per-redshift C(Delta) (universality check):")
    print(f"  {'z-bin':>14s}  {'D50':>7s}  {'w':>6s}  {'N':>7s}")
    z_bins = [(0.01, 0.08), (0.08, 0.15), (0.15, 0.25)]
    percurves = []
    for za, zb in z_bins:
        m = (z >= za) & (z < zb)
        if m.sum() < 200:
            continue
        cz, Cz, Ntz = completeness_curve(delta[m], det[m], edges)
        percurves.append((f"{za:.2f}-{zb:.2f}", cz, Cz))
        try:
            (d50z, wz), _ = curve_fit(erf_ramp, cz, Cz, p0=[0.0, 0.2], maxfev=10000)
            print(f"  {za:.2f}-{zb:.2f}   {d50z:+.3f}  {wz:.3f}  {int(m.sum()):7d}")
        except Exception:
            print(f"  {za:.2f}-{zb:.2f}   fit failed")

    print("\n  If per-z D50/w are all close to the global values -> C is universal")
    print("  in Delta, and the 1D ramp C(Delta) is enough (parametric erf).")
    print("  If they drift with z -> need C(m,z) 2D (tabulated).")

    # plot
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
        label=f"global erf (D50={d50:+.2f}, w={w:.2f})",
    )
    ax.scatter(cen, C, c="k", s=40, zorder=5, label="global measured")
    for name, cz, Cz in percurves:
        ax.plot(cz, Cz, "o-", ms=4, alpha=0.7, label=f"z {name}")
    ax.axvline(0, color="grey", ls=":", label="mlim")
    ax.set(
        xlabel=r"$\Delta = m_{\rm true} - m_{\rm lim}(z)$",
        ylabel="completeness C",
        xlim=(-1.5, 1.5),
        ylim=(-0.02, 1.02),
        title="Mock completeness ramp",
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("completeness_ramp.pdf")
    print("\n  saved completeness_ramp.pdf")


if __name__ == "__main__":
    main()
