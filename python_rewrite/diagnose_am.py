"""
============================================================
DIAGNOSTIC: do the abundance-matched masses actually follow the injected
MRP in the fit range?  Pure data inspection -- no fitting, no modelling.
============================================================

The AM assigns mass so that rank r -> mass where Vsurvey * int_M phi = r.
By construction the empirical n(>M) of ALL groups should equal the MRP
prediction EXACTLY. This checks it, in the mass range we actually fit, and
separately shows the detected subsample.

If empirical n(>M) matches the MRP for ALL groups  -> AM is fine; the bias
is in how we model the detected subsample (completeness).
If it does NOT match -> AM itself is the bug (e.g. rank clamp, mvir ties,
or a mock abundance that differs from the MRP normalisation).

Run:  python diagnose_am.py
============================================================
"""

import numpy as np
import recovery as R

ln10 = np.log(10.0)


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)

    m_all = np.sort(gv["log_mass_am"].values)  # all groups in volume
    m_det = np.sort(gv.loc[gv["detected"], "log_mass_am"].values)
    N_all = m_all.size

    # MRP cumulative prediction: Vsurvey * int_M^inf phi dm
    mg = np.arange(9.0, 16.5 + 1e-9, 0.001)
    phi = R.mrp_phi(mg, **R.TRUE)
    cum = np.cumsum((phi * 0.001)[::-1])[::-1] * Vsurvey  # predicted n(>M)

    def n_gt_emp(m, arr):  # empirical count above m
        return arr.size - np.searchsorted(arr, m)

    def n_gt_mrp(m):  # predicted count above m
        return np.interp(m, mg, cum)

    print(f"N_all groups in volume = {N_all}")
    print(f"N_detected             = {m_det.size}")
    print(f"Vsurvey                = {Vsurvey:.3e}")
    print(
        f"AM mass range: {m_all.min():.2f} .. {m_all.max():.2f}  "
        f"(median {np.median(m_all):.2f})"
    )
    print(
        f"MRP predicts above m=9.0: {n_gt_mrp(9.0):.0f}   "
        f"(clamp pile-up if N_all exceeds this)"
    )

    print(
        "\n--- n(>M): empirical (ALL) vs MRP prediction  [the AM construction check] ---"
    )
    print(
        f"  {'M':>6s}  {'emp_all':>10s}  {'MRP_pred':>10s}  {'ratio':>7s}  {'emp_det':>9s}"
    )
    for M in [12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5]:
        ea = n_gt_emp(M, m_all)
        em = n_gt_mrp(M)
        ed = n_gt_emp(M, m_det)
        ratio = ea / em if em > 0 else np.nan
        print(f"  {M:6.2f}  {ea:10d}  {em:10.0f}  {ratio:7.2f}  {ed:9d}")

    print("\n  ratio ~ 1.00 everywhere -> AM masses follow the MRP; bias is in")
    print("  the completeness/selection model, not AM.")
    print("  ratio drifting from 1 -> AM itself does not reproduce the MRP in")
    print("  the fit range; that is the upstream bug.")

    # extra: differential check -- counts in fit-range mass bins
    print("\n--- differential dN/dlogM: empirical (ALL) vs MRP, 0.25-dex bins ---")
    print(f"  {'bin':>14s}  {'emp':>8s}  {'MRP':>8s}  {'ratio':>7s}")
    edges = np.arange(12.0, 15.5 + 1e-9, 0.25)
    for a, b in zip(edges[:-1], edges[1:]):
        emp = int(((m_all >= a) & (m_all < b)).sum())
        pred = Vsurvey * np.trapezoid(
            R.mrp_phi(np.linspace(a, b, 50), **R.TRUE), np.linspace(a, b, 50)
        )
        r = emp / pred if pred > 0 else np.nan
        print(f"  [{a:.2f},{b:.2f})  {emp:8d}  {pred:8.0f}  {r:7.2f}")


if __name__ == "__main__":
    main()
