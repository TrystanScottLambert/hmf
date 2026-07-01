"""
============================================================
DIAGNOSTIC: is the marginalised phi* bias a normalisation problem?
============================================================

Reuses the exact recovery.py pipeline, rebuilds Lambda at the TRUE
parameters in numpy, and checks one decisive identity:

    logphi*(fit) - logphi*(truth) = log10( N_observed / Lambda(truth) )

(holding ms, al, be at truth -- phi* enters the likelihood as a pure
multiplicative constant, so its conditional MLE is analytic).

If Lambda(truth) >> N_observed, the model expects far more groups above
mlim than are actually detected, and the fit lowers phi* to compensate --
that IS the phi* bias, and it means the "complete above mlim" assumption
baked into Lambda is wrong.

Run:  python diagnose_norm.py
Needs: the WAVES parquet files (same as recovery.py) + numpy/scipy/pandas.
============================================================
"""

import numpy as np
from scipy.stats import norm
import recovery as R

ln10 = np.log(10.0)


def main():
    # ---- rebuild the exact GAMA-mock pipeline (as run_real_pipeline) ----
    print("Rebuilding the GAMA-mock pipeline ...")
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)

    rng = np.random.default_rng(42)
    z_obs, m_obs, sigma_obs, nfof = R.add_mass_errors(gv, rng)

    mlim_func, _, kind, _ = R.turnover_mlim(z_obs, m_obs)
    z_mids, V_sh = R.shell_volumes(sky_frac)
    mlim_sh = mlim_func(z_mids)
    sig_sh = R.sigma_eff_per_shell(z_obs, m_obs, sigma_obs, mlim_sh)

    above = m_obs > mlim_func(z_obs)
    N_fit = int(above.sum())
    print(f"  mlim(z) [{kind}]  N above mlim (observed) = {N_fit}")
    print(f"  sig_sh: min={sig_sh.min():.2f} max={sig_sh.max():.2f}")

    # ---- Lambda at truth: soft Phi boundary (the model) and sharp cut ----
    mgrid = np.linspace(mlim_sh.min() - 2.0, R.XHI, 4000)
    pg = R.mrp_phi(mgrid, **R.TRUE)
    Lam_soft = 0.0
    Lam_sharp = 0.0
    for j in range(len(V_sh)):
        w = norm.cdf((mgrid - mlim_sh[j]) / sig_sh[j])
        Lam_soft += V_sh[j] * np.trapezoid(pg * w, mgrid)
        hard = (mgrid >= mlim_sh[j]).astype(float)
        Lam_sharp += V_sh[j] * np.trapezoid(pg * hard, mgrid)

    print("\n--- Normalisation check (at TRUE parameters) ---")
    print(f"  N_observed          = {N_fit}")
    print(
        f"  Lambda_soft (model) = {Lam_soft:8.0f}   ratio Lam/N = {Lam_soft / N_fit:5.2f}"
    )
    print(
        f"  Lambda_sharp (cut)  = {Lam_sharp:8.0f}   ratio Lam/N = {Lam_sharp / N_fit:5.2f}"
    )
    print(f"\n  PREDICTED logphi* bias (soft) = {np.log10(N_fit / Lam_soft):+.2f} dex")
    print(f"  (observed fit bias was about -0.81 dex: -4.00 vs -3.19)")

    # ---- WHY: true completeness above mlim, from the mock itself ----
    # gv carries every group's AM (true) mass and detection flag.
    m_true_all = gv["log_mass_am"].values
    det_all = gv["detected"].values
    z_all = gv["zcos"].values
    mlim_all = mlim_func(z_all)

    # completeness as a function of (true mass - local mlim)
    print("\n--- True completeness vs true mass, relative to mlim(z) ---")
    print("  (fraction of ALL abundance-matched groups that are detected)")
    dm = m_true_all - mlim_all
    edges = np.array([-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    print(f"  {'m-mlim bin':>14s}  {'N_tot':>7s}  {'N_det':>7s}  {'complete':>9s}")
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (dm >= a) & (dm < b)
        nt = int(sel.sum())
        nd = int((sel & det_all).sum())
        if nt > 0:
            print(f"  [{a:+.2f},{b:+.2f})  {nt:7d}  {nd:7d}  {nd / nt:9.1%}")

    # completeness integrated above mlim (true mass > mlim)
    ab = m_true_all > mlim_all
    print(
        f"\n  Groups with TRUE mass above mlim: {int(ab.sum())} total, "
        f"{int((ab & det_all).sum())} detected "
        f"-> {(det_all[ab].mean()):.1%} complete"
    )
    print("\n  If completeness above mlim is well below 100%, Lambda (which assumes")
    print("  ~complete above mlim) over-predicts, and phi* absorbs it. That would")
    print("  mean the turnover sits too deep into the incomplete regime for the")
    print("  marginalised normalisation.")


if __name__ == "__main__":
    main()
