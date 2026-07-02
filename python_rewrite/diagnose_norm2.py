"""
============================================================
DIAGNOSTIC: is the stable logphi* bias a NORMALISATION/VOLUME bug, or the
beta-phi* degeneracy dragging phi* down?
============================================================

Key idea: phi* enters the likelihood as a pure multiplicative constant, so
with the SHAPE (ms, al, be) held at truth, the conditional MLE of logphi* is
analytic:

    lp_cond = lp_true + log10( N / Lambda(truth) )

  - If lp_cond ~ -3.19 (truth): Lambda(truth) ~ N, the normalisation is
    CORRECT, and the joint-fit logphi* bias (~-1.2) must come from the SHAPE
    moving -- i.e. the beta-phi* degeneracy (beta biased high pulls phi*
    down). That is a data/degeneracy issue, not a bug.
  - If lp_cond ~ -4.4: Lambda(truth) >> N, the normalisation genuinely
    over-counts -> a volume / z-range / V_sh bug to fix.

Also checks the volume plumbing directly: sum(V_sh) vs Vsurvey.

Run:  python diagnose_norm2.py
============================================================
"""

import numpy as np
import recovery as R

ln10 = np.log(10.0)
TV = np.array([13.51, -3.19, -1.27, 0.47])


def phi(m, th):
    ms, lp, al, be = th
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)
    z_mids, V_sh = R.shell_volumes(sky_frac)

    # --- volume plumbing ---
    print("--- volume consistency ---")
    print(f"  Vsurvey (used in AM)      = {Vsurvey:.4e}")
    print(f"  sum(V_sh) (used in Lambda)= {V_sh.sum():.4e}")
    print(f"  ratio sum(V_sh)/Vsurvey   = {V_sh.sum() / Vsurvey:.4f}   (want 1.0000)")

    m_true = gv["log_mass_am"].values
    det = gv["detected"].values.astype(bool)
    z_all = gv["zcos"].values
    sig_det = R.sigma_from_nfof(gv["n_gama"].values[det])

    def lam_sharp(th, mlim_sh):
        mg = np.linspace(mlim_sh.min() - 3, R.XHI, 2000)
        pg = phi(mg, th)
        return sum(
            V_sh[j] * np.trapezoid(pg * (mg >= mlim_sh[j]), mg)
            for j in range(len(V_sh))
        )

    cases = {"1_trunc": [], "2_+complete": [], "3_+eddington": []}
    for r in range(15):
        rng = np.random.default_rng(1000 + r)
        m_obs_det = m_true[det] + rng.normal(0, sig_det)
        z_det = z_all[det]
        try:
            mlim_func, _, _, _ = R.turnover_mlim(z_det, m_obs_det)
        except RuntimeError:
            continue
        mlim_sh = mlim_func(z_mids)
        Ltruth = lam_sharp(TV, mlim_sh)

        for key, mask, mm in [
            ("1_trunc", m_true > mlim_func(z_all), m_true),
            ("2_+complete", det & (m_true > mlim_func(z_all)), m_true),
            ("3_+eddington", None, None),
        ]:
            if key == "3_+eddington":
                mask3 = m_obs_det > mlim_func(z_det)
                N = int(mask3.sum())
            else:
                N = int(mask.sum())
            lp_cond = TV[1] + np.log10(N / Ltruth)
            cases[key].append((N, Ltruth, lp_cond))

    print("\n--- normalisation at truth, per sample definition ---")
    print(f"  logphi*_true = {TV[1]:.3f}")
    print(
        f"\n  {'case':>14s}  {'N':>6s}  {'Lam(truth)':>11s}  {'Lam/N':>7s}  "
        f"{'lp_cond':>9s}  {'joint_lp*':>10s}"
    )
    joint = {"1_trunc": -4.336, "2_+complete": -4.456, "3_+eddington": -4.548}
    for key in ["1_trunc", "2_+complete", "3_+eddington"]:
        arr = np.array(cases[key])
        N, L, lpc = arr.mean(axis=0)
        print(
            f"  {key:>14s}  {N:6.0f}  {L:11.0f}  {L / N:7.2f}  {lpc:9.3f}  {joint[key]:10.3f}"
        )
    print(
        "\n  (*joint_lp = the full-MLE logphi* from diagnose_systematic, for comparison)"
    )
    print("\n  Reading it:")
    print("  - lp_cond ~ -3.19 (Lam/N ~ 1) but joint_lp ~ -4.4  ->  normalisation")
    print("    is FINE; the phi* bias is the beta-phi* degeneracy (shape moving).")
    print("    Fix = constrain the cutoff (beta), not the volume.")
    print("  - lp_cond ~ -4.4 (Lam/N >> 1)  ->  real volume/normalisation bug.")


if __name__ == "__main__":
    main()
