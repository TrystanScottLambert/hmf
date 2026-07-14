"""
============================================================
CHECK: does the C>cmin floor break the Lambda normalisation?
============================================================

Coverage on marg_comp showed logphi* biased LOW by -0.25 dex and M* HIGH by
+0.15 dex, both tightly constrained. phi* is the normalisation, so the prime
suspect is Lambda vs the kept sample.

phi* enters the likelihood multiplicatively, so with the shape at truth the
conditional MLE of logphi* is analytic:

    lp_cond = lp_true + log10( N_kept / Lambda(truth) )

If Lambda(truth) != N_kept, that mismatch IS the phi* bias. Predicted bias of
-0.25 dex needs Lambda/N ~ 10^0.25 ~ 1.8.

Compares three Lambda constructions at truth:
  A) HARD-ZERO  (what the Stan model does):  int phi * (C>cmin ? C : 0)
  B) SMOOTH     (no floor at all):           int phi * C
  C) SHIFTED    (integrate phi*C only above the mass where C crosses cmin,
                 per shell -- the 'consistent with the data cut' version)

and reports Lambda/N_kept and the implied logphi* bias for each.

Run:  python diagnose_lambda_floor.py
============================================================
"""

import numpy as np
from scipy.special import erf
import recovery as R

ln10 = np.log(10.0)
TV = dict(ms=14.13, lp=-3.96, al=-1.68, be=0.63)


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def Cz(m, mlim, d50, w):
    return 0.5 * (1 + erf((m - mlim - d50) / (np.sqrt(2) * w)))


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, _ = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)
    z_obs, m_obs, sigma_obs, _ = R.add_mass_errors(gv, np.random.default_rng(42))

    mlim_func, _, _, _ = R.turnover_mlim(z_obs, m_obs, zmin=R.ZMIN, zmax=R.ZLIMIT)
    z_mids, V_sh = R.shell_volumes(sky_frac, zmin=R.ZMIN, zmax=R.ZLIMIT)
    mlim_sh = mlim_func(z_mids)
    mlim_obj = mlim_func(z_obs)

    d50_obj = np.interp(z_obs, R.COMP_Z_PTS, R.COMP_D50_PTS)
    w_obj = np.interp(z_obs, R.COMP_Z_PTS, R.COMP_W_PTS)
    d50_sh = np.interp(z_mids, R.COMP_Z_PTS, R.COMP_D50_PTS)
    w_sh = np.interp(z_mids, R.COMP_Z_PTS, R.COMP_W_PTS)

    cmin = R.CMIN
    C_obj = Cz(m_obs, mlim_obj, d50_obj, w_obj)
    keep = C_obj > cmin
    N_kept = int(keep.sum())

    mg = np.linspace(mlim_sh.min() - 3.0, R.XHI, 4000)
    pg = phi(mg, **TV)

    lamA = lamB = lamC = 0.0
    for j in range(len(V_sh)):
        C = Cz(mg, mlim_sh[j], d50_sh[j], w_sh[j])
        lamA += V_sh[j] * np.trapezoid(pg * np.where(C > cmin, C, 0.0), mg)  # hard-zero
        lamB += V_sh[j] * np.trapezoid(pg * C, mg)  # smooth, no floor
        # C: integrate phi*C only above the mass where C crosses cmin
        above = C > cmin
        lamC += V_sh[j] * np.trapezoid((pg * C)[above], mg[above])

    print(f"  N_kept (C > {cmin})   = {N_kept}")
    print(f"  N_detected (all)     = {m_obs.size}")
    print(
        f"\n  {'Lambda variant':>26s}  {'Lambda(truth)':>13s}  {'Lam/N':>7s}  {'lp bias':>9s}"
    )
    for name, lam in [
        ("A) hard-zero (Stan model)", lamA),
        ("B) smooth, no floor", lamB),
        ("C) integrate above cross", lamC),
    ]:
        bias = np.log10(N_kept / lam)
        print(f"  {name:>26s}  {lam:13.0f}  {lam / N_kept:7.3f}  {bias:+9.3f}")

    print(f"\n  observed coverage lp bias = -0.254 dex  (M* bias +0.147)")
    print("\n  Reading it:")
    print("  - a variant whose 'lp bias' ~ -0.25 explains the coverage bias:")
    print("    that Lambda is over-predicting the kept count -> phi* absorbs it.")
    print("  - if ALL variants give Lam/N ~ 1 (bias ~ 0), the normalisation is")
    print("    fine and the M*/phi* residual is the erf shape, not the floor.")


if __name__ == "__main__":
    main()
