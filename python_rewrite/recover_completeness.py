"""
============================================================
STEP 2: does the completeness forward-model recover the truth? (numpy, mock)
============================================================

Replaces the sharp mlim cut with the measured erf completeness ramp
C(Delta) = 0.5(1+erf((Delta-D50)/(sqrt2 w))),  Delta = m - mlim(z),
applied CONSISTENTLY in both the per-object term and Lambda:

  per-object:  L_i   = int phi(m) C(m-mlim_i) N(x_i|m,sig_i) dm   (ALL detected)
  normalise:   Lambda= sum_j V_sh[j] int phi(m) C(m-mlim_j) dm

and keeps ALL detected groups (no mlim cut). Compared against the current
SHARP-CUT model (cut at mlim, soft-Phi Lambda, no C in per-object).

Truth (injected) = 14.13 / -3.96 / -1.68 / 0.63. beta fixed at 0.63.
If COMPLETENESS recovers alpha (~-1.68) where SHARP does not, the fix works
and we port it to Stan. numpy only, ~1-2 min.

Run:  python recover_completeness.py
============================================================
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm
from scipy.special import erf
import recovery as R

ln10 = np.log(10.0)
TV = np.array([14.13, -3.96, -1.68, 0.63])
BE = 0.63
D50, W = -0.208, 0.257  # from measure_completeness.py (global erf)

# per-z ramp parameters (Step 1): z-bin centre -> (D50, w)
Z_PTS = np.array([0.045, 0.115, 0.20])
D50_PTS = np.array([-0.148, -0.193, -0.232])
W_PTS = np.array([0.326, 0.256, 0.227])


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def Cramp(delta):
    return 0.5 * (1 + erf((delta - D50) / (np.sqrt(2) * W)))


def Cramp_z(delta, zval):
    """z-dependent ramp: D50(z), w(z) linearly interpolated from Step-1 bins."""
    d50 = np.interp(zval, Z_PTS, D50_PTS)
    wz = np.interp(zval, Z_PTS, W_PTS)
    return 0.5 * (1 + erf((delta - d50) / (np.sqrt(2) * wz)))


def main():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, Vsurvey = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)
    z_obs, m_obs, sigma_obs, _ = R.add_mass_errors(gv, np.random.default_rng(42))

    mlim_func, _, _, _ = R.turnover_mlim(z_obs, m_obs, zmin=R.ZMIN, zmax=R.ZLIMIT)
    z_mids, V_sh = R.shell_volumes(sky_frac, zmin=R.ZMIN, zmax=R.ZLIMIT)
    mlim_sh = mlim_func(z_mids)
    mlim_obj = mlim_func(z_obs)

    mg = np.linspace(mlim_sh.min() - 2.5, R.XHI, 2000)

    # ---------- COMPLETENESS model: all detected, C in per-object AND Lambda ----------
    xC, sC, mlC = m_obs, sigma_obs, mlim_obj  # ALL detected

    def lam_comp(th):
        pg = phi(mg, *th)
        return sum(
            V_sh[j] * np.trapezoid(pg * Cramp(mg - mlim_sh[j]), mg)
            for j in range(len(V_sh))
        )

    def sumL_comp(th):
        g = np.linspace(-6, 6, 61)
        mt = xC[:, None] + sC[:, None] * g[None, :]
        integ = (
            phi(mt, *th)
            * Cramp(mt - mlC[:, None])
            * np.exp(-0.5 * g[None, :] ** 2)
            / (sC[:, None] * np.sqrt(2 * np.pi))
        )
        val = np.trapezoid(integ, mt, axis=1)
        return np.sum(np.log(np.maximum(val, 1e-300)))

    def nll_comp(t3):
        th = [t3[0], t3[1], t3[2], BE]
        if not (0.1 < BE < 2.0):
            return 1e12
        return -(-lam_comp(th) + sumL_comp(th))

    # ---------- SHARP model (current): cut at mlim, soft-Phi Lambda, no C per-object ----------
    above = m_obs > mlim_obj
    xS, sS = m_obs[above], sigma_obs[above]
    sig_sh = R.sigma_eff_per_shell(
        z_obs, m_obs, sigma_obs, mlim_sh, zmin=R.ZMIN, zmax=R.ZLIMIT
    )

    def lam_sharp(th):
        pg = phi(mg, *th)
        return sum(
            V_sh[j] * np.trapezoid(pg * norm.cdf((mg - mlim_sh[j]) / sig_sh[j]), mg)
            for j in range(len(V_sh))
        )

    def sumL_sharp(th):
        g = np.linspace(-5, 5, 41)
        mt = xS[:, None] + sS[:, None] * g[None, :]
        integ = (
            phi(mt, *th)
            * np.exp(-0.5 * g[None, :] ** 2)
            / (sS[:, None] * np.sqrt(2 * np.pi))
        )
        val = np.trapezoid(integ, mt, axis=1)
        return np.sum(np.log(np.maximum(val, 1e-300)))

    def nll_sharp(t3):
        th = [t3[0], t3[1], t3[2], BE]
        return -(-lam_sharp(th) + sumL_sharp(th))

    print(f"  N detected (completeness fit) = {xC.size}")
    print(f"  N above mlim (sharp fit)      = {xS.size}")
    print(f"  completeness ramp: D50={D50}, w={W}  (mlim at C={Cramp(0.0):.2f})")

    # ---------- COMPLETENESS with z-DEPENDENT ramp (per-z Step-1 params) ----------
    d50_obj = np.interp(z_obs, Z_PTS, D50_PTS)
    wz_obj = np.interp(z_obs, Z_PTS, W_PTS)
    d50_sh = np.interp(z_mids, Z_PTS, D50_PTS)
    wz_sh = np.interp(z_mids, Z_PTS, W_PTS)

    def lam_compz(th):
        pg = phi(mg, *th)
        tot = 0.0
        for j in range(len(V_sh)):
            Cj = 0.5 * (
                1 + erf((mg - mlim_sh[j] - d50_sh[j]) / (np.sqrt(2) * wz_sh[j]))
            )
            tot += V_sh[j] * np.trapezoid(pg * Cj, mg)
        return tot

    def sumL_compz(th):
        g = np.linspace(-6, 6, 61)
        mt = xC[:, None] + sC[:, None] * g[None, :]
        delta = mt - mlC[:, None]
        Cz = 0.5 * (
            1 + erf((delta - d50_obj[:, None]) / (np.sqrt(2) * wz_obj[:, None]))
        )
        integ = (
            phi(mt, *th)
            * Cz
            * np.exp(-0.5 * g[None, :] ** 2)
            / (sC[:, None] * np.sqrt(2 * np.pi))
        )
        val = np.trapezoid(integ, mt, axis=1)
        return np.sum(np.log(np.maximum(val, 1e-300)))

    def nll_compz(t3):
        th = [t3[0], t3[1], t3[2], BE]
        return -(-lam_compz(th) + sumL_compz(th))

    r_sharp = optimize.minimize(
        nll_sharp,
        [14.2, -4.0, -1.4],
        method="Nelder-Mead",
        options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
    )
    r_comp = optimize.minimize(
        nll_comp,
        [14.2, -4.0, -1.6],
        method="Nelder-Mead",
        options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
    )
    r_compz = optimize.minimize(
        nll_compz,
        [14.2, -4.0, -1.6],
        method="Nelder-Mead",
        options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
    )

    # ---------- COMPLETENESS z-dep + C-FLOOR (drop groups below C_min) ----------
    # low-C faint groups get divided by a tiny C and dominate/tilt the fit;
    # keep only groups whose OWN completeness exceeds C_min (still model the ramp
    # in Lambda and per-object above the floor).
    def fit_floored(cmin):
        Cobj = 0.5 * (1 + erf((m_obs - mlim_obj - d50_obj) / (np.sqrt(2) * wz_obj)))
        keep = Cobj > cmin
        xk, sk = m_obs[keep], sigma_obs[keep]
        d50k, wzk, mlk = d50_obj[keep], wz_obj[keep], mlim_obj[keep]

        # Lambda: integrate phi*C but only over the mass range where C > cmin,
        # per shell, to match the floored per-object sample.
        def lam_fl(th):
            pg = phi(mg, *th)
            tot = 0.0
            for j in range(len(V_sh)):
                Cj = 0.5 * (
                    1 + erf((mg - mlim_sh[j] - d50_sh[j]) / (np.sqrt(2) * wz_sh[j]))
                )
                tot += V_sh[j] * np.trapezoid(pg * np.where(Cj > cmin, Cj, 0.0), mg)
            return tot

        def sumL_fl(th):
            g = np.linspace(-6, 6, 61)
            mt = xk[:, None] + sk[:, None] * g[None, :]
            Cz = 0.5 * (
                1
                + erf((mt - mlk[:, None] - d50k[:, None]) / (np.sqrt(2) * wzk[:, None]))
            )
            integ = (
                phi(mt, *th)
                * Cz
                * np.exp(-0.5 * g[None, :] ** 2)
                / (sk[:, None] * np.sqrt(2 * np.pi))
            )
            return np.sum(np.log(np.maximum(np.trapezoid(integ, mt, axis=1), 1e-300)))

        def nll_fl(t3):
            th = [t3[0], t3[1], t3[2], BE]
            return -(-lam_fl(th) + sumL_fl(th))

        r = optimize.minimize(
            nll_fl,
            [14.2, -4.0, -1.6],
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return r.x, int(keep.sum())

    print(
        f"\n  truth              ms={TV[0]:.3f}  lp={TV[1]:.3f}  al={TV[2]:.3f}  be={TV[3]:.2f}"
    )
    print(
        f"  SHARP (cut)        ms={r_sharp.x[0]:.3f}  lp={r_sharp.x[1]:.3f}  al={r_sharp.x[2]:.3f}"
        f"   (d_al={r_sharp.x[2] - TV[2]:+.2f})"
    )
    print(
        f"  COMPLETE (global)  ms={r_comp.x[0]:.3f}  lp={r_comp.x[1]:.3f}  al={r_comp.x[2]:.3f}"
        f"   (d_al={r_comp.x[2] - TV[2]:+.2f})"
    )
    print(
        f"  COMPLETE (z-dep)   ms={r_compz.x[0]:.3f}  lp={r_compz.x[1]:.3f}  al={r_compz.x[2]:.3f}"
        f"   (d_al={r_compz.x[2] - TV[2]:+.2f})"
    )
    for cmin in (0.2, 0.3, 0.5):
        xx, nk = fit_floored(cmin)
        print(
            f"  z-dep + C>{cmin}       ms={xx[0]:.3f}  lp={xx[1]:.3f}  al={xx[2]:.3f}"
            f"   (d_al={xx[2] - TV[2]:+.2f})  N={nk}"
        )
    print("\n  Reading it:")
    print("  - a floor brings al toward -1.68 -> low-C noisy groups were tilting it;")
    print("    port the z-dependent ramp WITH that C floor to Stan.")
    print("  - floor doesn't help -> residual is the erf shape / true-vs-obs mass,")
    print("    and al ~ -1.48 is good enough given the Driver alpha prior; port as-is.")


if __name__ == "__main__":
    main()
