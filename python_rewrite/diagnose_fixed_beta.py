"""
============================================================
CHECK: does fixing beta = 0.47 remove the phi* bias?
============================================================

diagnose_norm2 showed the normalisation is correct and phi* is dragged down
purely by beta floating high. This pins beta at the literature value and
refits (ms, lp, al) -- if phi* snaps back to -3.19, the whole systematic was
the beta-phi* degeneracy, and "fix/constrain beta" is the cure.

Shows free-beta (4-param) vs beta-fixed (3-param) side by side.
Run:  python diagnose_betafix.py
============================================================
"""

import numpy as np
from scipy import optimize
import recovery as R

ln10 = np.log(10.0)
TV = np.array([13.51, -3.19, -1.27, 0.47])
BE_FIX = 0.47


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

    def negll_free(th, x, mlim_sh):
        if not (0.1 < th[3] < 2.0):
            return 1e12
        return -(
            -lam_sharp(th, mlim_sh) + np.sum(np.log(np.maximum(phi(x, th), 1e-300)))
        )

    def negll_fix(th3, x, mlim_sh):
        th = np.array([th3[0], th3[1], th3[2], BE_FIX])
        return -(
            -lam_sharp(th, mlim_sh) + np.sum(np.log(np.maximum(phi(x, th), 1e-300)))
        )

    def mle_free(x, mlim_sh):
        r = optimize.minimize(
            negll_free,
            TV,
            args=(x, mlim_sh),
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return r.x

    def mle_fix(x, mlim_sh):
        r = optimize.minimize(
            negll_fix,
            TV[:3],
            args=(x, mlim_sh),
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return np.array([r.x[0], r.x[1], r.x[2], BE_FIX])

    free = {"1_trunc": [], "2_+complete": [], "3_+eddington": []}
    fix = {"1_trunc": [], "2_+complete": [], "3_+eddington": []}
    for r in range(15):
        rng = np.random.default_rng(1000 + r)
        m_obs_det = m_true[det] + rng.normal(0, sig_det)
        z_det = z_all[det]
        try:
            mlim_func, _, _, _ = R.turnover_mlim(z_det, m_obs_det)
        except RuntimeError:
            continue
        mlim_sh = mlim_func(z_mids)
        masks = {
            "1_trunc": (m_true > mlim_func(z_all), m_true),
            "2_+complete": (det & (m_true > mlim_func(z_all)), m_true),
            "3_+eddington": (m_obs_det > mlim_func(z_det), m_obs_det),
        }
        for key, (mask, mm) in masks.items():
            x = mm[mask]
            free[key].append(mle_free(x, mlim_sh))
            fix[key].append(mle_fix(x, mlim_sh))

    print(f"\n  truth:  ms={TV[0]:.3f}  lp={TV[1]:.3f}  al={TV[2]:.3f}  be={TV[3]:.3f}")
    print(
        f"\n  {'case':>14s} {'fit':>10s}  {'ms':>7s} {'lp':>7s} {'al':>7s} {'be':>6s}"
        f"   (lp bias)"
    )
    for key in ["1_trunc", "2_+complete", "3_+eddington"]:
        mf = np.array(free[key]).mean(axis=0)
        mx = np.array(fix[key]).mean(axis=0)
        print(
            f"  {key:>14s} {'free-beta':>10s}  {mf[0]:7.3f} {mf[1]:7.3f} {mf[2]:7.3f} "
            f"{mf[3]:6.3f}   ({mf[1] - TV[1]:+.2f})"
        )
        print(
            f"  {'':>14s} {'beta=0.47':>10s}  {mx[0]:7.3f} {mx[1]:7.3f} {mx[2]:7.3f} "
            f"{mx[3]:6.3f}   ({mx[1] - TV[1]:+.2f})"
        )
    print(
        "\n  If beta=0.47 rows have lp ~ -3.19 and ms ~ 13.51  ->  fixing/constraining"
    )
    print(
        "  the cutoff cures the phi* bias; the systematic was the beta-phi* degeneracy."
    )


if __name__ == "__main__":
    main()
