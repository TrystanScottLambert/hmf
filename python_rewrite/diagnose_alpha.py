"""
============================================================
DIAGNOSTIC: is the alpha/phi* bias in the LIKELIHOOD, and is it the
sigma-treatment?  (GAMA-mock, priors OFF)
============================================================

Finds the unconstrained MLE (no priors) of both the SIMPLE and the
MARGINALISED likelihood on the real WAVES GAMA-mock, plus two controls.

What each line tells you:
  SIMPLE (obs=truth) : if its alpha ~ -1.27, the simple likelihood is fine
                       and the bias is MARGINALISATION-specific.
  MARG (real sigma)  : should reproduce the biased fit (alpha ~ -0.95).
  MARG (const sigma) : swaps the mass-dependent errors for one constant
                       sigma. If alpha recovers here, the mass-dependent
                       sigma / sig_sh treatment is the trigger.
  MARG (sigma -> 0)  : MUST converge to SIMPLE. If it does NOT, the marg
                       per-object term and/or Lambda_soft are inconsistent
                       with the simple ones in the no-error limit -> a bug,
                       not a modelling effect.

Run:  python diagnose_alpha.py
Needs: the WAVES parquet files + numpy/scipy/pandas.
============================================================
"""

import numpy as np
from scipy.stats import norm
from scipy import optimize
import recovery as R

ln10 = np.log(10.0)
TV = np.array([13.51, -3.19, -1.27, 0.47])


def build():
    groups, galaxies = R.load_catalogues(R.DATA_DIR)
    sky_frac = (
        R.sky_area_deg2(groups["ra"], groups["dec"]) * (np.pi / 180) ** 2 / (4 * np.pi)
    )
    gv, _ = R.abundance_match(groups, sky_frac)
    gv = R.gama_select(gv, galaxies)
    z_obs, m_obs, sigma_obs, _ = R.add_mass_errors(gv, np.random.default_rng(42))
    mlim_func, _, _, _ = R.turnover_mlim(z_obs, m_obs)
    z_mids, V_sh = R.shell_volumes(sky_frac)
    mlim_sh = mlim_func(z_mids)
    sig_sh = R.sigma_eff_per_shell(z_obs, m_obs, sigma_obs, mlim_sh)
    above = m_obs > mlim_func(z_obs)
    return dict(
        x=m_obs[above], sig=sigma_obs[above], mlim_sh=mlim_sh, V_sh=V_sh, sig_sh=sig_sh
    )


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def lam_soft(th, D):
    ms, lp, al, be = th
    mg = np.linspace(D["mlim_sh"].min() - 2, R.XHI, 1500)
    pg = phi(mg, ms, lp, al, be)
    return sum(
        D["V_sh"][j]
        * np.trapezoid(pg * norm.cdf((mg - D["mlim_sh"][j]) / D["sig_sh"][j]), mg)
        for j in range(len(D["V_sh"]))
    )


def lam_sharp(th, D):
    ms, lp, al, be = th
    mg = np.linspace(D["mlim_sh"].min() - 2, R.XHI, 1500)
    pg = phi(mg, ms, lp, al, be)
    return sum(
        D["V_sh"][j] * np.trapezoid(pg * (mg >= D["mlim_sh"][j]), mg)
        for j in range(len(D["V_sh"]))
    )


def sumL_marg(th, D, Nint=41):
    ms, lp, al, be = th
    x, sig = D["x"], D["sig"]
    g = np.linspace(-5, 5, Nint)
    mt = x[:, None] + sig[:, None] * g[None, :]
    u = mt - ms
    ph = be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))
    ga = np.exp(-0.5 * g[None, :] ** 2) / (sig[:, None] * np.sqrt(2 * np.pi))
    val = np.trapezoid(ph * ga, mt, axis=1)
    return np.sum(np.log(np.maximum(val, 1e-300)))


def sumL_simple(th, D):
    ms, lp, al, be = th
    return np.sum(np.log(np.maximum(phi(D["x"], ms, lp, al, be), 1e-300)))


def negll_marg(th, D):
    if not (0.1 < th[3] < 2.0):
        return 1e12
    return -(-lam_soft(th, D) + sumL_marg(th, D))


def negll_simple(th, D):
    if not (0.1 < th[3] < 2.0):
        return 1e12
    return -(-lam_sharp(th, D) + sumL_simple(th, D))


def mle(negll, D, label):
    r = optimize.minimize(
        negll,
        TV,
        args=(D,),
        method="Nelder-Mead",
        options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
    )
    dal = r.x[2] - TV[2]
    print(
        f"  {label:22s} ms={r.x[0]:7.3f}  lp={r.x[1]:7.3f}  "
        f"al={r.x[2]:7.3f}  be={r.x[3]:6.3f}   (d_al={dal:+.2f})"
    )
    return r.x


def main():
    D = build()
    print(
        f"N_fit={D['x'].size}  median sig={np.median(D['sig']):.3f}  "
        f"sig range {D['sig'].min():.2f}-{D['sig'].max():.2f}"
    )
    print(f"sig_sh: {np.round(D['sig_sh'], 2)}")
    print(
        "\n--- Unconstrained MLE, NO priors (where the likelihood actually peaks) ---"
    )
    print("  truth                  ms= 13.510  lp= -3.190  al= -1.270  be= 0.470")
    mle(negll_simple, D, "SIMPLE (obs=truth)")
    mle(negll_marg, D, "MARG (real sigma)")

    Dc = dict(D)
    s = float(np.median(D["sig"]))
    Dc["sig"] = np.full(D["x"].size, s)
    Dc["sig_sh"] = np.full(len(D["sig_sh"]), s)
    mle(negll_marg, Dc, "MARG (const sigma)")

    Dt = dict(D)
    Dt["sig"] = np.full(D["x"].size, 0.05)
    Dt["sig_sh"] = np.full(len(D["sig_sh"]), 0.05)
    mle(negll_marg, Dt, "MARG (sigma->0)")

    print("\n  Reading it:")
    print(
        "  - SIMPLE al ~ -1.27  &  MARG(real) al >> -1.27  -> bias is marginalisation-specific."
    )
    print(
        "  - MARG(const) recovers but MARG(real) doesn't    -> the mass-dependent sigma is the trigger."
    )
    print(
        "  - MARG(sigma->0) must match SIMPLE; if not, the marg terms are inconsistent (a bug)."
    )


if __name__ == "__main__":
    main()
