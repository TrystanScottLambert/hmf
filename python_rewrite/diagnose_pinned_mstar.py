"""
============================================================
Does pinning M* bias alpha?  (mock, known truth)
============================================================

On real data, anchoring M* with REFLEX gives:
    GAMA alone      M* = 14.60 (free)  ->  alpha = -1.87
    GAMA + REFLEX   M* = 14.14 (pinned) ->  alpha = -1.48
while Driver has alpha = -1.68 AT M* = 14.13. The 0.20 dex gap is
perpendicular to the M*-alpha degeneracy, so it is not a projection effect.

This runs the SAME experiment on the mock, where truth is known
(M* = 14.13, alpha = -1.68):

  A) M* FREE            -- what the completeness model does unaided
  B) M* PINNED to truth -- the REFLEX-anchored configuration
  C) M* pinned to 14.60 -- pinned at the value the real GAMA fit prefers

Reading it:
  - B gives alpha ~ -1.68  ->  ridge geometry is sound; the real-data -1.48
    is a genuine statement about the GAMA faint end.
  - B gives alpha ~ -1.48  ->  the completeness model shallows alpha whenever
    M* is pinned below where the group data wants it; the REFLEX-anchored
    alpha is an artefact and should not be quoted.

Uses the z-dependent completeness ramp + CMIN floor, i.e. the production
marg_comp likelihood. beta pinned at 0.63 as in production. numpy only.

Run:  python diagnose_pinned_mstar.py
============================================================
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm
import recovery as R

ln10 = np.log(10.0)
TRUE_MS, TRUE_LP, TRUE_AL, BE = 14.13, -3.96, -1.68, 0.63


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


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

    # z-dependent ramp, exactly as prep_comp / the Stan model
    d50_o = np.interp(z_obs, R.COMP_Z_PTS, R.COMP_D50_PTS)
    w_o = np.interp(z_obs, R.COMP_Z_PTS, R.COMP_W_PTS)
    d50_s = np.interp(z_mids, R.COMP_Z_PTS, R.COMP_D50_PTS)
    w_s = np.interp(z_mids, R.COMP_Z_PTS, R.COMP_W_PTS)

    C_obj = norm.cdf((m_obs - mlim_obj - d50_o) / w_o)
    keep = C_obj > R.CMIN
    x, sg = m_obs[keep], sigma_obs[keep]
    ml_k, d50_k, w_k = mlim_obj[keep], d50_o[keep], w_o[keep]
    print(f"  mock: kept {x.size} / {m_obs.size} (C > {R.CMIN})")
    print(
        f"  truth: M* = {TRUE_MS}, logphi* = {TRUE_LP}, alpha = {TRUE_AL}, beta = {BE}\n"
    )

    mg = np.linspace(mlim_sh.min() - 3.0, R.XHI, 1200)
    gq = np.linspace(-5, 5, 41)

    def lam(th):
        pg = phi(mg, *th)
        tot = 0.0
        for j in range(len(V_sh)):
            Cj = norm.cdf((mg - mlim_sh[j] - d50_s[j]) / w_s[j])
            tot += V_sh[j] * np.trapezoid(pg * np.where(Cj > R.CMIN, Cj, 0.0), mg)
        return tot

    def sumL(th):
        mt = x[:, None] + sg[:, None] * gq[None, :]
        Cm = norm.cdf((mt - ml_k[:, None] - d50_k[:, None]) / w_k[:, None])
        integ = (
            phi(mt, *th)
            * Cm
            * np.exp(-0.5 * gq[None, :] ** 2)
            / (sg[:, None] * np.sqrt(2 * np.pi))
        )
        return np.sum(np.log(np.maximum(np.trapezoid(integ, mt, axis=1), 1e-300)))

    def fit(ms_fixed=None):
        if ms_fixed is None:
            f = lambda t: -(-lam([t[0], t[1], t[2], BE]) + sumL([t[0], t[1], t[2], BE]))
            r = optimize.minimize(
                f,
                [14.2, -4.0, -1.6],
                method="Nelder-Mead",
                options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
            )
            return r.x[0], r.x[1], r.x[2]
        f = lambda t: (
            -(-lam([ms_fixed, t[0], t[1], BE]) + sumL([ms_fixed, t[0], t[1], BE]))
        )
        r = optimize.minimize(
            f,
            [-4.0, -1.6],
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return ms_fixed, r.x[0], r.x[1]

    print(
        f"  {'configuration':<28} {'M*':>8} {'logphi*':>9} {'alpha':>8} {'d_alpha':>9}"
    )
    for label, msf in [
        ("A) M* free", None),
        ("B) M* pinned to truth 14.13", TRUE_MS),
        ("C) M* pinned to 14.60", 14.60),
    ]:
        ms, lp, al = fit(msf)
        print(f"  {label:<28} {ms:8.3f} {lp:9.3f} {al:8.3f} {al - TRUE_AL:+9.3f}")

    print(
        "\n  B) ~ -1.68  -> geometry sound; real-data -1.48 is a real faint-end statement."
    )
    print("  B) ~ -1.48  -> pinning M* low shallows alpha; the REFLEX-anchored alpha")
    print("                 is a model artefact and should not be quoted.")


if __name__ == "__main__":
    main()
