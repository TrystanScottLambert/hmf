"""
============================================================
DIAGNOSTIC (option 3): localise the FIXED systematic seen in coverage.
SIMPLE model, numpy MLE, across error realisations -- no MCMC.
============================================================

The marg coverage showed a reproducible bias (ms~14.2, lp~-4.2, cov68=0%).
This asks WHERE it comes from by turning on one effect at a time, all in the
SIMPLE model (no marginalisation), so we can see if the bias is shared:

  1) ALL true masses above mlim, sharp Lambda
       -> pure TRUNCATION only. AM is exact, so these masses follow phi
          truncated at mlim, and Lambda_sharp models that truncation.
          If this is ALREADY biased -> the cut sitting ABOVE M* biases the
          MRP fit by itself (degeneracy / no knee). Nothing downstream can
          fix that; it's the sample.
  2) DETECTED true masses above mlim, sharp Lambda
       -> adds INCOMPLETENESS (sample is phi*C, Lambda assumes phi).
          bias(2) - bias(1) = the completeness effect.
  3) DETECTED observed masses above mlim, sharp Lambda   [= current SIMPLE]
       -> adds EDDINGTON (error scatter across the cut).
          bias(3) - bias(2) = the error/Eddington effect.

Compare the three mean biases to see which step introduces the systematic.
Run:  python diagnose_systematic.py
============================================================
"""

import numpy as np
from scipy import optimize
import recovery as R

ln10 = np.log(10.0)
TV = np.array([13.51, -3.19, -1.27, 0.47])
PN = ["ms", "lp", "al", "be"]


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
        mg = np.linspace(mlim_sh.min() - 3, R.XHI, 1500)
        pg = phi(mg, th)
        return sum(
            V_sh[j] * np.trapezoid(pg * (mg >= mlim_sh[j]), mg)
            for j in range(len(V_sh))
        )

    def negll(th, x, mlim_sh):
        if not (0.1 < th[3] < 2.0):
            return 1e12
        return -(
            -lam_sharp(th, mlim_sh) + np.sum(np.log(np.maximum(phi(x, th), 1e-300)))
        )

    def mle(x, mlim_sh):
        r = optimize.minimize(
            negll,
            TV,
            args=(x, mlim_sh),
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return r.x

    n_real = 15
    out = {"1_trunc": [], "2_+complete": [], "3_+eddington": []}
    Ns = {"1_trunc": [], "2_+complete": [], "3_+eddington": []}
    for r in range(n_real):
        rng = np.random.default_rng(1000 + r)
        m_obs_det = m_true[det] + rng.normal(0, sig_det)
        z_det = z_all[det]
        try:
            mlim_func, _, _, _ = R.turnover_mlim(z_det, m_obs_det)
        except RuntimeError:
            continue
        mlim_sh = mlim_func(z_mids)

        a1 = m_true > mlim_func(z_all)  # ALL true above mlim
        a2 = det & (m_true > mlim_func(z_all))  # DETECTED true above mlim
        a3 = m_obs_det > mlim_func(z_det)  # DETECTED obs above mlim

        out["1_trunc"].append(mle(m_true[a1], mlim_sh))
        Ns["1_trunc"].append(a1.sum())
        out["2_+complete"].append(mle(m_true[a2], mlim_sh))
        Ns["2_+complete"].append(a2.sum())
        out["3_+eddington"].append(mle(m_obs_det[a3], mlim_sh))
        Ns["3_+eddington"].append(int(a3.sum()))

    print(f"\n  SIMPLE-model MLE, mean over {n_real} realisations")
    print(f"  truth: ms={TV[0]:.3f}  lp={TV[1]:.3f}  al={TV[2]:.3f}  be={TV[3]:.3f}")
    print(
        f"\n  {'case':>14s}  {'N':>6s}  {'ms':>7s} {'lp':>7s} {'al':>7s} {'be':>7s}   (bias vs truth)"
    )
    for k in ["1_trunc", "2_+complete", "3_+eddington"]:
        arr = np.array(out[k])
        m = arr.mean(axis=0)
        b = m - TV
        print(
            f"  {k:>14s}  {int(np.mean(Ns[k])):6d}  "
            f"{m[0]:7.3f} {m[1]:7.3f} {m[2]:7.3f} {m[3]:7.3f}   "
            f"({b[0]:+.2f},{b[1]:+.2f},{b[2]:+.2f},{b[3]:+.2f})"
        )

    print("\n  Reading it:")
    print("  - case 1 already biased  -> the CUT ABOVE M* is the systematic by")
    print("    itself; completeness/errors/marginalisation are secondary. The")
    print("    fix has to change the SELECTION (reach below the knee) or the")
    print("    claim (fit only where the cut is), not the likelihood details.")
    print("  - case 1 clean, 2 biased -> completeness is the culprit.")
    print("  - cases 1-2 clean, 3 biased -> Eddington/errors (marg's job).")


if __name__ == "__main__":
    main()
