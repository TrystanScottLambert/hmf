"""
============================================================
DIAGNOSTIC: is the sharp turnover cut (which sits ABOVE M*) the source of
the shared M*/phi* bias?  SIMPLE model, true masses (no error confound).
============================================================

Two fits, head to head, on the same mock:

  A) SHARP CUT (current approach): keep only groups above mlim(z); normalise
     with a sharp cut,  Lambda = sum_j V_sh int_{mlim} phi dm.
     mlim ~ 13.6-14.4 sits ABOVE M*=13.51, so the knee is discarded.

  B) TRUE COMPLETENESS (candidate fix): keep ALL detected groups, including
     sub-mlim; normalise with the mock's real completeness curve,
     Lambda = sum_j V_sh int phi(m) C(m - mlim_j) dm,  C measured from the
     detection flags. The knee is retained and its incompleteness modelled.

If B recovers M*/phi*/alpha and A does not, the sharp cut above M* is the
shared bias, and the fix is completeness-weighted Lambda + keeping sub-mlim
data. If B ALSO fails, the problem is deeper (abundance matching / volume).

Run:  python diagnose_completeness.py
Needs: WAVES parquet files + numpy/scipy/pandas.
============================================================
"""

import numpy as np
from scipy import optimize
import recovery as R

ln10 = np.log(10.0)
TV = np.array([13.51, -3.19, -1.27, 0.47])


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


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

    m_true = gv["log_mass_am"].values  # true mass for ALL groups in volume
    det = gv["detected"].values.astype(bool)  # detection flag
    z_all = gv["zcos"].values
    mlim_all = mlim_func(z_all)
    return dict(
        m_true=m_true,
        det=det,
        z_all=z_all,
        mlim_all=mlim_all,
        mlim_sh=mlim_sh,
        V_sh=V_sh,
        sky_frac=sky_frac,
    )


def completeness_curve(D):
    """C(Delta) = detected fraction vs Delta = m_true - mlim(z), pooled over z."""
    Delta = D["m_true"] - D["mlim_all"]
    edges = np.arange(-2.0, 2.0 + 1e-9, 0.1)
    cen = 0.5 * (edges[:-1] + edges[1:])
    C = np.full(cen.size, np.nan)
    for k in range(cen.size):
        sel = (Delta >= edges[k]) & (Delta < edges[k + 1])
        if sel.sum() >= 10:
            C[k] = D["det"][sel].mean()
    # fill: below measured range -> 0, above -> 1, interpolate interior nans
    good = np.isfinite(C)
    C = np.interp(cen, cen[good], C[good])
    C = np.clip(C, 0.0, 1.0)
    print("  completeness curve C(m - mlim):")
    for dd in [-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]:
        print(f"    C({dd:+.2f}) = {np.interp(dd, cen, C):.2f}")
    return lambda d: np.interp(d, cen, C, left=0.0, right=1.0)


def lam_sharp(th, D):
    ms, lp, al, be = th
    mg = np.linspace(D["mlim_sh"].min() - 3, R.XHI, 2000)
    pg = phi(mg, ms, lp, al, be)
    return sum(
        D["V_sh"][j] * np.trapezoid(pg * (mg >= D["mlim_sh"][j]), mg)
        for j in range(len(D["V_sh"]))
    )


def lam_comp(th, D, Cfun):
    ms, lp, al, be = th
    mg = np.linspace(D["mlim_sh"].min() - 3, R.XHI, 2000)
    pg = phi(mg, ms, lp, al, be)
    return sum(
        D["V_sh"][j] * np.trapezoid(pg * Cfun(mg - D["mlim_sh"][j]), mg)
        for j in range(len(D["V_sh"]))
    )


def negll(th, x, lamfn):
    if not (0.1 < th[3] < 2.0):
        return 1e12
    ll = np.sum(np.log(np.maximum(phi(x, *th), 1e-300)))
    return -(-lamfn(th) + ll)


def mle(x, lamfn, label, D):
    r = optimize.minimize(
        negll,
        TV,
        args=(x, lamfn),
        method="Nelder-Mead",
        options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
    )
    lam_at_truth = lamfn(TV)
    print(
        f"  {label:26s} N={x.size:4d}  Lam(truth)={lam_at_truth:6.0f}  "
        f"ms={r.x[0]:7.3f} lp={r.x[1]:7.3f} al={r.x[2]:7.3f} be={r.x[3]:6.3f}"
    )
    return r.x


def main():
    D = build()
    Cfun = completeness_curve(D)

    # Variant A: sharp cut, only groups above mlim
    aboveA = D["det"] & (D["m_true"] > D["mlim_all"])
    xA = D["m_true"][aboveA]
    # Variant B: all detected groups (sub-mlim retained)
    xB = D["m_true"][D["det"]]

    print("\n--- SIMPLE model on TRUE masses (no error confound) ---")
    print(
        "  truth                       "
        "                       ms= 13.510 lp= -3.190 al= -1.270 be= 0.470"
    )
    mle(xA, lambda th: lam_sharp(th, D), "A) SHARP cut (>mlim)", D)
    mle(xB, lambda th: lam_comp(th, D, Cfun), "B) TRUE completeness", D)

    print("\n  Reading it:")
    print("  - B recovers, A doesn't -> the sharp cut above M* is the shared bias;")
    print("    fix = completeness-weighted Lambda + keep sub-mlim data (both models).")
    print("  - B also fails -> problem is upstream (abundance matching / volume).")
    print("  - Lam(truth) should ~ N for the correct normalisation of each variant.")


if __name__ == "__main__":
    main()
