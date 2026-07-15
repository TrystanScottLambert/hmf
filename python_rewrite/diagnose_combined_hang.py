"""
============================================================
Why does combined_comp hang? Evaluate the likelihood in numpy at the init point.
============================================================

Rebuilds the EXACT combined data (GAMA fixed ramp + SDSS adopted ramp), then
evaluates each survey's log-likelihood at the prior-centre parameters and
reports:
  - any NaN / inf / extreme sigma, mlim, V_sh
  - how many per-object integrals underflow to the -300 floor
  - the per-survey Lambda and log-sum, and whether the total is finite

A hang in Stan almost always means the init lands in a -300 wasteland (zero
gradient) or a NaN poisons the gradient. This finds which, in seconds.

Run:  python diagnose_combined_hang.py <gama_fits> <sdss_parquet>
============================================================
"""

import sys
import numpy as np
from scipy.stats import norm
import recovery as R

ln10 = np.log(10.0)
TH = dict(ms=14.13, lp=-3.96, al=-1.68, be=0.63)  # prior centre = init


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def audit(tag, x, sig, mlim_obj, d50, w, V_sh, mlim_sh, d50_sh, w_sh, sky_note):
    print(f"\n  === {tag} ===  N={x.size}  {sky_note}")
    for nm, v in [
        ("x_obs", x),
        ("sig", sig),
        ("mlim_obj", mlim_obj),
        ("V_sh", V_sh),
        ("mlim_sh", mlim_sh),
    ]:
        bad = (~np.isfinite(v)).sum()
        print(
            f"    {nm:9s} min={np.nanmin(v):.3g} max={np.nanmax(v):.3g} "
            f"med={np.nanmedian(v):.3g}  nonfinite={bad}"
        )
    if (sig <= 0).any():
        print(f"    !! {(sig <= 0).sum()} groups with sig<=0  (would break the grid)")
    if (sig > 1.0).any():
        print(f"    !! {(sig > 1.0).sum()} groups with sig>1.0 (huge grid span)")

    # per-object integral, exactly as Stan does it (Nint=31, +-5 sigma)
    g = np.linspace(-5, 5, 31)
    mt = x[:, None] + sig[:, None] * g[None, :]
    C = norm.cdf((mt - mlim_obj[:, None] - d50[:, None]) / w[:, None])
    integ = (
        phi(mt, **TH)
        * C
        * np.exp(-0.5 * g[None, :] ** 2)
        / (sig[:, None] * np.sqrt(2 * np.pi))
    )
    val = np.trapezoid(integ, mt, axis=1)
    floored = int((val <= 1e-300).sum())
    logsum = np.sum(np.log(np.maximum(val, 1e-300)))
    print(f"    per-object: {floored}/{x.size} underflow to -300 floor")
    print(
        f"    sum log L_i = {logsum:.1f}   ({'FINITE' if np.isfinite(logsum) else 'NOT FINITE'})"
    )

    # Lambda
    mg = np.linspace(mlim_sh.min() - 2.5, R.XHI, 1000)
    pg = phi(mg, **TH)
    Lam = 0.0
    for j in range(len(V_sh)):
        Cj = norm.cdf((mg - mlim_sh[j] - d50_sh[j]) / w_sh[j])
        Cj = np.where(Cj > R.CMIN, Cj, 0.0)
        Lam += V_sh[j] * np.trapezoid(pg * Cj, mg)
    print(
        f"    Lambda = {Lam:.1f}   Lambda/N = {Lam / x.size:.3f}   "
        f"({'FINITE' if np.isfinite(Lam) else 'NOT FINITE'})"
    )
    print(f"    survey loglik = {logsum - Lam:.1f}")
    return logsum - Lam


def main():
    gama_fits = sys.argv[1] if len(sys.argv) > 1 else "../data/G3CFoFGroupv10.fits"
    sdss_parq = sys.argv[2] if len(sys.argv) > 2 else "sdss_groups.parquet"
    sdss_frac = 0.2126803

    g_lm, g_sig, g_z, _ = R.load_real_gama(gama_fits)
    s_lm, s_sig, s_z, _ = R.load_sdss_groups(sdss_parq, zmin=0.01, zmax=0.08)

    gama_frac = 179.92 * (np.pi / 180) ** 2 / (4 * np.pi)

    g_mlim_func, _, _, _ = R.turnover_mlim(g_z, g_lm, zmin=R.ZMIN, zmax=R.ZLIMIT)
    gz_mids, gV_sh = R.shell_volumes(gama_frac, zmin=R.ZMIN, zmax=R.ZLIMIT)
    g_mlim_obj = g_mlim_func(g_z)
    g_d50 = np.interp(g_z, R.COMP_Z_PTS, R.COMP_D50_PTS)
    g_w = np.interp(g_z, R.COMP_Z_PTS, R.COMP_W_PTS)
    g_keep = norm.cdf((g_lm - g_mlim_obj - g_d50) / g_w) > R.CMIN

    s_mlim_func, _, sk, _ = R.turnover_mlim(s_z, s_lm, zmin=0.01, zmax=0.08)
    sz_mids, sV_sh = R.shell_volumes(sdss_frac, zmin=0.01, zmax=0.08)
    s_mlim_obj = s_mlim_func(s_z)
    s_d50 = np.interp(s_z, R.COMP_Z_PTS, R.COMP_D50_PTS)
    s_w = np.interp(s_z, R.COMP_Z_PTS, R.COMP_W_PTS)
    s_keep = norm.cdf((s_lm - s_mlim_obj - s_d50) / s_w) > R.CMIN

    print(
        f"  SDSS turnover kind = {sk}  mlim(0.01)={s_mlim_func(0.01):.2f} "
        f"mlim(0.08)={s_mlim_func(0.08):.2f}"
    )

    llA = audit(
        "GAMA",
        g_lm[g_keep],
        g_sig[g_keep],
        g_mlim_obj[g_keep],
        g_d50[g_keep],
        g_w[g_keep],
        gV_sh,
        g_mlim_func(gz_mids),
        np.interp(gz_mids, R.COMP_Z_PTS, R.COMP_D50_PTS),
        np.interp(gz_mids, R.COMP_Z_PTS, R.COMP_W_PTS),
        f"V={gV_sh.sum():.3e}",
    )
    llB = audit(
        "SDSS",
        s_lm[s_keep],
        s_sig[s_keep],
        s_mlim_obj[s_keep],
        s_d50[s_keep],
        s_w[s_keep],
        sV_sh,
        s_mlim_func(sz_mids),
        np.interp(sz_mids, R.COMP_Z_PTS, R.COMP_D50_PTS),
        np.interp(sz_mids, R.COMP_Z_PTS, R.COMP_W_PTS),
        f"V={sV_sh.sum():.3e}",
    )

    print(
        f"\n  TOTAL loglik at init = {llA + llB:.1f}   "
        f"({'FINITE -> Stan should move' if np.isfinite(llA + llB) else 'NOT FINITE -> found the hang'})"
    )
    print("\n  Reading it:")
    print("  - lots of -300 underflows or a NaN sigma/shell -> that's the hang;")
    print("    the init sits in a flat/poisoned region, zero gradient.")
    print("  - all finite, few underflows, Lambda/N ~ 1 -> the model is fine and")
    print("    the slowness is elsewhere (grid size / N too large).")


if __name__ == "__main__":
    main()
