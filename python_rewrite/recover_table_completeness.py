"""
============================================================
Does the TABULATED (measured) completeness remove the M*/phi* residual?
============================================================

diagnose_lambda_floor showed Lambda is fine (lp bias -0.03, not -0.25), so the
coverage residual (M* +0.15, lp -0.25) is REPRESENTATION ERROR: the erf is a
few percent wrong near/above mlim, which tilts M* and drags phi*.

This measures C(Delta) directly from the mock detection flags in z-bins,
builds a 2D interpolator C(Delta, z), and refits with:

  A) ERF     C (current marg_comp)
  B) TABLE   C (measured, interpolated -- no functional form)

both with the CMIN floor, all detected groups, C in per-object AND Lambda.
If TABLE recovers ms=14.13 / lp=-3.96 where ERF sits at +0.15 / -0.25, the erf
shape was the residual -> port the table to Stan.

Run:  python recover_table_completeness.py
============================================================
"""

import numpy as np
from scipy import optimize
from scipy.special import erf
import recovery as R

ln10 = np.log(10.0)
TV = np.array([14.13, -3.96, -1.68, 0.63])
BE = 0.63
CMIN = R.CMIN

D_EDGES = np.arange(-1.5, 1.5 + 1e-9, 0.1)
Z_EDGES = np.array([0.01, 0.06, 0.10, 0.14, 0.18, 0.25])


def phi(m, ms, lp, al, be):
    u = m - ms
    return be * ln10 * 10**lp * 10 ** ((al + 1) * u) * np.exp(-(10 ** (be * u)))


def erf_C(delta, z):
    d50 = np.interp(z, R.COMP_Z_PTS, R.COMP_D50_PTS)
    w = np.interp(z, R.COMP_Z_PTS, R.COMP_W_PTS)
    return 0.5 * (1 + erf((delta - d50) / (np.sqrt(2) * w)))


def build_table(m_true, det, z, mlim_func):
    """C(Delta, z) measured from detection flags; returns (dcen, zcen, Ctab)."""
    delta = m_true - mlim_func(z)
    dcen = 0.5 * (D_EDGES[:-1] + D_EDGES[1:])
    zcen = 0.5 * (Z_EDGES[:-1] + Z_EDGES[1:])
    Ctab = np.zeros((zcen.size, dcen.size))
    for i in range(zcen.size):
        zm = (z >= Z_EDGES[i]) & (z < Z_EDGES[i + 1])
        for k in range(dcen.size):
            dm = zm & (delta >= D_EDGES[k]) & (delta < D_EDGES[k + 1])
            n = int(dm.sum())
            Ctab[i, k] = det[dm].mean() if n >= 15 else np.nan
        # fill: below -> 0, above -> last good; interpolate interior
        row = Ctab[i]
        good = np.isfinite(row)
        if good.sum() >= 2:
            Ctab[i] = np.interp(
                dcen, dcen[good], row[good], left=0.0, right=row[good][-1]
            )
        else:
            Ctab[i] = erf_C(dcen, zcen[i])
    Ctab = np.clip(Ctab, 0.0, 1.0)
    return dcen, zcen, Ctab


def table_C(delta, z, dcen, zcen, Ctab):
    """bilinear-ish: interp in delta per z-row, then linear in z."""
    delta = np.atleast_1d(delta)
    z = np.atleast_1d(z)
    iz = np.clip(np.searchsorted(zcen, z) - 1, 0, zcen.size - 2)
    t = np.clip((z - zcen[iz]) / (zcen[iz + 1] - zcen[iz]), 0, 1)
    out = np.empty(delta.shape)
    flat_iz = (
        np.broadcast_to(iz, delta.shape) if iz.size > 1 else np.full(delta.shape, iz[0])
    )
    flat_t = (
        np.broadcast_to(t, delta.shape) if t.size > 1 else np.full(delta.shape, t[0])
    )
    for i in range(zcen.size - 1):
        m = flat_iz == i
        if not m.any():
            continue
        c0 = np.interp(delta[m], dcen, Ctab[i], left=0.0, right=Ctab[i][-1])
        c1 = np.interp(delta[m], dcen, Ctab[i + 1], left=0.0, right=Ctab[i + 1][-1])
        out[m] = (1 - flat_t[m]) * c0 + flat_t[m] * c1
    return out


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

    dcen, zcen, Ctab = build_table(
        gv["log_mass_am"].values,
        gv["detected"].values.astype(bool),
        gv["zcos"].values,
        mlim_func,
    )
    print("  measured C table (rows=z, cols=Delta), a few columns:")
    show = [np.argmin(abs(dcen - d)) for d in (-0.4, -0.2, 0.0, 0.2, 0.4)]
    print("        Delta:", "  ".join(f"{dcen[k]:+.1f}" for k in show))
    for i, zc in enumerate(zcen):
        print(f"    z={zc:.3f}: ", "  ".join(f"{Ctab[i, k]:.2f}" for k in show))
    print(
        "\n  erf at same points (z=0.115):",
        "  ".join(f"{erf_C(dcen[k], 0.115):.2f}" for k in show),
    )

    mg = np.linspace(mlim_sh.min() - 3.0, R.XHI, 1200)

    def make(kind):
        Cobj_at = (
            (lambda d, z: erf_C(d, z))
            if kind == "erf"
            else (lambda d, z: table_C(d, z, dcen, zcen, Ctab))
        )
        C_now = Cobj_at(m_obs - mlim_obj, z_obs)
        keep = C_now > CMIN
        xk, sk = m_obs[keep], sigma_obs[keep]
        mlk, zk = mlim_obj[keep], z_obs[keep]

        def lam(th):
            tot = 0.0
            pg = phi(mg, *th)
            for j in range(len(V_sh)):
                Cj = Cobj_at(mg - mlim_sh[j], np.full(mg.shape, z_mids[j]))
                tot += V_sh[j] * np.trapezoid(pg * np.where(Cj > CMIN, Cj, 0.0), mg)
            return tot

        def sumL(th):
            g = np.linspace(-6, 6, 61)
            mt = xk[:, None] + sk[:, None] * g[None, :]
            Cm = Cobj_at(mt - mlk[:, None], np.broadcast_to(zk[:, None], mt.shape))
            integ = (
                phi(mt, *th)
                * Cm
                * np.exp(-0.5 * g[None, :] ** 2)
                / (sk[:, None] * np.sqrt(2 * np.pi))
            )
            return np.sum(np.log(np.maximum(np.trapezoid(integ, mt, axis=1), 1e-300)))

        def nll(t3):
            th = [t3[0], t3[1], t3[2], BE]
            return -(-lam(th) + sumL(th))

        r = optimize.minimize(
            nll,
            [14.2, -4.0, -1.6],
            method="Nelder-Mead",
            options=dict(xatol=1e-4, fatol=1e-3, maxiter=6000),
        )
        return r.x, int(keep.sum())

    (xe, ne), (xt, nt) = make("erf"), make("table")
    print(f"\n  truth        ms={TV[0]:.3f}  lp={TV[1]:.3f}  al={TV[2]:.3f}")
    print(
        f"  A) ERF   C   ms={xe[0]:.3f}  lp={xe[1]:.3f}  al={xe[2]:.3f}   N={ne}"
        f"   (d_ms={xe[0] - TV[0]:+.2f}, d_lp={xe[1] - TV[1]:+.2f})"
    )
    print(
        f"  B) TABLE C   ms={xt[0]:.3f}  lp={xt[1]:.3f}  al={xt[2]:.3f}   N={nt}"
        f"   (d_ms={xt[0] - TV[0]:+.2f}, d_lp={xt[1] - TV[1]:+.2f})"
    )
    print("\n  Reading it:")
    print("  - TABLE d_ms/d_lp much smaller than ERF -> the erf shape WAS the")
    print("    residual; port the tabulated C to Stan (pass the grid as data).")
    print("  - TABLE no better -> residual is elsewhere (C measured vs TRUE mass")
    print("    but applied at observed-mass-convolved m; a deeper subtlety).")


if __name__ == "__main__":
    main()
