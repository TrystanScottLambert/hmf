#!/usr/bin/env python3
"""kappa(c) = R_g / sigma_sky for a truncated NFW halo, from first principles.

This is the derivation Tempel+2014 delegates to Bartelmann (1996) and Lokas &
Mamon (2001) rather than writing out.  An earlier attempt in this project got
it wrong and had to fall back on calibrating kappa empirically against his
published catalogue; ``--compare`` checks the derivation against that
calibration.

The step that was wrong before: sigma_sky is the second moment of the
**projected** profile,

    sigma_sky^2 = (1/2) * int Sigma(R) R^3 dR / int Sigma(R) R dR

matching his estimator sigma_sky^2 = (2N)^-1 sum R_i^2 (his eq. 18).  It is
NOT <r^2>/3 from a 3D sphere, and it is not a Rayleigh-converted half-mass
radius -- that construction applies only to the Hernquist profile, whose
projected second moment diverges.

Two conventions are genuinely ambiguous and are exposed as options, because
the paper and the released Fortran disagree:

* ``x_out``  -- the upper limit of the energy integral (his eq. 9).  The paper
  says R_200; ``groups_mass.f90`` integrates to ``par%nfw_r200*10.0`` while
  keeping M_200^2 in the numerator.  At c = 4.3 that is a factor 2.5 in R_g.
* ``trunc``  -- whether the projected density is truncated by subtracting
  f(c) at the boundary, as ``projected_density_nfw_rmax`` does.

Note the Fortran offers Correa2015 / Ragagnin2018 / Duffy2008 concentrations
and NOT Maccio+2008, which the 2014 paper quotes -- so that file is a later
version of the code than the DR10 catalogue we compare against.
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.integrate import quad


def m_of_x(x):
    """Enclosed-mass shape function, ln(1+x) - x/(1+x)."""
    return np.log1p(x) - x / (1.0 + x)


def P_of_y(y):
    """Energy integral int_0^y m(x)/(1+x)^2 dx, in closed form."""
    return 0.5 - np.log1p(y) / (1.0 + y) - 0.5 / (1.0 + y) ** 2


def f_proj(s):
    """Dimensionless projected NFW surface density; Sigma = 2 rho_s R_s f(s)."""
    s = np.asarray(s, dtype=float)
    out = np.empty_like(s)
    hi, lo = s > 1.0, s < 1.0
    eq = ~(hi | lo)
    with np.errstate(divide="ignore", invalid="ignore"):
        sh = s[hi]
        out[hi] = (1.0 - 2.0 / np.sqrt(sh ** 2 - 1.0)
                   * np.arctan(np.sqrt((sh - 1.0) / (sh + 1.0)))) / (sh ** 2 - 1.0)
        sl = s[lo]
        out[lo] = (1.0 - 2.0 / np.sqrt(1.0 - sl ** 2)
                   * np.arctanh(np.sqrt((1.0 - sl) / (1.0 + sl)))) / (sl ** 2 - 1.0)
    out[eq] = 1.0 / 3.0
    return out


def sigma_sky_over_rs(c, trunc=True):
    """sigma_sky / R_s for a halo of concentration c, projected inside R_200."""
    fc = float(f_proj(np.array([float(c)]))[0]) if trunc else 0.0
    g = lambda s: max(float(f_proj(np.array([s]))[0]) - fc, 0.0)
    num = quad(lambda s: g(s) * s ** 3, 0.0, c, limit=200)[0]
    den = quad(lambda s: g(s) * s, 0.0, c, limit=200)[0]
    return np.sqrt(0.5 * num / den)


def rg_over_rs(c, x_out_factor=1.0):
    """R_g / R_s = m(c)^2 / P(x_out), with x_out = x_out_factor * c."""
    return m_of_x(c) ** 2 / P_of_y(x_out_factor * c)


def kappa(c, x_out_factor=1.0, trunc=True):
    """The shape factor R_g / sigma_sky.  R_s cancels, so it depends on c only."""
    return rg_over_rs(c, x_out_factor) / sigma_sky_over_rs(c, trunc)


def c_maccio(m200_h):
    """Maccio+2008 eq. 12 as quoted by Tempel+2014, M in h^-1 Msun."""
    return 10.0 ** (0.83 - 0.098 * np.log10(np.asarray(m200_h) / 1e12))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compare", action="store_true",
                   help="check against the empirically calibrated table in "
                        "nessie_sdss_hmf.KAPPA_NFW_VAL")
    args = p.parse_args()

    print("  Reproducing the appendix's table 1 (x_out = c):")
    print(f"  {'c':>5s} {'m(c)':>8s} {'P(c)':>8s} {'Rg/Rs':>8s} {'Rg/R200':>8s}")
    for c in (3.0, 4.3, 6.0, 10.0):
        r = rg_over_rs(c)
        print(f"  {c:5.1f} {m_of_x(c):8.4f} {P_of_y(c):8.4f} {r:8.2f} {r / c:8.2f}")

    print("\n  kappa at c = 4.3 under the four convention choices:")
    for xf, lab in ((1.0, "x_out = c   (paper)"), (10.0, "x_out = 10c (Fortran)")):
        for tr in (False, True):
            k = kappa(4.3, xf, tr)
            print(f"    {lab:24s} {'truncated' if tr else 'untruncated':12s}"
                  f"  kappa = {k:6.3f}   M_Her/M_NFW = {4.582 / k:5.3f}")
    print("    Tempel+22 fig. 7 quotes M_Her/M_NFW = 1.55-1.75")

    if args.compare:
        import nessie_sdss_hmf as nsd
        print("\n  vs the empirical calibration (inverted from his published col15)")
        print(f"  {'logM':>6s} {'c':>6s} | " + " | ".join(
            f"{lab:>14s}" for lab in ("x=c untrunc", "x=c trunc", "x=10c untr")))
        print(f"  {'':6s} {'':6s} | " + " | ".join(
            f"{'kappa   ratio':>14s}" for _ in range(3)) + f"   {'empirical':>9s}")
        for lm, ke in zip(nsd.KAPPA_NFW_LOGM, nsd.KAPPA_NFW_VAL):
            c = c_maccio(10.0 ** lm)
            row = f"  {lm:6.2f} {c:6.3f} | "
            for xf, tr in ((1.0, False), (1.0, True), (10.0, False)):
                k = kappa(c, xf, tr)
                row += f"{k:6.3f} {k / ke:6.3f} | "
            print(row + f"  {ke:7.4f}")


if __name__ == "__main__":
    main()
