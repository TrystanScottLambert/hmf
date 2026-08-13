"""GAMA (Nessie DMU) halo mass function: hierarchical Bayes vs 1/Vmax.

Everything is drawn in DRIVER'S published units (A = 13.9, h = 0.6737,
Omega_M = 0.3147, M in Msun, phi in Mpc^-3 dex^-1).  That is deliberate.

The alternative -- drawing in recovery.py's (A = 10, h = 1) system -- needs FOUR
separate conversions on this one figure, each with a different rule:

    Driver GSR / GAMA-only  dynamical mass fit : h AND A
    Murray+21 LCDM          true-mass theory   : h only, no A
    1/Vmax binned points    mass and density   : h, A, and 3 log10 h on phi
    posterior                                  : identity (native)

Doing those by hand at the call site produced four bugs in a row.  Drawing in
Driver's units leaves exactly one conversion -- the posterior -- routed through
recovery.driver_to_recovery, and lets the three published curves and the binned
points be plotted natively with no arithmetic at all.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import driver_recovery as dr
import recovery as R

DRAWS = "gama_marg_tab_draws.csv"
OUT = "hmf_gama_bayes.pdf"

# Driver+22 published fits, in his own units -- plotted natively, no conversion.
DRIVER_GSR = (14.13, -3.96, -1.68, 0.63)  # combined GAMA+SDSS+REFLEX (abstract)
DRIVER_GAMA5 = (13.51, -3.19, -1.27, 0.47)  # his GAMA-only fit

# Closed-loop bias measured on the Shark/Nessie mock, in recovery units.
NESSIE_BIAS = dict(ms=-0.229, lp=+0.418, al=+0.011, be=-0.037)


def to_driver(ms, lp, al, be):
    """recovery units -> Driver's units (inverse of driver_to_recovery)."""
    dM, dPHI = R.driver_to_recovery("dynamical")
    return ms - dM, lp - dPHI, al, be


def band(x, draws, lo=16, hi=84):
    y = np.array([np.log10(R.mrp_phi(x, *p)) for p in draws])
    return np.percentile(y, lo, axis=0), np.percentile(y, hi, axis=0)


def main():
    d = np.genfromtxt(DRAWS, delimiter=",", names=True)
    raw = np.column_stack([d["ms"], d["lp"], d["al"], d["be"]])
    cor = raw - np.array([NESSIE_BIAS[k] for k in ("ms", "lp", "al", "be")])
    raw_d = np.array([to_driver(*p) for p in raw])
    cor_d = np.array([to_driver(*p) for p in cor])

    x = np.linspace(12.6, 15.5, 300)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    lx, ly, _ = dr.lcdm_curve_xy() if hasattr(dr, "lcdm_curve_xy") else (None,) * 3
    if lx is None:
        mx, my, fac, _ = dr.lcdm_curve()
        lx, ly = mx - 0.08, np.log10(my) - np.log10(fac) + 0.08
    ax.plot(lx, ly, "k--", lw=1.6, label=r"$\Lambda$CDM (Murray+21, $z=0.1$)")
    ax.plot(x, np.log10(R.mrp_phi(x, *DRIVER_GSR)), ":", c="#c8781e", lw=1.8,
            label="Driver+22 GSR")
    ax.plot(x, np.log10(R.mrp_phi(x, *DRIVER_GAMA5)), "--", c="#a03020", lw=1.4,
            label="Driver+22 GAMA-only")

    for dd, c, lab in ((raw_d, "#4878a8", "This work, hierarchical Bayes"),
                       (cor_d, "#5a9367", "  + closed-loop bias correction")):
        lo, hi = band(x, dd)
        ax.fill_between(x, lo, hi, color=c, alpha=0.35, lw=0, label=lab)

    # 1/Vmax on the SAME catalogue, built natively in Driver's units -- no
    # conversion applied or needed.  M/L cut as elsewhere (fixes the 14.2 bin).
    import new_gama_hmf as ng

    g3c, vlimit = ng.build_groups_new(verbose=False)
    g3c, _ = ng.apply_ml_cut(g3c, 1.0, verbose=False)
    b, _ = dr.bin_hmf(g3c, vlimit, seed=10, nmc=1001, verbose=False, nboot=200)
    ok = b.gamay > 0
    ax.errorbar(b.gamax[ok], np.log10(b.gamay[ok]),
                yerr=0.4343 * b.gamaf[ok], fmt="o", ms=5, mfc="w", mec="k",
                ecolor="k", lw=1, zorder=5, capsize=2,
                label="1/Vmax, same catalogue (not completeness corrected)")

    ax.set_xlim(12.6, 15.5)
    ax.set_ylim(-7.5, -2.0)
    ax.set_xlabel(r"$\log_{10}(M_{\rm halo}/M_\odot)$   [Driver+22 units: $A=13.9$, $h=0.6737$]")
    ax.set_ylabel(r"$\log_{10}\phi$  [Mpc$^{-3}$ dex$^{-1}$]")
    ax.set_title("GAMA (Nessie DMU) halo mass function")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT)
    print(f"wrote {OUT}")

    p = np.percentile(raw_d, [16, 50, 84], axis=0)
    print("\nposterior in Driver units (16/50/84):")
    for i, n in enumerate(("logM*", "logphi*", "alpha", "beta")):
        print(f"  {n:8s} {p[1, i]:8.3f}  +{p[2, i] - p[1, i]:.3f} -{p[1, i] - p[0, i]:.3f}")


if __name__ == "__main__":
    main()
