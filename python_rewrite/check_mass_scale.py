"""
============================================================
Are the real fit and the completeness table on the same mass scale?
============================================================

measure_completeness_nessie.py builds mock masses as

    MassA = mass_proxy * MASS_A          (MASS_A = 10, make_gama_dmu/config.py)

and mlim(z), Delta and therefore C(Delta,z) all follow from that.

load_real_gama() does NOT read MassA. It rebuilds the mass from VelDisp and
Rad50 and then divides by 10^masscorr(Nfof) -- Driver's old multiplicity
debiasing, which the new DMU does not apply (it offers MassAfunc instead).

If those two differ, Delta is computed on one scale in the mock and another in
the data, and the completeness correction is applied at the wrong masses.

Run:  python check_mass_scale.py [G3CFoFGroup.fits]
============================================================
"""

import sys
import numpy as np
import recovery as R


def main():
    path = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"
    )
    from astropy.io import fits as afits

    with afits.open(path) as h:
        t = h[1].data
    cols = set(t.columns.names)
    print(f"  {path}")
    print(
        f"  columns present: MassA={('MassA' in cols)}, "
        f"MassAfunc={('MassAfunc' in cols)}, MassProxy={('MassProxy' in cols)}"
    )

    Nfof = np.asarray(t["Nfof"], float)
    Zfof = np.asarray(t["Zfof"], float)
    sel = (Nfof >= R.MULTI) & (Zfof > R.ZMIN) & (Zfof < R.ZLIMIT)

    # what recovery.py actually fits
    lm_rec, sig, z, nf = R.load_real_gama(path)

    # the catalogue's own columns
    out = {}
    for c in ("MassA", "MassAfunc", "MassProxy"):
        if c in cols:
            v = np.asarray(t[c], float)[sel]
            v = v[np.isfinite(v) & (v > 0)]
            out[c] = np.log10(v)

    print(f"\n  {'definition':>26} {'median':>8} {'16%':>8} {'84%':>8} {'N':>6}")
    print(
        f"  {'recovery.py (rebuilt)':>26} {np.median(lm_rec):8.3f} "
        f"{np.percentile(lm_rec, 16):8.3f} {np.percentile(lm_rec, 84):8.3f} "
        f"{lm_rec.size:6d}"
    )
    for c, v in out.items():
        print(
            f"  {'log10 ' + c:>26} {np.median(v):8.3f} "
            f"{np.percentile(v, 16):8.3f} {np.percentile(v, 84):8.3f} {v.size:6d}"
        )

    if "MassA" in out and out["MassA"].size == lm_rec.size:
        d = lm_rec - out["MassA"]
        print(f"\n  recovery.py - log10(MassA):")
        print(
            f"    median {np.median(d):+.3f} dex   16-84% "
            f"[{np.percentile(d, 16):+.3f}, {np.percentile(d, 84):+.3f}]"
        )
        if abs(np.median(d)) < 0.005 and d.std() < 0.005:
            print("    -> identical: the two mass definitions agree.")
        else:
            print("    -> DIFFERENT. The completeness table was measured against")
            print("       MassA, so Delta = m - mlim(z) is not the same quantity")
            print("       in the mock and in the fit.")
            # is it the masscorr?
            mc = np.array(
                [
                    0.0,
                    0.0,
                    -2.672595e-01,
                    -1.513503e-01,
                    -1.259069e-01,
                    -9.006064e-02,
                    -5.466009e-02,
                    -6.666895e-02,
                    -1.988694e-02,
                    -2.439581e-02,
                    -2.067060e-02,
                    -1.812964e-02,
                    -1.556899e-02,
                    -1.313664e-02,
                    -1.743112e-02,
                    -7.965513e-03,
                    -1.257178e-02,
                    -7.064037e-03,
                    -3.963656e-03,
                    -1.271533e-02,
                    -2.664687e-03,
                    -1.691287e-03,
                ]
            )
            idx = np.clip(nf.astype(int) - 1, 0, mc.size - 1)
            pred = -mc[idx]
            print(
                f"\n    masscorr alone predicts a shift of median "
                f"{np.median(pred):+.3f} dex (range {pred.min():+.3f} to "
                f"{pred.max():+.3f})"
            )
            resid = d - pred
            print(
                f"    after removing masscorr, residual median "
                f"{np.median(resid):+.3f} dex, scatter {resid.std():.3f}"
            )
            if abs(np.median(resid)) < 0.02 and resid.std() < 0.02:
                print("    -> the whole difference is masscorr. Use MassA directly")
                print("       (--mass-col MassA) to match the completeness table.")
            else:
                print("    -> not just masscorr: the sigma^2 R reconstruction itself")
                print("       differs from the catalogue's mass_proxy. Check units,")
                print("       the radius column, and the h convention.")

    print(f"\n  A is a pure multiplicative constant: M = A sigma^2 R / G, so")
    print(f"  changing A shifts every log-mass by log10(A/A_old) and moves M* by")
    print(f"  exactly the same amount. alpha, beta and phi* are unaffected.")
    print(f"  The mock closed loop measures the resulting M* bias, so any A works")
    print(f"  provided the SAME A is used for the mock and the data.")


if __name__ == "__main__":
    main()
