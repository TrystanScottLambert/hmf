#!/usr/bin/env python3
"""Seed scan + penalty-range check for the Nessie SDSS leg.

CLAUDE.md's rule: never quote a fitted parameter from one Monte-Carlo seed.
This builds both group catalogues once, then re-bins and re-fits per seed.
``--fit-max`` sensitivity is checked in the same pass because Nessie reaches
higher mass than Tempel and ``make_massfn`` integrates its penalty over
``max(allx)``.
"""
from __future__ import annotations

import argparse

import numpy as np

import nessie_sdss_hmf as ns
import sdss_hmf
from driver_recovery import lcdm_curve
from sdss_hmf import bin_hmf

NAMES = ["log10(M*)", "log10(phi*)", "alpha", "beta"]


def run(grp, vol, phimrp, seeds, nmc, nboot, fit_max, fitter):
    out = []
    for s in seeds:
        b = bin_hmf(grp, vol, seed=s, nmc=nmc, verbose=False, nboot=nboot)
        if fitter is ns.fit:
            par, val, _, conv, n = fitter(b, vol, phimrp, fit_max=fit_max)
        else:
            par, val, _, conv, n = fitter(b, vol, phimrp)
        out.append([par[0], np.log10(abs(par[1])), par[2], par[3], val, conv, n])
    return np.array(out)


def summarise(a, label):
    print(f"\n  {label}  (n={len(a)})")
    for i, nm in enumerate(NAMES):
        print(f"    {nm:<13} median {np.median(a[:, i]):8.3f}   "
              f"scatter {a[:, i].std(ddof=1):6.3f}   "
              f"[{a[:, i].min():.3f}, {a[:, i].max():.3f}]")
    print(f"    chi2          median {np.median(a[:, 4]):8.3f}   "
          f"converged {int((a[:, 5] == 0).sum())}/{len(a)}   bins {int(a[0, 6])}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--nmc", type=int, default=1001)
    p.add_argument("--nboot", type=int, default=0)
    p.add_argument("--fit-max", type=float, default=None)
    args = p.parse_args()

    seeds = list(range(1, args.seeds + 1))
    _, _, _, phimrp = lcdm_curve()

    print("Building both catalogues ...")
    ngrp, nvol = ns.build_groups(verbose=False)
    ogrp, ovol = sdss_hmf.build_groups(verbose=False)
    print(f"  Nessie {len(ngrp)} groups, Tempel {len(ogrp)} groups")
    print(f"  seeds {seeds}, nmc={args.nmc}, nboot={args.nboot}, "
          f"fit_max={args.fit_max}")

    nes = run(ngrp, nvol, phimrp, seeds, args.nmc, args.nboot, args.fit_max,
              ns.fit)
    tem = run(ogrp, ovol, phimrp, seeds, args.nmc, args.nboot, None,
              sdss_hmf.fit)

    summarise(tem, "Tempel (Driver's SDSS)")
    summarise(nes, "Nessie SDSS")

    print("\n  difference of medians (Nessie - Tempel):")
    for i, nm in enumerate(NAMES):
        d = np.median(nes[:, i]) - np.median(tem[:, i])
        pooled = np.hypot(nes[:, i].std(ddof=1), tem[:, i].std(ddof=1))
        sig = f"{abs(d) / pooled:.1f}x scatter" if pooled > 0 else "-"
        print(f"    {nm:<13}{d:+9.3f}   ({sig})")


if __name__ == "__main__":
    main()
