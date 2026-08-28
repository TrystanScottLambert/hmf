#!/usr/bin/env python3
"""MRP best-fit parameters for every deliverable, in Driver+22 table 2 style.

Reads the saved chains (``chain_<name>.npz``) and reports the posterior median
with 16/84 credible intervals, which is how Driver+22 quote his table 2.  No
resampling: the chains are the record.

**Railed rows are marked, not tabulated as measurements.**  Several
configurations in this project have posteriors that pile against the prior
bound at ``log10 M* = 11.0`` (see mcmc_hmf.BOUNDS).  There the credible
interval is set by the prior and not by the data, so quoting it as an
uncertainty would be wrong -- the binned points and the plotted curve are still
sound, only the four MRP parameters are not.

Usage
-----
    uv run python mrp_table.py                 # text
    uv run python mrp_table.py --latex         # LaTeX for the paper
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from mcmc_hmf import BOUNDS

# name -> (chain file, row label as it should appear in the paper)
DELIVERABLES = [
    ("gama_nessie", "GAMA (Nessie)"),
    ("gama_sdss_nessie", "GAMA + SDSS (Nessie)"),
    ("combined_nessie_gama", "GAMA (Nessie) + SDSS (Tempel+14) + REFLEX II"),
    ("combined_nessie_all", "GAMA + SDSS (Nessie) + REFLEX II"),
]

# Driver+22 table 2, for comparison.
PUBLISHED = [
    ("GAMA5 (Driver+22)", (13.51, -3.19, -1.27, 0.47),
     [(1.37, 3.22), (0.67, 0.61), (0.19, 0.19), (0.08, 0.08)]),
    ("GSR (Driver+22)", (14.13, -3.96, -1.68, 0.63),
     [(0.43, 0.40), (0.55, 0.82), (0.21, 0.24), (0.25, 0.11)]),
]

PARAMS = [r"log$_{10}$(M$_*$/M$_\odot$)", r"log$_{10}(\phi_*)$",
          r"$\alpha$", r"$\beta$"]
PLAIN = ["log10(M*)", "log10(phi*)", "alpha", "beta"]


def railed(chain, tol=0.02):
    """True if any parameter's posterior piles against a prior bound.

    Tested on the 2nd percentile rather than the minimum: a handful of walkers
    touching the edge is not the same as the distribution being cut off there.
    """
    for i in range(chain.shape[1]):
        lo, hi = BOUNDS[i]
        span = hi - lo
        if (np.percentile(chain[:, i], 2) - lo) < tol * span:
            return True
        if (hi - np.percentile(chain[:, i], 98)) < tol * span:
            return True
    return False


def omega_total(d, ndraw=4000, seed=0):
    """Total Omega_M: median and the 16/84 offsets.

    Prefers ``omega_draws`` saved by ``combined_hmf`` -- those are the *same*
    draws the figure inset histograms, so the table and the figure cannot
    disagree.  Re-subsampling the chain here instead produced a table that
    differed from its own figure in the third decimal.

    Falls back to recomputing for chains written before that was saved.
    ``omega_matter`` takes **linear** phi; passing log10 phi returns a silently
    nonsense density.
    """
    if "omega_draws" in d.files:
        o = np.asarray(d["omega_draws"])
    else:
        import combined_hmf as ch      # lazy: importing it is not free
        chain = d["chain"]
        rng = np.random.default_rng(seed)
        i = rng.choice(len(chain), min(ndraw, len(chain)), replace=False)
        lin = np.column_stack([chain[i, 0], 10.0 ** chain[i, 1],
                               chain[i, 2], chain[i, 3]])
        o = np.array([ch.omega_matter(p) for p in lin])
    o = o[np.isfinite(o)]
    med = np.median(o)
    return med, med - np.percentile(o, 16), np.percentile(o, 84) - med


def summarise(chain):
    """median and the 16/84 offsets for each parameter."""
    med = np.median(chain, axis=0)
    lo = med - np.percentile(chain, 16, axis=0)
    hi = np.percentile(chain, 84, axis=0) - med
    return med, lo, hi


def rows(verbose=True):
    out = []
    for name, label in DELIVERABLES:
        path = f"chain_{name}.npz"
        if not os.path.exists(path):
            if verbose:
                print(f"  (missing {path} -- run combined_hmf.py --name {name})")
            continue
        d = np.load(path)
        med, lo, hi = summarise(d["chain"])
        r = railed(d["chain"])
        # A railed posterior's Omega_M is set by the prior bound, not the data,
        # which is how one of these ends up "21 sigma below Planck".  Do not
        # compute it at all rather than print a number nobody may use.
        om = (None, None, None) if r else omega_total(d)
        out.append(dict(label=label, med=med, lo=lo, hi=hi,
                        omega=om[0], omega_lo=om[1], omega_hi=om[2],
                        railed=r,
                        chi2=float(d["chi2_best"]) if "chi2_best" in d.files
                        else None,
                        nsamp=len(d["chain"])))
    return out


def as_text(rs):
    print("\nMRP fits -- posterior median with 16/84 credible intervals\n")
    head = f"  {'Sample':<44s}" + "".join(f"{p:>22s}" for p in PLAIN)
    print(head + f"{'chi2':>9s}{'Omega_M (total)':>22s}")
    print("  " + "-" * (len(head) + 7))
    for r in rs:
        cells = "".join(
            f"{r['med'][i]:>10.3f} +{r['hi'][i]:.2f}/-{r['lo'][i]:.2f} "
            for i in range(4))
        c = f"{r['chi2']:>9.2f}" if r["chi2"] is not None else f"{'--':>9s}"
        om = (f"{r['omega']:>10.3f} +{r['omega_hi']:.3f}/-{r['omega_lo']:.3f}"
              if r["omega"] is not None else f"{'--':>22s}")
        flag = "  <- RAILED, interval set by the prior" if r["railed"] else ""
        print(f"  {r['label']:<44s}{cells}{c}{om}{flag}")
    print("\n  Driver+22 published, for comparison")
    for label, v, err in PUBLISHED:
        cells = "".join(f"{v[i]:>10.3f} +{err[i][0]:.2f}/-{err[i][1]:.2f} "
                        for i in range(4))
        import combined_hmf as ch
        o = ch.omega_matter(np.array([v[0], 10.0 ** v[1], v[2], v[3]]))
        print(f"  {label:<44s}{cells}{'--':>9s}{o:>22.3f}")
    print(f"  {'Planck 2018':<44s}" + " " * (23 * 4)
          + f"{'--':>9s}{0.3147:>22.3f}")
    print("\n  Omega_M is the TOTAL, i.e. the fitted MRP integrated over all "
          "mass -- the\n  same quantity the figure inset plots.  39% of it lies "
          "below the fitted\n  range, so it is an extrapolation, and it scales "
          "as 10^delta under a\n  uniform mass-scale shift: the +/-0.15 dex "
          "systematic on the mass\n  calibration is x1.41, far larger than the "
          "quoted statistical error.")


def as_latex(rs):
    print(r"\begin{table*}")
    print(r"\centering")
    print(r"\caption{MRP function fits.  Values are the posterior median with "
          r"16th/84th percentile credible intervals.  Rows marked $\dagger$ "
          r"have a posterior that rails against the prior boundary at "
          r"$\log_{10}(M_*/{\rm M}_\odot) = 11$; there the interval reflects "
          r"the prior rather than the data, and the parameters should not be "
          r"quoted, although the binned measurements and the plotted curve "
          r"remain valid.  $\Omega_M$ is the total, i.e. the fitted MRP "
          r"integrated over all mass; 39 per cent of it lies below the fitted "
          r"range and it scales as $10^{\delta}$ under a uniform mass-scale "
          r"shift, so the $\pm0.15$~dex calibration systematic ($\times1.41$) "
          r"dominates the quoted statistical error.}")
    print(r"\begin{tabular}{lcccccc}")
    print(r"\hline")
    print("Sample & " + " & ".join(PARAMS)
          + r" & $\chi^2$ & $\Omega_M$ \\")
    print(r" & & Mpc$^{-3}$ dex$^{-1}$ & & & & (total) \\")
    print(r"\hline")
    for r in rs:
        cells = " & ".join(
            f"${r['med'][i]:.2f}^{{+{r['hi'][i]:.2f}}}_{{-{r['lo'][i]:.2f}}}$"
            for i in range(4))
        c = f"{r['chi2']:.1f}" if r["chi2"] is not None else "--"
        om = (f"${r['omega']:.3f}^{{+{r['omega_hi']:.3f}}}"
              f"_{{-{r['omega_lo']:.3f}}}$" if r["omega"] is not None else "--")
        print(f"{r['label']}{r'$\dagger$' if r['railed'] else ''} & {cells} "
              f"& {c} & {om} " + r"\\")
    print(r"\hline")
    for label, v, err in PUBLISHED:
        cells = " & ".join(
            f"${v[i]:.2f}^{{+{err[i][0]:.2f}}}_{{-{err[i][1]:.2f}}}$"
            for i in range(4))
        import combined_hmf as ch
        o = ch.omega_matter(np.array([v[0], 10.0 ** v[1], v[2], v[3]]))
        print(f"{label} & {cells} & -- & {o:.3f} " + r"\\")
    print(r"\hline")
    print(r"\end{tabular}")
    print(r"\label{tab:mrp}")
    print(r"\end{table*}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--latex", action="store_true")
    args = p.parse_args()
    rs = rows(verbose=not args.latex)
    if not rs:
        raise SystemExit("no chains found")
    (as_latex if args.latex else as_text)(rs)


if __name__ == "__main__":
    main()
