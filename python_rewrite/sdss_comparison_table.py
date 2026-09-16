#!/usr/bin/env python3
"""Group multiplicities in redshift slices: Tempel+14 SDSS DR10 vs Nessie.

Writes ``sdss_comparison_statistics.tex``: the ``tabular`` block only, ready to
\\input into the paper, in the same style as ``group_comparison_statistics.tex``.

Unlike the GAMA comparison this one is like-for-like by construction: the Nessie
galaxy file *is* Tempel's table 1 (all 584,447 rows match on RA/Dec with
identical ``rmag``), so the footprint, the flux limit and the input redshifts are
the same and the only thing that differs is the group finder.  There is
therefore no field axis to split on, and the rows are redshift slices instead --
which is where the finder difference lives.

Redshifts are each catalogue's own group redshift, as the HMF legs use them:
Tempel's ``col10`` (``sdss_hmf.GRP_COLS['zcl']``) and Nessie's
``median_redshift`` (``nessie_sdss_hmf`` line 249).  The 0.015 < z < 0.08 row is
the HMF selection range, and its N >= 5 entries are the 4,847 and 4,824 groups
the two SDSS legs fit.

Catalogue statistics otherwise -- no mass cut and no multiplicity floor beyond
the column itself.  Thresholds are inclusive (``--exclusive`` for N > threshold),
as in ``group_summary_table.py``.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

from group_summary_table import SDSS_GROUPS, THRESHOLDS, counts
from sdss_hmf import ZLIMIT, ZMIN

TEMPEL_GROUPS = "../data/sdssdr10table2.fits"

SLICES = [(0.00, 0.02), (0.02, 0.04), (0.04, 0.06), (0.06, 0.08)]


def tempel_multiplicity():
    """nrich and group redshift for every Tempel+14 group."""
    t = Table.read(TEMPEL_GROUPS)
    return (np.asarray(t["col2"]).astype(int), np.asarray(t["col10"]))


def nessie_multiplicity():
    """multiplicity and group redshift for every Nessie SDSS group."""
    d = pd.read_parquet(SDSS_GROUPS)
    return (d["multiplicity"].values.astype(int), d["median_redshift"].values)


def rows_for_table(exclusive=False):
    """(label, Tempel counts, Nessie counts, rule_after)."""
    n_old, z_old = tempel_multiplicity()
    n_new, z_new = nessie_multiplicity()

    def cnt(n, z, lo, hi):
        m = (z >= lo) & (z < hi)
        return counts(n[m], THRESHOLDS, exclusive)

    out = [(rf"${lo:.2f} \leq z < {hi:.2f}$",
            cnt(n_old, z_old, lo, hi), cnt(n_new, z_new, lo, hi), False)
           for lo, hi in SLICES]
    out[-1] = out[-1][:3] + (True,)
    out.append((rf"${ZMIN} < z < {ZLIMIT}$",
                cnt(n_old, z_old, ZMIN, ZLIMIT),
                cnt(n_new, z_new, ZMIN, ZLIMIT), False))
    out.append(("All groups",
                counts(n_old, THRESHOLDS, exclusive),
                counts(n_new, THRESHOLDS, exclusive), False))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exclusive", action="store_true",
                   help="count N > threshold rather than N >= threshold")
    p.add_argument("--out", default="sdss_comparison_statistics.tex")
    args = p.parse_args()

    rows = rows_for_table(args.exclusive)
    op = ">" if args.exclusive else r"\geq"
    k = len(THRESHOLDS)

    print(f"  multiplicity counts ({'N >' if args.exclusive else 'N >='} threshold)"
          "   Tempel+14 | Nessie")
    for name, old, new, _ in rows:
        fmt = lambda c: "".join(f"{v:>8d}" for v in c)
        label = name.replace("$", "").replace(r"\leq", "<=")
        print(f"  {label:<18s}{fmt(old)}  |{fmt(new)}")

    def cells(c):
        return [f"{v:,}".replace(",", r"\,") for v in c]

    thr = " & ".join(rf"$N {op} {t}$" for t in THRESHOLDS)
    lines = [
        r"\begin{tabular}{l" + "r" * k + "r" * k + "}",
        r"\hline",
        rf" & \multicolumn{{{k}}}{{c}}{{Tempel et al. (2014)}} & "
        rf"\multicolumn{{{k}}}{{c}}{{This work}} \\",
        f"Redshift & {thr} & {thr} " + r"\\",
        r"\hline",
    ]
    for name, old, new, rule in rows:
        lines.append(f"{name} & " + " & ".join(cells(old) + cells(new)) + r" \\")
        if rule:
            lines.append(r"\hline")
    lines += [r"\hline", r"\end{tabular}"]

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  wrote {args.out}")


if __name__ == "__main__":
    main()
