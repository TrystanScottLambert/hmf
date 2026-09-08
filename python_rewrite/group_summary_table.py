#!/usr/bin/env python3
"""Multiplicity summary of the two Nessie group catalogues -> LaTeX.

Writes ``group_summary_statistics.tex``: the ``tabular`` block only, no
preamble and no ``document`` wrapper, ready to \\input into the paper.

Both catalogues are read in full -- no redshift cut, no mass cut, no
multiplicity floor beyond the column itself.  These are catalogue statistics,
not the HMF sample (which additionally requires N >= 5, 0.015 < z < 0.25 and
MassAfunc > 10 for GAMA, and N >= 5, 0.015 < z < 0.08 for SDSS).

**Thresholds are inclusive**: the "N >= 3" column counts groups with three or
more members.  Driver+22 writes the same cut as "N > 2"; pass ``--exclusive``
to use the strict reading instead, which shifts every column by one.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from astropy.table import Table

# Same files the catalogue figures use (plot_catalogues.py).
GAMA_GROUPS = "/Users/00115372/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"
SDSS_GROUPS = ("/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS"
               "/sdss_groups.parquet")

THRESHOLDS = [2, 3, 5, 20]


GAMA_REGIONS = ["g09", "g12", "g15", "g23"]


def gama_multiplicity():
    """Nfof and region for every group in the DMU."""
    t = Table.read(GAMA_GROUPS)
    n = np.asarray(t["Nfof"]).astype(int)
    reg = np.asarray(t["GAMARegion"])
    if reg.dtype.kind == "S":                 # FITS gives bytes
        reg = np.char.decode(reg)
    return n, np.asarray(reg).astype(str)


def sdss_multiplicity():
    return pd.read_parquet(SDSS_GROUPS)["multiplicity"].values.astype(int)


def counts(n, thresholds, exclusive=False):
    return [int((n > t).sum()) if exclusive else int((n >= t).sum())
            for t in thresholds]


def rows_for_table():
    """(label, multiplicities, rule_after) in the order they appear."""
    n_g, reg = gama_multiplicity()
    out = [(r.upper(), n_g[reg == r], False) for r in GAMA_REGIONS]
    out.append(("GAMA (all)", n_g, True))      # total after its fields
    out.append(("SDSS", sdss_multiplicity(), False))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exclusive", action="store_true",
                   help="count N > threshold rather than N >= threshold")
    p.add_argument("--out", default="group_summary_statistics.tex")
    args = p.parse_args()

    rows = rows_for_table()
    op = ">" if args.exclusive else r"\geq"

    print(f"  multiplicity counts ({'N >' if args.exclusive else 'N >='} threshold)")
    hdr = "  " + f"{'Survey':<12s}" + "".join(f"{t:>10d}" for t in THRESHOLDS)
    print(hdr + f"{'total':>10s}")
    body = []
    for name, n, rule in rows:
        c = counts(n, THRESHOLDS, args.exclusive)
        print("  " + f"{name:<12s}" + "".join(f"{v:>10d}" for v in c)
              + f"{len(n):>10d}")
        body.append((name, c, rule))

    col = "l" + "r" * len(THRESHOLDS)
    lines = [
        r"\begin{tabular}{" + col + "}",
        r"\hline",
        "Survey & " + " & ".join(rf"$N {op} {t}$" for t in THRESHOLDS) + r" \\",
        r"\hline",
    ]
    for name, c, rule in body:
        lines.append(f"{name} & " + " & ".join(f"{v:,}".replace(",", r"\,")
                                               for v in c) + r" \\")
        if rule:                       # rule under the GAMA total
            lines.append(r"\hline")
    lines += [r"\hline", r"\end{tabular}"]

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  wrote {args.out}")


if __name__ == "__main__":
    main()
