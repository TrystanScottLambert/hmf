#!/usr/bin/env python3
"""Per-field group multiplicities: G3C v10 (Robotham+11 / Driver+22) vs Nessie.

Writes ``group_comparison_statistics.tex``: the ``tabular`` block only, ready to
\\input into the paper, in the same style as ``group_summary_statistics.tex``.

The two catalogues do not cover the same fields.  v10 has G02 (GroupID prefix
4) and no G23; the new DMU has G23 and no G02.  So the like-for-like comparison
is the G09+G12+G15 subtotal; each catalogue's full total is given after it.
v10 carries no region column, so its field comes from the GroupID prefix
(1 = G09, 2 = G12, 3 = G15, 4 = G02), verified against the group RA/Dec.

Catalogue statistics, not the HMF sample -- no redshift, mass or multiplicity
cut beyond the column itself.  Thresholds are inclusive (``--exclusive`` for
N > threshold), as in ``group_summary_table.py``.
"""

from __future__ import annotations

import argparse

import numpy as np
from astropy.table import Table

from group_summary_table import THRESHOLDS, counts, gama_multiplicity

V10_GROUPS = "../data/G3CFoFGroupv10.fits"
V10_PREFIX = {1: "g09", 2: "g12", 3: "g15", 4: "g02"}

FIELDS = ["g02", "g09", "g12", "g15", "g23"]
OVERLAP = ["g09", "g12", "g15"]


def v10_multiplicity():
    """Nfof and region for every group in G3CFoFGroupv10."""
    t = Table.read(V10_GROUPS)
    n = np.asarray(t["Nfof"]).astype(int)
    prefix = np.asarray(t["GroupID"]).astype(int) // 100000
    return n, np.array([V10_PREFIX[p] for p in prefix])


def rows_for_table(exclusive=False):
    """(label, v10 counts or None, Nessie counts or None, rule_after)."""
    n_old, reg_old = v10_multiplicity()
    n_new, reg_new = gama_multiplicity()

    def cnt(n, reg, fields):
        m = np.isin(reg, fields)
        return counts(n[m], THRESHOLDS, exclusive) if m.any() else None

    out = [(f.upper(), cnt(n_old, reg_old, [f]), cnt(n_new, reg_new, [f]), False)
           for f in FIELDS]
    out[-1] = out[-1][:3] + (True,)
    out.append(("G09+G12+G15", cnt(n_old, reg_old, OVERLAP),
                cnt(n_new, reg_new, OVERLAP), False))
    out.append(("All fields", counts(n_old, THRESHOLDS, exclusive),
                counts(n_new, THRESHOLDS, exclusive), False))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exclusive", action="store_true",
                   help="count N > threshold rather than N >= threshold")
    p.add_argument("--out", default="group_comparison_statistics.tex")
    args = p.parse_args()

    rows = rows_for_table(args.exclusive)
    op = ">" if args.exclusive else r"\geq"
    k = len(THRESHOLDS)

    print(f"  multiplicity counts ({'N >' if args.exclusive else 'N >='} threshold)"
          "   v10 | Nessie")
    for name, old, new, _ in rows:
        fmt = lambda c: "".join(f"{v:>8d}" for v in c) if c else f"{'--':>8s}" * k
        print(f"  {name:<12s}{fmt(old)}  |{fmt(new)}")

    def cells(c):
        if c is None:
            return ["--"] * k
        return [f"{v:,}".replace(",", r"\,") for v in c]

    thr = " & ".join(rf"$N {op} {t}$" for t in THRESHOLDS)
    lines = [
        r"\begin{tabular}{l" + "r" * k + "r" * k + "}",
        r"\hline",
        rf" & \multicolumn{{{k}}}{{c}}{{G$^3$C v10}} & "
        rf"\multicolumn{{{k}}}{{c}}{{This work}} \\",
        f"Field & {thr} & {thr} " + r"\\",
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
