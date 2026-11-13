# regen/ — the FIGURES.md deliverable set, re-run on the rebuilt catalogues (2026-11-13)

After the 2026-11 data loss (see CLAUDE.md, "READ FIRST"), the Nessie catalogues
had to be rebuilt from unseeded randoms. Everything here was produced by
`run_all.sh` with the exact FIGURES.md commands, only `--outdir regen` and the
mass-error paths changed. Nothing in the parent directory was overwritten.

Inputs:
* Nessie GAMA DMU rebuilt with the new randoms and `OVERFACTOR = 100`: 1807 groups (was 1833).
* Nessie SDSS parquets regenerated with nessie-py 1.1.1 (sqrt(3) fix): 4823 groups (was 4824).
* `gama_masserr.csv`, `sdss_masserr.csv` — remeasured on the rebuilt catalogues.
  The SDSS curve agrees with the committed one to <= 0.004 dex.

## Two directories, because two things changed

| | REFLEX conversion | compare with |
|---|---|---|
| `regen/` | **corrected** (current default since 2026-09-07) | nothing committed — new result |
| `regen/reflex_as_published/` | Driver's `allhmf.r` line 307 (+0.638 dex) | the committed `chain_combined_*.npz` (2026-09-01), like-for-like |

The committed combined chains PREDATE the REFLEX correction. Comparing them with
`regen/` mixes the rebuild with the REFLEX change; use `reflex_as_published/`
to isolate the rebuild.

### The rebuild alone (like-for-like)

| | logM* | alpha | beta | chi2 | Omega_M |
|---|---|---|---|---|---|
| GSR, Nessie GAMA + Tempel SDSS, committed | 13.983 | -1.618 | 0.583 | 230.7 | 0.199 |
| GSR, Nessie GAMA + Tempel SDSS, **rebuilt** | 13.874 | -1.573 | 0.552 | 228.3 | 0.189 |
| GSR, all Nessie, committed | 12.792 | -0.986 | 0.390 | 215.7 | 0.131 |
| GSR, all Nessie, **rebuilt** | 12.522 | -0.881 | 0.365 | 211.2 | 0.128 |

Every shift is well inside the 16/84 intervals (~0.3 sigma). The rebuild does
not change the science.

### With the corrected REFLEX (the current default)

Every native-units combined row RAILS (`mrp_table.txt`): anchored GSR chi2
525, all-Nessie 654. This is the calibration spread CLAUDE.md predicted ("Two
errors that cancelled"): corrected REFLEX sits ~0.6 dex below the dynamical
legs, and one MRP cannot pass through both. It is a real result, not a
regression — and it is why outstanding item 0 (a fitted per-leg mass offset)
is now the main open problem.

### Shape-only fits (already like-for-like: committed chains are 2026-09-14)

Anchored row: logM* 12.52 vs 12.61, alpha -0.95 vs -1.03, chi2 197.2 vs 174.4.
GAMA-only moved most (alpha -0.40 vs -0.80) but rails in both and was never
quotable.

## Files

`hmf_*`, `corner_*`, `chain_*` — the seven deliverables; `mrp_table*.txt/.tex`
— tables built from these chains (`mrp_table.py --chain-dir regen`);
`logs/` — full stdout of each run; `vuvuzela_sdss.png` — the remeasured SDSS
mass-error curve.

Re-run: `regen/run_all.sh` (~3 min, 7 chains in parallel) and
`regen/run_all.sh --reflex-as-published` (~2 min). emcee is unseeded, so a
re-run moves medians by a few hundredths.

## `old_conventions/` — what the committed figures actually looked like

The committed figures (2026-09-01) were drawn with BOTH of Driver's conventions:
his REFLEX normalisation (`--reflex-as-published`) and his Murray "LCDM" curve
(`--lcdm driver`). Both became non-default on 2026-09-07. `regen/` therefore
looks very different from the committed figures for reasons that have nothing
to do with the rebuilt data: REFLEX drops 0.638 dex and the black line drops
~0.2-0.5 dex, leaving every dynamical leg visibly above it.

`run_all.sh --old-conventions` re-runs the four native-units figures with both
flags. They reproduce the committed figures (anchored GSR 13.888 / -1.579 /
chi2 228.3 / Omega_M 0.190, against 13.983 / -1.618 / 230.7 / 0.199).
Cosmetic differences only: the black line's legend now reads "MRP normalised to
Omega_M", and there is an extra "same fit using Driver+22 GAMA" curve, which
was added to the code after 2026-09-01.

## `mlcut08/` — the paper's Figs 5-7 with the M/L cut at 0.8 dex

On the rebuilt DMU the 205509 pathology reappears as **group 205415** (N = 6,
z = 0.020, sigma = 378 km/s, RA 181.4 Dec 1.58, Vmax 0.27% of the survey): the
same loose z ~ 0.02 association as the old 205509/205469, split differently by
the new randoms. Its M/L excess is 0.83 dex, so the 1.0 cut misses it, and it
carries ~1/3 of the 14.0 bin (-3.673; -4.019 without it; the lost catalogue had
-4.070). Because 14.0 is the edge of the 14-15 matching window it inflates
GAMA's normalisation integral and drags every rescaled GAMA point down -- that
is why `regen/hmf_shape_gama.png` misses the Murray curve above 14.

With `--ml-cut 0.8` the shape fits return to the paper's Table 4:

| | paper Table 4 | regen, cut 1.0 | regen, cut 0.8 |
|---|---|---|---|
| GAMA | 12.52 / -0.78 / 13.4 | 11.90 / -0.40 / 19.0 | 12.13 / -0.53 / 14.2 |
| GAMA + SDSS | 11.68 / -0.97 / 108.2 | 11.59 / -0.90 / 127.4 | 11.70 / -0.96 / 112.2 |
| GAMA + SDSS + REFLEX | 12.60 / -1.02 / 185.9 | 12.52 / -0.95 / 197.2 | 12.70 / -1.03 / 184.7 |

(logM* / alpha / chi2.) The 0.8 threshold was chosen after seeing the figure;
the paper must say so, and quote the number of groups removed.
