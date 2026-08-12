# Driver+22's halo mass function: reproduction, and the same method on new data

## Status

**Step 1 (reproduce `gamahmf.r`) is done and verified bit-exactly against R.**
The binned HMF matches published Table 1; the fit matches R's `optim` to 10
significant figures; the combined multi-survey fit matches published Table 2 to
better than 0.035 in every parameter.

**Step 2 (new GAMA DMU) is done for GAMA.** The Nessie catalogue has been run
through the identical method, both GAMA-only and in the combined fit.

**Step 3 (Nessie SDSS) is now built** — `nessie_sdss_hmf.py`, 4824 groups. All
three blockers turned out to be answerable from the data. Like GAMA-only, the
SDSS-only *fit* is ill-posed; see "The Nessie SDSS leg".

**Step 4 (Nessie SDSS in the combined fit) is done** — two separate figures,
`--sdss auto` and `--sdss nessie`. Everything else is built and validated.

**Outstanding:** why the Nessie SDSS leg fails to anchor the combined fit (it
does not; Tempel's does). See "In the combined fit".

Do not "fix" the things the old version of this file listed as problems. Most of
them were either resolved or were never problems. Read "The central finding"
before changing anything.

---

## The central finding

`gamahmf.r` fits with `optim(..., maxit=500)`. That call returns
`convergence = 1` — the Nelder-Mead simplex is **still marching when the budget
runs out**. It is not a converged minimum, and this governs how every GAMA-only
number in the paper should be read.

Released from `maxit=500`, the GAMA-only fit does not find a better answer, it
leaves physics entirely:

| | maxit=500 | converged | bounded |
|---|---|---|---|
| Driver v10, Poisson | 13.493, α=−1.27 | **9.910, α=−0.15** | 12.500 (rail) |
| Driver v10, bootstrap | 13.415, α=−1.29 | **8.286, α=+0.16** | 12.500 (rail) |
| Nessie, Poisson | 12.240, α=−0.56 | **9.577, α=+0.42** | 12.500 (rail) |
| Nessie, bootstrap | 12.990, α=−1.21 | **10.138, α=−0.05** | 12.500 (rail) |

Every converged answer has a **lower** χ² than Driver's (16.79 vs 20.78 for
Driver/Poisson). The χ² surface decreases monotonically toward low M\* across the
whole physical range — there is no interior minimum to find. Bounded, it rails at
whatever lower limit you impose.

So **logM\* = 13.51 is the Murray+21 start point decaying, not a measurement.**
`maxit=500` is not a bug hiding a good answer; it is an accidental regulariser,
and it is the only thing producing a physical-looking number at all.

**Adding REFLEX II fixes it.** The combined fit has a genuine interior minimum:

```
GR  (GAMA+REFLEX)  maxit=500/1000/5000/50000 -> logM* = 13.3839 every time, conv=0
GSR (Driver's)     maxit=500  -> 14.150, conv=1
                   maxit=2000 -> 14.142, conv=0, 533 fevals
```

The x-ray data anchors the exponential cutoff, which pins β and M\*. Truncation
then costs 0.008 dex instead of 5.

**Consequence for the paper:** quote the combined fit. Treat GAMA-only MRP
parameters — Driver's and ours — as not measurable from GAMA alone. The GAMA-only
*binned points* are perfectly good; it is only the four-parameter fit to them that
is ill-posed.

This is also why the GAMA-only old/new comparison cannot settle whether the two
catalogues' faint-end slopes genuinely differ (see "The catalogues do not agree at
the faint end"): the seed-to-seed scatter is a substantial fraction of the
difference being measured.

---

## The Python pipeline

| file | role |
|---|---|
| `driver_recovery.py` | Fig. 4 reproduction, and the shared R-compatibility layer everything else imports |
| `new_gama_hmf.py` | Nessie GAMA DMU through the identical method; old-vs-new GAMA-only figure; `compute_zmax` |
| `sdss_hmf.py` | Port of `sdsshmf.r`; generates `sdsshmf5.csv` |
| `nessie_sdss_hmf.py` | Nessie SDSS through the identical method; generates `sdsshmfNessie5.csv` |
| `scan_nessie_sdss.py` | Seed and penalty-range scans for the SDSS legs |
| `combined_hmf.py` | Port of `allhmf.r`; the GSR fit and Ω_M inset |
| `robust_hmf.py` | Bootstrap errors, Vmax floor scan, flagged-group table |

`driver_recovery.py` holds the pieces that make faithfulness possible:

* `RRandom`, `r_qnorm` — R's `set.seed` + Mersenne-Twister + inversion normal
  (Wichura AS241), bit-for-bit. numpy's normals give a different draw, and the
  fit is sensitive to it.
* `r_optim_nm` — verbatim port of `nmmin` from R's `src/appl/optim.c`, including
  the initial simplex construction and the `funcount <= maxit` do-while.
* `r_maghist` / `r_weighted_hist` — the two histograms are **not** the same
  function. `maghist` trims to the break range and is right-closed;
  `weighted.hist` is left-closed, extends the last break, and has no `na.rm`, so
  one NA weight turns a whole bin NA.
* `bootstrap_error`, `bin_index`, `bin_hmf`, `fit_mrp`, `fit_selection`,
  `co_dist` / `co_vol`.

**Changing `driver_recovery.py`'s defaults breaks the verified reproductions.**
New behaviour goes behind an opt-in argument (`nboot`, `mmax`, `vmax_floor_frac`
were all added this way).

---

## Validation numbers — treat these as regression tests

| check | expected |
|---|---|
| Binned HMF vs published Table 1 | raw counts identical in all 21 bins; log10 φ identical to 3 dp |
| `driver_recovery.py --seed 1` | logM\* = 13.5822 (R gives 13.5821771438), 501 fevals |
| `driver_recovery.py --seed 10` | 13.493 / −3.184 / −1.268 / 0.463 vs published 13.51 / −3.19 / −1.27 / 0.47 |
| MC error band, seed 10 | logM\* 12.04–13.78 vs Driver's published 13.51 +0.26/−1.51 |
| `combined_hmf.py` GSR, Driver's GAMA | 14.150 / −3.995 / −1.695 / 0.640 vs Table 2's 14.13 / −3.96 / −1.68 / 0.63 |
| `combined_hmf.py --myoption GR` | 13.741 vs published 13.72 |
| `sdss_hmf.py` | 13.213 vs published SDSS5 13.38 |

The GSR check is the strong one: it exercises the GAMA HMF, the SDSS HMF, the
REFLEX conversions, the combined objective, the penalty and the optimiser
simultaneously, and they all agree to ~0.02.

**R is installed on this machine** (`/usr/local/bin/Rscript`) with `celestial`,
`plotrix`, `magicaxis`, `Rfits` and `data.table`. If a number is ever in doubt,
run Driver's script directly and diff. That is how the port was verified.

> **Never quote a GAMA-only fitted parameter from a single seed.** The MC seed
> alone moves α by 0.1–0.2 and logM\* by more than a dex, so any single-seed
> comparison between catalogues is noise. Run several seeds (or
> `driver_recovery.py --seed-scan N`) and report the median and scatter. This
> already caused one wrong conclusion — a claimed α agreement of 0.011 that was
> really a 0.45 disagreement. Single-seed numbers elsewhere in this file are
> illustrative only; the entries in the table above are reproducible checks
> against R and the paper, which is different.

---

## Two pathological groups

The 1/Vmax estimator is dominated by a handful of nearby groups at the
multiplicity threshold. `zmax` is set by the group's **5th brightest** member —
by construction, since the group leaves the sample when that galaxy does — and at
low redshift the 5th member is an intrinsically faint dwarf that disappears
almost immediately. Result: Vmax ≈ 0.1% of the survey and a weight several
hundred times typical.

**205509** (g12, Nfof=5, z=0.0192, σ=482 km/s, logM 14.26, vmax/vlimit=2.4e−3):
carries **38% of the 14.2 bin** on its own. Its members run r = 14.05, 15.34,
16.65, 18.06, **18.36** — the last one, M_r = −15.34, sets zmax = 0.033. Also note
σ = 482 km/s from five redshifts; a real 10^14.3 halo at z = 0.02 would have far
more than five GAMA members, so the mass is probably inflated too.
*Fix:* bootstrap errors (`--nboot`), which raise its fractional error from 0.301
to 0.501 and drop its χ² leverage to **0.36×**. That is the right treatment
because 1/Vmax is unbiased — the honest response is to state the uncertainty
properly, not to delete the object.

But be precise about what it does and does not do. It biases **that bin**
(−3.46 measured, against −4.20 for Driver's 14.2 bin), and it does **not** affect
the fitted parameters: removing group 205509 entirely moves the median α from
−0.870 to −0.871.

**300223** (g15, Nfof=258, z=0.139, σ=782, logM 15.58): a *good* cluster, not a
pathology — but it is alone in the 15.6 bin, and `max(allx)` sets the range the
penalty term integrates over. It moves the GAMA-only answer by **0.7–1.9 dex**
while contributing almost nothing to χ². This one *is* a large effect on the fit.
*Fix:* `--fit-max 15.5`, matching the range Driver's catalogue reaches.

### The catalogues do not agree at the faint end

Over an identical 12.8–15.4 range, with bootstrap errors, **8 MC seeds**:

| α | Driver v10 | Nessie as-is | Nessie −205509 |
|---|---|---|---|
| median | **−1.321** | **−0.870** | **−0.871** |
| scatter | 0.106 | 0.217 | 0.165 |

Nessie's GAMA-only faint-end slope is shallower than Driver's by about **0.45**,
with roughly twice the seed-to-seed scatter. Neither outlier explains it.

An earlier version of this file claimed the two agreed to 0.011. That was seed 10
alone, and it does not survive a seed scan — see the warning under "Validation
numbers". The residual difference is real but **not quotable**: the GAMA-only fit
is the ill-posed one, so this is a hint that something differs between the
catalogues at the faint end, not a measurement of it.

`flagged_groups.csv` lists 42 new + 37 old offenders (criteria: Vmax below the
floor, or >10% of a bin holding ≥10 groups). **The same pathology is already in
Driver's catalogue** — one group is 89% of his 12.0 bin — so this is a statement
about the method, not about Nessie.

### The Vmax floor is available but deliberately not used

`vmax_floor_frac` is Driver's own `vlimitmin` (line 229), set at 1/1000 so it
never binds. Raising it does suppress the spikes, but it is a **biased**
estimator and the bias is mass-dependent: it bites hardest in the low-mass bins
where nearby low-Vmax groups are legitimate, flattening α from about −1.0 to
−0.4. That is exactly the quantity you want to report. Bootstrap errors achieve
the same protection with no bias. Keep the floor as a robustness check only.

---

## Files

### Driver's scripts (the specification) — all three now in `python_rewrite/`

| file | what it is |
|---|---|
| `gamahmf.r` | The binned GAMA 1/Vmax HMF and the MRP fit. Lines 224-320 and 355-420 hold nearly everything. |
| `allhmf.r` | The combined multi-survey figure and fit. Came from `~/Downloads`. |
| `sdsshmf.r` | The SDSS pipeline that produces `sdsshmf5.csv`. Came from `~/Downloads`. |

Read them end to end before changing the corresponding Python. The first attempt
at this read `gamahmf.r` in fragments and guessed the gaps; four guesses were
wrong.

### Data

| file | location | notes |
|---|---|---|
| `G3CFoFGroupv10.fits` | `../data/` | Old GAMA groups, 3 equatorial fields, 179.92 deg² |
| `GAMAGalsInGroups.csv` | `../data/` | Members with per-galaxy `zmax_19p8` |
| `reflex.csv` | `../data/` | REFLEX II, 43 rows, columns `x`, `Curve1` |
| `tpigg.dat` | `../data/` | 2PIGG, 4 columns. Plotted, never fitted. |
| `elmo.csv` | `../data/` | Tempel+14 SDSS DR10 **curve**, 200 rows × 5 cols. Plotted, never fitted — this is *not* the SDSS HMF that gets fitted. |
| `sdssdr10table1/2.fits` | `../data/` | Tempel+14 galaxies (584,449) and groups (88,662) |
| `G3CFoFGroup.fits`, `G3CGal.fits` | `/Users/00115372/Desktop/my_tools/make_gama_dmu/` | New GAMA DMU, 4 regions, 238.11 deg² |
| `sdss_groups.parquet`, `sdss_galaxies.parquet` | `/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS/` | Nessie SDSS |

### Generated (Driver's V1..V8 column order)

`gamahmfGAMA5.csv`, `gamahmfNessie5.csv`, `sdsshmf5.csv`, `flagged_groups.csv`.

V1 bin centre, V2 raw N, V3 log10(weighted counts), V4 log10(φ_corr),
V5 Poisson, V6 Monte-Carlo, V7 cosmic variance, V8 combined fractional error.
`allhmf.r` uses only V1, V4, V8.

### The paper is a data source

**Published Table 1 *is* Driver's binned GAMA table** (`gamahmfGAMA5.csv`), and
Table 2 is his fits for every sample combination. Extract with
`pdftotext -layout driver2022.pdf`. Both are hardcoded in the code as
`PUBLISHED_TABLE1` / `PUBLISHED_TABLE2` for the regression checks.

### Figures

`driver_fig4.pdf` (reproduction), `hmf_gama_only_nessie.pdf` (old vs new,
GAMA-only, with MC draws), `hmf_combined_nessie.pdf` (combined, with Driver's fit
faint behind), `hmf_robust*.pdf`.

---

## Discrepancies inside Driver's own scripts

Easy to trip on. Each is reproduced where it belongs rather than harmonised.

| quantity | `gamahmf.r` | `allhmf.r` | `sdsshmf.r` |
|---|---|---|---|
| GAMA area | 179.92 deg² | 175 deg² | — |
| `parscale` | `c(1,1,1,0.5)` | `c(1,1,1,0.1)` | `c(1,1,1,0.1)` |
| `logbin` | 0.2 | 0.2 | **0.1** |
| `fitbinwid` | 0.01 | 0.001 | 0.01 |
| `rootnerr` | `1/sqrt(N)` → Inf at N=0 | — | `sqrt(N)/N` → **NaN** at N=0 |
| `volumesdss` | — | no zmin subtraction | subtracts zmin |
| mlimit | 12.7 | 12.7 GAMA / 12.9 SDSS | 12.9 |

Other reproduced quirks:

* **The `weightszlimit` overwrite bug** (`gamahmf.r` 307-308, `sdsshmf.r` 304-305):
  the second assignment overwrites the first and its else-branch is `vmax`, so the
  upper clip is discarded. Harmless (line 303 already caps zmax) but reproduced.
  Worth mentioning to Driver.
* `allhmf.r` line 354 draws `rnorm(length(gama$V4))` for REFLEX, so 14 GAMA-length
  deviates get recycled across 43 REFLEX points.
* `allhmf.r` applies cosmic variance in **log** space (`cv = V4 * rnorm(...)`,
  where V4 is log10 φ) but the measurement error in linear space.
* `sdsshmf.r` line 300 uses `1.1 * zcl` where `gamahmf.r` line 302 uses `Zfof`.
* Hand-fixes: `gamahmf.r` line 310 sets `GroupID==100622` to 1e9; `sdsshmf.r`
  line 306 sets `idcl==81455` to `volumesdss`.

---

## Sentinels — check every file

**`GAMAGalsInGroups.csv`: `GroupID == 0` is ungrouped** (109,052 of 184,081).
Harmless for lookups of real groups, but any `groupby` builds a spurious
109,052-member group. Drop it first.

**New GAMA DMU: the ungrouped sentinel is `0`, not `-1`.** 30,039 galaxies in
`G3CGal.fits`; there are no `-1` and no `99999`. `make_gama_dmu` already applies
the offset and remap, and IDs run 100001–505225. An earlier version of this file
warned about `-1`; that check is a no-op.

**`GroupID == 100622`** is a known-bad object in the **old** catalogue only. The
same ID exists in the new one (g09 offset is 1e5) but is an entirely different,
good group — do not apply the hand-fix there.

**Shark mock magnitudes: `-999`.** `np.isfinite()` does not catch it. Filter on a
physical range: `~np.isfinite(m) | (m < -40) | (m > 0)`.

**Group masses:** 431 new groups have `VelDisp == 0` → `MassAfunc == 0`. Driver's
`MassAfunc > 1E1` cut removes them; keep that cut and keep it *before* the mass
rebuild.

**The new DMU emits no `zmax`.** `new_gama_hmf.compute_zmax` inverts the DMU's own
selection — `AbsoluteMagR + distmod(z) + ke(z) = 19.65` — using
`make_gama_dmu`'s k+e polynomial and its `FlatCosmology(1.0, 0.25)`. Verified by
recomputing `AbsoluteMagR` from `ApparentMagR`: agrees to **3e-6 mag**. Note the
old `zmax_19p8` column is at r < 19.8 and would give every new group too much
volume. Control run: using the polynomial on the *old* data shifts logM\* by only
−0.035, so the zmax method is not driving any old/new difference.

---

## Cosmology

`ho = 67.37`, `omegam = 0.3147` **everywhere** — volumes and masses both. The
`(100/ho)` factor is already inside `mymass`; do not apply any further h
conversion. Sanity check: `vlimit` = 2.082936e+07 Mpc³ for the old GAMA area
(≈3.27× the h=1 value), `cosvariance` = 0.066122.

A previous attempt reassigned a module-level `H0` at runtime, which changed the
masses but not the volumes because the distance function had captured the old
value. Masses in one cosmology and volumes in another is worse than either being
consistently wrong. The current code sets it once at module level and derives
everything from `co_dist` / `co_vol`; keep it that way.

Reference values for the new data: area 238.11 deg², `vlimit` = 2.756602e+07,
`cosvariance` = 0.060440, 1833 groups (g09 403, g12 501, g15 500, g23 429).
Old data: 1733 groups, **0 with NA zmax** — the `Nfof >= 5` cut guarantees every
group has at least 5 members in the file, so the NA-weight question that an
earlier version of this file raised does not arise.

---

## Method reference

Unchanged from Driver; recorded here so the code can be checked against it.

**Constants** (`gamahmf.r` 213-253): `zmin = 0.015` (not 0.01), `zlimit = 0.25`,
`multi = 5`, `mlimit = 12.7`, `magica = 13.9`, `myoption = "GAMA"`,
`logbin = 0.2`, `massx = seq(10.3, 16.1, 0.2)`.

**Masses** — `myoption="GAMA"` rebuilds from the velocity dispersion, it does not
read `MassAfunc`:

```
mymass    = magica * (VelDisp*1000)^2 * Rad50 * parsec*1e6 / (G*msol) * (100/ho)
MassAfunc = mymass / 10^masscorr[Nfof]
```

**Vmax** comes from the members, not from inverting a mass limit:

```
zmax = sort(members' zmax_19p8, decreasing=TRUE)[multi]   # [2] if Nfof==2
zmax = max(zmax, Zfof); zmax = min(zmax, zlimit)
vmax = V(zmax) - V(zmin)
```

Every group is kept; nothing is cut on mass.

**Eddington correction** is a Monte-Carlo *ratio of counts*, not a deconvolution:
smear the masses 1001 times by the per-group `log10MassErr`, take
`edb = mean(mock counts) / measured counts`, and divide. `gamay = counts/(logbin*edb)`.

**Errors:** `mcerr` from the 0.66 quantile of the MC scatter, `rootnerr` from the
**raw** counts, `gamaf = sqrt(mcerr² + rootnerr²)`, then NA→0.9999, Inf→0.0,
≥1→0.9999.

**The fit:** `phi` is **linear** in the parameter vector, not log10. Start from
the Murray+21 LCDM point (`mstarmrp = 14.42947`, `alphamrp = -1.864908`,
`betamrp = 0.7097976`, `phimrp = A/factor` with `A = 1.727006e-19`). Only bins
with `gamay > 0` and `gamax > mlimit` enter.

---

## The Nessie SDSS leg — built

`nessie_sdss_hmf.py`, with `scan_nessie_sdss.py` for the seed/range scans. Only
group construction differs from `sdss_hmf.py`; the binning, table and fit are
**imported** from it so the two legs cannot drift apart. 4824 groups against
Tempel's 4847.

The three questions the old version of this file raised were all answerable from
the data, and none needed a convention chosen:

1. **No k-correction had to be invented.** The Nessie galaxy file *is* Tempel's
   table 1 — all 584,447 rows match `sdssdr10table1.fits` on (RA, Dec) with
   **identical** `rmag` (max |Δ| = 0.0) — so `absmag_r` comes across by join and
   `sdsshmf.r` line 269's analytic `dmax` is used verbatim. Tempel's `col9`
   redshift is used inside `dmax`, not Nessie's `zobs`: the two differ by up to
   0.0012 (frame convention) and the formula must stay self-consistent with the
   absolute magnitude it is paired with.
2. **The two ID columns are two different group finders.** `GroupID` is Tempel's
   (max 88662, matching `sdssdr10table2.fits`); `group_id` is Nessie's (max
   84890, matching `sdss_groups.parquet`). Both use **`-1`** for ungrouped —
   297,202 and 300,068 galaxies respectively. Only `group_id` joins the parquet
   pair. They agree for just 49% of rows, so joining on the wrong one silently
   mixes the two catalogues.
3. **The mass range was much less alarming than it looked.** 10^16.1 is real but
   comes from pairs and triples at high z. Within the selection actually used
   (multiplicity ≥ 5, z < 0.08) `estimated_mass` runs to 10^15.64 against
   Tempel's 10^15.08, with 19 groups above 10^15 against Tempel's 4, and a
   median only 0.09 dex above his. **`mass_proxy` is not the analogue** — it
   sits ~0.9 dex low and is uncalibrated. Use `estimated_mass`.

Other facts pinned down while building:

* **Nessie's cosmology is H0 = 70, Ω_M = 0.30**, recovered from its own
  `co_dist` column to 4e−8. Masses therefore carry `(70/ho)`, mirroring
  Driver's `(67.8/ho)` for Tempel.
* **Tempel's `rank` is by absolute magnitude**, ascending — verified exactly
  (1.000) on his own groups, against 0.545 for apparent. `nessie_sdss_hmf`
  reproduces that ranking inside each Nessie group to pick the rank-5 member.
* `sdsshmf.r`'s hand-fix of `idcl == 81455` is **not** applied to Nessie. Same
  lesson as GAMA's `GroupID 100622`: the ID exists in the new catalogue but is a
  different object.
* Column decode for `sdssdr10table1.fits`: `col4` idcl, `col6` rank, `col9` z,
  `col13` RA, `col14` Dec, `col26` r, `col31` M_r.
* `sdsshmf.r` was missing from `python_rewrite/` despite this file claiming
  otherwise; restored from `~/Downloads`.

### The SDSS-only fit is ill-posed too

Same disease as GAMA-only, and if anything worse. Nessie SDSS reaches logM 15.65
against Tempel's 15.05, and `make_massfn` integrates its penalty over
`max(allx) + 1..10 bins`, so **where you cut sets the answer** (seed 10):

| `--fit-max` | bins | logM\* | α | β | χ² | conv |
|---|---|---|---|---|---|---|
| none (15.65) | 28 | 13.306 | −1.170 | 0.442 | 126.7 | 1 |
| 15.4 | 25 | **14.128** | −1.440 | 0.731 | 103.3 | 1 |
| 15.2 | 23 | **12.395** | −0.537 | 0.388 | 76.7 | 1 |
| 15.0 | 21 | 12.849 | −0.579 | 0.467 | 67.4 | 0 |
| 14.8 | 19 | 13.700 | −0.943 | 0.731 | 71.5 | 0 |

logM\* swings 1.7 dex and α swings 0.9 on the range choice alone — far larger
than the seed scatter. **Do not quote a Nessie SDSS-only fitted parameter.**
`--fit-max` is provided for exactly this check.

Over 8 MC seeds, full range, Poisson errors — **0/8 converged for both**, the
`maxit=500` story again:

| | Tempel | Nessie |
|---|---|---|
| logM\* | 13.793 ± 0.168 | 13.762 ± 0.428 |
| log φ\* | −3.570 ± 0.251 | −3.402 ± 0.512 |
| α | **−1.746 ± 0.078** | **−1.431 ± 0.181** |
| β | 0.557 ± 0.043 | 0.507 ± 0.101 |
| χ² | 59.4 | 133.8 |
| bins | 22 | 28 |

Nessie's α is shallower by **+0.315** (1.6× the pooled scatter), with ~2.3×
Tempel's seed-to-seed scatter. Same direction as the GAMA discrepancy (+0.45) —
two independent surveys, same sign.

**But matching the range does not reconcile them, it amplifies the gap** — which
is what identifies the problem. 8 seeds, both capped at logM ≤ 15.0, 21 bins
each:

| α | Tempel | Nessie | diff |
|---|---|---|---|
| Poisson | −1.657 ± 0.087 | −0.585 ± 0.079 | **+1.072** (9.1×) |
| bootstrap | −1.863 ± 0.084 | −0.606 ± 0.083 | **+1.258** (10.6×) |

So the measured old/new difference is itself a function of the range chosen:
+0.32 at full range, +1.07 matched at 15.0. Read the two scans together and the
diagnosis is clean:

* **At a fixed range, the seeds are tight** — α scatters by only ±0.08.
* **Across range choices, α moves by ~0.9.** The systematic is an order of
  magnitude larger than the statistical error, and it is not reduced by matching
  ranges, bootstrapping, or averaging seeds.

That is the signature of an ill-posed fit, not of a measured catalogue
difference. The GAMA-only lesson transfers intact: **the SDSS-only binned points
are fine; the four-parameter fit to them is not.** Quote the combined fit.

Note seed 10 alone gives α = −1.170, outside the entire seed 1–8 range. The
"never quote one seed" rule earned its place again — but note it is the *range*,
not the seed, that does the real damage here.

Driver's own SDSS (Tempel+14) remains the `S` leg: 4847 groups, `volumesdss` =
3.100035e+07 Mpc³, area 7221 deg² (Nessie shares the footprint exactly — same
galaxies), z < 0.08, r < 17.77, `logbin = 0.1`, masses `mass * 1e12 * (67.8/ho)`
with no A factor and no multiplicity debiasing.

### In the combined fit — and the finding that matters

`combined_hmf.py --sdss` now takes `nessie` alongside `auto` (Tempel+14). The
two are kept as **separate figures**, never overlaid, and the default `--out`
follows `--sdss` so they cannot overwrite each other:

| figure | GAMA | SDSS |
|---|---|---|
| `hmf_combined_nessie.pdf` | Nessie | Tempel+14 |
| `hmf_combined_nessie_sdss.pdf` | Nessie | Nessie |

Swapping only the SDSS leg (GSR, Nessie GAMA, seed 10, Driver's `maxit=500`)
moves logM\* by −0.40 and α by +0.29:

| | logM\* | log φ\* | α | β | χ² |
|---|---|---|---|---|---|
| SDSS = Tempel | 14.362 | −4.368 | −1.807 | 0.738 | 244.2 |
| SDSS = Nessie | 13.967 | −3.697 | −1.521 | 0.597 | 274.8 |

**But only Driver's Tempel leg actually anchors the fit.** Released from the
budget, the two behave completely differently:

| SDSS leg | maxit=500 | maxit=5000 | fevals | χ² released |
|---|---|---|---|---|
| Tempel | 14.362, conv=**0** | **14.362, identical** | 247 | 244.2 |
| Nessie (full, 15.65) | 13.967, conv=1 | **12.431**, α=−0.879 | 895 | 241.1 |
| Nessie capped 15.05 | 14.064, conv=1 | **9.844**, α=+0.278 | 2769 | 186.6 |
| Nessie capped 15.0 | 13.729, conv=1 | **12.647**, α=−0.822 | 749 | 186.0 |

Tempel's SDSS gives a **genuine interior minimum** — converged in 247 fevals,
bit-identical at maxit 500 and 5000. Every Nessie SDSS variant stops on the
budget at 500 and, released, **runs away to low M\* with a lower χ²** — the
exact signature "The central finding" describes for GAMA-only. Capping the range
does not rescue it; capped at 15.05 it goes further, to logM\* = 9.84.

So **the fitted curve in `hmf_combined_nessie_sdss.pdf` is the `maxit=500`
accident, not a measurement** — the Murray+21 start point decaying, same as
Driver's GAMA-only number. Both figures use identical settings, so they are
honestly comparable as *figures*; but only the Tempel one carries a quotable
fit.

Read together with the SDSS-only scans, the pattern is consistent: the Nessie
SDSS **binned points are usable and interesting** — above 14.5 they track REFLEX
noticeably better than Tempel's, which fall below it, and they extend to 15.65 —
but the Nessie SDSS catalogue does not constrain the MRP cutoff well enough to
pin M\* and β. Driver's Tempel leg does. **Quote the Tempel combined fit; show
the Nessie SDSS points.**

`--maxit` (added for this check) and `--sdss-max` are both wired through.
Anything reporting `convergence = 1` should be re-run at maxit 2000-5000 before
it is believed.

### GS — dropping REFLEX II entirely

`--myoption GS` fits GAMA + SDSS only. REFLEX and 2PIGG are still *plotted*, as
Driver plots things he does not fit; the legend reads "Best fit MRP function to
GS".

| figure | GAMA | SDSS fitted |
|---|---|---|
| `hmf_combined_nessie_GS_sdss.pdf` | Nessie | Nessie |
| `hmf_combined_nessie_GS.pdf` | Nessie | Tempel+14 |

At Driver's `maxit=500`, seed 10:

| | logM\* | log φ\* | α | β | χ² | pts |
|---|---|---|---|---|---|---|
| GS, SDSS = Nessie | 13.735 | −3.383 | −1.339 | 0.511 | 188.3 | 43 |
| GS, SDSS = Tempel | 14.076 | −3.950 | −1.863 | 0.558 | 152.2 | 37 |

The configuration itself is sound — run on Driver's GAMA + his SDSS it
reproduces published GS (14.35, −4.38, −1.96, 0.60) as 14.265, −4.225, −1.915,
0.576, worst parameter 0.155.

**But without REFLEX nothing anchors the fit, for either SDSS catalogue.** This
is the sharpest demonstration of "The central finding" in the file:

| config | maxit=500 | maxit=5000 | fevals | conv |
|---|---|---|---|---|
| GS, SDSS = Tempel | 14.076, conv=1 | **10.541**, α=−0.841, χ² 133 | 1497 | 0 |
| GS, SDSS = Nessie | 13.735, conv=1 | **7.044**, α=+0.739, χ² 150 | 5001 | **1** |
| GSR, SDSS = Tempel | 14.362, conv=**0** | **14.362, identical** | 247 | 0 |
| GSR, SDSS = Nessie | 13.967, conv=1 | 12.431, α=−0.879 | 895 | 0 |

Both GS fits stop on the budget and, released, collapse to nonsense with a lower
χ² — Tempel's to logM\* 10.5, Nessie's to 7.0, where it has still not converged
after 5001 evaluations. **`GSR` with Driver's Tempel SDSS is the only
configuration anywhere in this project with a genuine interior minimum.**

So the two GS figures are honest *figures* — identical pipeline, correct points,
correct errors — but their fitted curves are `maxit=500` artefacts and neither is
quotable. The x-ray anchor is not optional; it is the whole reason the combined
fit works.

### Clamping the noisy high-mass tail — `--fit-max`

The high-mass ends of both survey legs are single-group bins whose fractional
errors have hit the 0.9999 ceiling, so they add no constraint while still
setting `max(allx)` and therefore the penalty range:

| leg | bins above 15.0 | raw N in them |
|---|---|---|
| GAMA (Nessie) | 15.60, 15.40 | 1, 1 |
| SDSS (Nessie) | 15.65, 15.55, 15.45, 15.35, 15.25, 15.15, 15.05 | 1, 1, 2, 1, 3, 2, 10 |
| SDSS (Tempel) | 15.05 | 4 |

`--fit-max` clamps **both survey legs** at a given logM. REFLEX is deliberately
left uncapped — it is the anchor and its high-mass points are the constraint,
not noise. `--sdss-max` still clamps SDSS alone; when both are given the tighter
wins. Clamped runs get a `_clamp<value>` filename suffix.

A clamp of **15.0 is the principled choice**: it is the largest cut that leaves
every retained bin with N ≥ 4, in all three legs. Do not pick the clamp by which
value converges.

**Clamping lowers χ² but does not, in general, restore an interior minimum:**

| SDSS | clamp | maxit=500 | maxit=5000 | conv |
|---|---|---|---|---|
| Tempel | none | 14.076 | 10.541 | 0 |
| Tempel | 15.2 | 14.046 | 9.122 | 0 |
| Tempel | 15.0 | 13.570 | 11.418 | 0 |
| Tempel | 14.8 | 13.914 | 12.016 | 0 |
| Nessie | none | 13.735 | 7.044 | 1 (!) |
| Nessie | 15.2 | 12.712 | 10.343 | 0 |
| Nessie | 15.0 | 12.895 | 10.874 | 0 |
| Nessie | **14.8** | **13.565** | **13.565, identical** | **0, 305 fevals** |

Every Tempel GS variant still runs away no matter where it is clamped. Exactly
one configuration converges — Nessie SDSS clamped at 14.8, giving logM\* =
13.565, α = −0.861, β = 0.659 in 305 evaluations, bit-identical at 500 and 5000.

**Treat that single convergence with suspicion rather than relief.** It appears
at 14.8 but not at 15.0 or 15.2, and not at all for Tempel; an answer that
depends that sharply on where you truncate is the ill-posedness showing itself,
not a measurement emerging. Quoting it would be picking the clamp that produces
a converged number, which is the error this file has warned about twice already.

The useful conclusion is the honest one: the noisy tail *is* worth clamping — it
drops GS χ² from 188 to 108 for Nessie — but clamping does not make GS
well-posed. Only `GSR` with Tempel SDSS is.

**Still to do:** why the Nessie SDSS leg fails to anchor. The excess at
13.0-14.3 relative to Tempel and the extra 0.6 dex of high-mass reach are the
obvious suspects; the `--sdss-max` scan shows it is not the reach alone.

---

## Deliverables

1. **GAMA-only plot** — binned points with errors, the fitted MRP with its
   Monte-Carlo band, Driver's result for comparison. `hmf_gama_only_nessie.pdf`.
2. **Combined plot** — GAMA + SDSS + REFLEX II fitted, 2PIGG and Tempel shown but
   not fitted, following `allhmf.r`. Two of them, same pipeline, differing only
   in the SDSS leg: `hmf_combined_nessie.pdf` (SDSS = Tempel+14) and
   `hmf_combined_nessie_sdss.pdf` (SDSS = Nessie). Only the first carries a
   quotable fit — see "In the combined fit".
3. **GAMA + SDSS only, no REFLEX** (`--myoption GS`) —
   `hmf_combined_nessie_GS_sdss.pdf` (SDSS = Nessie) and
   `hmf_combined_nessie_GS.pdf` (SDSS = Tempel+14). Points are sound; neither
   fit is quotable, because dropping the x-ray anchor makes the fit ill-posed
   for both catalogues. See "GS — dropping REFLEX II entirely".

Keep them as separate figures.
