# Driver+22's halo mass function: reproduction, and the same method on new data

## Status

**Step 1 (reproduce `gamahmf.r`) is done and verified bit-exactly against R.**
The binned HMF matches published Table 1; the fit matches R's `optim` to 10
significant figures; the combined multi-survey fit matches published Table 2 to
better than 0.035 in every parameter.

**Step 2 (new GAMA DMU) is done for GAMA.** The Nessie catalogue has been run
through the identical method, both GAMA-only and in the combined fit.

**Step 3 (Nessie SDSS) is built** — `nessie_sdss_hmf.py`, 4824 groups. All
three blockers were answerable from the data. Like GAMA-only, the SDSS-only
*fit* is ill-posed; see "The Nessie SDSS leg".

**Step 4 (Nessie SDSS in the combined fit) is done** — separate figures for
`--sdss auto` and `--sdss nessie`, plus GS (no REFLEX) and clamped variants.

**Step 5 (the 14.2 bin) is solved** — a mass-to-light consistency cut,
`--ml-cut 1.0`, drops 5 of 1833 groups and moves the bin 0.53 dex while every
other bin moves ≤ 0.05. See "The 14.2 bin — solved".

### What this session established about the mass scale

This is the main result, and it is largely negative:

* The Nessie/Tempel SDSS mass offset **decomposes exactly** (per-group residual
  +0.0004 dex) into four terms: the Hernquist-vs-NFW profile choice, Nessie's
  `cbrt(3)`-for-`sqrt(3)` bug, the gapper-vs-rms dispersion estimator, and
  sigma_sky. On a like-for-like NFW-to-NFW basis the offset is **+0.153 dex**,
  of which the dispersion estimator is ~64%.
* **The estimator difference is not the velocity-error term.** Removing it
  changes sigma by 0.0001 dex. Tempel uses the plain rms (his eq. 3), Nessie the
  gapper. `compare_dispersion_estimators.py` isolates this for a referee.
* **Forcing both legs onto one estimator makes things worse, in both
  directions** — see "Forcing a common mass estimator: tried, rejected". Each
  catalogue's native calibration agrees with REFLEX better than any common
  recipe. The mass scale is a calibration choice, not a formula bug to fix.
* **No configuration anywhere in this project is a global minimum.** Multi-start
  shows even Driver's own GSR has a lower-chi^2 solution at low M\*. `conv = 0`
  means the simplex stopped, not that it found the best answer. See
  "Multi-start".

An earlier version of this header claimed that correcting the mass shift
"restores the combined fit's interior minimum". **That was wrong** — it came
from a uniform-shift approximation landing in a different basin, and does not
survive the mass-dependent calculation or a multi-start check.

**Outstanding:** the sigma_sky units question (h^-1 Mpc vs Mpc, ~0.17 dex),
which also bears on Driver's own SDSS normalisation; and Nessie's `sqrt(3)` bug,
which should go upstream regardless.

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
| `vuvuzela.py` | SDSS-specific version of Driver's fig. 3 mass-error calibration; also `--mass-audit` |
| `compare_dispersion_estimators.py` | gapper vs Tempel's rms at fixed everything else; referee figure |
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

**Every deliverable figure carries Driver+22's published best fit in the
background**, in tan (`PUBLISHED_C = "#c8781e"`, lw 3, alpha 0.55, lowest
zorder). It is his *printed* answer from table 2 — `PUBLISHED_TABLE2[myoption]`
for the combined figures, `dr.PUBLISHED_FIT` (GAMA5) for the GAMA-only one — and
is a **different line** from the existing dash-dot "same fit using Driver+22
GAMA", which is our refit of his catalogue. Both are labelled so they cannot be
confused.

`combined_hmf.plot_combined` lays its legend out by hand, so the row spacing is
computed from the number of entries (`n_rows`, `y0 = -4.95`, `y_floor = -7.62`).
Adding the published entry as a fixed-`dy` row pushed the last line off the axes
and into the tick labels. If you add another legend row, that arithmetic is what
keeps it on the canvas — do not hardcode `dy` again.

`driver_fig4.pdf` deliberately does **not** get the extra line: it is the
reproduction of Driver's own figure, and our fit already matches his to 10
significant figures, so a second curve would be indistinguishable noise.


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

Tempel's SDSS gives an interior minimum that is stable under the budget —
converged in 247 fevals, bit-identical at maxit 500 and 5000. **It is not the
global minimum**; see "Multi-start: none of these are global minima".

Every Nessie SDSS variant stops on the budget at 500 and, released, **runs away to low M\* with a lower χ²** — the
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
after 5001 evaluations. **`GSR` with Driver's Tempel SDSS is the only configuration that is stable
under the optimiser budget** — but multi-start shows even it is only a local
minimum (see below).

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

## Why the Nessie SDSS masses are high

`estimated_mass` is Nessie's `m_total` (`fof/src/group_properties.rs`):

```
sky_disp = sqrt( sum r_proj^2 / (2N (1+z)^2) )
grav_rad = 4.582 * sky_disp
M        = 2.325e12 * grav_rad * (3^(1/3) * sigma_gapper / 100)^2
```

The last line is **Tempel+2014 equation 8** — the same formula Tempel himself
uses. Only `(70/ho)` is applied on top, for the h difference. (`mass_proxy` is
the Robotham+2011 raw mass `R50*sigma^2/G` with no A factor, which is why it
sits ~0.9 dex low. It is not the analogue.)

Tested on **7258 groups whose membership is byte-identical in both
catalogues** — same galaxies, so the group finder cannot be the cause:

| term | dex |
|---|---|
| **observed** log10(M_Nessie / M_Tempel) | **+0.285** (scatter 0.432) |
| from sigma — Nessie's gapper runs 14% high, and M ∝ sigma² | +0.114 |
| from the R_g constant — Nessie **4.582** vs Tempel's effective **3.932** | +0.066 |
| residual — the sky-dispersion measure itself | +0.105 |

Not an h slip: log10(1/0.7) = +0.155, which matches nothing here, and the h
factor is already applied. Tempel's 3.932 was recovered by inverting his own
Eq 8 against his published sigma (`col12`) and sky dispersion (`col13`); the
ratio scatter is 0.045, so it is a genuine constant, not a fit.

**A +0.285 dex shift on a slope this steep is what produces both the elevated
points at 13.0-14.3 and the tail to 15.65.** So the Nessie SDSS HMF sits high
because of the mass estimator, not the grouping. This also reframes the earlier
note that Nessie tracks REFLEX better above 14.5 — that may be the offset, not
better recovery. Reproduce with `python vuvuzela.py --mass-audit`.

---

## The SDSS vuvuzela — `vuvuzela.py`

Driver applies the **GAMA-derived** mass-error curve to SDSS unchanged
(`sdsshmf.r` line 281 = `gamahmf.r` line 270, same hardcoded `xx`/`yy`). The
derivation code was never published — grepping the machine finds only copies of
the resulting arrays. Nessie's Rust has the right idea
(`calculate_error_function` documents the exact procedure) but
`create_mass_error_track` is still `todo!()`.

`vuvuzela.py` measures it directly for SDSS. Method, from the figure 3 caption:
groups with N > 20, remove members one at a time, recompute the mass at each
multiplicity, take the 16/50/84 quantiles. Members are removed **faintest-first
in apparent magnitude**, since that is what happens to a group receding through
a flux limit. 677 tracks.

Nessie's estimator is ported to Python and **validated against Nessie's own
output**: the gapper sigma is exact (max |Δ| = 0.00000) and, using Nessie's
stored centre, the mass agrees to 0.0006 dex. The only discrepancy is a
tie-break in `calculate_iterative_center_idx` affecting ~1% of groups.

Results (`sdss_masserr.csv`, `vuvuzela_sdss.pdf`):

| | ratio SDSS / GAMA |
|---|---|
| sigma vs Driver's `yy` | **0.995** |
| bias (q50) vs Driver's `masscorr` | 1.119 |

**Driver's reuse of the GAMA curve for SDSS is justified** — now empirically,
rather than by assumption. The measured SDSS scatter is within half a per cent
of his GAMA curve across N = 3-22, and the bias curve matches `masscorr` in sign
and shape. Robotham+2011 is far too wide (0.866 at N=3 against our 0.575).

**The one thing a single sigma cannot capture is the skew.** The tracks are
strongly asymmetric — median |q16|/|q84| = **1.72**, e.g. at N=3 the quantiles
are −0.852 / −0.171 / +0.298. Losing members drives mass *down* far more than up.
Driver's Eddington Monte-Carlo smears with a symmetric Gaussian, so it models
neither the skew nor the bias (he handles the bias separately, via `masscorr`).
Whether that matters to the HMF is untested — it would need the MC redone with
the empirical asymmetric kernel.

Nothing in the HMF pipeline uses these numbers yet; `sdss_masserr.csv` is
produced and left for a deliberate opt-in, per the rule that new behaviour does
not change verified defaults.

---

## The mass estimator is what breaks the anchor — confirmed

Shifting the Nessie SDSS masses onto Tempel's scale **restores the interior
minimum**. GSR, Nessie GAMA + REFLEX, seed 10:

| SDSS leg | maxit=500 | maxit=5000 | fevals | conv |
|---|---|---|---|---|
| Tempel (reference) | 14.362 | **14.362 identical** | 247 | 0 |
| Nessie as-is | 13.967 | 12.431 | 895 | runs away |
| Nessie shifted −0.102 | 14.116 | 11.002 | 2169 | runs away |
| **Nessie shifted −0.285** | 14.121 | **14.112** | 549 | **0** |
| Nessie `robotham` | 11.961 | 10.989, χ²=709 | 909 | catastrophic |

Only the **−0.285** shift works — the offset measured on identical-membership
groups, i.e. the true *estimator* offset. The −0.102 sample-median offset does
not. That distinction is the evidence: it is the estimator, not a normalisation
mismatch, and correcting it recovers logM\* = 14.11 against Tempel's 14.36 and
Driver's published 14.13.

χ² is 352.6 there against Tempel's 244.2, so it is a *converged* fit, not a
better one. Do not read −0.285 as a calibration to adopt; it is a diagnostic
that localises the problem to the mass scale.

### The two Nessie legs do not use the same estimator

They never did, and this was not deliberate:

| leg | estimator | debiasing |
|---|---|---|
| Nessie GAMA (`new_gama_hmf.py` line 104) | Robotham+2011, `13.9 * R50 * sigma^2 / G` | ÷ `10^masscorr` |
| Nessie SDSS (`nessie_sdss_hmf.py`) | Tempel+2014 eq. 8, `4.582 * sky_disp` | none |

Driver mixes them the same way — Robotham for GAMA, Tempel's own for SDSS — so
each leg is faithful to *its* parent catalogue. But the two Nessie legs are then
not on a common scale, which matters as soon as they are fitted together.

**Making them consistent makes it worse, not better.** `--mass-mode robotham`
puts the SDSS leg on the GAMA estimator (`mass_proxy * 13.9 / 10^masscorr`) and
gives masses **+0.391 dex** above Tempel, max logM 15.85, χ² 709, and a fit that
collapses to logM\* 10.99. So the inconsistency is not the thing to fix by
harmonising upward.

**This implicates the GAMA leg.** The Nessie GAMA leg uses exactly the estimator
that, applied to SDSS, lands 0.39 dex high. If the same inflation is present in
the GAMA masses it would push the binned points up and flatten the fitted slope
— which is the direction of the unexplained GAMA-only α discrepancy (+0.45,
shallower than Driver's). Untested, and the obvious next thing to test.

### The two Nessie catalogues are built in different cosmologies

* GAMA DMU: `FlatCosmology(1.0, 0.25)` — h = 1.0, Ω_M = 0.25
  (`make_gama_dmu/config.py` line 9)
* Nessie SDSS: h = 0.7, Ω_M = 0.30 (recovered from `co_dist` to 4e−8)

**The h scaling is handled correctly and every leg lands in the same
ho = 67.37 system:**

| leg | radius units | conversion | why |
|---|---|---|---|
| Nessie GAMA | `Rad50` in Mpc/h (DMU built at H0 = 100) | `(100/ho)` inside `mymass` | M ∝ R, and R[Mpc/h] → R[Mpc] is `×100/ho` |
| Nessie SDSS | `r50`, `sky_disp` in physical Mpc at h = 0.7 | `(70/ho)` | D ∝ 1/H0 at fixed z and angle |
| Tempel | Tempel's own, h = 0.678 | `(67.8/ho)` | Driver's, unchanged |

**The Ω_M difference is not corrected, but it is negligible** — an earlier
version of this file called it "worth straightening out", which overstated it.
Mass goes as radius, radius as distance, so the error is `d log10 D` at the
group redshift:

| leg | Ω_M used | at z | error |
|---|---|---|---|
| Nessie GAMA | 0.25 vs 0.3147 | 0.05 → 0.25 | +0.001 → **+0.005 dex** |
| Nessie SDSS | 0.30 vs 0.3147 | 0.02 → 0.08 | +0.0001 → **+0.0004 dex** |

At most 0.005 dex, against a 0.285 dex estimator offset — two orders of
magnitude smaller than the effect that actually matters. Not worth fixing.

**Volumes are all in one cosmology.** Every `vmax`, `vlimit` and `volumesdss`
comes from `co_dist`/`co_vol` at ho = 67.37, Ω_M = 0.3147, for every leg. The
new GAMA `zmax` is derived in the DMU's own cosmology, but that is correct
rather than inconsistent: it inverts the DMU's own selection, and a redshift at
which a galaxy drops out of a flux-limited sample is not cosmology-dependent
(verified to 3e-6 mag).

Also note `make_gama_dmu/config.py` sets `MASS_A = 10`, where Driver uses
A = 13.9. The DMU's own `MassA`/`MassAfunc` columns therefore differ from what
this pipeline computes — which is harmless only because `new_gama_hmf` rebuilds
the mass from `Rad50` and `VelDisp` rather than reading those columns.

---

## Why the Nessie masses run high — the full picture

Two independent effects, and they compound. Measured on the **7258 SDSS groups
whose membership is byte-identical** in both catalogues, so the group finder is
excluded by construction. The decomposition is exact per group (residual median
−0.00000):

| N | log10(σ_N/σ_T) | log10(M_N/M_T) |
|---|---|---|
| 5-6 | +0.073 | **+0.327** |
| 7-9 | +0.050 | +0.260 |
| 10-14 | +0.035 | +0.200 |
| 15-24 | +0.020 | +0.148 |
| 25+ | **+0.007** | **+0.101** |

**1. A constant ~+0.10 dex floor.** It survives at N >= 25 where the dispersions
agree to 1.6%, so it is not a sigma effect. `tempel2014.pdf` settles what it is,
and it is **two** errors that partly cancel — see "Hernquist vs NFW" below.



**2. A small-N bias in the gapper σ**, adding up to +0.22 dex more at N = 5-6.
M ∝ σ², so a 17% σ overestimate is 0.15 dex of mass. Nessie's
`velocity_dispersion_gapper` subtracts a *constant* velocity-error term
(`sigma_err_squared`, which is why `velocity_dispersion_gap_err` is 7.071068 =
sqrt(50) for every group) and applies no small-N correction.

Because the SDSS sample is dominated by N = 5-6 groups, the population median
lands at +0.27 to +0.285.

### The GAMA leg looks the same, though the test is weaker

Nessie's GAMA dispersions also run above Driver's v10, at every multiplicity:

| Nfof | v10 median σ | Nessie median σ | Δ log10 σ | Δ log10 M |
|---|---|---|---|---|
| 5 | 201.1 | 243.0 | +0.082 | **+0.176** |
| 6 | 225.3 | 252.2 | +0.049 | +0.074 |
| 7-9 | 252.9 | 286.6 | +0.054 | +0.070 |
| 10-14 | 284.4 | 334.8 | +0.071 | +0.083 |
| 15+ | 446.9 | 551.3 | +0.091 | +0.061 |

**Caveat: the group finders differ, so memberships differ — this is not the
clean estimator test the SDSS comparison is.** Different membership legitimately
changes σ. But the direction matches, and the largest offset is at Nfof = 5,
which is exactly where the pathological groups live. Treat as a strong hint, not
a measurement.

Note the GAMA leg does **not** use `estimated_mass`; it rebuilds Robotham+A from
`VelDisp` and `Rad50`, the same formula Driver uses. So effect 1 above (the
4.582 constant) does **not** apply to GAMA — only the dispersion difference
does.

### The Vmax floor is now reachable from the pipeline

`--vmax-floor` on `combined_hmf.py` (default 1e-3 = Driver's `vlimitmin`, which
never binds), threaded through the GAMA and both SDSS legs. Raising it caps
1/Vmax. Effect on the GAMA 14.2 bin, seed 10:

| floor | groups capped | 14.2 log φ | 205509's share |
|---|---|---|---|
| 0.001 (Driver's) | 3 | −3.464 | **37.7%** |
| 0.005 | 14 | −3.624 | 22.2% |
| 0.010 | 23 | −3.752 | 13.3% |
| 0.020 | 39 | −3.863 | 7.9% |
| 0.050 | 95 | −3.945 | 3.7% |

**The floor does not close the gap with Driver.** His 14.2 bin is −4.20. Even at
a floor of 0.05, which crushes 205509 to a 3.7% share and clips 95 groups, the
bin only reaches −3.945 — still 0.25 dex high. So the 14.2 excess is **not one
group**; it is consistent with the broader mass offset above, which raises the
whole distribution and moves groups up into the bin.

That reframes the fix: the floor is cosmetic for this bin. The mass scale is the
real lever.

Do **not** read fitted parameters off a floor scan — the GAMA-only fit is
ill-posed and single-seed, and in the scan above α swings −0.56 → +0.33
non-monotonically, which is noise, not a trend.

---

## Hernquist vs NFW — the constant offset, settled from `tempel2014.pdf`

Two independent errors, running opposite ways, which is why the naive "same
formula" comparison looked nearly right while being wrong twice over.

### 1. Nessie's 4.582 is correct — for Hernquist

Section 4.2, verbatim: *"Equating this with the Hernquist Re above, we obtain
Re = 1.386 sigma_sky = 1.8153a, hence a = 0.764 sigma_sky and **Rg = 6a =
4.582 sigma_sky**."* Nessie implements the coefficient correctly.

### 2. But Nessie has a real bug in the sigma factor

Eq. 8 is `Mtot = 2.325e12 (Rg/Mpc) (sigma_v / 100 km/s)^2`, and section 4 says
*"Assuming dynamical symmetry, the real (3D) velocity dispersion in groups would
thus be **sigma_v = sqrt(3) sigma_v1D**."*

`fof/src/group_properties.rs` codes:

```rust
2.325e12 * gravitational_radius
        * ((3_f64.powf(1. / 3.)) * los_velocity_dispersion / 100.).powi(2)
```

`3f64.powf(1./3.)` is the **cube** root, 1.4422. It must be the **square** root,
1.7321. Since M goes as the square of that factor, Nessie's masses are low by
`3^(2/3)/3 = 3^(-1/3)`, i.e. **−0.159 dex**. This is a genuine bug and worth
reporting upstream.

### 3. The published column Driver uses is NFW, not Hernquist

The paper releases both — appendix column list `[11] mass_nfw` and
`[12] mass_her`. Three independent lines of evidence say `col15` of
`sdssdr10table2.fits` is `mass_nfw`:

* Inverting eq. 8 on his own `col12` (sigma_v) and `col13` (sigma_sky) with the
  correct sqrt(3) gives **kappa = Rg/sigma_sky = 2.755**, nowhere near the
  Hernquist 4.582.
* Hernquist kappa is a **constant**; NFW kappa is solved iteratively per group
  and so is **mass-dependent**. The implied kappa trends 2.802 → 2.671 from
  logM 12.2 to 15.0, i.e. not constant.
* `4.582 / kappa` runs **1.636 → 1.715 rising with mass**, matching Fig. 7's
  *"the Hernquist masses are 1.55−1.75 more massive, depending on the system
  mass"* including the direction of the trend.

Scanning all 25 FITS columns, `col15` is the only one whose implied kappa is
even approximately constant (relative scatter 0.016). **The extract we have does
not contain `mass_her` at all** — worth requesting, since it would make the
comparison direct.

### The reconciliation

On 2500 identical-membership groups with N >= 5:

| term | dex |
|---|---|
| observed log10(M_nessie / M_tempel col15) | **+0.277** |
| correcting cbrt(3) → sqrt(3) *raises* Nessie by | +0.159 |
| converting Hernquist → NFW *lowers* Nessie by | −0.221 |
| net of the two | **−0.062** |
| so a fully corrected Nessie sits at | +0.215 |
| measured sigma + sigma_sky differences account for | +0.158 |

The remainder is the small-N gapper bias and the mass dependence of kappa.

**Fixing only the sqrt(3) bug would make the disagreement worse**, +0.277 →
+0.436. The two must be handled together.

`--mass-mode tempel_nfw` applies both corrections, putting the Nessie SDSS leg
on the same footing as Driver's `col15`. It is opt-in; `tempel_eq8` remains the
default so nothing already validated moves.

### Computing a proper NFW mass for the Nessie groups

`--mass-mode tempel_nfw` now does this properly, mass-dependent rather than by a
single ratio:

```
M_hern = estimated_mass * 3**(1/3)          # sqrt(3) where Nessie codes cbrt(3)
M_nfw  = M_hern * kappa_NFW(M_nfw) / 4.582  # iterated; kappa is weakly mass dep.
```

**The first-principles derivation does not reproduce Tempel's catalogue**, and
that is worth knowing before anyone tries again. Following section 4.1 — NFW
truncated at R200, `c200` from Maccio+2008 eq. 12, `Rg = G M^2/|U|` from eq. 9,
and sigma_sky^2 = <R^2>/2 = <r^2>/3 for a projected sphere — gives kappa
*rising* with mass (3.02 at logM 12.2 to 3.27 at 15.0) where his falls (2.80 to
2.67). Note `Rs` cancels, so kappa is a function of `c200` alone; the mismatch
is therefore in the sigma_sky → Rg step, which the paper delegates to
Bartelmann (1996) and Lokas & Mamon (2001) rather than writing out.

So kappa is instead **calibrated against his published catalogue**, by inverting
eq. 8 on his own `col12` and `col13` over 37 365 groups. That is exact by
construction: feeding his sigma and sigma_sky back through the calibrated
relation recovers his `col15` to **median −0.00003 dex, 16/84 ±0.0007, with 0%
of groups off by more than 0.05**.

| logM | kappa | | logM | kappa |
|---|---|---|---|---|
| 11.1 | 2.8672 | | 13.5 | 2.7343 |
| 12.1 | 2.8119 | | 14.1 | 2.7009 |
| 12.9 | 2.7671 | | 14.9 | 2.6594 |

kappa changes by only 0.0086 per dex, which is why the h convention used to
evaluate it does not matter.

Effect on identical-membership groups (N >= 5, 3000 groups):

| | median log10(M_nessie/M_tempel) |
|---|---|
| Nessie as-is (buggy Hernquist) | +0.296 |
| Nessie → NFW, both fixes | **+0.228** |

The remaining +0.23 is the sigma and sigma_sky estimator difference — the
small-N gapper bias — which is a separate problem and still open.

### The NFW-mass combined figures

`combined_hmf.py --sdss nessie --mass-mode tempel_nfw`. Affects the Nessie SDSS
leg only; the Tempel runs are untouched.

| figure | option |
|---|---|
| `hmf_combined_nessie_sdss_nfw.pdf` | GSR |
| `hmf_combined_nessie_GS_sdss_nfw.pdf` | GS |
| `hmf_combined_nessie_GS_sdss_nfw_clamp15.pdf` | GS, clamped 15.0 |

**The NFW masses genuinely improve the data-model agreement.** GSR chi^2 at
Driver's `maxit=500`:

| SDSS leg | logM\* | α | χ² |
|---|---|---|---|
| Nessie as-is (buggy Hernquist) | 13.967 | −1.521 | 274.8 |
| **Nessie → NFW** | **13.900** | **−1.461** | **234.6** |
| Tempel col15 | 14.362 | −1.807 | 244.2 |

So on the corrected mass scale the Nessie SDSS leg now fits *better than
Tempel's own*, and its binned points track REFLEX visibly more closely from 13.5
to 15.5. That is the real result here.

**It does not fix the ill-posedness, and an earlier claim in this file that it
did was wrong.** That claim came from a crude constant-ratio version of the
conversion; with the correct mass-dependent kappa the fit behaves like all the
others:

| start | maxit | logM\* | α | χ² | conv |
|---|---|---|---|---|---|
| Murray+21 | 500 | 13.900 | −1.461 | 234.6 | **1** |
| Murray+21 | 5000 | **12.345** | −0.852 | **218.2** | 0 |
| 13.0 | 5000 | 14.378 | −1.917 | 294.9 | 0 |
| 12.0 | 5000 | **10.945** | −0.324 | **214.8** | 0 |

Released from the budget it still runs away, and the lowest chi^2 is still at
low M\*. The 14.287 / conv=0 result reported earlier does not survive the
correct calculation — it was an artefact of the uniform-shift approximation
landing in a different basin.

**Bottom line:** use `tempel_nfw` because it is the right mass scale for
comparison with Driver's leg and it demonstrably fits better. Do not use it as
evidence that the combined fit is well-posed. All three NFW figures report
`convergence = 1`, so their fitted curves are `maxit=500` artefacts like every
other GAMA-or-Nessie configuration in this project.

### The sigma difference is gapper vs rms, not the velocity error

Tempel does **not** use the gapper. His eq. 3 is the plain rms:

```
sigma_v^2 = 1/[(1+z_m)^2 (n-1)] * sum (v_i - v_mean)^2
```

Recomputed from the members on 2000 identical-membership groups (N >= 5):

| estimator | median log10(est / col12) |
|---|---|
| gapper minus 50, i.e. Nessie's | +0.0533 |
| gapper with **no** error subtraction | +0.0534 |
| **Tempel eq. 3 rms, recomputed** | **+0.0003** |

Two conclusions:

* **The velocity-error term is irrelevant.** Removing it entirely moves sigma by
  0.0001 dex — 50 (km/s)^2 against a typical sigma^2 of 40 000 is 0.06%. The
  constant `velocity_dispersion_gap_err = 7.071068` is a red herring.
* **Eq. 3 reproduces his column exactly** (+0.0003 at every multiplicity), so
  the identification of `col12` is certain and the whole difference is
  gapper vs rms. The gap is N-dependent (+0.071 at N = 5-6, +0.018 at N >= 15),
  which is the two estimators diverging where the statistics are poor.

That accounts for ~+0.107 dex of the mass offset, since M goes as sigma^2.

`compare_dispersion_estimators.py` isolates exactly this for a referee:
`dispersion_estimator_comparison.pdf` puts NFW mass from the gapper against NFW
mass from eq. 3, on the **same groups, same members, same sigma_sky, same
kappa(M)** -- the only thing changing is sigma_v. 4000 identical-membership
groups, coloured by multiplicity:

| N | median Δlog10 M |
|---|---|
| 5-6 | **+0.141** |
| 7-9 | +0.095 |
| 10-14 | +0.065 |
| 15-24 | +0.038 |
| 25+ | +0.026 |
| **all** | **+0.108** |

The N-dependence is the whole story, and it is monotonic: the two estimators
agree for rich groups and diverge where the statistics are poor. Since the SDSS
sample is dominated by N = 5-6, the population offset lands at +0.108.

### An unresolved units problem in sigma_sky

`col13` is in **h^-1 Mpc** (eq. 4, comoving; the paper uses H0 = 100h,
Omega_m = 0.27). Nessie's `sky_disp` is in **physical Mpc** at h = 0.7. That is
a 0.169 dex unit difference that nothing in our chain accounts for.

Recomputing eq. 4 in his units with Nessie's centre gives **−0.111 dex** against
`col13`, so h alone does not explain it either — there is a residual difference
in the measure itself, most likely the group-centre definition, since he uses
his own FoF centre. The two effects net to the +0.054 measured between Nessie's
`sky_disp` and `col13`. **This one is not fully isolated.**

**This touches Driver's leg, not just ours.** If `col15` is in 1e12 h^-1 M_sun
then the physical mass is `col15 * 1e12 / h`, but `sdsshmf.r` line 278 does
`mass * 1e12 * (67.8/ho)` — treating it as already physical and applying only a
1.006 rescale. That is a potential **0.17 dex on the absolute mass scale of the
entire SDSS leg**, his published one included. His combined fit does reproduce
his own table 2, so he is internally consistent; but the SDSS normalisation
rests on that convention being the right reading. Worth settling before any
absolute mass scale is quoted.

### The offset that actually matters: NFW vs NFW

Everything above quoted the raw +0.28 dex between Nessie's `estimated_mass` and
Tempel's `col15`. **That number is not meaningful** — it compares a Hernquist
mass with an NFW mass. Once both legs are NFW (`--mass-mode tempel_nfw`), on
2500 identical-membership groups with N >= 5:

| comparison | dex |
|---|---|
| **Nessie NFW vs Tempel col15** | **+0.153** |
| of which, gapper → rms | +0.098 |
| of which, Nessie sigma_sky → his sigma_sky | +0.047 |
| both swapped, i.e. the closure check | **+0.0008** |

The closure line is the important one: swap both ingredients and Nessie's
pipeline reproduces his published `col15` to 0.0008 dex. **Nothing else is
unaccounted for.** The dispersion estimator is ~64% of the residual and
sigma_sky the rest.

Strongly multiplicity dependent, as expected:

| N | dex |
|---|---|
| 5-6 | **+0.192** |
| 7-9 | +0.133 |
| 10-14 | +0.092 |
| 15+ | +0.043 |

**Is +0.153 dex a lot?** Per group, no — it is half the random error at N = 5
(sigma_log10M = 0.326). As a systematic on 4824 groups, yes: the statistical
error on the mean is 0.008 dex, so it is ~19 sigma. And where the HMF is steep
it dominates:

| logM | dlogphi/dlogM | phi shift from +0.153 dex |
|---|---|---|
| 13.0-13.5 | −0.24 | 0.04 dex |
| 13.5-14.0 | −0.36 | 0.06 dex |
| **14.0-14.5** | **−1.81** | **0.28 dex (x1.9)** |
| **14.5-15.0** | **−1.89** | **0.29 dex (x1.9)** |

Below the knee it is negligible; above logM 14 it is a factor of ~1.9 in number
density, in exactly the range REFLEX anchors the fit. That is why the mass scale
kept surfacing as an anchoring problem.

---

## Multi-start: none of these are global minima

`conv = 0` means the simplex stopped, not that it found the best answer. Every
convergence claim in this file was made from Driver's single Murray+21 start
point. Restarting elsewhere (GSR, seed 10, maxit 5000):

| SDSS leg | start M\* | logM\* | α | χ² |
|---|---|---|---|---|
| Tempel col15 | Murray+21 | 14.362 | −1.807 | 244.17 |
| Tempel col15 | 13.0 | 14.385 | −1.922 | 257.65 |
| **Tempel col15** | **12.0** | **12.535** | −1.077 | **243.61** |
| Nessie → NFW | Murray+21 | 14.287 | −1.683 | 237.71 |
| Nessie → NFW | 13.0 | 14.373 | −1.914 | 283.01 |
| **Nessie → NFW** | **12.0** | **11.023** | −0.325 | **201.29** |

**Even Driver's own GSR configuration has a lower-χ² solution at low M\***
(243.61 against 244.17). The margin is only 0.56, so the two basins are close to
degenerate — which is why it looked stable. For the corrected Nessie leg the
low-M\* basin is preferred by 36 in χ², decisively.

The same knife-edge shows in a uniform mass shift. At maxit 5000 every shift
"converges", but bimodally and non-monotonically:

| shift | logM\* | χ² | | shift | logM\* | χ² |
|---|---|---|---|---|---|---|
| 0.000 | 12.431 | 241.1 | | −0.080 | 11.009 | 230.8 |
| −0.020 | 9.710 | 214.4 | | −0.100 | 11.350 | 223.9 |
| −0.040 | 12.073 | 230.6 | | −0.150 | **14.290** | 269.4 |
| −0.062 | **14.287** | 237.0 | | −0.200 | **14.296** | 275.9 |

Landing at ~14.29 versus ~11 depends on where you start and how far the mass
scale moved, not on which fits better — the low-M\* answers consistently have
the *lower* χ².

**So the mass corrections do not fix the ill-posedness.** What
`--mass-mode tempel_nfw` does is make the high-M\* basin reachable and stable
from Driver's standard start point, giving logM\* = 14.287, α = −1.683 — which
sits essentially on his published 14.13, −1.68. That is a meaningful result
about *consistency with Driver*, but it is not evidence that 14.29 is the best
fit to the data.

This is the same story as "The central finding", now shown to apply to the
combined fit as well as the GAMA-only one: `maxit = 500` from Murray+21 is an
accidental regulariser that keeps the answer physical. Driver's published
numbers rest on it, and so does any number we quote by the same route.

Below start M\* ≈ 11 the penalty term overflows and χ² returns ~1e35 after a
handful of evaluations. Numerical, not physical, but it means naive multi-start
over a wide range produces garbage rows rather than errors.

---

## The 14.2 bin — solved

The bin sat at log phi = −3.464 against Driver's −4.204, with **one group
carrying 38% of it**. The cause is now identified and the fix is a one-line
cut.

### What group 205509 actually is

N_fof = 5 at z = 0.019, sigma = 482 km/s, logM = 14.26, vmax/vlimit = 2.4e−3.
Its five members:

| M_r | dv (km/s) |
|---|---|
| −19.54 | −559 |
| −15.34 | −300 |
| −17.17 | 0 |
| −15.88 | +327 |
| −18.61 | +344 |

One bright galaxy and four dwarfs spread over 900 km/s. **At z = 0.019 GAMA
reaches M_r ~ −15**, so a genuine 10^14.3 halo here would have *hundreds* of
members, not five. The dispersion is being set by two galaxies at the extremes
of a loose association. 205469 (N = 8, dv −627 to +557, all dwarfs) is the same
story.

So the mass is not real — and because the group sits at z = 0.019 its Vmax is
0.2% of the survey, giving it a colossal 1/Vmax weight. A wrong mass *and* an
extreme weight.

### The diagnostic: mass-to-light

M/L is the one quantity that says whether a *mass* is credible independently of
its Vmax. Real systems sit near log10(M/L_r) ~ 2.5-3.5 and the ratio rises
smoothly with mass, so the residual against a running median is a clean, nearly
mass-independent outlier statistic.

| group | share of the 14.2 bin | log10(M/L) | median at that mass |
|---|---|---|---|
| 205509 | **37.7%** | **4.38** | 3.22 |
| 205469 | **13.5%** | **4.70** | 3.22 |
| 301298 | 5.9% | 4.17 | 3.22 |

Those three are 57% of the bin, at 10-30x the M/L of comparable systems.

### The cut

`ml_excess()` and `apply_ml_cut()` in `new_gama_hmf.py`, exposed as `--ml-cut`
on `new_gama_hmf.py` and `combined_hmf.py`. **1.0 dex is the recommended
value** — well beyond the 97.5th percentile of the excess distribution (+0.60),
so it removes only extreme outliers.

| cut | groups dropped | 14.2 phi | 14.0 phi | 13.0 phi |
|---|---|---|---|---|
| none | 0 | −3.464 | −4.118 | −3.549 |
| **1.0** | **5 of 1833** | **−3.993** | −4.070 | −3.510 |
| 0.8 | 11 | −4.081 | −4.075 | −3.498 |
| 0.6 | 47 | −4.104 | −4.098 | −3.493 |

Surgical: 5 groups (0.27%) move the 14.2 bin by 0.53 dex while every other bin
moves by <= 0.05.

**It is not forcing agreement.** The same cut on Driver's v10 removes 1 group
and changes his 14.2 bin by 0.001 dex — his catalogue does not have this
pathology, so the cut is removing something specific to Nessie. Apply it to
both regardless, so the comparison stays fair.

### Why this rather than a Vmax floor

The floor caps the *weight* of a group whose mass is still believed, and it is a
biased estimator whose bias is mass-dependent (see "The Vmax floor"). It also
does not work here: even at a floor of 0.05, which clips 95 groups, the bin only
reaches −3.945. The M/L cut removes objects whose mass is not credible in the
first place, which is the actual problem.

### Effect on the fit

GSR with Driver's Tempel SDSS, seed 10:

| | logM\* | α | χ² | fevals | conv |
|---|---|---|---|---|---|
| no cut | 14.362 | −1.807 | 244.17 | 247 | 0 |
| **M/L cut 1.0** | **14.361** | **−1.800** | **227.41** | 235 | 0 |

**The answer does not move (0.001 in logM\*, 0.007 in α) but χ² falls by 17.**
That is the signature of a good outlier cut: it removes points that were
genuinely discrepant rather than reshaping the fit.

It does **not** fix the multi-modality. Starting from M\* = 12 still finds
χ² = 225.86 against 227.41, so the low-M\* basin is still marginally preferred —
by 1.55 now, against 0.56 before the cut. See "Multi-start".

### Figures with the cut applied

| figure | contents |
|---|---|
| `hmf_old_vs_new_mlcut1.pdf` | GAMA-only, old vs new |
| `hmf_combined_nessie_mlcut1.pdf` | GSR, SDSS = Tempel |
| `hmf_combined_nessie_sdss_mlcut1_nfw.pdf` | GSR, SDSS = Nessie on the NFW scale |
| `hmf_combined_nessie_GS_sdss_mlcut1_nfw.pdf` | GS, SDSS = Nessie on the NFW scale |

`flagged_ml_outliers.csv` lists the dropped groups.

### What "drop" means, precisely

The group is removed from the sample before binning: it contributes no counts
and no 1/Vmax weight, and it is excluded from the Eddington Monte-Carlo. **The
survey volume is unchanged**, so phi genuinely falls in the bin it occupied.

It is **not** relocated to a lower mass. Re-added at the mass implied by the
median M/L for its luminosity, 205509 and 205469 would land near logM 12.70 — a
1.56 dex drop — and the affected low-mass bins would move by under 0.01 dex,
since they hold hundreds of groups. Not doing so is deliberate: assigning a mass
*from* luminosity and then measuring the halo mass function would be circular.
Dropping keeps the estimator purely dynamical and claims only that the mass is
not believed.

### A mass shift is NOT the residual explanation

After the cut, Nessie's 14.0 bin sits 0.215 dex *below* Driver's and its 14.2
bin 0.212 dex *above*, which looks like the function displaced along the mass
axis. **Tested and refuted.** Scanning a uniform shift applied to the Nessie
masses, agreement is *best with no shift*:

| shift | rms diff, 12.8-15.0 |
|---|---|
| **0.00** | **0.150** |
| −0.05 | 0.239 |
| −0.10 | 0.319 |
| −0.20 | 0.356 |

The two binned HMFs agree to 0.150 dex rms with medians of −0.058 (13-14) and
+0.082 (14-15) — small and of *opposite* sign, i.e. scatter, not a systematic
offset. The earlier note that Nessie GAMA masses run +0.176 dex high at
N_fof = 5 was a population median across two different group finders; it does
not show up as a coherent shift in the HMF and should not be treated as one.

---

## Forcing a common mass estimator: tried, rejected

Both directions were built, and both are worse than leaving each catalogue on
its own calibration. **The code and figures have been removed**; recover them
from commits `2eb8b8c` and `f170229` if ever needed. Recorded here so nobody
spends the time again.

| configuration | GSR chi^2 |
|---|---|
| GAMA Robotham + SDSS Nessie NFW *(the mix Driver effectively uses)* | **234.6** |
| GAMA Robotham + SDSS Tempel col15 | 244.2 |
| GAMA NFW + SDSS Nessie NFW | 327.1 |
| GAMA NFW + SDSS Nessie NFW + M/L cut | 409.6 |
| GAMA Robotham + SDSS Nessie Robotham | 719.7 |

* **SDSS put on GAMA's Robotham+A**: masses land **+0.39 dex** above Tempel's,
  lifting the SDSS points clear of the REFLEX sequence. Without REFLEX (GS) the
  chi^2 is 255.0 against 719.7 with it, so it is specifically the x-ray anchor
  that disagrees.
* **GAMA put on Tempel's NFW**: masses drop **-0.41 dex** (-0.46 at N_fof = 5,
  -0.35 above 30) and the GAMA points fall *below* REFLEX, SDSS and 2PIGG above
  logM 14.

The two estimators differ by ~0.4 dex in a consistent sense, and each
catalogue's **native** calibration agrees with REFLEX better than either forced
common recipe. Driver's A = 13.9 was tuned so GAMA masses match external mass
scales; Tempel's NFW was tuned for his own survey. Imposing either on the other
breaks that agreement, in both directions.

**So the mass scale is not something a change of formula fixes.** The offsets
are calibration choices, each internally coherent with the survey it was built
for. Keep each leg on its own calibration, state which was used, and treat the
~0.4 dex between them as a systematic to be quoted rather than removed.

One caveat that limited the GAMA-NFW test and is still unresolved: Tempel's
kappa was calibrated on *his* sigma_sky, which is in h^-1 Mpc comoving, while
the GAMA sigma_sky was computed in Mpc. M is linear in sigma_sky, so ~0.17 dex
of that -0.41 is uncertain. See "An unresolved units problem in sigma_sky".

## The reference line is fixed, and does not follow --myoption

`DRIVER_ABSTRACT_FIT = (14.13, -3.96, -1.68, 0.63)` in both `combined_hmf.py`
and `new_gama_hmf.py` — Driver's headline GSR result as quoted in his abstract
and table 2. **Every figure draws this same curve**, whatever `--myoption` is
set to, labelled "Driver+22 published GSR fit".

It previously used `PUBLISHED_TABLE2[myoption]`, so the reference moved from
plot to plot (GS was drawing 14.35, −4.38, −1.96, 0.60). That made figures
non-comparable at a glance.

Note the abstract states log10(phi\*) = **−3.96**, not −3.95.

The *validation printout* still uses `PUBLISHED_TABLE2[myoption]`, which is a
different question: there we check the port against his fit to the same sample
combination. Do not unify those two.

---

## `recovery.py` — audit of the hierarchical Bayesian route

Command audited:

```
uv run python recovery.py --realgama --gama-model marg --mass-col MassA \
  --mlim-form linear \
  --gama-fits .../make_gama_dmu/G3CFoFGroup.fits --gama-area 238.11
```

The model is a Poisson point process,
`target = -Lambda + sum_i log(integral)`, with
`Lambda = sum_j V_sh[j] * int_{mlim_sh[j]}^inf phi(m) dm`. **Everything hinges
on `mlim(z)`** — it sets the normalisation and the per-object integration
bounds.

### 1. `mlim(z)` is the histogram *mode*, not a limit — 54% of groups fall below it

`turnover_mlim()` takes the **mode of the observed mass histogram** in each
z-bin. For a smooth distribution that sits near the middle, not at the faint
edge. On this catalogue:

| | |
|---|---|
| groups below their own `mlim(z)` | **987 / 1833 (53.8%)** |
| more than 1 sigma below | 567 (30.9%) |
| more than 2 sigma below | 238 (13.0%) |
| more than 3 sigma below | 82 (4.5%) |
| median `(m - mlim)/sigma` | −0.15 |

`marg` integrates each object from `lo_i = mlim[i]` upward, so a group below
`mlim` has its Gaussian peak **outside** the integration range. Beyond ~3 sigma
the integrand underflows and Stan takes the `target += -100` branch — a constant
carrying no information about the HMF, for ~4.5% of the catalogue, with 31%
severely truncated. That is the failure mode.

### 2. `--gama-model marg` is the sharp-cut model — the wrong one here

The argparse help says it outright: *"'marg' (sharp cut), 'marg_comp'
(completeness forward-model)"*. `MARG_CODE`'s own docstring describes
`marg_comp` as replacing *"the sharp mlim cut with the measured completeness
ramp ... all detected groups above the floor (no mlim cut). This is the boundary
fix."* **`marg_comp` is the model built for exactly this problem and it is not
being used.**

### 3. `--mass-col MassA` silently drops 0.315 dex

Reading the column bypasses the rebuild, and the DMU's `MassA` is not on
Driver's scale:

| term | dex |
|---|---|
| `MASS_A = 10` in `make_gama_dmu/config.py` vs Driver's 13.9 | +0.143 |
| the column is never h-scaled; the rebuild applies `(100/H0)` | +0.172 |
| **total** | **+0.315** (measured: +0.315) |

`recovery.py` sets `A_SCALE = 10.0` and `H0 = 100.0`, so it is *internally*
consistent with the DMU convention. But the prior `ms ~ normal(14.0, 1.5)` is
written for Driver-scale masses, so with masses 0.3 dex low the prior is
off-centre. Note the run's own printout ("MassA median 13.446, rebuilt would
give 13.521") understates the gap, because that rebuild also uses A = 10.

### 4. The `mlim` functional form is unstable

Forcing `linear` vs `quad` moves `mlim(0.01)` from **13.03 to 12.40** — 0.63 dex
— and AIC prefers quad (−69.3 vs −66.0). The code comments already flag this:
*"a 0.075 dex shift in the masses flipped it linear -> quad and moved mlim by
0.69 dex at low z"*. `--mlim-form linear` suppresses the symptom, not the cause.

### 5. The structural problem, and the fix this project already has

**GAMA's selection is not a mass limit.** A group enters when it has 5 members
brighter than the flux limit, which depends on how its luminosity function is
sampled, not on its mass. No smooth `mlim(z)` can represent that, which is why
the mode-based estimate lands in the middle of the data.

The 1/Vmax pipeline solves this exactly and per object: `zmax` from the **5th
brightest member's own magnitude**, then `vmax = V(zmax) - V(zmin)` — see
`new_gama_hmf._vmax_from_members` and "Method reference". That is the correct
selection function, already computed, already validated.

The hierarchical model should consume it: replace `V_sh` / `mlim_sh` with a
per-object effective volume `V_i = vmax_i`, so `Lambda = int phi(m) <V(m)> dm`
with the per-group volumes rather than a shell decomposition behind a fitted
mass limit. That removes `turnover_mlim` entirely along with items 1, 2 and 4.

### The run itself is healthy — the problem is bias, not sampling

Run to completion (`python -u`; without `-u` the output buffers and it looks
hung, which it is not):

```
[marg] Rhat=1.005  min ESS=976  divergences=0  treedepth>=10: 0%
```

The sampler is fine. What is wrong is the answer:

| param | "true" | median | sd | bias (sd) |
|---|---|---|---|---|
| `ms` | 13.958 | 14.286 | 0.276 | **+1.19** |
| `lp` | −3.445 | −4.168 | 0.447 | **−1.62** |
| `al` | −1.680 | −1.819 | 0.095 | **−1.46** |
| `be` | 0.630 | 0.553 | 0.140 | −0.55 |

M\* high, phi\* low, alpha steep — a *coherent* 1.2-1.6 sigma pull, not noise.
That is what a model sees when it believes the sample is truncated at `mlim`
while half of it actually lies below: the objects below the limit stop
contributing, so the fit compensates by steepening the slope and dropping the
normalisation.

The script prints the warning itself — `N above mlim: 846 / 1833 (46.2%)` —
right before it fits.

### Suggested order

1. **Fix the mass scale.** Drop `--mass-col MassA`, or add +0.315 dex, so the
   masses and the `ms ~ normal(14.0, 1.5)` prior sit on the same scale.
   *(minutes)*
2. **Replace `mlim(z)` with the per-group Vmax.** This is the real fix and it
   should now be done *first among the substantive changes*, not last.
   `new_gama_hmf._vmax_from_members` already computes the exact selection per
   object — `zmax` from the 5th brightest member, `vmax = V(zmax) − V(zmin)`.
   Feed `V_i` into Stan in place of `V_sh` / `mlim_sh`, so
   `Lambda = sum_i` over per-group volumes rather than shells behind a fitted
   mass limit. That deletes `turnover_mlim` and removes problems 1, 2 and 4
   together — **and it should be fast**, because the selection becomes a
   precomputed per-object number instead of an integral over a soft boundary.
3. **`marg_comp` is not a practical route.** It is the model written for the
   boundary problem, but measured on this catalogue it runs at **48-67 seconds
   per iteration** against ~0.4 s for `marg` — a ~125x slowdown. It reached
   1000 of 3000 iterations in two hours; finishing would take **40+ hours per
   chain**. Iterations are hardcoded at 1500/1500 (line 3662) with no CLI flag.
   Do not plan around it unless the per-object integral is made much cheaper.
4. **The completeness ramp** (`COMP_D50_PTS` / `COMP_W_PTS`) was measured on a
   *mock* at r < 19.65 and is tabulated against `Delta = m - mlim(z)`, so it
   inherits whatever `mlim` does. If step 2 lands it is not needed at all.

Step 1 is minutes. Step 2 is the one that makes the method sound, and on the
evidence above it is also the only one that is computationally viable.


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
