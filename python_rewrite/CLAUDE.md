# Driver+22's halo mass function: reproduction, and the same method on new data

## Status

**Step 1 (reproduce `gamahmf.r`) is done and verified bit-exactly against R.**
The binned HMF matches published Table 1; the fit matches R's `optim` to 10
significant figures; the combined multi-survey fit matches published Table 2 to
better than 0.035 in every parameter.

**Step 2 (new GAMA DMU) is done for GAMA.** The Nessie catalogue has been run
through the identical method, both GAMA-only and in the combined fit.

**Outstanding:** the Nessie SDSS leg (see "Remaining work"). Everything else is
built and validated.

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

## Remaining work: the Nessie SDSS leg

`sdss_groups.parquet` (84,890 groups) and `sdss_galaxies.parquet` (584,447
galaxies) in `/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS/`.
Groups carry `multiplicity`, `median_redshift`, `r50`,
`velocity_dispersion_gap`, `mass_proxy`, `estimated_mass`.

Resolve before building:

1. **No absolute magnitude, no k+e.** Galaxies have `rmag` and `zobs` only, so
   the zmax inversion needs a k-correction decided on — unlike the GAMA DMU,
   where the correction was recoverable from `make_gama_dmu` and self-consistent
   by construction.
2. **Two ID columns**, `GroupID` and `group_id` (Tempel's original vs Nessie's
   presumably). Disambiguate before joining.
3. **Masses reach logM 15.7**, roughly 1.4 dex above REFLEX at the same mass. In
   that volume you would expect ~0.01 haloes above 10^15.7. Given what group
   300223 did to the GAMA fit through the penalty range, check this before it
   goes anywhere near `max(allx)`. Consider capping at logM < 15.

Driver's own SDSS (Tempel+14) stays as the `S` leg meanwhile: 4847 groups,
`volumesdss` = 3.100035e+07 Mpc³, area 7221 deg², z < 0.08, r < 17.77,
`logbin = 0.1`, masses `mass * 1e12 * (67.8/ho)` with no A factor and no
multiplicity debiasing.

---

## Deliverables

1. **GAMA-only plot** — binned points with errors, the fitted MRP with its
   Monte-Carlo band, Driver's result for comparison. `hmf_gama_only_nessie.pdf`.
2. **Combined plot** — GAMA + SDSS + REFLEX II fitted, 2PIGG and Tempel shown but
   not fitted, following `allhmf.r`. `hmf_combined_nessie.pdf`.

Keep them as separate figures.
