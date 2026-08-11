# Reproducing Driver+22's GAMA halo mass function, then applying it to new data

## The goal

Two steps, in order:

1. **Reproduce `gamahmf.r` exactly** on the old GAMA data, recovering the
   published GAMA-only MRP parameters:

   | parameter | target |
   |---|---|
   | log10 M\* | 13.51 |
   | log10 phi\* | -3.19 |
   | alpha | -1.27 |
   | beta | 0.47 |

   These are in Driver's own units: **ho = 67.37, omegam = 0.3147, A = 13.9**.
   Do not attempt step 2 until the binned points and the fit both match. If the
   points differ, everything downstream is meaningless.

2. **Apply the identical method to the new data** — the new GAMA DMU (four
   regions including G23, ProFound photometry, r < 19.65, Nessie group finder)
   and the Nessie SDSS catalogue. The comparison "same method, new data" is what
   isolates the effect of the data from the effect of the method.

The previous session got alpha to within 0.035 of the target but never matched
M\* or phi\*, and the last change made it worse. The known problems are listed
under "Where the previous attempt failed".

---

## Files

### Driver's scripts (the specification)

| file | what it is |
|---|---|
| `gamahmf.r` | **The one that matters.** Builds the binned GAMA 1/Vmax HMF and fits the MRP. Everything in step 1 comes from here. |
| `allhmf.r` | The multi-survey figure: reads pre-binned HMFs for GAMA/SDSS/REFLEX/2PIGG, does the h-conversions, fits the combined MRP. Use for the combined plot and for the comparison-data conversions. |
| `plot.R`, `run.R` | The Stan pipeline, unrelated to the 1/Vmax method. Ignore for this exercise. |

Read `gamahmf.r` **end to end before writing any code.** The previous attempt
read it in fragments and guessed the gaps; four of those guesses were wrong and
each one shifted the answer. Lines 224-320 and 355-420 hold nearly everything.

### Data

| file | location | notes |
|---|---|---|
| `G3CFoFGroupv10.fits` | `../data/` | Old GAMA groups (3 equatorial fields, 179.92 deg^2). **Use this for step 1.** |
| `GAMAGalsInGroups.csv` | `../data/` | Member galaxies with the per-galaxy `zmax_19p8`. Needed for Vmax. |
| `reflex.csv` | `../data/` | REFLEX II binned HMF. Columns `x`, `Curve1`. |
| `tpigg.dat` | `../data/` | 2PIGG binned HMF, 4 columns. |
| `elmo.csv` | `../data/` | Tempel+14 SDSS band, 4 columns, headerless. |
| `hmfparams_gsr.csv` | `../data/` | Driver's GSR MCMC chains, 10001 rows, `logMstar,logphi,alpha,beta`. |
| `G3CFoFGroup.fits` | `/Users/00115372/Desktop/my_tools/make_gama_dmu/` | **New** GAMA DMU, 4 regions, 238.11 deg^2. For step 2. |
| `sdss_groups.parquet` | `/Users/00115372/Desktop/my_tools/nessie_tutorials/python/SDSS/` | Nessie SDSS groups. |

`gamahmfGAMA5.csv` and `sdsshmf5.csv` — Driver's own binned output tables — are
**not available**. Without them you cannot check the binned points directly
against his, only the fitted parameters. If they can be obtained from him, get
them: it would turn a four-parameter comparison into a point-by-point one and
make debugging far faster.

---

## Sentinel values — check every file

Each of these silently poisons the result if missed. Check all of them.

**`GAMAGalsInGroups.csv`: `GroupID == 0` is the ungrouped sentinel.**
109052 of 184081 galaxies. It does not corrupt lookups of real groups, but any
`groupby` will build a spurious 109052-member "group". Drop it first.

```python
gig = gig[gig["GroupID"] != 0]
```

**Nessie group IDs: `-1` means ungrouped.**
`make_gama_dmu` adds `GROUP_ID_OFFSET` (1e5) and then maps 99999 -> 0, so the
raw value is -1. A previous attempt filtered on 99999 instead, which lumped
every ungrouped galaxy into one 78754-member group and looked exactly like the
friends-of-friends having percolated.

**Shark mock magnitudes: `-999`.**
`np.isfinite()` does **not** catch it, and a luminosity of `10^(0.4*1003)`
overflows silently. Filter on a physical range instead:

```python
bad = ~np.isfinite(m) | (m < -40) | (m > 0)
```

**Group masses:** check for zeros and negatives before `log10`.
`gamahmf.r` line ~310 also removes one group by hand:
```r
g3c$MassAfunc[g3c$GroupID==100622] = 1E9
```
That is a known-bad object. Reproduce it.

**Empty histogram bins** give `1/sqrt(0) = inf` in the Poisson error term, which
then overflows when squared. Mask them rather than letting `nan_to_num` guess.

---

## Step 1: reproducing `gamahmf.r`

### Constants (all from the file — do not guess)

```
zmin    = 0.015          (line 242)   NOT 0.01
zlimit  = 0.25
multi   = 5              (line 239)
mlimit  = 12.7           (line 240)
magica  = 13.9           (line 238)
myoption= "GAMA"         (line 241)
logbin  = 0.2
massx   = seq(10.3, 16.1, 0.2)
ho      = 67.37
omegam  = 0.3147
area    = 179.92 deg^2
```

### Masses

`myoption = "GAMA"` means the mass is **rebuilt from the velocity dispersion**,
not read from the catalogue's `MassAfunc` column:

```
mymass   = magica * (VelDisp*1000)^2 * Rad50 * parsec*1e6 / (G*msol) * (100/ho)
MassAfunc = mymass / 10^masscorr[Nfof]
```

`masscorr` is the hardcoded multiplicity-debiasing array (index 1-based on
`Nfof`, values for Nfof = 3..22). The `(100/ho)` factor is already inside
`mymass` — do not apply any further h conversion to the masses.

### Vmax — from the member galaxies

This is the piece most likely to be got wrong. Vmax does **not** come from
inverting a mass limit. It comes from the members:

```r
if (Nfof == 2) zmax = sort(members' zmax_19p8, decreasing=TRUE)[2]
else           zmax = sort(members' zmax_19p8, decreasing=TRUE)[multi]

zmax = ifelse(zmax < Zfof,  Zfof,  zmax)      # line 302
zmax = ifelse(zmax > zlimit, zlimit, zmax)    # line 303
vmax = V(zmax) - V(zmin)                      # lines 304-305
```

i.e. the redshift at which the 5th-brightest member drops below the magnitude
limit, so the group would fall below N >= 5. Every group is kept; nothing is cut
on mass.

**Note the bug at lines 307-308:**

```r
g3c$weightszlimit = ifelse(g3c$vmax>vlimit,    vlimit,    g3c$vmax)
g3c$weightszlimit = ifelse(g3c$vmax<vlimitmin, vlimitmin, g3c$vmax)
```

The second line overwrites the first, and its else-branch is `g3c$vmax` rather
than `weightszlimit`, so the upper clip is discarded. It happens not to matter
(line 303 already caps `zmax`, so `vmax <= vlimit`), but reproduce the behaviour
and mention it to Driver.

`vlimitmin = vlimit/1000` (line 229).

**Groups with fewer than `multi` members in the file**: `sort()[5]` returns `NA`
in R. In the old data 206 of 1939 groups are affected — 11 per cent of the
sample. What R's `weighted.hist` does with `NA` weights determines whether those
groups vanish silently or take whole bins with them. Establish this, because it
is a plausible source of the remaining mismatch. Do not simply drop them without
checking what R does.

### Binning and the Eddington correction

```r
gamahmf  = maghist(log10(MassAfunc), breaks=massx)                    # raw counts
gamahmf2 = weighted.hist(log10(MassAfunc), w=1/weightszlimit, breaks=massx)

for i in 1:1001:
    mockmass    = log10(MassAfunc) + rnorm(n, 0, log10MassErr)
    mockcounts[i,] = weighted.hist(mockmass, w=1/weightszlimit, breaks=massx)$counts

meancounts = colMeans(mockcounts)
edb        = meancounts / gamahmf2$counts        # inf/NA -> 1.0
gamay      = gamahmf2$counts / (logbin * edb)
```

The correction is a **Monte-Carlo ratio of counts**, not a deconvolution of a
model. Smearing scatters objects between bins; dividing by `edb` removes the
resulting excess. No MRP is involved.

`log10MassErr` is the multiplicity -> sigma lookup (the "vuvuzela"), floored at
0.1 dex, NA -> 0.03.

**Watch `edb` in sparse bins.** With one or two groups it can reach 10 or more
and crushes or inflates that bin by a dex. Those bins also get `f -> 1`. In the
old data the bins above logM 15.0 hold 1-4 groups and behave wildly.

### Errors

```r
mcerr[i]  = quantile((meancounts[i] - mockcounts[,i])^2, 0.66)^0.5 / meancounts[i]
rootnerr  = 1/sqrt(gamahmf$counts)               # RAW counts, not weighted
gamaf     = sqrt(mcerr^2 + rootnerr^2)
gamaf[is.na]   = 0.9999
gamaf[is.inf]  = 0.0
gamaf          = ifelse(gamaf >= 1, 0.9999, gamaf)
```

### The fit

```r
massfn = chi^2 in log10 space, sigma_log = allf/ln(10)
         + penalty 2*vlimit*sum(phi(m > max(allx)))*logbin

gamafit = optim(par = c(mstarmrp, phimrp, alphamrp, betamrp),
                fn  = massfn,
                control = list(maxit=500, reltol=1e-8, parscale=c(1,1,1,0.5)))
```

Four things matter here:

* **`phi` is LINEAR in the parameter vector**, not log10.
* **The start is the Murray+21 LCDM point**: `mstarmrp = 14.42947`,
  `alphamrp = -1.864908`, `betamrp = 0.7097976`, `phimrp = A/factor` with
  `A = 1.727006e-19` (lines 246-253).
* **`maxit = 500`, one start, no restarts.** This is not incidental. With more
  iterations or multiple starts the optimiser finds a lower-chi^2 but unphysical
  branch at M\* ~ 11.5, where the whole fitted range sits in the exponential
  tail. R's optimiser never reaches it from the Murray start. Reproducing his
  settings is the point.
* **`parscale = c(1,1,1,0.5)`** rescales the parameters internally. scipy has no
  equivalent, and without one you are optimising `mstar ~ 14` alongside
  `phi ~ 1e-4`. Fit `log10(phi)` instead and convert — same effect, and it is
  what makes the fit stable.

Only bins with `gamay > 0` and `gamax > mlimit` enter the fit.

### Parameter errors

Perturb and refit:

```r
cv         = gamay * rnorm(n, 0, cosvariance)
mockgamay  = gamay + gamay*rnorm(n, 0, gamaf) + cv
```
with `cosvariance = cosvar(vlimit/3, 3)` and

```r
cosvar(V,N) = ((219.7 - 52.4*log10(V) + 3.21*(log10(V))^2) / sqrt(N)) / 100
```

---

## Cosmology — the trap that broke the last attempt

`gamahmf.r` uses `ho = 67.37, omegam = 0.3147` **everywhere**: lines 228, 304,
305 all call `cosdist(..., H0=ho)`. Volumes and masses must both use it.

Consequences of getting this wrong:

* Running volumes at h = 1 makes them `0.6737^3 = 0.306x` his, so
  `phi = counts/vmax` comes out **+0.515 dex high**.
* Masses carry `(100/ho)`, so his are **+0.172 dex above** an h = 1 build.

**The specific failure last time:** `recovery.py`'s module-level `H0` was
reassigned at runtime, which changed the masses but **not** the volumes, because
`comoving_distance` caches or captures `H0` rather than reading it live. The
result was masses in one cosmology and volumes in another — worse than either
being consistently wrong. Verify with a printed `vlimit` that actually changes
when you change the cosmology:

```python
# vlimit at ho=67.37 should be ~3.27x the h=1 value
```

Decide the cosmology once, at the top, and confirm every derived quantity
follows it before doing anything else.

---

## Where the previous attempt failed

Best result reached, against the published 13.51 / -3.19 / -1.27 / 0.47:

| parameter | got | target | diff |
|---|---|---|---|
| log M\* | 12.80 (at a bound) | 13.51 | -0.71 |
| log phi\* | -2.14 | -3.19 | +1.05 |
| alpha | -1.235 | -1.27 | **+0.035** |
| beta | 0.381 | 0.47 | -0.089 |

alpha and beta were close; M\* railed against a bound and phi\* followed it up
the M\*-phi\* degeneracy. Known outstanding problems, in the order worth
attacking:

1. **Cosmology inconsistency** (above). Fix first — it alone accounts for
   roughly +0.51 dex of the phi\* offset.
2. **The binned points zigzag.** Values alternated by 0.1-0.4 dex bin to bin
   (-2.677, -2.756, -3.102, -3.011, -3.168) while the raw counts were perfectly
   smooth (97, 120, 198, 198, 251). That means `vmax` is not a smooth function
   of mass, which should not happen. Suspect the member join.
3. **1/Vmax is barely correcting anything.** `zmax_19p8` has a median of 0.298
   against `zlimit = 0.25`, so most groups saturate at `vlimit` and
   `median vmax/vlimit = 0.69`. If the weighting is near-uniform, phi traces the
   raw counts, which turn over — and the corrected HMF should not. Check whether
   this is expected or a symptom of the join.
4. **206 groups with fewer than 5 members in the file.** Establish what R does
   with the resulting `NA` weights.
5. **M\* railing.** Once 1-3 are fixed this may resolve itself; if not, the fit
   is being dragged by the sparse high-mass bins.

**Diagnostic that will save the most time:** print the full binned table —
`logM, N, log10(weighted counts), log10(gamay), rootnerr, mcerr, gamaf` — which
is exactly what `gamahmf.r` line 449 writes out. Comparing that table against
his (if the file can be obtained) localises the problem immediately, rather than
inferring it from four fitted numbers.

---

## Step 2: the new data

Only once step 1 matches.

### New GAMA

`G3CFoFGroup.fits` from `make_gama_dmu`. Columns match the old file plus
`GAMARegion` and `VelDispErr`.

* **Area 238.11 deg^2** (fractional 0.005771988), four regions g09/g12/g15/g23.
* **Do not apply the old `IterCenDec > -3.5` cut** — it selects the equatorial
  fields and deletes G23 entirely. Select on `GAMARegion` instead.
* `MassA = mass_proxy * 10` — A = 10 is already applied. `MassAfunc` adds the
  functional correction, ~0.11 dex above `MassA`.
* Nessie runs at **h = 1.0, omegam = 0.25** (Robotham+11), so the catalogue
  masses are in Msun/h. Driver's script assumes ho = 67.37. Reconcile
  explicitly; do not let the two mix.
* **The magnitude limit changed to r < 19.65** (was 19.8). `zmax_19p8` in the old
  member file is therefore wrong for the new catalogue — it would give every
  group too large a volume. A new member file with per-galaxy zmax at 19.65 and
  the new Nessie GroupIDs is needed. If `make_gama_dmu` does not emit one, zmax
  is computable per galaxy from its absolute magnitude and the k+e correction.

### SDSS

`sdss_groups.parquet`, fractional area 0.2126803, **z < 0.08** (Driver's cut —
beyond it the sample becomes massive-only and the turnover runs away).

**Known problem:** the masses reach logM 15.7 and sit ~1.4 dex above REFLEX at
the same mass. In that volume (1.17e7 Mpc^3 at h = 1) you would expect ~0.5
halos above 10^15.3 and 0.01 above 10^15.7. Those masses are not credible and
they force alpha to rail when SDSS is included. Investigate before using, or cap
at logM < 15.

### Comparison data

`allhmf.r` has the exact h-conversions for each set. REFLEX fractional errors are
`1/sqrt(20)`, with `1/sqrt(3)` at the two endpoints.

---

## Deliverables

1. **GAMA-only plot**: binned points with errors, the fitted MRP with its
   Monte-Carlo band, and Driver's published GAMA-only curve. On the old data
   these should overlap; on the new data the difference is the result.
2. **Combined plot**: GAMA + SDSS + REFLEX II fitted, with 2PIGG shown but not
   fitted, following `allhmf.r`.

Keep them as separate figures.
