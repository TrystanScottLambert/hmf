# HMF figure set — what each file is, and whether its fit may be quoted

Every figure here uses the **Nessie GAMA** catalogue with its **own measured
mass-error curve** (`gama_masserr.csv`, from `vuvuzela_gama.py`) rather than
Driver+22's hardcoded arrays, the mass-to-light cut at 1.0 dex, and the MCMC
posterior rather than Nelder-Mead. Nessie SDSS legs likewise carry
`sdss_masserr.csv` and use **`--mass-mode tempel_rms`**: the NFW profile with
Tempel's own rms dispersion (his eq. 3) in place of Nessie's gapper, and the
`cbrt(3)` -> `sqrt(3)` bug fixed.

The GAMA leg keeps GAMA's own estimator -- Robotham+2011 with A = 13.9,
rebuilt from `VelDisp` and `Rad50` -- because that is what the GAMA catalogue
is calibrated on. Each leg stays on its parent survey's calibration; see
CLAUDE.md, "Forcing a common mass estimator: tried, rejected".

**Why `tempel_rms` and not `tempel_nfw`.** With the gapper dispersion the
chi2 profile in log10 M* is monotonic to the prior bound (dchi2 = 16.2 at
logM* = 13.5, i.e. the data actively drive M* down). Swapping in the rms
dispersion turns that into a genuine interior minimum near logM* = 12.0-12.3.
Reproduce either with `combined_hmf.py --profile-mstar`.

Driver's own legs — his v10 GAMA comparison points and the Tempel+14 SDSS leg —
keep **his** error curve, because it was measured on his catalogue.

| file | contents | width |
|---|---|---|
| `hmf_gama_nessie.pdf` | GAMA (Nessie) alone, with Driver+22's published **GAMA-only** fit behind | 1 column |
| `hmf_gama_sdss_nessie.pdf` | GAMA + SDSS, both Nessie only -- no REFLEX, 2PIGG or Driver+22 points | 1 column |
| `hmf_combined_nessie_gama.pdf` | GAMA (Nessie) + SDSS (Tempel+14) + REFLEX II | full page |
| `hmf_combined_nessie_all.pdf` | GAMA + SDSS (both Nessie) + REFLEX II | full page |

Each has a matching `corner_<name>.pdf` and `chain_<name>.npz`. The parameter
table is produced by `uv run python mrp_table.py [--latex]`, which reads the
chains.

## Exact commands

```bash
S="--mcmc --mcmc-steps 25000 --ml-cut 1.0 --mass-err gama_masserr.csv --fit-min 12.9"
N="--sdss nessie --mass-mode tempel_rms --sdss-mass-err sdss_masserr.csv"

uv run python combined_hmf.py $S --myoption G   --sdss none --no-extras \
    --pub-fit gama --name gama_nessie
uv run python combined_hmf.py $S --myoption GS  $N --no-extras --no-driver-gama \
    --name gama_sdss_nessie
uv run python combined_hmf.py $S --myoption GSR --sdss auto --full-page \
    --name combined_nessie_gama
uv run python combined_hmf.py $S --myoption GSR $N --full-page \
    --name combined_nessie_all
```

## Which fits are quotable

`mrp_table.py` marks a row **RAILED** when its posterior piles against the
prior bound (`mcmc_hmf.BOUNDS`, log10 M* = 11.0). A railed row's credible
interval is set by the prior, not the data, so its four MRP parameters must not
be quoted — **but its binned points and plotted curve are still sound**, and so
is `Omega_M(> 12.7)`, which is an integral over the observed range. See
CLAUDE.md, "The parameters are meaningless but THE CURVES ARE FINE".

Historically only GSR-with-an-anchor is well posed: GAMA-only always rails, and
so does GS, because without REFLEX II nothing pins the exponential cutoff.
Check the table rather than assuming.

## Deliberately still PNG / Nelder-Mead

`driver_fig4.pdf` (the bit-exact reproduction of Driver's figure 4) and the
published-table-2 validation printout. Both are regression tests against R;
converting them would destroy what they test.

## The low-mass cut: 12.9, not Driver's 12.7

`allhmf.r` fits GAMA above 12.7 and SDSS above 12.9, and we reproduce both.
But 12.7 was calibrated on **his** catalogue; Nessie's GAMA completeness limit
is higher, and the 12.80 bin -- the first one his cut admits -- is visibly
incomplete:

| logM | log phi | raw N |
|---|---|---|
| **12.80** | **-4.210** | 36 |
| 13.00 | -3.411 | 68 |
| 13.20 | -3.091 | 131 |

0.8 dex below its neighbour, in the direction the HMF should be rising, on half
the counts. `--fit-min 12.9` drops it. **chi2 falls by ~23 for that one point**
-- it was a ~4.8 sigma outlier on its own -- and alpha *steepens* from -1.42 to
-1.62 rather than flattening.

The cut is justified by the completeness argument, not by the improvement: the
12.80 bin was identified as the anomalous one before any of these fits were
run. CLAUDE.md's warning that `--fit-min` flattens alpha was measured at 13.1
and 13.3, where genuinely complete bins start being deleted; it does not apply
at 12.9.

The SDSS side needs no change -- the Nessie SDSS bins straddling its 12.9 cut
run -3.58, -3.58, -3.46, -3.00, -3.35, -3.39 on 150-320 counts each, with no
cliff.

Bins below the cut are not plotted. Driver shows his as open symbols, but ours
are noise rather than a legible turnover (12.60 -> -4.16 on N=23, 12.40 ->
-2.96 on N=16, 12.20 -> -3.88 on N=11), so showing them would confuse rather
than justify the cut.

## Colours

Each fitted curve carries the colour of the data it was fitted to: **red** for
the Nessie GAMA points and their fit, **cornflower** for Driver+22's points and
the same fit run on his catalogue, **tan** for his published fit. Until
2026-08-28 our fit was drawn in exactly the cornflower of his data points,
which paired the wrong line with the wrong sample.

## Omega_M

The quoted `Omega_M` is the **total** -- the fitted MRP integrated over all mass
-- which is what the figure inset plots and what `mrp_table.py` tabulates. The
anchored fit gives **0.199 +0.045/-0.025**, reproducing Driver+22's own 0.205.

Three caveats belong with it every time:

* **39% of it is extrapolated.** Only 61% comes from the fitted range
  (12.9-15.6); 36% comes from 10-12.9 and 3% from below 10. It is not a
  measurement of the matter in haloes.
* **It scales as exactly 10^delta** under a uniform mass-scale shift. The
  +/-0.15 dex calibration systematic is a factor 1.41, dwarfing the ~15%
  statistical error; a +0.20 dex shift reaches Planck's 0.3147 exactly. So the
  ~2.6 sigma deficit is **not** evidence of missing matter.
* **The reference is not an independent prediction.** Murray's MRP is
  normalised so integrating it to zero mass returns Omega_M by construction, so
  falling below 0.3147 means our fitted *shape* differs from LCDM's -- most of
  it our alpha = -1.62 against his -1.865, extrapolated over ~3 decades.

`Omega_M(>12.9) = 0.121 +/- 0.004` is the part that is 100% inside the data and
is worth keeping alongside. Railed rows carry no Omega_M at all: their tight
intervals come from the prior bound, which is how one of them ends up "21 sigma
below Planck".

## The LCDM reference curve and REFLEX

Since 2026-09-07 the dashed black line is **Tinker+2008 M200c at Planck 2018**
(`lcdm_tinker.py`, cached in `lcdm_tinker.csv`), not `allhmf.r`'s Murray MRP
renormalised so all matter sits in haloes -- that convention is +0.22 dex high.
`--lcdm driver` restores it; the choice is plotting only and changes no fitted
number.

REFLEX II now carries the corrected `dn/dM -> dn/dlog10M` conversion by default
(`allhmf.r` line 307's `+1` should be `log10(ln 10)`, a 0.638 dex error).
`--reflex-as-published` restores Driver's version and is REQUIRED for the
published-table-2 validation, which is otherwise suppressed.

See `REFLEX_ERROR.md`.

## Diagnostics

`combined_hmf.py --profile-mstar` prints the chi2 profile in log10 M*, fixing
it and minimising the other three parameters at each point. This is what
distinguishes "the data prefer a low M*" from "the data do not constrain M*",
which a railed corner plot alone cannot tell you.

## The shape-only fit — `normalised_shape_fit.py`

Each leg's density is multiplied by a constant so its integral over 10^14-10^15
matches the anchor's (Nessie SDSS), the absolute normalisations are discarded,
and one MRP is fitted to the rescaled points.  The legs sit on different mass
calibrations, so their amplitudes are not mutually comparable; this asks what
shape they collectively prefer without reconciling the calibrations first.

| figure | command |
|---|---|
| `hmf_shape_all.png`, `corner_shape_all.png` | `uv run python normalised_shape_fit.py --mcmc --name shape_all` |
| `hmf_shape_gt14.png`, `corner_shape_gt14.png` | `uv run python normalised_shape_fit.py --mcmc --fit-lo 14 --name shape_gt14` |

Full page width (`plotting.FULL_PAGE`), PNG.  The dashed black line is
**Murray+21 LCDM rescaled the same way** -- not the Tinker curve the other
figures use, because on a rescaled axis only the shape is comparable.  The grey
band is the 10^14-10^15 matching range.  Tempel+14 SDSS is plotted but **not
fitted** (the fit is GAMA + Nessie SDSS + REFLEX).

**`--from-chain` redraws from `chain_<name>.npz` without re-sampling.**  Use it
for any restyling; a fresh 25k-step run per figure is not needed to change a
line colour.

### `--fit-lo 14` is a robustness check, not the primary

| | all (>12.9) | >14 only |
|---|---|---|
| log10 M* | 12.60 -0.73/+0.59 | 12.45 -0.95/+1.05 |
| alpha | **-1.02 +/-0.24** | **-1.08 +/-0.50** |
| beta | 0.367 +/-0.06 | 0.349 +/-0.08 |

Same answer, intervals roughly doubled, and alpha becomes unmeasured -- the
faint-end slope lives below the knee by definition.  There is also a
circularity: the amplitudes are matched over 14-15 and then fitted over >14, so
the agreement in that band is partly manufactured.  Report it as *restricting to
the matched range leaves the shape unchanged within errors*, which is a useful
null result, and lead with the full-range fit.

### phi* is a nuisance parameter -- marginalise, do not fix

After rescaling, phi* is the anchor's arbitrary amplitude and carries no
physics.  It is still *constrained* by the fit and degenerate with everything
else (rho(M*, phi*) = -0.85, rho(M*, alpha) = -0.97, rho(M*, beta) = +0.98), so
fixing it at its median halves the M* interval and cuts beta's by 60% while the
medians barely move.  That tightening is manufactured.  Sample it, marginalise,
and footnote it as set by the anchor.

### The parameters look nothing like Driver's; the curves agree

Ours 12.60 / -1.02 / 0.367 against his GSR 14.13 / -1.68 / 0.63.  Put his fit
on the same rescaled axis and the two curves differ by **0.059 dex rms over
13-15** (0.007 dex at logM 14), with chi2 191.6 against 196.8 on the same 84
points -- Delta chi2 = 5.2 for 4 parameters, i.e. just outside 1 sigma.  The
parameter difference is motion along the degeneracy ridge, not disagreement.
Both fits are poor in absolute terms (chi2/dof ~ 2.4): the legs still disagree
by ~0.5 dex below logM 14 after rescaling, which is the calibration spread and
the actual open problem.
