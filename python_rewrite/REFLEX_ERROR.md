# The REFLEX II leg: one arithmetic error, and one thing that is not an error

Recorded 2026-09-02. Two separate issues, and the distinction matters:

* **A genuine arithmetic error** in `allhmf.r` line 307 — a `+1` where
  `log10(ln 10) = 0.362` is required — which puts the entire REFLEX leg
  **0.638 dex too high**. This is a bug, it is in Driver+2022 as published, and
  we reproduce it faithfully.
* **A cosmology difference** between Böhringer's mass function and the
  Planck-normalised LCDM curve everything else is compared against. This is
  **not** a conversion error and cannot be fixed by rescaling. It was hidden by
  the arithmetic error, and it is the more consequential of the two.

Reproduce everything below with `combined_hmf.py --reflex-norm` (corrects the
arithmetic) and `--reflex-fix` (a separate, smaller error-bar issue, see the
last section). Both default off, so nothing already validated has moved.

---

## 1. The arithmetic error

### The line

`allhmf.r` line 307:

```r
reflex$x       = reflex$x + log10(70/ho)                              # line 306
reflex$Curve1  = reflex$Curve1 + 4.0*log10(ho/70) - 14.0 + reflex$x + 1
#                                                                      ^^^
```

Ported at `combined_hmf.py:157` (`load_reflex`).

### What the input is

`reflex.csv` is a digitisation of **Böhringer, Chon & Fukugita 2017 (A&A 608,
A65), Fig. 2** — the column name `Curve1` is WebPlotDigitizer's default. Its
axes read:

| axis | label |
|---|---|
| y | `density [Mpc^-3 h70^4 (10^14 M_sun)^-1]` |
| x | `cluster mass [10^14 h70^-1 M_sun]` |

So `Curve1` is `log10(dn/dM)` **per 10^14 M_sun**, not per dex.

### The correct conversion

Three of the four terms are right:

| term | required | line 307 | |
|---|---|---|---|
| h scaling of the density | `h70^4` -> `4.0*log10(ho/70)` | same | correct |
| density is per 10^14 Msun | `-14.0` | same | correct |
| mass is in h70^-1 Msun | `+log10(70/ho)` (line 306) | same | correct |
| **dn/dM -> dn/dlog10M** | `+ x + log10(ln 10)` | **`+ x + 1`** | **wrong** |

Converting a differential in mass to a differential in dex is
`dn/dlog10M = ln(10) * M * dn/dM`, i.e. `+ log10(M) + log10(ln 10)`.
`log10(ln 10) = 0.3622`, so the `+1` over-normalises by **0.6378 dex**, a factor
4.34.

### Evidence, three independent ways

Böhringer's §5.1 gives an analytic fit to the cumulative mass function (his
eq. 1, parameters in his Table 1), so the digitised points can be checked
against his own published function:

**(a) Normalisation.** Measured offset of Driver's converted points from
Böhringer's eq. 1: **+0.584 dex** with the `+1`, **-0.054 dex** with
`log10(ln 10)`. The residual -0.05 has the right sign and size for his fit
being z = 0 while the data are z = 0.102.

**(b) Shape — this is the decisive one.** The two possible readings of `Curve1`
differ by a factor `M ln 10`, which swings 1.6 dex across the plotted range.
Comparing raw `Curve1` against eq. 1 with **no unit conversion at all**:

| reading of `Curve1` | median resid | scatter | tilt over 2 dex |
|---|---|---|---|
| `dn/dM` per 10^14 Msun | -0.054 | 0.129 | **-0.03** |
| `dn/dlog10M` (per dex) | -0.975 | 0.438 | **-2.00** |

The per-mass reading is flat; the per-dex reading tilts by 2 dex. This does not
depend on the axis label, on any h convention, or on the digitisation quality.

**(c) My implementation of eq. 1 is validated against his own text.**
Integrating it reproduces his §5.2 numbers: **13.09%** of matter in haloes above
10^13 h70^-1 Msun against his stated "14 +/- 1%", and **4.21%** above 10^14
against his "4.4 +/- 0.4%". Checked two ways — differentiating eq. 1, and
integrating by parts using `n(>M)` only — agreeing to four significant figures.

### Consequence

REFLEX carries **54% of the statistical weight** of the combined fit (43 points
at a flat 1/sqrt(20), against 22 SDSS and 14 GAMA points). Being 0.638 dex high,
it drags the fitted curve up with it, which is why the error is invisible in the
published figure: the fit follows REFLEX and the other two legs are pulled along.

| | REFLEX - SDSS | REFLEX - GAMA |
|---|---|---|
| as published | +0.36 dex | -0.17 dex |
| corrected | -0.28 dex | -0.81 dex |

---

## 2. What is NOT the error: the cosmology

**The conversion from Böhringer's cosmology to Planck is not where the mistake
is, and no such conversion is possible.**

Driver converts the **h convention** (h70 -> ho = 67.37) and that part is
correct — it is a units change, Böhringer expresses masses in h70^-1 Msun and
densities in h70^4 Mpc^-3.

But Böhringer's mass function is derived at his **best-fit cosmology,
Omega_m = 0.285 +/- 0.04 and sigma_8 = 0.776 +/- 0.07** (his §2), while
`lcdm_curve` and the rest of the analysis use **Planck 2018, Omega_m = 0.3147,
sigma_8 ~ 0.811**. That is a difference of physical model, not of units. You
cannot rescale a mass function from one sigma_8 to another with a multiplicative
factor — the cluster mass function is exponentially sensitive to sigma_8 and the
correction is mass-dependent.

Its size, from Böhringer's own published integrals:

| fraction of matter in haloes | Böhringer §5.2 | our LCDM curve |
|---|---|---|
| above 10^13 h70^-1 Msun | **14 +/- 1%** | **30.2%** |
| above 10^14 h70^-1 Msun | **4.4 +/- 0.4%** | **10.3%** |

A factor 2.2-2.4, or **+0.37 dex at logM = 14**. This is his prose against our
curve — it involves neither `reflex.csv` nor line 307.

So correctly normalised, his points sit ~0.4 dex below a Planck-normalised
LCDM curve. **That is not evidence that the points are wrong.** His
sigma_8 = 0.776 +/- 0.07 is consistent with Planck's 0.811 at well under 1
sigma, and cluster-count determinations preferring a slightly lower sigma_8 than
the CMB is a long-standing and much-discussed result. The offset is the expected
consequence of his fitted cosmology.

**The methodological problem this exposes.** GAMA and SDSS are direct 1/Vmax
counts; the REFLEX points are an X-ray luminosity function converted through an
LX-M scaling relation calibrated within a specific cosmology, and compared in
the source paper against a Tinker+2008 model at that cosmology. Combining the
three in a single chi^2 and calling the result a measurement is therefore not
clean, and the ~0.3-0.4 dex offset between them is a cosmology/calibration
difference rather than a discrepancy to be fitted away.

The arithmetic error concealed this: `+1` moved REFLEX up by 0.638 dex and
landed it, coincidentally, on a Planck-LCDM curve it has no reason to lie on.

---

## 3. A third, smaller issue: the top error bar

Böhringer §3: *"the binned mass function of the REFLEX sample with 20 clusters
per bin, except for the bin at **the lowest masses**, which has only three
clusters."*

Driver+22 §4.3: *"contain 20 clusters per bin except for **the highest mass
bin**, which contains 3 clusters."*

Böhringer is right by construction — equal-occupancy bins run out at the
low-mass end. `allhmf.r` hedges and inflates **both** ends:

```r
reflexf    = 1/20^0.5 + reflexx*0.0
reflexf[43] = 1/3^0.5   # highest mass -- follows the erroneous text
reflexf[1]  = 1/3^0.5   # lowest mass  -- correct per Böhringer
```

So the highest-mass point (logM = 15.38) carries a 0.577 fractional error
instead of 0.224 and is **down-weighted 6.7x** — the single point that most
constrains the exponential cutoff. `--reflex-fix` removes it.

Two further details from Böhringer that Driver's treatment discards: his own
uncertainty is **18.1%** from marginalising over cosmological and scaling
parameters, not 22.4% from sqrt(20) — but it is a *coherent* band, so
Driver's sqrt(n) is arguably the better choice for a per-point chi^2 even though
it is not what Böhringer quotes. And he attributes the two or three lowest-mass
points to **cosmic variance from the local southern underdensity**; those enter
our fit at full weight, unflagged.

---

## 4. Effect on the fits

GSR, Nessie GAMA, `--ml-cut 1.0 --fit-min 12.9`:

| | logM* | alpha | beta | chi2 |
|---|---|---|---|---|
| Tempel SDSS, as published | 13.98 +0.26/-0.36 | -1.62 | 0.58 | 232 |
| Tempel SDSS, `--reflex-norm --reflex-fix` | 13.47 +0.37/-0.68 | -1.73 | 0.47 | **669** |
| Nessie SDSS, as published | 12.79 +0.53/-0.66 | -0.99 | 0.39 | 231 |
| Nessie SDSS, `--reflex-norm --reflex-fix` | **railed at 11.01** | -0.67 | 0.28 | **804** |

chi^2 triples and the chi^2 profile in logM* becomes monotonic to the prior
bound. Figures: `hmf_reflexnorm_gama.pdf`, `hmf_reflexnorm_all.pdf`.

**This is why `--reflex-norm` is not adopted as the default.** The arithmetic is
right, but applying it does not produce a better fit — it exposes a real ~0.3
dex offset between the X-ray and dynamical mass scales that the combined chi^2
has no way to represent. Fitting a single MRP across that offset is not
meaningful in either configuration.

## 5. Before acting on this

1. **Check Driver's published Fig. 11** against our uncorrected points, to
   confirm the error is in the paper as printed and not only in the script.
2. **Check the corrected normalisation against an independent X-ray mass
   function** — Vikhlinin+2009 appears in Böhringer's own Fig. 4.
3. **Verify our LCDM curve independently.** It is the other half of the factor
   2.4, and CLAUDE.md records that `lcdm_curve`'s normalisation is imposed by
   construction (its `factor` forces Omega_M = 0.3147), so it has not been
   checked against an external calculation.
4. **Consider fitting a relative mass-scale offset** between the X-ray and
   dynamical legs rather than forcing them onto one curve — that turns the
   tension into a measured quantity instead of a chi^2 penalty.
5. **Tell Driver.** Item 1 of this document is a straightforward error with a
   large effect on his headline result.
