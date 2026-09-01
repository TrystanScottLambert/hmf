# Omega_M: why the inset was removed, and what replaces it

Recorded 2026-08-31 so the decision can be revisited. `combined_hmf.py
--omega-inset` restores `allhmf.r`'s inset exactly as it was; nothing has been
deleted.

## What the inset plotted

`omega_matter(par)` integrates the **fitted MRP curve** (not the data) over
log10 M from 0 to 18 and divides by the critical density:

```
Omega_M = sum_over_logM [ phi(logM) * 10^logM ] * dlogM * Msun / Mpc^3 / rho_crit
```

i.e. the first moment of the mass function -- number of haloes at each mass
times the mass of each -- as a fraction of rho_crit = 1.26e11 Msun/Mpc^3. The
histogram was that calculation repeated over 1000 posterior draws; the shading
was its 16/84 percentiles, so the width is parameter uncertainty and nothing
else.

Physically: **the fraction of the critical density in collapsed, virialized
haloes**, as implied by extrapolating our fitted mass function over all masses.

## The four reasons it was dropped

**1. 80% of its information is off the panel.** Splitting the difference
between our two anchored fits by mass range:

| range | GSR/Tempel | GSR/Nessie | difference |
|---|---|---|---|
| below 10 (unseen) | 0.0066 | 0.0000 | +0.0065 |
| 10-12.9 (unseen) | 0.0710 | 0.0202 | **+0.0508** |
| 12.9-15.6 (plotted) | 0.1210 | 0.1063 | +0.0146 |
| **total** | **0.1985** | **0.1266** | +0.0719 |

Only 20% of the gap between the two fits is visible in the panel. So the
inset's unique content is the **extrapolation**, not the measurement -- which
is worse than being redundant with the curve.

**2. 39% of the number itself is extrapolated.** Of the 0.199 for the anchored
fit, 61% comes from the fitted range 12.9-15.6, 36% from 10-12.9 and 3% from
below 10. It is not a measurement of the matter in haloes.

**3. The reference line is guaranteed by construction.** Murray+21's MRP is
normalised so that integrating it to zero mass returns Omega_M exactly
(`lcdm_curve`'s `factor`). Plotting our value against a dotted line at Planck's
0.3147 therefore compares our fitted *shape* with LCDM's, not our matter
census with the universe's. It invites a "missing matter" reading that the
measurement cannot support.

**4. The plotted error bar omits the term that dominates.** The total scales as
exactly 10^delta under a uniform mass-scale shift -- no boundary term, unlike
the restricted version -- so:

| shift | Omega_M |
|---|---|
| 0 | 0.198 |
| +0.10 dex | 0.249 |
| +0.15 dex | 0.280 |
| **+0.20 dex** | **0.314  (= Planck, exactly)** |

Measured mass-scale systematics in this project: **+0.153 dex** (Nessie vs
Tempel SDSS, NFW-to-NFW), **~0.17 dex** (unresolved sigma_sky h convention),
**~0.4 dex** (GAMA Robotham vs SDSS NFW estimators). Any one of them exceeds
the +/-0.045 statistical error the shading showed. The 2.6 sigma deficit is
absorbed entirely by calibration.

## What is kept

* **In the table** (`mrp_table.py`): the total Omega_M with its 16/84
  interval, read from the same `omega_draws` the inset used, so the numbers
  cannot drift apart. Railed rows carry `--`.
* **The honest quotation**:

      Omega_M = 0.199 +0.045/-0.024 (stat) x 10^(+/-0.15) (syst, mass scale)

* **The measured companion**: `Omega_M(>12.9) = 0.121 +/- 0.004`, entirely
  inside the data, no extrapolation.
* **The corroboration that matters**: Driver+22's own published GSR fit gives
  **0.205** through the identical integral, so the ~35% deficit belongs to the
  method, not to the Nessie catalogue.

## What would make the total quotable again

Settling the sigma_sky h^-1 Mpc vs Mpc convention (outstanding item 3). It is
~0.17 dex, the same size as the whole deficit, and until it is resolved the
total cannot distinguish "less matter in haloes" from "our masses are 0.2 dex
low".
