# HMF figure set — what each file is, and whether its fit may be quoted

Every figure here uses the **Nessie GAMA** catalogue with its **own measured
mass-error curve** (`gama_masserr.csv`, from `vuvuzela_gama.py`) rather than
Driver+22's hardcoded arrays, the mass-to-light cut at 1.0 dex, and the MCMC
posterior rather than Nelder-Mead. Nessie SDSS legs likewise carry
`sdss_masserr.csv` and use the NFW mass scale (`--mass-mode tempel_nfw`).

Driver's own legs — his v10 GAMA comparison points and the Tempel+14 SDSS leg —
keep **his** error curve, because it was measured on his catalogue.

| file | contents | width |
|---|---|---|
| `hmf_gama_nessie.pdf` | GAMA (Nessie) alone, with Driver+22's published **GAMA-only** fit behind | 1 column |
| `hmf_gama_sdss_nessie.pdf` | GAMA + SDSS, both Nessie, no REFLEX or 2PIGG | 1 column |
| `hmf_combined_nessie_gama.pdf` | GAMA (Nessie) + SDSS (Tempel+14) + REFLEX II | full page |
| `hmf_combined_nessie_all.pdf` | GAMA + SDSS (both Nessie) + REFLEX II | full page |

Each has a matching `corner_<name>.pdf` and `chain_<name>.npz`. The parameter
table is produced by `uv run python mrp_table.py [--latex]`, which reads the
chains.

## Exact commands

```bash
S="--mcmc --mcmc-steps 25000 --ml-cut 1.0 --mass-err gama_masserr.csv"
N="--sdss nessie --mass-mode tempel_nfw --sdss-mass-err sdss_masserr.csv"

uv run python combined_hmf.py $S --myoption G   --sdss none --no-extras \
    --pub-fit gama --name gama_nessie
uv run python combined_hmf.py $S --myoption GS  $N --no-extras \
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
