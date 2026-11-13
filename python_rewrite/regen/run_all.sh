#!/bin/bash
# Regenerate the FIGURES.md deliverable set into regen/ after the 2026-11 data loss.
# Inputs: rebuilt GAMA DMU (OVERFACTOR=100, 1807 groups), regenerated Nessie SDSS
# parquets (nessie-py 1.1.1, sqrt(3) fix in estimated_mass), regen/*_masserr.csv.
set -u
cd "$(dirname "$0")/.."
PY=.venv/bin/python
O=regen
S="--mcmc --mcmc-steps 25000 --ml-cut 1.0 --mass-err $O/gama_masserr.csv --fit-min 12.9 --outdir $O"
N="--sdss nessie --mass-mode tempel_rms --sdss-mass-err $O/sdss_masserr.csv"
L="--mcmc --mass-err $O/gama_masserr.csv --sdss-mass-err $O/sdss_masserr.csv --outdir $O"
mkdir -p $O/logs
if [ -z "${1:-}" ]; then
$PY combined_hmf.py $S --myoption G   --sdss none --no-extras --pub-fit gama --name gama_nessie        > $O/logs/gama_nessie.log 2>&1 &
$PY combined_hmf.py $S --myoption GS  $N --no-extras --no-driver-gama     --name gama_sdss_nessie   > $O/logs/gama_sdss_nessie.log 2>&1 &
$PY combined_hmf.py $S --myoption GSR --sdss auto --full-page             --name combined_nessie_gama > $O/logs/combined_nessie_gama.log 2>&1 &
$PY combined_hmf.py $S --myoption GSR $N --full-page                      --name combined_nessie_all > $O/logs/combined_nessie_all.log 2>&1 &
$PY normalised_shape_fit.py $L --name shape_all > $O/logs/shape_all.log 2>&1 &
$PY normalised_shape_fit.py $L --single --legs gama_nessie --show gama_nessie,gama_driver --name shape_gama > $O/logs/shape_gama.log 2>&1 &
$PY normalised_shape_fit.py $L --single --legs gama_nessie,sdss_nessie --show gama_nessie,sdss_nessie --name shape_gama_sdss > $O/logs/shape_gama_sdss.log 2>&1 &
wait
$PY mrp_table.py --chain-dir $O         > $O/mrp_table.txt 2>&1
$PY mrp_table.py --chain-dir $O --latex > $O/mrp_table.tex 2>&1
$PY mrp_table.py --chain-dir $O --shape         > $O/mrp_table_shape.txt 2>&1
$PY mrp_table.py --chain-dir $O --shape --latex > $O/mrp_table_shape.tex 2>&1
echo ALL DONE
fi

# Like-for-like with the committed GSR chains (2026-09-01), which predate the
# REFLEX dn/dM -> dn/dlogM correction becoming the default (2026-09-07).
# Isolates the effect of the catalogue rebuild from the REFLEX change.
if [ "${1:-}" = "--reflex-as-published" ]; then
  R=$O/reflex_as_published; mkdir -p $R/logs
  S2="${S/--outdir $O/--outdir $R} --reflex-as-published"
  $PY combined_hmf.py $S2 --myoption GSR --sdss auto --full-page --name combined_nessie_gama > $R/logs/combined_nessie_gama.log 2>&1 &
  $PY combined_hmf.py $S2 --myoption GSR $N --full-page          --name combined_nessie_all  > $R/logs/combined_nessie_all.log 2>&1 &
  wait
  $PY mrp_table.py --chain-dir $R > $R/mrp_table.txt 2>&1
  echo REFLEX-AS-PUBLISHED DONE
fi

# The conventions the committed (2026-09-01) figures were drawn with: Driver's
# REFLEX normalisation AND Driver's Murray "LCDM" curve.  Both were replaced by
# default on 2026-09-07 ("Two errors that cancelled" in CLAUDE.md).
if [ "${1:-}" = "--old-conventions" ]; then
  R=$O/old_conventions; mkdir -p $R/logs
  S3="${S/--outdir $O/--outdir $R} --reflex-as-published --lcdm driver"
  $PY combined_hmf.py $S3 --myoption G   --sdss none --no-extras --pub-fit gama --name gama_nessie        > $R/logs/gama_nessie.log 2>&1 &
  $PY combined_hmf.py $S3 --myoption GS  $N --no-extras --no-driver-gama     --name gama_sdss_nessie   > $R/logs/gama_sdss_nessie.log 2>&1 &
  $PY combined_hmf.py $S3 --myoption GSR --sdss auto --full-page             --name combined_nessie_gama > $R/logs/combined_nessie_gama.log 2>&1 &
  $PY combined_hmf.py $S3 --myoption GSR $N --full-page                      --name combined_nessie_all > $R/logs/combined_nessie_all.log 2>&1 &
  wait
  $PY mrp_table.py --chain-dir $R         > $R/mrp_table.txt 2>&1
  $PY mrp_table.py --chain-dir $R --latex > $R/mrp_table.tex 2>&1
  echo OLD-CONVENTIONS DONE
fi
