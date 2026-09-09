#!/usr/bin/env bash
# Copy the current figure set into the paper directory, as PNG.
#
# Everything lands as PNG.  The four native-unit deliverables and their corner
# plots are written as PDF by combined_hmf.py and are RASTERISED here rather
# than regenerated: emcee is not seeded, so re-running those chains would move
# mrp_table.py's numbers away from the figures they belong to.  Regenerate them
# deliberately (and rebuild the table) if the fits themselves need to change.
set -euo pipefail

DEST=${1:-~/Desktop/papers/hmf/plots}
DEST=$(eval echo "$DEST")
DPI=300
mkdir -p "$DEST"
echo "destination: $DEST"

# --- native-unit deliverables (PDF source) --------------------------------
for f in hmf_gama_nessie hmf_gama_sdss_nessie \
         hmf_combined_nessie_gama hmf_combined_nessie_all \
         corner_gama_nessie corner_gama_sdss_nessie \
         corner_combined_nessie_gama corner_combined_nessie_all
do
    pdftoppm -png -r "$DPI" -singlefile "$f.pdf" "$DEST/$f"
    echo "  rasterised  $f.pdf -> $f.png"
done

# --- shape-normalised fits (normalised_shape_fit.py) ----------------------
# hmf_shape_gama and hmf_shape_gama_sdss are single column; hmf_shape_all is
# full page.  shape_gt14 is deliberately NOT here -- see FIGURES.md.
for f in hmf_shape_gama hmf_shape_gama_sdss hmf_shape_all \
         corner_shape_gama corner_shape_gama_sdss corner_shape_all
do
    cp "$f.png" "$DEST/"
    echo "  copied      $f.png"
done

# --- catalogue and diagnostic figures -------------------------------------
for f in cone_g09 cone_g12 cone_g15 cone_g23 sdss_puck sky_polar \
         vuvuzela_gama vuvuzela_sdss completeness_ridge
do
    cp "$f.png" "$DEST/"
    echo "  copied      $f.png"
done

# --- warn about leftovers -------------------------------------------------
# \includegraphics{name} with no extension prefers PDF, so a stale PDF would
# silently win over the PNG this script just wrote.
shopt -s nullglob
stale=("$DEST"/*.pdf)
if (( ${#stale[@]} )); then
    echo
    echo "  NOTE: ${#stale[@]} PDF(s) still in $DEST.  \\includegraphics{name}"
    echo "        without an extension will prefer these over the new PNGs:"
    printf '          %s\n' "${stale[@]##*/}"
fi
