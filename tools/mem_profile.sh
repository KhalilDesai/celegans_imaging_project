#!/bin/sh
# Run the whole pipeline one step per process, recording each step's wall time and peak RSS.
# Same steps and order as `python -m embd.pipeline reproduce` / `all`.
# usage: tools/mem_profile.sh <out_dir>      -> <out_dir>/mem_profile.tsv (step, seconds, peak RSS in MB)
set -eu
OUT="$(mkdir -p "$1" && cd "$1" && pwd)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="$ROOT/.venv/bin/python"
export PYTHONPATH="$ROOT/src" EMBD_OUT_DIR="$OUT/" EMBD_SQLITE="$OUT/embd.sqlite" MPLBACKEND=Agg
case "$(uname)" in Darwin) TIMEFLAG=-l; SCALE=1048576 ;; *) TIMEFLAG=-v; SCALE=1024 ;; esac
REPORT="$OUT/mem_profile.tsv"
printf 'step\tseconds\tpeak_rss_mb\n' > "$REPORT"

run() {  # run <label> [VAR=value ...] -- <pipeline args>
  label="$1"; shift
  envs=""; while [ "$1" != "--" ]; do envs="$envs $1"; shift; done; shift
  log="$OUT/mem_$label.log"
  t0=$(date +%s)
  env $envs /usr/bin/time $TIMEFLAG "$PY" -m embd.pipeline "$@" > "$log" 2>&1 || { tail -20 "$log"; exit 1; }
  t1=$(date +%s)
  # macOS: "<bytes>  maximum resident set size"; GNU time: "Maximum resident set size (kbytes): <kb>"
  rss=$(grep -i 'maximum resident set size' "$log" | grep -oE '[0-9]+' | head -1)
  printf '%s\t%s\t%s\n' "$label" $((t1 - t0)) $((rss / SCALE)) | tee -a "$REPORT"
}

mkdir -p "$OUT/canonical"
rm -f "$EMBD_SQLITE"
"$PY" -c "from embd_data import sheet as s; s.convert_genes_table('$ROOT/data/All_genes_mySQL_table.xlsx', '$OUT/canonical/genes_table.csv'); s.write_gene_key('$ROOT/data/Index_EMBD_ID_Gene_names.xlsx', '$OUT/canonical/gene_key.csv')"
run init-db -- init-db --genes-csv "$OUT/canonical/genes_table.csv" --gene-key "$OUT/canonical/gene_key.csv"
for s in opt gls ms all; do run "phenodist_$s" EMBD_PARAMSET=$s -- phenodist $s; done
run bestmatches -- bestmatches
run padrank -- padrank
run derived -- derived
for s in opt gls ms all; do
  run "vectors_$s" EMBD_PARAMSET=$s -- vectors $s
  run "matrices_$s" EMBD_PARAMSET=$s -- matrices $s
  run "umap_$s" -- umap $s
done
