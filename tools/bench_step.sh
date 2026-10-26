#!/bin/sh
# Re-run selected pipeline steps into a fresh dir and compare outputs byte-for-byte with a reference run.
# usage: tools/bench_step.sh <out_dir> <reference_run_dir> <step> [<paramset>]
set -eu
OUT="$1"; REF="$2"; STEP="$3"; PS="${4:-}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$OUT"
[ -f "$OUT/embd.sqlite" ] || cp "$REF/embd.sqlite" "$OUT/embd.sqlite"
for f in gene_key.csv phenoDist_rank_04142022.csv phenoDist_rank_GLS_ONLY05102022.csv phenoDist_rank_MS_ONLY05102022.csv phenoDist_rank_ALL_PARAMS05102022.csv; do
  [ -f "$OUT/$f" ] || [ "$STEP" = phenodist ] || cp "$REF/$f" "$OUT/$f"
done
cp "$REF/gene_key.csv" "$OUT/gene_key.csv"
env EMBD_OUT_DIR="$OUT" EMBD_SQLITE="$OUT/embd.sqlite" PYTHONPATH="$ROOT/src" ${PS:+EMBD_PARAMSET=$PS} \
  "$ROOT/.venv/bin/python" -m embd.pipeline $STEP $PS > "$OUT/$STEP$PS.log" 2>&1
tail -1 "$OUT/timings_port.jsonl"
