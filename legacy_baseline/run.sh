#!/bin/sh
# Run one legacy-baseline step in the container. Repo is mounted read-only (original code
# cannot be modified); outputs go to build/legacy_work (the "Z:" folder lives inside it).
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec docker run --rm --platform linux/amd64 \
  -v "$ROOT":/repo:ro \
  -v "$ROOT/build/canonical":/canon:ro \
  -v "$ROOT/build/legacy_work":/work \
  -e EMBD_GENES_CSV=/canon/genes_table.csv -e EMBD_GENE_KEY_CSV=/canon/gene_key.csv \
  ${EMBD_DOCKER_ARGS:-} embd-legacy:py27 "$@"
