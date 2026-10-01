#!/usr/bin/env bash
#
# chunk_fasta.sh — split a FASTA file into multiple smaller FASTA files.
#
# Usage:
#   ./chunk_fasta.sh <input.fasta> <chunk_size>
#
#   input.fasta   Path to input FASTA file (multi-line sequences supported)
#   chunk_size    Max number of sequences per output chunk
#
# Chunked FASTA files are written to a new folder named
# <input_basename>_chunked, created alongside the input file, containing
# <input_basename>_chunk001.fasta, _chunk002.fasta, etc.

set -euo pipefail

usage() {
  echo "Usage: $0 <input.fasta> <chunk_size>" >&2
  exit 1
}

if [[ $# -ne 2 ]]; then
  usage
fi

input="$1"
chunk_size="$2"

if [[ ! -f "$input" ]]; then
  echo "Error: input file '$input' not found" >&2
  exit 1
fi

if ! [[ "$chunk_size" =~ ^[0-9]+$ ]] || [[ "$chunk_size" -lt 1 ]]; then
  echo "Error: chunk_size must be a positive integer" >&2
  exit 1
fi

indir="$(cd "$(dirname "$input")" && pwd)"
base="$(basename "$input")"
base="${base%.*}"
outdir="${indir}/${base}_chunked"

mkdir -p "$outdir"

awk -v chunk_size="$chunk_size" -v outdir="$outdir" -v base="$base" '
  /^>/ {
    if (seq_count % chunk_size == 0) {
      if (out) close(out)
      chunk_num = int(seq_count / chunk_size) + 1
      out = sprintf("%s/%s_chunk%03d.fasta", outdir, base, chunk_num)
    }
    seq_count++
  }
  { print > out }
  END {
    if (out) close(out)
    n_chunks = (seq_count == 0) ? 0 : int((seq_count - 1) / chunk_size) + 1
    printf "Wrote %d sequence(s) into %d chunk(s) in %s\n", seq_count, n_chunks, outdir
  }
' "$input"
