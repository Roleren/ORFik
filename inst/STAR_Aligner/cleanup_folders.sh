#!/bin/bash
# Organize STAR logs only after successful processing.
set -eo pipefail
out_dir=$1
[[ -d "$out_dir" ]] || { echo "Output directory does not exist: $out_dir" >&2; exit 1; }
for step in aligned ncRNA_depletion phix_depletion rRNA_depletion tRNA_depletion contaminants_depletion; do
  dir="$out_dir/$step"
  [[ -d "$dir" ]] || continue
  mkdir -p "$dir/LOGS"
  shopt -s nullglob
  logs=("$dir/"*_Log.* "$dir/"*.out.tab)
  if (( ${#logs[@]} )); then mv -- "${logs[@]}" "$dir/LOGS/"; fi
  # STAR normally removes its own temporary directories on success.
  for tmp in "$dir/"*_STARtmp; do
    [[ ! -d "$tmp" ]] || rm -r -- "$tmp"
  done
done
