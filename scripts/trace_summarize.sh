#!/usr/bin/env bash
set -euo pipefail

OUTDIR=/tmp/trace_selected
SUMMARY="$OUTDIR/summary.txt"

if [ ! -f "$SUMMARY" ]; then
  echo "summary not found: $SUMMARY" >&2
  exit 1
fi

echo "Using summary: $SUMMARY"

# Iterate rows; summary is tab-separated with header
tail -n +2 "$SUMMARY" | while IFS=$'\t' read -r index diff_count mmap_matches stream_matches out_path seq; do
  # ensure numeric defaults
  mmap_matches=${mmap_matches:-0}
  stream_matches=${stream_matches:-0}
  if [ "$mmap_matches" -eq 0 ] && [ "$stream_matches" -eq 0 ]; then
    continue
  fi

  printf "\n===== seq #%s\t diff=%s\t mmap=%s\t stream=%s\npath=%s\nseq=%s\n" "$index" "$diff_count" "$mmap_matches" "$stream_matches" "$out_path" "$seq"

  for mode in MMAP STREAM; do
    echo "--- $mode counts ---"
    for marker in RAW_G4 FAMILY MERGED_G4 STREAM_HIT; do
      cnt=$(grep -F "[$mode]" "$out_path" | grep -c "$marker" || true)
      printf "%s: %s\n" "$marker" "$cnt"
    done

    echo
    echo "--- $mode samples (first 8) ---"
    grep -F "[$mode]" "$out_path" | head -n8 || true
    echo
  done

  echo "---- diff patch context ----"
  grep -nF "$seq" /tmp/trace_diff_max45.patch || echo "(none in patch)"
  echo "--------------------------------------------------"

done
