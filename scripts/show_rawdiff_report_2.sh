#!/usr/bin/env bash
set -euo pipefail

REPORT=/tmp/trace_selected/raw_diff_reports/rawdiff_2_gggaccgagtacgggaccgagtacgggaccagtacggg.report.txt
if [ ! -f "$REPORT" ]; then
  echo "report not found: $REPORT" >&2
  exit 1
fi

echo "=== Showing report: $REPORT ==="
# print header
sed -n '1,60p' "$REPORT"

echo
echo "--- stream_only examples (first 20) ---"
awk '/^--- stream_only/,/^$/' "$REPORT" | sed '1d' | sed -n '1,20p'

echo
echo "--- mmap_only examples (first 20) ---"
awk '/^--- mmap_only/,/^$/' "$REPORT" | sed '1d' | sed -n '1,20p'

echo
echo "--- Full counts ---"
grep -nE 'mmap_total_raw|stream_total_raw|stream_only_count|mmap_only_count' "$REPORT" || true
