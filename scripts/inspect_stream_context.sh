#!/usr/bin/env bash
set -euo pipefail

REPORT=/tmp/trace_selected/raw_diff_reports/rawdiff_2_gggaccgagtacgggaccgagtacgggaccagtacggg.report.txt
STREAM_LOG=/tmp/trace_stream/run.log
OUT=/tmp/trace_selected/raw_diff_reports/stream_contexts.txt

if [ ! -f "$REPORT" ]; then echo "report missing: $REPORT" >&2; exit 1; fi
if [ ! -f "$STREAM_LOG" ]; then echo "stream log missing: $STREAM_LOG" >&2; exit 1; fi

> "$OUT"

echo "Inspecting stream-only RAW candidates from report: $REPORT" >> "$OUT"

# extract start/end pairs more robustly using regex across the stream_only section
awk '/^--- stream_only/,/^$/' "$REPORT" | sed '1d' | sed '/^$/q' > /tmp/__stream_only_lines.txt

# find all start=NN end=MM occurrences
grep -oE 'start=[0-9]+[[:space:]]+end=[0-9]+' /tmp/__stream_only_lines.txt | sort -u > /tmp/__stream_only_pairs.txt

while IFS= read -r pair; do
  start=$(echo "$pair" | sed -E 's/start=([0-9]+)\s+end=.*/\1/')
  end=$(echo "$pair" | sed -E 's/.*end=([0-9]+)/\1/')
  if [ -z "${start:-}" ] || [ -z "${end:-}" ]; then
    echo "Could not parse start/end from pair: $pair" >> "$OUT"
    continue
  fi
  echo "\n=== Candidate start=$start end=$end ===" >> "$OUT"
  # find matching lines in stream log and print 8 lines of context around each match
  grep -n "RAW_G4: start=$start end=$end" "$STREAM_LOG" | while IFS= read -r match; do
    lineno=$(echo "$match" | cut -d: -f1)
    echo "-- found at line $lineno --" >> "$OUT"
    startln=$((lineno-6)); if [ $startln -lt 1 ]; then startln=1; fi
    endln=$((lineno+6))
    sed -n "${startln},${endln}p" "$STREAM_LOG" >> "$OUT"
    echo "----------------" >> "$OUT"
  done
  # if no matches found, record that
  if ! grep -q "RAW_G4: start=$start end=$end" "$STREAM_LOG"; then
    echo "(no RAW_G4 match in stream log for start=$start end=$end)" >> "$OUT"
  fi

done < /tmp/__stream_only_lines.txt

echo "Wrote context to $OUT"
