#!/usr/bin/env bash
set -euo pipefail

SEQ='gggaccgagtacgggaccgagtacgggaccagtacggg'
STREAM_LOG=/tmp/trace_stream/run.log
OUT=/tmp/trace_selected/raw_diff_reports/stream_contexts_by_seq.txt

if [ ! -f "$STREAM_LOG" ]; then echo "stream log missing: $STREAM_LOG" >&2; exit 1; fi

> "$OUT"

echo "Searching stream log for sequence: $SEQ" >> "$OUT"

# find lines with the sequence and RAW_G4, capture start/end
grep -n "RAW_G4: .*seq=$SEQ" "$STREAM_LOG" | while IFS= read -r line; do
  lineno=$(echo "$line" | cut -d: -f1)
  text=$(echo "$line" | cut -d: -f2-)
  echo "-- match at line $lineno --" >> "$OUT"
  startln=$((lineno-6)); if [ $startln -lt 1 ]; then startln=1; fi
  endln=$((lineno+6))
  sed -n "${startln},${endln}p" "$STREAM_LOG" >> "$OUT"
  echo "----------------" >> "$OUT"
done

# list unique start:end pairs found
echo "\nUnique start:end pairs found for sequence:" >> "$OUT"
grep "RAW_G4: .*seq=$SEQ" "$STREAM_LOG" | sed -n 's/.*start=\([0-9]\+\).*end=\([0-9]\+\).*/start=\1 end=\2/p' | sort -u >> "$OUT"

echo "Wrote $OUT"
