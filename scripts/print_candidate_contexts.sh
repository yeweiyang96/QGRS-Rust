#!/usr/bin/env bash
set -euo pipefail
STREAM_LOG=/tmp/trace_stream/run.log
OUT=/tmp/trace_selected/raw_diff_reports/candidate_contexts.txt
pairs=("49 86" "63 100" "69 104")
> "$OUT"
for p in "${pairs[@]}"; do
  start=$(echo "$p" | awk '{print $1}')
  end=$(echo "$p" | awk '{print $2}')
  echo "\n=== Candidate start=$start end=$end ===" >> "$OUT"
  grep -n "start=$start end=$end" "$STREAM_LOG" | while IFS= read -r m; do
    lineno=$(echo "$m" | cut -d: -f1)
    echo "-- at line $lineno --" >> "$OUT"
    s=$((lineno-8)); if [ $s -lt 1 ]; then s=1; fi
    e=$((lineno+8))
    sed -n "${s},${e}p" "$STREAM_LOG" >> "$OUT"
    echo "----------------" >> "$OUT"
  done
  if ! grep -q "start=$start end=$end" "$STREAM_LOG"; then
    echo "(no match in stream log)" >> "$OUT"
  fi
done

echo "Wrote $OUT"
