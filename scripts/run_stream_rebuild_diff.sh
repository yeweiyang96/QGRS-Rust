#!/usr/bin/env bash
set -euo pipefail

BINARY=./target/release/qgrs
INPUT=./dme.fa
MMAP_OUT=/tmp/trace_mmap
STREAM_OUT=/tmp/trace_stream

# ensure mmap results exist (we won't re-run mmap)
if [ ! -d "$MMAP_OUT" ]; then echo "mmap output dir missing: $MMAP_OUT" >&2; exit 1; fi

# prepare stream out dir
rm -rf "$STREAM_OUT"
mkdir -p "$STREAM_OUT"

# run stream mode and capture log
"$BINARY" --file "$INPUT" --output-dir "$STREAM_OUT" --mode stream --max-g4-length 45 2>&1 | tee "$STREAM_OUT/run.log"

# merge CSVs
for d in "$MMAP_OUT" "$STREAM_OUT"; do
  out="$d/combined_sorted.csv"
  files=( "$d"/*.csv )
  if [ ${#files[@]} -eq 0 ]; then echo "no csv in $d"; continue; fi
  header="$(head -n1 "${files[0]}")"
  {
    printf "%s\n" "$header"
    for f in "${files[@]}"; do tail -n +2 "$f"; done | sort -t, -k1,1n -k2,2n

  } > "$out"
  echo "wrote $out"
done

# compute shasum and unified diff
shasum -a 256 "$MMAP_OUT/combined_sorted.csv" "$STREAM_OUT/combined_sorted.csv" > /tmp/trace_md5s_after_dedup.txt
 diff -u "$MMAP_OUT/combined_sorted.csv" "$STREAM_OUT/combined_sorted.csv" > /tmp/trace_diff_max45_after_dedup.patch || true

ls -l "$MMAP_OUT/combined_sorted.csv" "$STREAM_OUT/combined_sorted.csv" /tmp/trace_diff_max45_after_dedup.patch /tmp/trace_md5s_after_dedup.txt
wc -l "$MMAP_OUT/combined_sorted.csv" "$STREAM_OUT/combined_sorted.csv" /tmp/trace_diff_max45_after_dedup.patch || true

echo "done"
