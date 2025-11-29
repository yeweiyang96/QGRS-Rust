#!/usr/bin/env python3
"""
Trace extractor for QGRS run logs

This script parses a unified diff (CSV patch) and the combined_sorted CSV header
to discover which CSV column contains the sequence string, then selects the top-N
sequences with most diff entries and extracts matching pipeline debug lines from
both mmap and stream run logs. It writes per-sequence trace files and a summary.

Usage example (do NOT run automatically; user requested generation only):

  python3 scripts/trace_extractor.py \
    --diff /tmp/trace_diff_max45.patch \
    --mmap-csv /tmp/trace_mmap/combined_sorted.csv \
    --stream-csv /tmp/trace_stream/combined_sorted.csv \
    --mmap-log /tmp/trace_mmap/run.log \
    --stream-log /tmp/trace_stream/run.log \
    --top-n 10 \
    --out-dir /tmp/trace_selected

Output:
 - per-sequence trace files in the `--out-dir` directory, named `seq_<i>.txt`.
 - `summary.txt` in the same dir with sequence, diff-count and file path.

Note: This script intentionally avoids heavy memory usage and reads logs line-by-line.
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict


def find_sequence_column(csv_path):
    """Return the 0-based index of the sequence column.
    Heuristics: look for header names containing 'seq' or 'sequence', fallback to last column.
    """
    with open(csv_path, 'r', encoding='utf-8', errors='replace') as f:
        header = f.readline().strip()
        if not header:
            raise RuntimeError(f"Empty or missing header in {csv_path}")
        # CSV may be comma-separated; use csv.reader to parse header robustly
        reader = csv.reader([header])
        cols = next(reader)
        # normalize
        lower = [c.lower() for c in cols]
        for name in ('seq', 'sequence', 'seq_str', 'g4_sequence'):
            if name in lower:
                return lower.index(name)
        # fallback to last column
        return len(cols) - 1


def parse_diff_for_sequences(diff_path, seq_col_index, sample_csv_for_schema=None):
    """Parse unified diff and return Counter of sequence -> occurences.
    We scan lines that start with '+' or '-' and are not the diff metadata (i.e. skip '+++'/'---'/'@@').
    We attempt to CSV-parse the line content and extract the sequence column by index.
    """
    counts = Counter()
    if not os.path.isfile(diff_path):
        raise FileNotFoundError(diff_path)
    with open(diff_path, 'r', encoding='utf-8', errors='replace') as f:
        for raw in f:
            if not raw:
                continue
            if raw.startswith('+++') or raw.startswith('---') or raw.startswith('@@'):
                continue
            if raw[0] not in ('+', '-'):
                continue
            line = raw[1:].rstrip('\n')
            # some diff lines may include context or be empty
            if not line.strip():
                continue
            # CSV-parse the line
            try:
                # csv.reader expects an iterable of lines
                row = next(csv.reader([line]))
            except Exception:
                # fallback: treat whole line as the sequence
                seq = line.strip()
                counts[seq] += 1
                continue
            if seq_col_index < len(row):
                seq = row[seq_col_index].strip()
            else:
                # fallback: use last column
                seq = row[-1].strip() if row else ''
            counts[seq] += 1
    return counts


def scan_logs_for_sequence(sequence, mmap_log_path, stream_log_path, out_dir, max_lines_per_section=10000):
    """Extract matching DEBUG lines for `sequence` from both logs and write to a file in out_dir.
    Returns the output filepath.
    """
    os.makedirs(out_dir, exist_ok=True)
    safe_name = sequence.replace('/', '_').replace('\n', '_')[:120]
    out_path = os.path.join(out_dir, f'seq_{safe_name}.txt')

    # markers to capture (include both RAW/FAMILY/MERGED/STREAM_HIT)
    markers = ('RAW_G4', 'FAMILY', 'MERGED_G4', 'STREAM_HIT')

    def scan_file(path, mode_label, writer):
        if not path or not os.path.isfile(path):
            writer.write(f"--- {mode_label}: MISSING log file: {path}\n")
            return 0
        matches = 0
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if sequence in line and any(m in line for m in markers):
                    writer.write(f"[{mode_label}] "+line)
                    matches += 1
                    if matches >= max_lines_per_section:
                        writer.write(f"... truncated after {max_lines_per_section} matches in {mode_label}\n")
                        break
        return matches

    with open(out_path, 'w', encoding='utf-8') as out:
        out.write(f"Sequence: {sequence}\n\n")
        out.write("--- mmap log matches ---\n")
        mcount = scan_file(mmap_log_path, 'MMAP', out)
        out.write(f"--- end mmap ({mcount} matches) ---\n\n")
        out.write("--- stream log matches ---\n")
        scount = scan_file(stream_log_path, 'STREAM', out)
        out.write(f"--- end stream ({scount} matches) ---\n")
    return out_path, mcount, scount


def main():
    parser = argparse.ArgumentParser(description='Extract pipeline logs for top-N differing sequences')
    parser.add_argument('--diff', required=True, help='Unified diff patch file (CSV diff)')
    parser.add_argument('--mmap-csv', required=True, help='Combined mmap CSV to learn header/sequence column')
    parser.add_argument('--stream-csv', required=True, help='Combined stream CSV (not required for schema but accepted)')
    parser.add_argument('--mmap-log', required=True, help='MMAP run log (e.g. /tmp/trace_mmap/run.log)')
    parser.add_argument('--stream-log', required=True, help='STREAM run log (e.g. /tmp/trace_stream/run.log)')
    parser.add_argument('--top-n', type=int, default=10, help='Top N differing sequences to extract')
    parser.add_argument('--out-dir', default='/tmp/trace_selected', help='Output dir for per-sequence traces')
    args = parser.parse_args()

    print(f"Reading CSV header from {args.mmap_csv} to determine sequence column...", file=sys.stderr)
    seq_col = find_sequence_column(args.mmap_csv)
    print(f"Using sequence column index: {seq_col}", file=sys.stderr)

    print(f"Parsing diff {args.diff} to count differing sequences...", file=sys.stderr)
    counts = parse_diff_for_sequences(args.diff, seq_col)
    if not counts:
        print("No differing sequences found in the diff.", file=sys.stderr)
        sys.exit(0)

    top = counts.most_common(args.top_n)
    print(f"Top {len(top)} differing sequences:" , file=sys.stderr)
    for i, (seq, cnt) in enumerate(top, 1):
        print(f"{i}. count={cnt} seq={seq[:80]}...", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    summary_lines = []
    for i, (seq, cnt) in enumerate(top, 1):
        out_path, mcount, scount = scan_logs_for_sequence(seq, args.mmap_log, args.stream_log, args.out_dir)
        summary_lines.append((i, seq, cnt, mcount, scount, out_path))

    summary_file = os.path.join(args.out_dir, 'summary.txt')
    with open(summary_file, 'w', encoding='utf-8') as s:
        s.write('index\tdiff_count\tmmap_matches\tstream_matches\tout_path\tsequence\n')
        for rec in summary_lines:
            i, seq, cnt, mcount, scount, out_path = rec
            s.write(f"{i}\t{cnt}\t{mcount}\t{scount}\t{out_path}\t{seq}\n")

    print(f"Per-sequence traces written under {args.out_dir}. Summary: {summary_file}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
