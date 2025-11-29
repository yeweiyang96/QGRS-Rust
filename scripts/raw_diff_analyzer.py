#!/usr/bin/env python3
"""
Compare RAW_G4 candidate sets between MMAP and STREAM for sequences listed in summary.txt.

Produces per-sequence reports under the specified output directory. Uses the per-sequence
trace files generated earlier in `/tmp/trace_selected/seq_*.txt` produced by
`scripts/trace_extractor.py`.

Usage (invoked by wrapper script):
  python3 scripts/raw_diff_analyzer.py --summary /tmp/trace_selected/summary.txt --out-dir /tmp/trace_selected/raw_diff_reports

"""
import argparse
import os
import re
from collections import namedtuple

Candidate = namedtuple('Candidate', ['start', 'end', 'gscore', 'seq'])

RAW_RE = re.compile(r'start=(\d+)\s+end=(\d+)\s+gscore=(\d+)\s+seq=(.+)$')


def parse_seq_file(path):
    """Return two sets: mmap_set, stream_set of Candidate parsed from the file"""
    mmap_set = set()
    stream_set = set()
    if not os.path.isfile(path):
        return mmap_set, stream_set
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if 'DEBUG RAW_G4' not in line:
                continue
            mode = None
            if line.startswith('[MMAP]'):
                mode = 'MMAP'
                content = line[len('[MMAP] '):]
            elif line.startswith('[STREAM]'):
                mode = 'STREAM'
                content = line[len('[STREAM] '):]
            else:
                # sometimes the per-seq file contains lines without the bracket prefix,
                # try to detect presence of markers
                if 'MMAP' in line:
                    mode = 'MMAP'
                    content = line
                elif 'STREAM' in line:
                    mode = 'STREAM'
                    content = line
                else:
                    content = line
            m = RAW_RE.search(content)
            if not m:
                # try a more permissive fallback: find start= and end= and gscore and seq=
                try:
                    start = int(re.search(r'start=(\d+)', content).group(1))
                    end = int(re.search(r'end=(\d+)', content).group(1))
                    gscore = int(re.search(r'gscore=(\d+)', content).group(1))
                    seq = re.search(r'seq=(.+)$', content).group(1).strip()
                except Exception:
                    continue
            else:
                start = int(m.group(1))
                end = int(m.group(2))
                gscore = int(m.group(3))
                seq = m.group(4).strip()
            cand = Candidate(start, end, gscore, seq)
            if mode == 'MMAP':
                mmap_set.add(cand)
            elif mode == 'STREAM':
                stream_set.add(cand)
            else:
                # unknown mode: skip
                pass
    return mmap_set, stream_set


def safe_name(seq):
    s = seq.replace('/', '_').replace('\n', '_')[:120]
    return s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True, help='summary.txt produced earlier')
    parser.add_argument('--out-dir', required=True, help='output reports dir')
    parser.add_argument('--max-samples', type=int, default=50, help='max examples per diff side')
    args = parser.parse_args()

    out_base = args.out_dir
    os.makedirs(out_base, exist_ok=True)

    if not os.path.isfile(args.summary):
        print('summary not found:', args.summary)
        return

    reports = []
    with open(args.summary, 'r', encoding='utf-8', errors='replace') as s:
        header = s.readline()
        for line in s:
            line = line.rstrip('\n')
            if not line:
                continue
            # fields: index\tdiff_count\tmmap_matches\tstream_matches\tout_path\tsequence
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            idx = parts[0].strip()
            diff_count = int(parts[1])
            mmap_matches = int(parts[2])
            stream_matches = int(parts[3])
            out_path = parts[4]
            seq = parts[5]
            if mmap_matches == 0 and stream_matches == 0:
                continue
            mmap_set, stream_set = parse_seq_file(out_path)
            stream_only = sorted(stream_set - mmap_set, key=lambda c: (c.start, c.end, -c.gscore))
            mmap_only = sorted(mmap_set - stream_set, key=lambda c: (c.start, c.end, -c.gscore))

            report_path = os.path.join(out_base, f'rawdiff_{idx}_{safe_name(seq)}.report.txt')
            with open(report_path, 'w', encoding='utf-8') as out:
                out.write(f'Sequence: {seq}\n')
                out.write(f'Index: {idx}\n')
                out.write(f'diff_count: {diff_count}\n')
                out.write(f'mmap_matches: {mmap_matches}\n')
                out.write(f'stream_matches: {stream_matches}\n')
                out.write('\n')
                out.write(f'mmap_total_raw: {len(mmap_set)}\n')
                out.write(f'stream_total_raw: {len(stream_set)}\n')
                out.write(f'stream_only_count: {len(stream_only)}\n')
                out.write(f'mmap_only_count: {len(mmap_only)}\n')
                out.write('\n')
                out.write('--- stream_only (examples) ---\n')
                for c in stream_only[:args.max_samples]:
                    out.write(f'start={c.start} end={c.end} gscore={c.gscore} seq={c.seq}\n')
                if len(stream_only) > args.max_samples:
                    out.write(f'... ({len(stream_only)-args.max_samples} more)\n')
                out.write('\n')
                out.write('--- mmap_only (examples) ---\n')
                for c in mmap_only[:args.max_samples]:
                    out.write(f'start={c.start} end={c.end} gscore={c.gscore} seq={c.seq}\n')
                if len(mmap_only) > args.max_samples:
                    out.write(f'... ({len(mmap_only)-args.max_samples} more)\n')
            reports.append((idx, seq, report_path, len(stream_only), len(mmap_only)))

    # Summary across sequences
    summary_path = os.path.join(out_base, 'rawdiff_summary.txt')
    with open(summary_path, 'w', encoding='utf-8') as s:
        s.write('idx\tseq\tstream_only\tmmap_only\treport\n')
        for rec in reports:
            s.write(f'{rec[0]}\t{rec[1]}\t{rec[3]}\t{rec[4]}\t{rec[2]}\n')

    print(f'Wrote {len(reports)} reports to {out_base}. Summary: {summary_path}')

if __name__ == '__main__':
    main()
