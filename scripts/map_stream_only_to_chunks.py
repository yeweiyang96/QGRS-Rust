#!/usr/bin/env python3
import re
from collections import defaultdict
from pathlib import Path

stream_csv=Path('/tmp/trace_stream/combined_sorted.csv')
mmap_csv=Path('/tmp/trace_mmap/combined_sorted.csv')
logfile=Path('/tmp/trace_stream/run.log')
out=Path('/tmp/stream_only_chunk_mapping.csv')

if not stream_csv.exists() or not mmap_csv.exists():
    print('Error: missing one of the combined CSVs:')
    print(stream_csv, stream_csv.exists())
    print(mmap_csv, mmap_csv.exists())
    raise SystemExit(1)
if not logfile.exists():
    print('Error: missing run.log at', logfile)
    raise SystemExit(1)

# read CSVs skipping header
with stream_csv.open() as f:
    s_lines=[l.strip() for l in f if l.strip()]
with mmap_csv.open() as f:
    m_lines=[l.strip() for l in f if l.strip()]

s_hdr=s_lines[0]
m_hdr=m_lines[0]
stream_rows=s_lines[1:]
mmap_set=set(m_lines[1:])

stream_only=[r for r in stream_rows if r not in mmap_set]
print(f'stream_total={len(stream_rows)} mmap_total={len(m_lines)-1} stream_only_count={len(stream_only)}')

# parse run.log into lines
with logfile.open() as f:
    raw_lines=[l.rstrip('\n') for l in f]
log_lines=list(enumerate(raw_lines))

# find indices of STREAM_CHUNK lines
chunk_regex=re.compile(r'DEBUG STREAM_CHUNK \(offset=(\d+), len=(\d+)\): contains target; snippet="([^"]*)"')
merged_regex=re.compile(r'DEBUG MERGED_G4(?: \(offset=(\d+)\))?: start=(\d+) end=(\d+) gscore=(\d+) seq=([a-z]+)')
streamhit_regex=re.compile(r'DEBUG STREAM_HIT \(offset=(\d+)\): start=(\d+) end=(\d+) gscore=(\d+) seq=([a-z]+)')

chunk_positions=[]
for idx,line in log_lines:
    m=chunk_regex.search(line)
    if m:
        off=int(m.group(1)); ln=int(m.group(2)); snip=m.group(3)
        chunk_positions.append((idx,off,ln,snip))

# helper to find nearest previous chunk for a given log index
def find_prev_chunk_for_idx(log_idx):
    for i in range(len(chunk_positions)-1,-1,-1):
        if chunk_positions[i][0] <= log_idx:
            return chunk_positions[i]
    return None

results=[]
not_found=[]
# limit for detailed sample rows - produce full mapping but CSV limited to first N for brevity
SAMPLE_LIMIT = 10000

for row in stream_only:
    parts=row.split(',')
    if len(parts) < 9:
        continue
    start=parts[0]; end=parts[1]; seq=parts[8]
    matches=[]
    # search for exact MERGED_G4 or STREAM_HIT lines that mention start & end & seq
    for idx,line in log_lines:
        if 'MERGED_G4' in line or 'STREAM_HIT' in line:
            if f'start={start} end={end}' in line and seq in line:
                matches.append((idx,line))
    if not matches:
        # fallback: find MERGED_G4/STREAM_HIT lines that include the sequence
        for idx,line in log_lines:
            if ('MERGED_G4' in line or 'STREAM_HIT' in line) and seq in line:
                matches.append((idx,line))
    if not matches:
        not_found.append((row,[]))
        continue
    entry_matches=[]
    for idx,line in matches[:8]:
        prev_chunk=find_prev_chunk_for_idx(idx)
        if prev_chunk:
            entry_matches.append({'log_idx':idx,'log_line':line,'chunk_idx':prev_chunk[0],'chunk_offset':prev_chunk[1],'chunk_len':prev_chunk[2],'chunk_snippet':prev_chunk[3]})
        else:
            entry_matches.append({'log_idx':idx,'log_line':line,'chunk_idx':None,'chunk_offset':None,'chunk_len':None,'chunk_snippet':None})
    results.append((row,entry_matches))

# Summarize
by_chunk=defaultdict(int)
ambiguous=0
for row,ems in results:
    if not ems:
        continue
    offsets=set(e['chunk_offset'] for e in ems if e['chunk_offset'] is not None)
    if len(offsets)>1:
        ambiguous+=1
    for o in offsets:
        by_chunk[o]+=1

print('\nmapped_count=',len(results),'not_found_count=',len(not_found),'ambiguous_chunk_offsets=',ambiguous)
print('\nTop chunk offsets (count):')
for k,v in sorted(by_chunk.items(),key=lambda kv:(-kv[1],kv[0]))[:30]:
    print(k, v)

# write CSV of sample mappings
with out.open('w') as f:
    f.write('row_csv,start,end,sequence,matched_entries_count,chunk_offsets_samples,chunk_snippets_samples,sample_log_lines\n')
    for row,ems in results[:SAMPLE_LIMIT]:
        parts=row.split(',')
        start, end, seq = parts[0], parts[1], parts[8]
        offs=';'.join(str(e['chunk_offset']) for e in ems if e['chunk_offset'] is not None)
        snips='||'.join((e['chunk_snippet'] or '').replace('"','""') for e in ems[:3])
        logs='||'.join(e['log_line'].replace('"','""') for e in ems[:3])
        f.write('"'+row.replace('"','""')+'",%s,%s,%s,%d,"%s","%s","%s"\n' % (start,end,seq,len(ems),offs,snips,logs))

print('\nWrote detailed sample report to', out)
print('Done')
