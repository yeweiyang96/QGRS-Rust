#!/usr/bin/env python3
"""Export NCBI region metadata from genomic.gff and split genomic.fna by seqid."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
from urllib.parse import unquote


CSV_FIELDS = [
    "assembly_id",
    "seqid",
    "region_start",
    "region_end",
    "region_length",
    "source",
    "strand",
    "genome",
    "mol_type",
    "is_circular",
    "name",
    "strain",
    "plasmid_name",
    "taxon_id",
    "fna_header",
    "fna_length",
    "length_match",
    "status",
    "attributes_json",
]

MISMATCH_STATUSES = {"missing_in_fna", "length_mismatch"}


@dataclass
class RegionRow:
    assembly_id: str
    seqid: str
    region_start: int
    region_end: int
    region_length: int
    source: str
    strand: str
    genome: str
    mol_type: str
    is_circular: str
    name: str
    strain: str
    plasmid_name: str
    taxon_id: str
    attributes_json: str
    fna_header: str = ""
    fna_length: Optional[int] = None
    length_match: str = ""
    status: str = "missing_in_fna"

    def as_csv_row(self) -> Dict[str, str]:
        return {
            "assembly_id": self.assembly_id,
            "seqid": self.seqid,
            "region_start": str(self.region_start),
            "region_end": str(self.region_end),
            "region_length": str(self.region_length),
            "source": self.source,
            "strand": self.strand,
            "genome": self.genome,
            "mol_type": self.mol_type,
            "is_circular": self.is_circular,
            "name": self.name,
            "strain": self.strain,
            "plasmid_name": self.plasmid_name,
            "taxon_id": self.taxon_id,
            "fna_header": self.fna_header,
            "fna_length": "" if self.fna_length is None else str(self.fna_length),
            "length_match": self.length_match,
            "status": self.status,
            "attributes_json": self.attributes_json,
        }


def info(message: str) -> None:
    print(f"[INFO] {message}", file=sys.stderr)


def warn(message: str) -> None:
    print(f"[WARN] {message}", file=sys.stderr)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def parse_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in raw.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = unquote(value)
        else:
            attrs[part] = ""
    return attrs


def extract_taxon_id(dbxref: str) -> str:
    for token in dbxref.split(","):
        token = token.strip()
        if token.startswith("taxon:"):
            return token.split(":", 1)[1]
    return ""


def parse_sequence_region(line: str) -> Optional[Tuple[str, int]]:
    parts = line.split()
    if len(parts) != 4:
        return None
    seqid = parts[1]
    try:
        start = int(parts[2])
        end = int(parts[3])
    except ValueError:
        return None
    return seqid, (end - start + 1)


def parse_gff_regions(gff_path: Path, assembly_id: str) -> List[RegionRow]:
    seq_region_lengths: Dict[str, int] = {}
    rows: List[RegionRow] = []

    with open_text(gff_path) as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue

            if line.startswith("##sequence-region "):
                parsed = parse_sequence_region(line)
                if parsed is None:
                    warn(
                        f"{assembly_id}: malformed ##sequence-region at {gff_path.name}:{line_no}"
                    )
                else:
                    seqid, length = parsed
                    seq_region_lengths[seqid] = length
                continue

            if line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                warn(f"{assembly_id}: malformed GFF row at {gff_path.name}:{line_no}")
                continue

            seqid, source, feature_type, start_raw, end_raw, _, strand, _, attrs_raw = cols
            if feature_type != "region":
                continue

            try:
                start = int(start_raw)
                end = int(end_raw)
            except ValueError:
                warn(
                    f"{assembly_id}: invalid region coordinates at {gff_path.name}:{line_no}"
                )
                continue

            region_length = end - start + 1
            attrs = parse_attributes(attrs_raw)

            seq_region_length = seq_region_lengths.get(seqid)
            if seq_region_length is not None and seq_region_length != region_length:
                warn(
                    f"{assembly_id}: region length mismatch between ##sequence-region "
                    f"({seq_region_length}) and feature row ({region_length}) for {seqid}"
                )

            row = RegionRow(
                assembly_id=assembly_id,
                seqid=seqid,
                region_start=start,
                region_end=end,
                region_length=region_length,
                source=source,
                strand=strand,
                genome=attrs.get("genome", ""),
                mol_type=attrs.get("mol_type", ""),
                is_circular=attrs.get("Is_circular", ""),
                name=attrs.get("Name", ""),
                strain=attrs.get("strain", ""),
                plasmid_name=attrs.get("plasmid-name", ""),
                taxon_id=extract_taxon_id(attrs.get("Dbxref", "")),
                attributes_json=json.dumps(attrs, ensure_ascii=False, sort_keys=True),
            )
            rows.append(row)

    return rows


def iter_fasta_records(fna_path: Path) -> Iterator[Tuple[str, str, List[str], int]]:
    header = ""
    seqid = ""
    seq_lines: List[str] = []
    seq_len = 0

    with open_text(fna_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line:
                continue

            if line.startswith(">"):
                if header:
                    yield seqid, header, seq_lines, seq_len
                header = line[1:].strip()
                seqid = header.split()[0] if header else ""
                seq_lines = []
                seq_len = 0
                continue

            seq_lines.append(line)
            seq_len += len(line.strip())

    if header:
        yield seqid, header, seq_lines, seq_len


def write_fasta(path: Path, header: str, seq_lines: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{header}\n")
        for line in seq_lines:
            handle.write(f"{line}\n")


def write_regions_csv(path: Path, rows: List[RegionRow]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.as_csv_row())


def append_summary_rows(writer: csv.DictWriter, rows: List[RegionRow]) -> None:
    for row in rows:
        writer.writerow(row.as_csv_row())


def collect_matches(assembly_dir: Path, pattern: str) -> List[Path]:
    return sorted(assembly_dir.glob(pattern))


def find_single_assembly_file(assembly_dir: Path, suffix: str) -> Optional[Path]:
    matches = collect_matches(assembly_dir, f"*_genomic.{suffix}")
    matches.extend(collect_matches(assembly_dir, f"*_genomic.{suffix}.gz"))
    matches = sorted(matches)
    if len(matches) == 1:
        return matches[0]
    return None


def process_assembly(
    assembly_dir: Path,
    output_root: Path,
) -> Tuple[Optional[List[RegionRow]], int]:
    assembly_id = assembly_dir.name
    gff_path = find_single_assembly_file(assembly_dir, "gff")
    fna_path = find_single_assembly_file(assembly_dir, "fna")

    if gff_path is None or fna_path is None:
        gff_candidates = collect_matches(assembly_dir, "*_genomic.gff") + collect_matches(
            assembly_dir, "*_genomic.gff.gz"
        )
        fna_candidates = collect_matches(assembly_dir, "*_genomic.fna") + collect_matches(
            assembly_dir, "*_genomic.fna.gz"
        )
        warn(
            f"{assembly_id}: skipped (expect exactly one *_genomic.gff(.gz) and "
            f"one *_genomic.fna(.gz), got gff={len(gff_candidates)}, fna={len(fna_candidates)})"
        )
        return None, 0

    rows = parse_gff_regions(gff_path, assembly_id)
    out_assembly_dir = output_root / assembly_id
    out_assembly_dir.mkdir(parents=True, exist_ok=True)
    split_dir = out_assembly_dir / "fna_split"

    by_seqid: Dict[str, List[RegionRow]] = {}
    for row in rows:
        by_seqid.setdefault(row.seqid, []).append(row)

    seen_seqids: set[str] = set()
    written_split_files = 0
    for seqid, header, seq_lines, seq_len in iter_fasta_records(fna_path):
        target_rows = by_seqid.get(seqid)
        if not target_rows:
            continue

        seen_seqids.add(seqid)
        should_write = False
        for row in target_rows:
            row.fna_header = header
            row.fna_length = seq_len
            if seq_len == row.region_length:
                row.length_match = "true"
                row.status = "ok"
                should_write = True
            else:
                row.length_match = "false"
                row.status = "length_mismatch"

        if should_write:
            split_dir.mkdir(parents=True, exist_ok=True)
            write_fasta(split_dir / f"{seqid}.fna", header, seq_lines)
            written_split_files += 1

    mismatch_count = 0
    warned: set[Tuple[str, str]] = set()
    for row in rows:
        if row.seqid not in seen_seqids:
            row.status = "missing_in_fna"
            row.fna_header = ""
            row.fna_length = None
            row.length_match = ""

        if row.status in MISMATCH_STATUSES:
            mismatch_count += 1
            key = (row.seqid, row.status)
            if key not in warned:
                warned.add(key)
                if row.status == "missing_in_fna":
                    warn(f"{assembly_id}: seqid {row.seqid} missing in FNA")
                elif row.status == "length_mismatch":
                    warn(
                        f"{assembly_id}: seqid {row.seqid} length mismatch "
                        f"(region={row.region_length}, fna={row.fna_length})"
                    )

    write_regions_csv(out_assembly_dir / "regions.csv", rows)
    info(
        f"{assembly_id}: regions={len(rows)} mismatches={mismatch_count} "
        f"split_files={written_split_files}"
    )
    return rows, mismatch_count


def discover_assembly_dirs(input_root: Path) -> List[Path]:
    subdirs = sorted(p for p in input_root.iterdir() if p.is_dir())
    if subdirs:
        return subdirs

    # Fallback: treat the root itself as one assembly if no subdirectories exist.
    return [input_root]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Root directory containing assembly subdirectories.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Directory for assembly CSV files and split FASTA outputs.",
    )
    parser.add_argument(
        "--on-mismatch",
        choices=["warn_skip", "strict"],
        default="warn_skip",
        help="Mismatch policy for missing seqid or length mismatch (default: warn_skip).",
    )
    parser.add_argument(
        "--summary-name",
        default="regions.all.csv",
        help="Filename for summary CSV under output-root (default: regions.all.csv).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_root: Path = args.input_root
    output_root: Path = args.output_root

    if not input_root.exists() or not input_root.is_dir():
        print(f"[ERROR] input-root is not a directory: {input_root}", file=sys.stderr)
        return 2

    output_root.mkdir(parents=True, exist_ok=True)

    assembly_dirs = discover_assembly_dirs(input_root)
    if not assembly_dirs:
        warn(f"No assembly directories found under: {input_root}")
        return 0

    summary_path = output_root / args.summary_name
    processed_assemblies = 0
    skipped_assemblies = 0
    total_rows = 0
    total_mismatches = 0

    with summary_path.open("w", encoding="utf-8", newline="") as summary_handle:
        summary_writer = csv.DictWriter(summary_handle, fieldnames=CSV_FIELDS)
        summary_writer.writeheader()

        for assembly_dir in assembly_dirs:
            rows, mismatch_count = process_assembly(assembly_dir, output_root)
            if rows is None:
                skipped_assemblies += 1
                continue

            append_summary_rows(summary_writer, rows)
            processed_assemblies += 1
            total_rows += len(rows)
            total_mismatches += mismatch_count

            if args.on_mismatch == "strict" and mismatch_count > 0:
                print(
                    f"[ERROR] strict mode: mismatches found in {assembly_dir.name}",
                    file=sys.stderr,
                )
                return 1

    info(
        f"Done. processed={processed_assemblies} skipped={skipped_assemblies} "
        f"rows={total_rows} mismatches={total_mismatches} summary={summary_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
