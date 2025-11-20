#!/usr/bin/env python3
"""Compare two CSV hit files and validate their FASTA-backed sequences."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Sequence


@dataclass(frozen=True)
class Hit:
    start: int
    end: int
    length: int
    tetrads: int
    y1: int
    y2: int
    y3: int
    gscore: int
    sequence: str

    @classmethod
    def from_row(cls, row: dict[str, str]) -> "Hit":
        try:
            return cls(
                start=int(row["start"]),
                end=int(row["end"]),
                length=int(row["length"]),
                tetrads=int(row["tetrads"]),
                y1=int(row["y1"]),
                y2=int(row["y2"]),
                y3=int(row["y3"]),
                gscore=int(row["gscore"]),
                sequence=row["sequence"].strip().lower(),
            )
        except KeyError as exc:  # pragma: no cover - defensive parsing
            missing = exc.args[0]
            raise ValueError(f"Missing column '{missing}' in CSV row: {row}") from exc


@dataclass
class ValidationResult:
    ok: bool
    message: str | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_a", type=Path, help="First CSV file (treated as reference A)")
    parser.add_argument("csv_b", type=Path, help="Second CSV file (reference B)")
    parser.add_argument("fasta", type=Path, help="FASTA file containing chromosome/reference sequence")
    parser.add_argument(
        "--chromosome",
        required=True,
        help="FASTA header (without >) to use when validating sequences",
    )
    parser.add_argument(
        "--report-limit",
        type=int,
        default=20,
        help="Maximum number of per-row details to print for each difference direction (default: 20).",
    )
    parser.add_argument(
        "--case-sensitive",
        action="store_true",
        help="Enable case-sensitive FASTA vs CSV sequence comparison (disabled by default).",
    )
    return parser.parse_args()


def read_hits(csv_path: Path) -> List[Hit]:
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return [Hit.from_row(row) for row in reader]


def counter_difference(lhs: Counter[Hit], rhs: Counter[Hit]) -> List[Hit]:
    unique: List[Hit] = []
    for hit, count in lhs.items():
        remainder = count - rhs.get(hit, 0)
        if remainder > 0:
            unique.extend([hit] * remainder)
    unique.sort(key=lambda h: (h.start, h.end, h.sequence))
    return unique


def load_fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, List[str]] = {}
    current_name: str | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                sequences.setdefault(current_name, [])
            else:
                if current_name is None:
                    raise ValueError("FASTA content appears before any header")
                sequences[current_name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def validate_hit(hit: Hit, reference: str, *, case_sensitive: bool) -> ValidationResult:
    if hit.start < 0 or hit.end < 0:
        return ValidationResult(False, "start/end must be non-negative")
    if hit.end < hit.start:
        return ValidationResult(False, "end precedes start")
    if hit.end > len(reference):
        return ValidationResult(False, "end exceeds reference length")
    slice_ = reference[hit.start : hit.end]
    if len(slice_) != hit.length:
        return ValidationResult(False, "length mismatch with FASTA slice")
    left = slice_ if case_sensitive else slice_.lower()
    right = hit.sequence if case_sensitive else hit.sequence.lower()
    if left != right:
        return ValidationResult(False, "sequence mismatch vs FASTA")
    return ValidationResult(True, None)


def report_differences(
    *,
    label: str,
    hits: Sequence[Hit],
    reference: str,
    case_sensitive: bool,
    limit: int,
) -> int:
    print(f"\n{label}: {len(hits)} unique row(s)")
    failures = 0
    capped = hits if limit < 0 else hits[:limit]
    for hit in capped:
        validation = validate_hit(hit, reference, case_sensitive=case_sensitive)
        if not validation.ok:
            failures += 1
        seq_preview = hit.sequence[:30] + ("…" if len(hit.sequence) > 30 else "")
        status = "OK" if validation.ok else f"FAIL ({validation.message})"
        print(
            f"  start={hit.start:>6} end={hit.end:>6} len={hit.length:>3} "
            f"gscore={hit.gscore:>2} status={status} seq={seq_preview}"
        )
    remaining = len(hits) - len(capped)
    if remaining > 0:
        print(f"  … {remaining} more row(s) omitted (raise --report-limit to see all)")
    if limit < 0:  # show all rows
        failures = sum(
            not validate_hit(hit, reference, case_sensitive=case_sensitive).ok for hit in hits
        )
    else:
        failures += sum(
            not validate_hit(hit, reference, case_sensitive=case_sensitive).ok
            for hit in hits[len(capped) :]
        )
    return failures


def main() -> None:
    args = parse_args()
    hits_a = read_hits(args.csv_a)
    hits_b = read_hits(args.csv_b)
    counter_a = Counter(hits_a)
    counter_b = Counter(hits_b)

    seqs = load_fasta_sequences(args.fasta)
    if args.chromosome not in seqs:
        available = ", ".join(sorted(seqs)) or "<none>"
        raise SystemExit(
            f"Chromosome '{args.chromosome}' not found in FASTA. Available headers: {available}"
        )
    reference = seqs[args.chromosome]

    a_minus_b = counter_difference(counter_a, counter_b)
    b_minus_a = counter_difference(counter_b, counter_a)

    print(
        f"Comparing {args.csv_a} (rows={len(hits_a)}) vs {args.csv_b} (rows={len(hits_b)})"
    )
    failures = 0
    failures += report_differences(
        label="Rows only in csv_a",
        hits=a_minus_b,
        reference=reference,
        case_sensitive=args.case_sensitive,
        limit=args.report_limit,
    )
    failures += report_differences(
        label="Rows only in csv_b",
        hits=b_minus_a,
        reference=reference,
        case_sensitive=args.case_sensitive,
        limit=args.report_limit,
    )

    if failures:
        print(f"\nValidation completed with {failures} failure(s).", file=sys.stderr)
        raise SystemExit(1)
    print("\nValidation completed with no sequence mismatches.")


if __name__ == "__main__":
    main()
