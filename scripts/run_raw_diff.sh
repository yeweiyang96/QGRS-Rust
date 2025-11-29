#!/usr/bin/env bash
set -euo pipefail

# Wrapper to run the raw-diff analyzer; keeps the invocation short per user request.
PY=python3
SCRIPT=./scripts/raw_diff_analyzer.py
SUMMARY=/tmp/trace_selected/summary.txt
OUTDIR=/tmp/trace_selected/raw_diff_reports

mkdir -p "$OUTDIR"
$PY "$SCRIPT" --summary "$SUMMARY" --out-dir "$OUTDIR"
