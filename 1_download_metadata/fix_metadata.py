#!/usr/bin/env python3
"""
fix_metadata.py — Validate and fix iModulonMiner metadata TSV files.

Addresses every issue encountered when feeding metadata into the
iModulonMiner 2_process_data Nextflow pipeline (with fastq-dl):

  1. Column-count mismatches   – rows with fewer tabs than the header
  2. Empty trailing rows       – blank lines that cause assertion errors
  3. Missing/invalid fields    – LibraryLayout, Platform, Experiment
  4. Whitespace contamination  – leading/trailing spaces in fields
  5. Duplicate Experiment IDs  – the pipeline asserts uniqueness

Produces a fixed TSV and a human-readable log of everything it changed.

Usage:
    python fix_metadata.py INPUT.tsv [-o OUTPUT.tsv]

If -o is omitted the output is written to INPUT_fixed.tsv.
"""

import argparse
import csv
import os
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# Validation constants
# ---------------------------------------------------------------------------
VALID_LAYOUTS = {"SINGLE", "PAIRED"}
VALID_PLATFORMS = {"ILLUMINA", "ION_TORRENT", "ABI_SOLID", "LS454",
                   "BGISEQ", "CAPILLARY", "COMPLETE_GENOMICS",
                   "HELICOS", "OXFORD_NANOPORE", "PACBIO_SMRT"}

# Columns the pipeline reads directly — must exist in the header
REQUIRED_COLUMNS = {"Experiment", "LibraryLayout", "Platform", "Run"}
# Optional columns referenced during routing
OPTIONAL_COLUMNS = {"R1", "R2"}


def fix_metadata(input_path: str, output_path: str) -> None:
    """Read, validate, fix, and write the metadata TSV."""

    issues: list[str] = []      # human-readable log of every change
    dropped_rows: list[str] = []

    # ------------------------------------------------------------------
    # 1.  Read raw lines so we can detect column-count mismatches
    # ------------------------------------------------------------------
    with open(input_path, "r", newline="") as fh:
        raw_lines = fh.readlines()

    if not raw_lines:
        print("ERROR: input file is empty.", file=sys.stderr)
        sys.exit(1)

    # Parse header
    header_line = raw_lines[0].rstrip("\n\r")
    headers = header_line.split("\t")
    n_cols = len(headers)

    # Strip whitespace from header names
    cleaned_headers = [h.strip() for h in headers]
    if cleaned_headers != headers:
        issues.append(f"Header: stripped whitespace from column names")
        headers = cleaned_headers

    # Check for required columns
    header_set = set(headers)
    missing_required = REQUIRED_COLUMNS - header_set
    if missing_required:
        print(f"ERROR: missing required columns: {missing_required}",
              file=sys.stderr)
        sys.exit(1)

    # Ensure R1 and R2 columns exist (pipeline expects them)
    for col in ("R1", "R2"):
        if col not in header_set:
            headers.append(col)
            n_cols += 1
            issues.append(f"Header: added missing '{col}' column")

    # Column index lookup
    col_idx = {name: i for i, name in enumerate(headers)}

    # ------------------------------------------------------------------
    # 2.  Process data rows
    # ------------------------------------------------------------------
    fixed_rows: list[list[str]] = []
    seen_experiments: dict[str, int] = {}   # experiment -> first row num

    for row_num, raw_line in enumerate(raw_lines[1:], start=2):
        line = raw_line.rstrip("\n\r")

        # Skip completely empty lines
        if not line or line.replace("\t", "").strip() == "":
            issues.append(f"Row {row_num}: dropped empty row")
            dropped_rows.append(f"Row {row_num}")
            continue

        fields = line.split("\t")

        # Strip whitespace from every field
        stripped = [f.strip() for f in fields]
        if stripped != fields:
            issues.append(f"Row {row_num}: stripped whitespace from fields")
            fields = stripped

        # ----------------------------------------------------------
        # 2a. Fix column count
        # ----------------------------------------------------------
        if len(fields) < n_cols:
            deficit = n_cols - len(fields)
            issues.append(
                f"Row {row_num}: had {len(fields)} columns "
                f"(expected {n_cols}), padded {deficit} empty field(s)")
            fields.extend([""] * deficit)
        elif len(fields) > n_cols:
            # Extra trailing empty fields — trim them
            extra = fields[n_cols:]
            if all(f == "" for f in extra):
                issues.append(
                    f"Row {row_num}: had {len(fields)} columns "
                    f"(expected {n_cols}), trimmed {len(extra)} trailing "
                    f"empty field(s)")
                fields = fields[:n_cols]
            else:
                issues.append(
                    f"Row {row_num}: WARNING — has {len(fields)} columns "
                    f"with non-empty extra fields; keeping first {n_cols}")
                fields = fields[:n_cols]

        # ----------------------------------------------------------
        # 2b. Validate Experiment ID
        # ----------------------------------------------------------
        exp_id = fields[col_idx["Experiment"]]
        if not exp_id:
            issues.append(f"Row {row_num}: dropped row with empty Experiment")
            dropped_rows.append(f"Row {row_num} (empty Experiment)")
            continue

        # ----------------------------------------------------------
        # 2c. Check for duplicate Experiment IDs
        # ----------------------------------------------------------
        if exp_id in seen_experiments:
            issues.append(
                f"Row {row_num}: DUPLICATE Experiment '{exp_id}' "
                f"(first seen row {seen_experiments[exp_id]}); dropped")
            dropped_rows.append(f"Row {row_num} (dup {exp_id})")
            continue
        seen_experiments[exp_id] = row_num

        # ----------------------------------------------------------
        # 2d. Validate LibraryLayout
        # ----------------------------------------------------------
        layout = fields[col_idx["LibraryLayout"]]
        if layout.upper() in VALID_LAYOUTS and layout != layout.upper():
            issues.append(
                f"Row {row_num} ({exp_id}): normalised LibraryLayout "
                f"'{layout}' -> '{layout.upper()}'")
            fields[col_idx["LibraryLayout"]] = layout.upper()
        elif layout.upper() not in VALID_LAYOUTS:
            issues.append(
                f"Row {row_num} ({exp_id}): WARNING — invalid "
                f"LibraryLayout '{layout}' (expected SINGLE or PAIRED)")

        # ----------------------------------------------------------
        # 2e. Validate Platform
        # ----------------------------------------------------------
        platform = fields[col_idx["Platform"]]
        if platform.upper() in VALID_PLATFORMS and platform != platform.upper():
            issues.append(
                f"Row {row_num} ({exp_id}): normalised Platform "
                f"'{platform}' -> '{platform.upper()}'")
            fields[col_idx["Platform"]] = platform.upper()
        elif platform.upper() not in VALID_PLATFORMS:
            issues.append(
                f"Row {row_num} ({exp_id}): WARNING — unrecognised "
                f"Platform '{platform}'")

        # ----------------------------------------------------------
        # 2f. Validate Run field is not empty
        # ----------------------------------------------------------
        run_val = fields[col_idx["Run"]]
        if not run_val:
            issues.append(
                f"Row {row_num} ({exp_id}): WARNING — empty Run field")

        fixed_rows.append(fields)

    # ------------------------------------------------------------------
    # 3.  Write fixed output
    # ------------------------------------------------------------------
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(fixed_rows)

    # ------------------------------------------------------------------
    # 4.  Summary
    # ------------------------------------------------------------------
    n_warnings = sum(1 for i in issues if "WARNING" in i)
    n_fixes = len(issues) - n_warnings

    print(f"\n{'='*60}")
    print(f"  iModulonMiner Metadata Fixer")
    print(f"{'='*60}")
    print(f"  Input:       {input_path}")
    print(f"  Output:      {output_path}")
    print(f"  Header cols: {n_cols}")
    print(f"  Input rows:  {len(raw_lines) - 1}  (excluding header)")
    print(f"  Output rows: {len(fixed_rows)}")
    print(f"  Dropped:     {len(dropped_rows)}")
    print(f"  Fixes:       {n_fixes}")
    print(f"  Warnings:    {n_warnings}")
    print(f"{'='*60}")

    if issues:
        print("\nDetails:")
        for issue in issues:
            tag = "  ⚠ " if "WARNING" in issue else "  ✓ "
            print(f"{tag}{issue}")

    # Platform / Layout summary
    layout_counts = Counter(
        row[col_idx["LibraryLayout"]] for row in fixed_rows)
    platform_counts = Counter(
        row[col_idx["Platform"]] for row in fixed_rows)

    print(f"\nLibraryLayout distribution:")
    for k, v in layout_counts.most_common():
        print(f"  {k:20s} {v}")

    print(f"\nPlatform distribution:")
    for k, v in platform_counts.most_common():
        print(f"  {k:20s} {v}")

    if not issues:
        print("\n  No issues found — metadata is clean.")

    print()


# -----------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate and fix iModulonMiner metadata TSV files.")
    parser.add_argument("input", help="Path to the input metadata TSV")
    parser.add_argument(
        "-o", "--output",
        help="Path for the fixed output TSV "
             "(default: <input>_fixed.tsv)")
    args = parser.parse_args()

    if args.output:
        out = args.output
    else:
        base, ext = os.path.splitext(args.input)
        out = f"{base}_fixed{ext}"

    fix_metadata(args.input, out)
