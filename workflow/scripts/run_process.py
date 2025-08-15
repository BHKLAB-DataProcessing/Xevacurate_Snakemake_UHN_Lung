#!/usr/bin/env python3
import argparse
from pathlib import Path
from xevacurate_new.process_function import (
    process_all_files,
    save_processed_data,
    save_qc_report,
)

def main():
    p = argparse.ArgumentParser(
        description="Process scanned Excel files into a unified TSV and QC report."
    )
    p.add_argument(
        "--manifest",
        required=True,
        help="Path to results/scan/all_file_scan.json (the scan manifest).",
    )
    p.add_argument(
        "--base-dir",
        required=True,
        help="Root directory used to resolve relative file paths from the manifest "
             "(e.g., data/input/drug_screen).",
    )
    p.add_argument(
        "--output-dir",
        default="results/process",
        help="Directory to write outputs (default: results/process).",
    )
    p.add_argument(
        "--tsv-name",
        default="processed_data.tsv",
        help="Output TSV filename (default: processed_data.tsv).",
    )
    p.add_argument(
        "--qc-name",
        default="qc_report.json",
        help="QC report filename (default: qc_report.json).",
    )

    args = p.parse_args()

    manifest = Path(args.manifest)
    base_dir = Path(args.base_dir)
    out_dir = Path(args.output_dir)

    if not manifest.exists():
        p.error(f"--manifest not found: {manifest}")
    if not base_dir.exists():
        p.error(f"--base-dir not found: {base_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Process all files
    all_data, qc_issues = process_all_files(str(manifest), str(base_dir))

    # Save outputs
    processed_tsv = out_dir / args.tsv_name
    qc_report = out_dir / args.qc_name

    save_processed_data(all_data, str(processed_tsv))
    save_qc_report(qc_issues, str(qc_report))

    print(f"Wrote {processed_tsv}")
    print(f"Wrote {qc_report}")

if __name__ == "__main__":
    main()
