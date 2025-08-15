#!/usr/bin/env python3
import argparse
from pathlib import Path
from xevacurate_new.scan_function import collect_metadata, save_metadata

def main():
    parser = argparse.ArgumentParser(
        description="Scan drug screening Excel files and collect metadata."
    )
    parser.add_argument(
        "--drug-screen-dir",
        required=True,
        help="Path to raw drug screening data directory"
    )
    parser.add_argument(
        "--output-dir",
        default="results/scan",
        help="Directory to write outputs (default: results/scan)"
    )
    parser.add_argument(
        "--remove-list",
        required=False,
        help="Optional path to Excel file with models to remove"
    )
    args = parser.parse_args()

    raw_dir = Path(args.drug_screen_dir)
    out_dir = Path(args.output_dir)
    remove_list = Path(args.remove_list) if args.remove_list else None

    if not raw_dir.exists():
        parser.error(f"--drug-screen-dir not found: {raw_dir}")

    if remove_list and not remove_list.exists():
        print(f"Warning: --remove-list not found: {remove_list}; proceeding without it.")
        remove_list = None

    out_dir.mkdir(parents=True, exist_ok=True)

    # Collect metadata
    metadata = collect_metadata(str(raw_dir), str(remove_list) if remove_list else None)

    # Save JSON
    out_json = out_dir / "all_file_scan.json"
    save_metadata(metadata, str(out_json))
    print(f"Wrote {out_json}")

if __name__ == "__main__":
    main()
