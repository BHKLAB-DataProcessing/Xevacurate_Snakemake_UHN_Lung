#!/usr/bin/env python3
import argparse
from pathlib import Path
from xevacurate_new.update_function import (
    map_model,
    update_drug_data,
    update_batch_data,
    save_updated_tsv,
)

def main():
    p = argparse.ArgumentParser(
        description="Update model/drug/batch data from processed TSV and mapping files."
    )
    p.add_argument(
        "--processed-tsv", required=True,
        help="Path to results/process/processed_data.tsv"
    )
    p.add_argument(
        "--model-map-file", default="data/input/additional/all_model_clinical.csv",
        help="CSV with model metadata (default: data/input/additional/all_model_clinical.csv)"
    )
    p.add_argument(
        "--drug-map-file", default="data/input/additional/all_drug_mapping.xlsx",
        help="XLSX with drug mapping (default: data/input/additional/all_drug_mapping.xlsx)"
    )
    p.add_argument(
        "--out-dir", default="results/update",
        help="Directory to write outputs (default: results/update)"
    )
    # Output filenames (override if you like)
    p.add_argument("--model-out",    default="update_model_data.tsv")
    p.add_argument("--drug-out",     default="update_model_drug_data.tsv")
    p.add_argument("--filtered-out", default="drug_filtered_out.tsv")
    p.add_argument("--batch-out",    default="update_model_drug_batch_data.tsv")

    args = p.parse_args()

    processed   = Path(args.processed_tsv)
    model_map   = Path(args.model_map_file)
    drug_map    = Path(args.drug_map_file)
    out_dir     = Path(args.out_dir)

    if not processed.exists():
        p.error(f"--processed-tsv not found: {processed}")
    if not model_map.exists():
        p.error(f"--model-map-file not found: {model_map}")
    if not drug_map.exists():
        p.error(f"--drug-map-file not found: {drug_map}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) model mapping
    model_df = map_model(str(processed), str(model_map))
    model_out = out_dir / args.model_out
    save_updated_tsv(model_df, str(model_out))

    # 2) drug mapping (returns updated + filtered)
    drug_df, filtered_df = update_drug_data(str(model_out), str(drug_map))
    drug_out = out_dir / args.drug_out
    save_updated_tsv(drug_df, str(drug_out))

    filtered_out = out_dir / args.filtered_out
    save_updated_tsv(filtered_df, str(filtered_out))

    # 3) batch updates
    batch_df = update_batch_data(str(drug_out))
    batch_out = out_dir / args.batch_out
    save_updated_tsv(batch_df, str(batch_out))

    print(f"Wrote {model_out}")
    print(f"Wrote {drug_out}")
    print(f"Wrote {filtered_out}")
    print(f"Wrote {batch_out}")

if __name__ == "__main__":
    main()
