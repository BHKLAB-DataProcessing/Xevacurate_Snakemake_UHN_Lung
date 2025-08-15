#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from xevacurate_new.build_function import (
    build_model,
    build_drug,
    build_experiment,
    build_expDesign,
    build_modToBiobaseMap_pdata,
    save_csv,
)

def _norm_key(k: str) -> str:
    """Normalize pdata keys to stable filenames used in your repo."""
    kl = k.lower()
    if kl in {"rnaseq", "rna_seq", "rna-seq"}:
        return "RNASeq"
    if kl in {"cnv"}:
        return "CNV"
    if kl in {"mut", "mutation", "wes_mutation"}:
        return "mutation"
    # default: keep original, but strip spaces
    return k.replace(" ", "_")

def main():
    p = argparse.ArgumentParser(
        description="Build model/drug/experiment/expDesign/modToBiobaseMap and *_pdata.csv from updated batch TSV."
    )
    p.add_argument(
        "--in", dest="input_tsv", required=True,
        help="Path to results/update/update_model_drug_batch_data.tsv"
    )
    p.add_argument(
        "--out-dir", default="results/build",
        help="Directory to write outputs (default: results/build)"
    )
    p.add_argument(
        "--omics-dir", default="data/input/omics",
        help="Directory with omics inputs for *_pdata.csv (default: data/input/omics)"
    )
    args = p.parse_args()

    input_tsv = Path(args.input_tsv)
    out_dir   = Path(args.out_dir)
    omics_dir = Path(args.omics_dir)

    if not input_tsv.exists():
        p.error(f"--in not found: {input_tsv}")
    if not omics_dir.exists():
        p.error(f"--omics-dir not found: {omics_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Load the updated batch file
    df = pd.read_csv(input_tsv, sep="\t", low_memory=False)

    # 1) model.csv
    model_df = build_model(df)
    model_path = out_dir / "model.csv"
    save_csv(model_df, str(model_path))

    # 2) drug.csv
    drug_df = build_drug(df)
    drug_path = out_dir / "drug.csv"
    save_csv(drug_df, str(drug_path))

    # 3) experiment.csv
    experiment_df = build_experiment(df)
    experiment_path = out_dir / "experiment.csv"
    save_csv(experiment_df, str(experiment_path))

    # 4) expDesign.csv
    expdesign_df = build_expDesign(df)
    expdesign_path = out_dir / "expDesign.csv"
    save_csv(expdesign_df, str(expdesign_path))

    # 5) modToBiobaseMap.csv + *_pdata.csv
    mod2bio_df, pdata_dict = build_modToBiobaseMap_pdata(df, str(omics_dir))

    mod2bio_path = out_dir / "modToBiobaseMap.csv"
    save_csv(mod2bio_df, str(mod2bio_path))

    written = [model_path, drug_path, experiment_path, expdesign_path, mod2bio_path]

    for k, pdata in pdata_dict.items():
        key = _norm_key(str(k))
        out_path = out_dir / f"{key}_pdata.csv"
        save_csv(pdata, str(out_path))
        written.append(out_path)

    for pth in written:
        print(f"Wrote {pth}")

if __name__ == "__main__":
    main()
