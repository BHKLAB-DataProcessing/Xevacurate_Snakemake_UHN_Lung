# Xevacurate Tsao 2022 — PDX Drug-Response Curation

This repository curates [Tsao Lab](https://wwwlabs.uhnresearch.ca/labs/tsao/) lung cancer drug-screen data into a reproducible **XevaSet** object using a Snakemake workflow, Python utilities, and an R step executed inside Docker (Bioconductor 3.20).

The pipeline scans & validates raw inputs, processes and updates annotations, builds Xeva input CSVs, calculate drug response statistics, and finally creates an `.rds` XevaSet.

---

## Table of Contents

- [What You Get](#what-you-get)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Inputs Expected](#inputs-expected)
- [Configuration](#configuration)
- [Workflow Outline](#workflow-outline)
- [Outputs](#outputs)
- [Load the XevaSet in R](#load-the-xevaset-in-r)
- [Troubleshooting](#troubleshooting)
- [Repo Layout](#repo-layout)
- [Citation](#citation)
- [License](#license)

---

## What You Get

- A fully reproducible Snakemake workflow orchestrated via **pixi** (Python env manager).
- Clean CSV inputs for Xeva: `model.csv`, `drug.csv`, `experiment.csv`, `expDesign.csv`, `modToBiobaseMap.csv`, and per-omics `*_pdata.csv`.
- An **XevaSet** saved as `results/xevaset/<NAME>.rds` (name configurable).
- Logs at each step under `logs/`.

---

## Requirements

- **Docker** (Desktop) : <https://www.docker.com/>
- **pixi** (to provide Snakemake + Python deps): <https://pixi.sh/>
- **Git** (to clone this repo)
- Internet access during the R step (biomaRt Ensembl gene mapping & PubChem calls if enabled).

> **Windows users:** Docker Desktop with WSL2 backend works fine with this workflow.

---

## Quick Start

```bash
# 1) Clone
git clone git@github.com:BHKLAB-DataProcessing/Xevacurate_Snakemake_UHN_Lung.git
cd Xevacurate_Snakemake_UHN_Lung

# 2) Dry run
PYTHONPATH=src pixi run snakemake -n

# 3) Execute (adjust cores to taste)
PYTHONPATH=src pixi run snakemake --cores 4 --printshellcmds
```

The final XevaSet appears at:
```
results/xevaset/UHN_Tsao_Lung_DrugResponse_2022_v1.rds   # default, configurable
```

## Inputs Expected
Place the following under data/input/:
```
data/
└── input/
    ├── drug_screen/                # raw drug-response inputs (time/volume etc.)
    ├── additional/
    │   ├── models_to_remove.xlsx
    │   ├── all_model_clinical.csv
    │   └── all_drug_mapping.xlsx
    └── omics/
        ├── data_expression_mRNA.txt   # RNA TPM/CPM matrix (rows: genes, cols: samples)
        ├── Tsao_Lung_WES_mutation_2025.tsv  # WES mutation matrix
        └── data_CNA.tsv               # CNV matrix
```
**Notes**

 - Gene IDs in omics matrices are expected as **HGNC symbols** (rows). The R script will annotate via **biomaRt**.

 - Sample (column) IDs must match what `*_pdata.csv` files (generated during build) will contain. The R script ensures ordering & sanity checks.

## Configuration
All knobs live in `config.yaml`. Minimal example:
```
docker:
  xeva_r_image: docker.io/gfeng2023/xeva-build-r:bioc-3.20  # public image with Xeva deps

xevaset:
  name: Tsao_Lung_2022
  rds_name: UHN_Tsao_Lung_DrugResponse_2022_v1.rds
```

## Workflow Outline

1. **scan** → `workflow/scripts/run_scan.py`  
   - Scans `data/input/drug_screen` and writes a manifest + QC info.  
   - **Output:** `results/scan/all_file_scan.json`

2. **process** → `workflow/scripts/run_process.py`  
   - Validates/normalizes raw data and emits a clean TSV + QC.  
   - **Outputs:** `results/process/processed_data.tsv`, `results/process/qc_report.json`

3. **update** → `workflow/scripts/run_update.py`  
   - Maps models/drugs, applies removal lists, updates batches.  
   - **Outputs:**  
     - `results/update/update_model_data.tsv`  
     - `results/update/update_model_drug_data.tsv`  
     - `results/update/drug_filtered_out.tsv`  
     - `results/update/update_model_drug_batch_data.tsv`

4. **build** → `workflow/scripts/run_build.py`  
   - Produces Xeva input CSVs:  
     - `results/build/model.csv`  
     - `results/build/drug.csv`  
     - `results/build/experiment.csv`  
     - `results/build/expDesign.csv`  
     - `results/build/modToBiobaseMap.csv`  
     - `results/build/CNV_pdata.csv`  
     - `results/build/mutation_pdata.csv`  
     - `results/build/RNASeq_pdata.csv`  
   - Plus feature annotations created by the R step:  
     - `results/build/rna_fdata.csv`  
     - `results/build/mutation_fdata.csv`  
     - `results/build/cnv_fdata.csv`

5. **create_xevaset (Dockerized R)** → `workflow/scripts/create_xevaset.R`  
   - Runs inside the configured Docker image (Bioconductor 3.20).  
   - Builds `ExpressionSet`s, merges with experiment & design, computes responses, and writes the final xevaset as `.rds` file.

## Outputs
```
results/
├── scan/
│   └── all_file_scan.json
├── process/
│   ├── processed_data.tsv
│   └── qc_report.json
├── update/
│   ├── update_model_data.tsv
│   ├── update_model_drug_data.tsv
│   ├── drug_filtered_out.tsv
│   └── update_model_drug_batch_data.tsv
├── build/
│   ├── model.csv
│   ├── drug.csv
│   ├── experiment.csv
│   ├── expDesign.csv
│   ├── modToBiobaseMap.csv
│   ├── CNV_pdata.csv
│   ├── mutation_pdata.csv
│   ├── RNASeq_pdata.csv
│   ├── rna_fdata.csv
│   ├── mutation_fdata.csv
│   └── cnv_fdata.csv
└── xevaset/
    └── UHN_Tsao_Lung_DrugResponse_2022_v1.rds   # configurable
```

## Load the XevaSet in R
```
# From the repo root or any R session with Xeva installed:
xset <- readRDS("results/xevaset/UHN_Tsao_Lung_DrugResponse_2022_v1.rds")

# Inspect
xset
responseType(xset)       # Which response measures were computed (e.g., slope, AUC, angle, abc, mRECIST)
modelInfo(xset)[1:5,]
head(summarizeResponse(xset, by="model"))
```

## Troubleshooting

- `docker: command not found`

  Install **Docker Desktop** (and enable **WSL2** on Windows).

  Restart your shell/terminal.

- **Bioconductor/CRAN packages fail to install**
  
  The default image `gfeng2023/xeva-build-r:bioc-3.20` already contains the required R/Bioc packages (e.g., **Xeva**, **Biobase**, **biomaRt**).

  If you build your own image, use **Bioconductor 3.20** (or newer) and avoid packages unavailable in older Bioc releases (e.g., `AnnotationGx` is not on 3.19).

- **Omics sample mismatches**
  
  The R script performs sanity checks. Ensure **column names** of omics matrices match the `*_pdata.csv` **IDs** produced in the build step.

- **Windows path quirks**
  
  Always run from **inside WSL2 (Ubuntu)** and keep the repo on the WSL filesystem (e.g., under `~/`).

  Bind-mounting Windows paths can be slower and occasionally tricky.

## Repo Layout
```
.
├── config.yaml
├── src/                      # Python package(s): xevacurate_new/...
├── workflow/
│   ├── rules/                # Snakemake rules (*.smk)
│   │   └── curation.smk
│   ├── scripts/
│   │   ├── run_scan.py
│   │   ├── run_process.py
│   │   ├── run_update.py
│   │   ├── run_build.py
│   │   └── create_xevaset.R
│   └── Snakefile            # includes rules/curation.smk
├── data/
│   └── input/               # place your inputs here (see Inputs Expected)
├── results/                 # pipeline outputs
└── logs/                    # stepwise logs
```

## Citation

If you use this workflow or the produced **XevaSet** in a publication, please cite:

- **Xeva (Bioconductor)**
  
  *Mer A, Haibe-Kains B (2025). Xeva: Analysis of patient-derived xenograft (PDX) data. [doi:10.18129/B9.bioc.Xeva](https://doi.org/10.18129/B9.bioc.Xeva).*

  *R package version 1.24.0. https://bioconductor.org/packages/Xeva.*

- **The Original Dataset/Source**

  *Weiss J, Pham NA, Pintilie M, Li M, Liu G et al (2022) Optimizing Drug Response Study Design in Patient-Derived Tumor Xenografts. [doi.org/10.1177/11769351221136056](https://doi.org/10.1177/11769351221136056).*


## License

Add your license here (e.g., **MIT**).  
If you include third-party data, ensure you have the rights to redistribute it and include any required attribution or data-use terms.


## Contributions

Contributions are welcome!  
