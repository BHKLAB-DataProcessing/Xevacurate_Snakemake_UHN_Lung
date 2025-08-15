from xevacurate_new.scan_function import collect_metadata, save_metadata
from xevacurate_new.process_function import process_all_files, save_processed_data, save_qc_report
from xevacurate_new.update_function import process_model_data, process_drug_data, process_batch_data, save_updated_tsv
from xevacurate_new.build_function import build_model, build_drug, build_experiment, build_expDesign, build_modToBiobaseMap_pdata, build_pdata, save_csv #generate_drug, generate_experiment, generate_expDesign, generate_modToBiobaseMap, 
import os
import pandas as pd

BASE_DIR = "/Users/guanqiaofeng/Documents/BHK/Xeva/xevacurate_new"
DRUG_SCREEN_DIR = "/Users/guanqiaofeng/Documents/BHK/Xeva/xevacurate_new/data/input/drug_screen/Completed_Jan2024-Dec2024"
OUTPUT_DIR = "/Users/guanqiaofeng/Documents/BHK/Xeva/xevacurate_new/data/out"

# Step 1.1. scan_xlsx.py
# input: raw drug screening data directory
# output: XEVA_OFFICIAL_SHARING.json

metadata = collect_metadata(DRUG_SCREEN_DIR, BASE_DIR)
save_metadata(metadata, os.path.join(OUTPUT_DIR, "XEVA_OFFICIAL_SHARING.json"))

# Step 1.2 process_xlsx.py
# input: raw drug screening data directory & XEVA_OFFICIAL_SHARING.json
# output: processed_data.tsv
#         qc_report.json
METADATA_FILE = os.path.join(OUTPUT_DIR, "XEVA_OFFICIAL_SHARING.json")
PROCESSED_TSV_FILE = os.path.join(OUTPUT_DIR, "processed_data.tsv")
QC_REPORT_FILE = os.path.join(OUTPUT_DIR, "qc_report.json")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)
    
# Process all files using the metadata
all_data, qc_issues = process_all_files(METADATA_FILE, BASE_DIR)
    
# Save processed data and QC report
save_processed_data(all_data, PROCESSED_TSV_FILE)
save_qc_report(qc_issues, QC_REPORT_FILE)

# Step2.1: update_model.py update modelID & sampleID
# input: processed_data.tsv
# output: update_model_data.tsv
updated_df = process_model_data(PROCESSED_TSV_FILE)
updated_data = os.path.join(OUTPUT_DIR, "update_model_data.tsv")
save_updated_tsv(updated_df, updated_data)

# Step2.2: update_drug.py update drug information
# input: update_model_data.tsv
#       drug_dose_map.tsv
# output: update_model_data.tsv
DRUG_MAP_FILE = os.path.join(BASE_DIR, "data/input/additional/drug_dose_map.tsv")
updated_df = process_drug_data(updated_data, DRUG_MAP_FILE)
updated_data = os.path.join(OUTPUT_DIR, "update_model_drug_data.tsv")
save_updated_tsv(updated_df, updated_data)

# Step2.3: update_batch.py update batch information
# input: update_model_data.tsv
# output: update_model_drug_batch_data.tsv
updated_df = process_batch_data(updated_data)
updated_data = os.path.join(OUTPUT_DIR, "update_model_drug_batch_data.tsv")
save_updated_tsv(updated_df, updated_data)

# STEP3: generate_inputs.py Generate Input files for XevaSet Creation 
# input: update_model_drug_batch_data.tsv
df = pd.read_csv(updated_data, sep="\t")

# Step3.1: generate model.csv
model_data = os.path.join(OUTPUT_DIR, "model.csv")
model_df = build_model(df)
save_csv(model_df, model_data)

# Step3.2: generate drug.csv
drug_data = os.path.join(OUTPUT_DIR, "drug.csv")
drug_df = build_drug(df)
save_csv(drug_df, drug_data)

# Step3.3: generate experiment.csv
experiment_data = os.path.join(OUTPUT_DIR, "experiment.csv")
experiment_df = build_experiment(df)
save_csv(experiment_df, experiment_data)

# Step3.4:generate expDesign.csv
expDesign_data = os.path.join(OUTPUT_DIR, "expDesign.csv")
expDesign_df = build_expDesign(df)
save_csv(expDesign_df, expDesign_data)

# Step3.5: generate modToBiobaseMap.csv and {WES|RNAseq|CNV|ATACseq}_pdata.csv
omics_data = os.path.join(BASE_DIR, "data/input/omics") # input omics data directory
modToBiobaseMap_data = os.path.join(OUTPUT_DIR, "modToBiobaseMap.csv")
modToBiobaseMap_df, pdata_dict = build_modToBiobaseMap_pdata(df, omics_data)
save_csv(modToBiobaseMap_df, modToBiobaseMap_data)
for type, pdata_df in pdata_dict.items():
    pdata = os.path.join(OUTPUT_DIR, f"{type}_pdata.csv")
    save_csv(pdata_df, pdata)
