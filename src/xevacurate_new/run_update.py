import os
import argparse
from xevacurate_new.update_function import update_batch_data, save_updated_tsv, map_model, update_drug_data
from xevacurate_new.config import DATA_ADDITIONAL_DIR, RESULTS_DIR

# Ensure the output directory exists
OUTPUT_DIR = os.path.join(RESULTS_DIR, "update/")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# # Update model information
# PROCESSED_TSV_FILE = os.path.join(RESULTS_DIR, "process/processed_data.tsv")

# ## Step2.3: update_batch.py update batch information
# # input: update_model_data.tsv
# # output: update_model_drug_batch_data.tsv

# # Update batach information
# batch_updated_data = update_batch_data(PROCESSED_TSV_FILE)

# # Save updated data
# batch_updated_file = os.path.join(OUTPUT_DIR, "update_model_drug_batch_data.tsv")
# save_updated_tsv(batch_updated_data, batch_updated_file)

## Step2.1.1: update_model.py update modelID & sampleID
# input: processed_data.tsv
# output: update_model_data.tsv

# Update model information
PROCESSED_TSV_FILE = os.path.join(RESULTS_DIR, "process/processed_data.tsv")
MODEL_MAP_FILE = os.path.join(DATA_ADDITIONAL_DIR, "all_model_clinical.csv")
map_model_data = map_model(PROCESSED_TSV_FILE, MODEL_MAP_FILE)

# Save updated data
map_model_file = os.path.join(OUTPUT_DIR, "update_model_data.tsv")
save_updated_tsv(map_model_data, map_model_file)

## Step2.2: update_drug.py update drug information
# input: update_model_data.tsv
#       drug_dose_map.tsv
# output: update_model_data.tsv

# Update drug information
#DRUG_MAP_FILE = os.path.join(DATA_ADDITIONAL_DIR, "drug_dose_map.tsv")
DRUG_MAP_FILE = os.path.join(DATA_ADDITIONAL_DIR, "all_drug_mapping.xlsx")
drug_updated_data, filtered_data = update_drug_data(map_model_file, DRUG_MAP_FILE)

# Save updated data
drug_updated_file = os.path.join(OUTPUT_DIR, "update_model_drug_data.tsv")
save_updated_tsv(drug_updated_data, drug_updated_file)

filtered_file = os.path.join(OUTPUT_DIR, "drug_filtered_out.tsv")
save_updated_tsv(filtered_data, filtered_file)

## Step2.3: update_batch.py update batch information
# input: update_model_data.tsv
# output: update_model_drug_batch_data.tsv

# Update batach information
batch_updated_data = update_batch_data(drug_updated_file)

# Save updated data
batch_updated_file = os.path.join(OUTPUT_DIR, "update_model_drug_batch_data.tsv")
save_updated_tsv(batch_updated_data, batch_updated_file)

