import os
from xevacurate_new.scan_function import collect_metadata, save_metadata
from xevacurate_new.config import DRUG_SCREEN_DIR, RESULTS_DIR, DATA_ADDITIONAL_DIR

# Step 1.1. scan_xlsx.py
# input: raw drug screening data directory
# output: XEVA_OFFICIAL_SHARING.json

# Collect metadata
remove_list_path = os.path.join(DATA_ADDITIONAL_DIR, "models_to_remove.xlsx")
remove_list = remove_list_path if os.path.exists(remove_list_path) else None
metadata = collect_metadata(DRUG_SCREEN_DIR, remove_list)

# Ensure the output directory exists
OUTPUT_DIR = os.path.join(RESULTS_DIR, "scan/")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Save metadata to the specified JSON file
output_json = os.path.join(OUTPUT_DIR, "all_file_scan.json")
output_excel = os.path.join(OUTPUT_DIR, "all_file_scan.xlsx")
save_metadata(metadata, output_json, output_excel)
