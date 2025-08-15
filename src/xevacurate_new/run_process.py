import os
from xevacurate_new.process_function import process_all_files, save_processed_data, save_qc_report
from xevacurate_new.config import BASE_DIR, RESULTS_DIR

# Step 1.2 process_xlsx.py
# input: raw drug screening data directory & XEVA_OFFICIAL_SHARING.json
# output: processed_data.tsv
#         qc_report.json

# Process all files using the metadata
metadata_file = os.path.join(RESULTS_DIR, "scan/all_file_scan.json")
all_data, qc_issues = process_all_files(metadata_file, BASE_DIR)

# Ensure the output directory exists
OUTPUT_DIR = os.path.join(RESULTS_DIR, "process/")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Save processed data and QC report
processed_tsv_file = os.path.join(OUTPUT_DIR, "processed_data.tsv")
qc_report_file = os.path.join(OUTPUT_DIR, "qc_report.json")
save_processed_data(all_data, processed_tsv_file)
save_qc_report(qc_issues, qc_report_file)