import pandas as pd
import os
import json
from datetime import datetime

def process_file(file_entry: dict, base_dir: str) -> list:
    """
    Process a single Excel file based on its metadata.

    Args:
        file_entry (dict): A dictionary with keys 'file_name' and 'file_path'.
        base_dir (str): The base directory for relative path calculation.
        
    Returns:
        list: A list of DataFrames created from processing the file.
    """
    processed_dfs = []

    file_path = file_entry["file_path"]
    file_name = file_entry["file_name"]


    # Read the Excel file (assuming no header)
    df = pd.read_excel(file_path, header=None)
    #print(f"\nReading file: {file_path}")
    #print(df)
    
    # Extract metadata first
    metadata = {}
    mouse_row_idx = None
    
    for idx, row in df.iterrows():
        if pd.notna(row[0]):
            label = str(row[0]).strip()
            value = row[1] if pd.notna(row[1]) else None
            
            if "PHLC Sample" in label:
                metadata['phlc_sample'] = str(value).strip() if value else None
            elif "Passage Number" in label:
                metadata['passage'] = str(value).strip() if value else None
            elif "Start Date" in label:
                try:
                    metadata['start_date'] = pd.to_datetime(value).strftime('%Y-%m-%d') if pd.notna(value) else None
                    #print(f"Found Start Date: {metadata['start_date']}")
                except:
                    metadata['start_date'] = None
            elif "Start Treatment Date" in label:
                try:
                    metadata['start_treatment_date'] = pd.to_datetime(value).strftime('%Y-%m-%d') if pd.notna(value) else None
                    #print(f"Found Start Treatment Date: {metadata['start_treatment_date']}")
                except:
                    metadata['start_treatment_date'] = None
            elif "End Date" in label:
                try:
                    metadata['end_date'] = pd.to_datetime(value).strftime('%Y-%m-%d') if pd.notna(value) else None
                    #print(f"Found End Date: {metadata['end_date']}")
                except:
                    metadata['end_date'] = None
            elif "Mouse Number" in label:
                mouse_row_idx = idx
                break

    #print(f"Dates found: Start={metadata.get('start_date')}, Treatment={metadata.get('start_treatment_date')}, End={metadata.get('end_date')}")

    if mouse_row_idx is None:
        raise ValueError("Could not find mouse number row")

    # Get mouse numbers (from row after "Mouse Number")
    mouse_numbers = []
    mouse_number_row = df.iloc[mouse_row_idx]  # Get the row AFTER Mouse Number
    # Start from column 1 and process all columns
    for col in range(1, len(mouse_number_row)):
        val = mouse_number_row[col]
        if pd.notna(val) and str(val).strip() != '':
            try:
                mouse_num = int(float(str(val).strip()))
                mouse_numbers.append((col, mouse_num))  # Store column index with mouse number
            except (ValueError, TypeError):
                continue

    # print(f"Found mouse numbers: {[num for _, num in mouse_numbers]}")

    # Get treatments (from the correct row - 2 rows after Mouse Number)
    treatments = []
    treatment_row = df.iloc[mouse_row_idx + 2]  # Changed from mouse_row_idx + 1 to mouse_row_idx + 2
    for col_idx, _ in mouse_numbers:
        val = treatment_row[col_idx]
        if pd.notna(val) and str(val).strip() != '':
            treatments.append(str(val).strip())
        else:
            treatments.append('Control')

    # print(f"Found treatments: {treatments}")

    # Create ModelID from PHLC Sample
    if metadata.get('phlc_sample'):
        model_id = f"PHLC{metadata['phlc_sample']}"
    else:
        raise ValueError("Could not find PHLC Sample number")

    # Process measurements
    data_start_row = mouse_row_idx + 2

    for idx in range(data_start_row, df.shape[0]):
        day = df.iloc[idx, 0]
        if pd.notna(day) and isinstance(day, (int, float)):
            for i, (col_idx, mouse_number) in enumerate(mouse_numbers):
                measurement = df.iloc[idx, col_idx]
                if pd.notna(measurement) and isinstance(measurement, (int, float)):
                    sample_id = f"{model_id}_P{mouse_number}"
                    
                    new_df = pd.DataFrame([{
                        'file_name': file_name,
                        'file_checksum': file_entry.get('file_checksum', ''),
                        'file_last_modified_date': file_entry.get('file_last_modified_date', ''),
                        'ModelID': model_id,
                        'SampleID': sample_id,
                        'passage': metadata.get('passage'),
                        'tissue': 'lung',
                        'Drug': treatments[i] if i < len(treatments) else 'Control',
                        'start_date': metadata.get('start_date'),
                        'start_treatment_date': metadata.get('start_treatment_date'),
                        'end_date': metadata.get('end_date'),
                        'time': float(day),
                        'volume': float(measurement)
                    }])
                    processed_dfs.append(new_df)
                    #print(new_df)

    return processed_dfs

#######
def process_all_files(metadata_file: str, base_dir: str) -> (list, list):
    """
    Process all Excel files listed in a metadata JSON file.

    Args:
        metadata_file (str): Path to the metadata JSON file.
        base_dir (str): Base directory to resolve relative paths.
        
    Returns:
        tuple: (List of processed DataFrames, List of QC issues)
    """
    with open(metadata_file, "r") as f:
        file_list = json.load(f)
    
    all_data = []
    qc_issues = []
    
    for file_entry in file_list:
        file_path = os.path.join(base_dir, file_entry.get("file_path", ""))
        # file_name = file_entry.get("file_name","")
        tag = file_entry.get("tag", "")
        if tag:
            qc_issues.append({"file": file_path, "issue": tag})
        else:
            try:
                dfs = process_file(file_entry, base_dir)
                if dfs:
                    all_data.extend(dfs)
                    qc_issues.append({"file": file_path, "issue": "file processed successfully"})
                else:
                    qc_issues.append({"file": file_path, "issue": "processed dataframe is empty"})
            except Exception as e:
                qc_issues.append({"file": file_path, "issue": str(e)})

    return all_data, qc_issues

def save_processed_data(all_data: list, processed_tsv: str) -> None:
    """
    Save processed data to TSV and pickle files.

    Args:
        all_data (list): List of processed DataFrames.
        processed_tsv (str): Path to the TSV output file.
        processed_pkl (str): Path to the pickle output file.
    """
    if all_data:
        processed_data = pd.concat(all_data, ignore_index=True)
        # assumes `all_data` is a pandas DataFrame
        processed_data = processed_data.sort_values(
            by=["file_name", "SampleID"],
            ascending=[True, True],
            na_position="last",      # push NaNs to the end
            kind="mergesort"         # stable sort
        ).reset_index(drop=True)
        processed_data.to_csv(processed_tsv, index=False, sep="\t")
        print(f"Processed data saved to: {processed_tsv}")
    else:
        print("No data processed to save.")

def save_qc_report(qc_issues: list, qc_report_file: str) -> None:
    """
    Save quality control (QC) issues to a JSON file.

    Args:
        qc_issues (list): List of QC issues.
        qc_report_file (str): Path to the QC report file.
    """
    with open(qc_report_file, "w") as f:
        json.dump(qc_issues, f, indent=4)
    print(f"QC report saved to: {qc_report_file}")
