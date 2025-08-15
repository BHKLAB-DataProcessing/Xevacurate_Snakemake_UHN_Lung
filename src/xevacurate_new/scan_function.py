import os
import re
import json
import hashlib
import logging
import pandas as pd
from datetime import datetime
from typing import List, Dict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_checksum(file_path: str) -> str:
    """
    Calculate the SHA256 checksum of a file.
    
    Args:
        file_path (str): The path to the file.
        
    Returns:
        str: The SHA256 checksum as a hexadecimal string.
    """
    sha256_hash = hashlib.sha256()
    try:
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
    except Exception as e:
        logging.error(f"Failed to read file {file_path}: {e}")
        return ""
    
    return sha256_hash.hexdigest()

def collect_metadata(data_dir: str, remove_file: str) -> List[Dict[str, str]]:
    """
    Traverse a directory and collect metadata of all .xlsx files.
    
    Args:
        data_dir (str): The directory to scan for Excel files.
        base_dir (str): The base directory for relative path calculation.
        
    Returns:
        List[Dict[str, str]]: A list of metadata dictionaries.
    """
    if remove_file is not None:
        remove = pd.read_excel(remove_file)
        remove_list = remove[remove['leave.out']=='yes']['file.name'].tolist()
    else:
        remove_list = []

    seen_cores = set() # use this to verify when a file shows up the 2+ times
    metadata_list = []
    for root, _, files in os.walk(data_dir):
        folder_name = os.path.basename(root)
        for file in files:
            if file.endswith(".xlsx") and not file.startswith("~$"): 
                file_path = os.path.join(root, file)
                relative_file_path = os.path.relpath(file_path, data_dir)
                file_core = re.sub(r"(?i)^x?copy of ", "", file)  # (?i) = case-insensitive, ^ = start of string
                
                if "res" in file_core.lower():
                    model_type = "resistant"
                else:
                    model_type = "parental"
                
                tags = []
                if file_core in seen_cores:
                    tags.append("duplicated_file")
                if "mouse" in file_core.lower():
                    tags.append("mouse")
                if file in remove_list:
                    tags.append("in_remove_list")
                tag = ", ".join(tags)

                seen_cores.add(file_core)
                
                try:
                    modified_time = os.path.getmtime(file_path)
                    modified_date = datetime.utcfromtimestamp(modified_time).strftime('%Y-%m-%d %H:%M:%S')
                    checksum = calculate_checksum(file_path)
                    file_size_bytes = os.path.getsize(file_path)
                    file_size_kb = round(file_size_bytes / 1024, 2)
                    metadata_list.append({
                        "folder_name": folder_name,
                        "file_name": file,
                        "file_name_core": file_core,
                        "file_path": file_path,
                        "model_type": model_type,
                        "file_size_kb": file_size_kb,
                        "last_modified_date": modified_date,
                        "checksum": checksum,
                        "tag":tag
                    })
                except Exception as e:
                    logging.error(f"Error processing file {file_path}: {e}")
    
    return metadata_list

def save_metadata(metadata_list: List[Dict[str, str]], output_json: str) -> None:
    """
    Save the metadata to a JSON file.
    
    Args:
        metadata_list (List[Dict[str, str]]): The list of metadata to save.
        output_json (str): The path to the JSON output file.
    """
    try:
        with open(output_json, "w") as json_file:
            json.dump(metadata_list, json_file, indent=4)
        logging.info(f"Metadata Json saved: {output_json}")
                
        # # Save as Excel
        # df = pd.DataFrame(metadata_list)
        # df.to_excel(output_excel, index=False)
        # logging.info(f"Metadata Excel saved: {output_excel}")

    except Exception as e:
        logging.error(f"Failed to save JSON file: {e}")
