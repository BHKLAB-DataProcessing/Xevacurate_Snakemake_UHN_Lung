import pandas as pd
import os
import re

def update_model_id(samplename: str, filename: str):
    """
    Extracts ModelID, SampleID, and Model_Type from samplename or filename.

    Tries to parse `samplename` first. If that fails, falls back to `filename`.

    Args:
        samplename (str): Clean model name string (e.g., 'MODEL_ID P2' or 'MODEL_ID P2 DRUG RES P3')
        filename (str): Original filename of the Excel file

    Returns:
        tuple: (model_id, sample_id, model_type)

    Raises:
        ValueError: If both samplename and filename fail to parse
    """

    def _extract_from_parts(parts, label):
        passage_numbers = [p for p in parts if re.match(r"^P\d+$", p)]

        if "RES" in parts:
            model_type = "Resistant"
            if len(passage_numbers) > 1: # new format would all follow this type, with 2 Pn
                first_passage = passage_numbers[0]
                resistant_passage = passage_numbers[-1]
                try:
                    index_second_last = parts.index(passage_numbers[-2])
                    model_part = parts[:index_second_last]
                    resistant_drug = parts[index_second_last + 1]
                except Exception as e:
                    raise ValueError(f"Error extracting resistant model info from {label}: {e}")

                model_id = "_".join(model_part).replace(" ", "").upper() + f"_{first_passage}_{resistant_drug}_RES"
                sample_id = f"{model_id}_{resistant_passage}"

            elif len(passage_numbers) == 1: # older format, some does not contain the first Pn, only contains the last Pn
                resistant_passage = passage_numbers[0]
                # remove the first occurrence of that passage number
                parts.remove(resistant_passage)
                index_second_last = parts.index('RES')
                resistant_drug = parts[index_second_last - 1]
                model_part = parts[:index_second_last - 1]

                if model_part and resistant_drug: # both fields are non-empty
                    model_id = "_".join(model_part).replace(" ", "").upper() + f"_{resistant_drug}_RES"
                elif resistant_drug: # only one filed is non-empty, that will be model_id, empty drug always means 945
                    model_id = "_".join(resistant_drug) + "_945" + "_RES"
                else:
                    model_id = ""
                sample_id = f"{model_id}_{resistant_passage}"
            
            else: # older format, some does not contain any Pn
                if parts[-1] != "RES":
                    parts.append(parts.pop(parts.index("RES")))
                model_id = "_".join(parts).replace(" ", "").upper()
                sample_id = model_id + "_Pn"

        else:
            model_type = "Parental"
            if not passage_numbers:
                model_id = "_".join(parts).replace(" ", "").upper()
                sample_id = model_id + "_Pn"
            else:
                model_part = parts[:parts.index(passage_numbers[-1])]
                model_id = "_".join(model_part).replace(" ", "").upper()
                sample_id = f"{model_id}_{passage_numbers[-1]}"

        return model_id, sample_id, model_type

    # First try extract information from samplename
    samplename = samplename.upper()
    _pattern1 = re.compile(r'P[XM]\s*\+?\s*(\d+)', re.IGNORECASE)
    samplename = _pattern1.sub(r'P\1', samplename)
    samplename = (samplename.replace("#2- COHORT", "")
                  .replace("- OCT 12/16", "")
                  .replace(" (+METGEL)", "")
                  .replace("(+METGEL FPI)", "")
                  .replace("TNBC", "")
                  .replace("AXILLA", "")
                  .replace("(66820)", "")
                  .replace("PNA", "")
                  .replace("P6-", "P6 ")
                  .replace("BPTO.", "BPTO")
                  .replace("BXTO.", "BXTO")
                  .replace("PARENTAL", "")
                  .replace("INS B 0", "INSB0")
                  .replace("INSB014 ST", "INSB014")
                  .replace("INSB019 CT2", "INSB019")
                  .replace("NOTCH ", "NOTCH")
                  .replace("NOTCHB02P5-", "NOTCHB02 P5 ")
                  .replace("REF-OV-", "REFOV")
                  .replace("REF OV-", "REFOV")
                  .replace("REF S ", "REFS")
                  .replace("REF OV-", "REFOV")
                  .replace("REF OV ", "REFOV")
                  .replace("REV OV ", "REFOV ")
                  .replace("REF B ", "REFB")
                  .replace("REF ", "REF")
                  .replace("B BL", "")
                  .replace("BXBL", "")
                  .replace("BX BL", "")
                  .replace("BL", "")
                  .replace("REF06", "REF006")
                  .replace("(JULY 12, 2016)", "")
                  .replace("(JULY 12/16)", "")
                  .replace("BT-", "BT")
                  .replace("BVR-0-", "BVR0")
                  .replace("BVR-M-", "BVRM")
                  .replace("REF010 P3 945", "REF010 P3 945 RES")
                  .replace("RES 945/BMN/TAX", "945/BMN/TAX RES")
                  .replace("RES 945/BMN", "945/BMN RES")
                  .replace("RES 945+ABT80", "945/ABT80 RES")
                  .replace("RES 945LD", "945LD RES")
                  .replace("RES 945W", "945W RES")
                  .replace("RES 945", "945 RES")
                  .replace("945 WEEKLY RES", "945W RES")
                  .replace("945 WKLY RES", "945W RES")
                  .replace("RES CDX","CDX RES")
                  .replace("ERIBULIN", "ERIB")
                  .replace("RES ERIB LD", "ERIBLD RES")
                  .replace("RES ERIB", "ERIB RES")
                  .replace("ERIB LD RES", "ERIBLD RES")
                  .replace("ERIB LD RE", "ERIBLD RES")
                  .replace("ERIB ELD", "ERIBELD")
                  .replace("RES ERIB RES", "ERIB RES")
                  .replace("RES TAX", "TAX RES")
                  .replace("RES BMN RES", "BMN RES")
                  .replace(" TAXOL", " TAX")
                  .replace("/TAXOL", "TAX")
                  .replace("RES TTK/TAX", "TTK/TAX RES")
                  .replace("TAX 15MG", "TAXLD")
                  .replace("TAX 15", "TAXLD")
                  .replace("CARBORES", "CARBO RES")
                  .replace("NOTCH06TAXRESPM+6","NOTCH06 TAX RES PM+6")
                  .replace("DOCETAXOLRES", "DOCE RES")
                  .replace("DOCETAXOL", "DOCE")
                  .replace("DOCETAXEL", "DOCE")
                  .replace("RES CX", "CX RES")
                  .replace("RES CARBO", "CARBO RES")
                  .replace("RES BMN", "BMN RES")
                  .replace("RED010", "REF010")
                  .replace("P2- NOV 3/16", "P2")
                  .replace("P2- NOV 9/16", "P2")
                  .replace("RES TTK8/TAX20", "TTK8/TAX20 RES")
                  .replace("TTK+TAX", "TTK/TAX")
                  .replace("REF001CARBORESP3", "REF001 CARBO RES P3")
                  .replace("RES TAX15", "TAXLD RES")
                  .replace(" TAX15 ", " TAXLD ")
                  .replace("RES TTK", "TTK RES")
                )
    parts = re.split(r"\s+", samplename.strip().upper())
    (model_id, sample_id, model_type) = _extract_from_parts(parts, f"samplename '{samplename}'")
    
    #return _extract_from_parts(parts, f"samplename '{samplename}'")

    # Second, try extract from file name
    name_wo_ext = os.path.splitext(filename)[0].upper()
    # older files names usually start with "COPY OF", remove it before extract information
    if name_wo_ext.startswith("XCOPY OF "):
        name_wo_ext = name_wo_ext[len("XCOPY OF "):]
    if name_wo_ext.startswith("COPY OF "):
        name_wo_ext = name_wo_ext[len("COPY OF "):]
    name_wo_ext = (name_wo_ext.replace("REF B ", "REFB")
                       .replace("REF ", "REF")
                       .replace("REF S ", "REFS")
                       .replace("BL", "")
                       .replace("TTK+TAX", "TTK/TAX")
                       .replace("TNBC", "")
                       .replace("(66820)", "")
                    )
        #parts = name_wo_ext.strip().split(" ")
    parts = re.split(r"\s+", name_wo_ext.strip())
    (fmodel_id, fsample_id, fmodel_type) = _extract_from_parts(parts, f"samplename '{samplename}'")

    # return the longer names extracted from filename or samplename, 
    if len(samplename) < 4 or "00:00:00" in samplename:
        return (fmodel_id, fsample_id, fmodel_type)
    elif model_type == "Parental":
        return (model_id, sample_id, model_type)
    elif model_type != fmodel_type: # model_type Resistant; fmodel_type Parental
        return (model_id, sample_id, model_type)
    elif len(fsample_id) > len(sample_id):
        return (fmodel_id, fsample_id, fmodel_type)
    else:
        return (model_id, sample_id, model_type)

def dose_number_unit_split(dose_series: pd.Series) -> pd.DataFrame:
    """
    Splits a pandas Series of dose strings into two columns: numeric and unit.

    Args:
        dose_series (pd.Series): Series of dose strings (e.g., "10mg/kg", "500mcg", "8mg/kg+20mcg").
    
    Returns:
        pd.DataFrame: A DataFrame with two columns: 'dose_numeric' (as float) and 'dose_unit' (str).
    """
    
    # Function to process each dose string
    def split_composite_dose(dose_str):
        if pd.isna(dose_str):
            return '', ''
        # Find all number+unit groups
        parts = re.findall(r'([0-9]*\.?[0-9]+)\s*([a-zA-Z/µ%]+)', dose_str)
        if not parts:
            return '', ''
        doses = '+'.join(p[0] for p in parts)
        units = '+'.join(p[1] for p in parts)
        return doses, units

    # Apply to the series
    extracted = dose_series.apply(split_composite_dose)
    result = pd.DataFrame(extracted.tolist(), columns=['dose', 'dose.unit'])
    
    return result

def assign_mouse_ids(group: pd.DataFrame) -> pd.DataFrame:
    """
    Assigns mouse IDs for a group of rows based on the 'time' column.
    A new mouse ID is assigned whenever a row with time == 0 is encountered.

    Args:
        group (pd.DataFrame): A DataFrame group with the same 'batch' value.
    
    Returns:
        pd.DataFrame: The updated group with a new column 'model.id'.
    """
    counter = 0
    current_mouse_id = None
    mouse_ids = []
    time = 0
    
    # Iterate through rows in original order
    for _, row in group.iterrows():
        if row["time"] < time:
            counter += 1
            current_mouse_id = f"{row['batch']}.m{counter}"
            time = row["time"]
        mouse_ids.append(current_mouse_id)
    
    group["model.id"] = mouse_ids
    return group

def update_model_data(raw_data: str) -> pd.DataFrame:
    """
    Reads raw model data from a TSV file, updates the ModelID and SampleID.

    Args:
        raw_data (str): Path to the input TSV file with raw data.
    
    Returns:
        pd.DataFrame: The updated DataFrame.
    """
    # Read the raw data from the TSV file
    df = pd.read_csv(raw_data, sep="\t", low_memory=False)

    # Apply the update_model_id function to the SampleID and filename columns
    # Apply the function row-wise and expand into three new columns
    df[["modelID", "sample_id.update", "modelType"]] = df.apply(
        lambda row: pd.Series(update_model_id(row["sampleID"], row["file.name"])),
        axis=1
    )

    # Rename columns to preserve original values and use updated ones
    df = df.rename(columns={
        'sampleID': 'sampleID.original',
        'sample_id.update': 'sampleID'
    })

    # Apply the update_strain_location_source to the mouse_strain, implantation_location, implantation_tissue_source column
    df['mouse.strain'] = (
        df['mouse.strain']
            .str.removeprefix("Strain: ")
            .where(df['mouse.strain'].str.startswith("Strain:"), "")
    )

    df['implantation.location'] = (
        df['implantation.location']
            .str.removeprefix("Location: ")
            .where(df['implantation.location'].str.startswith("Location:"), "")
    )

    allowed = ["Fresh Passage", "Cell Injection", "Frozen Viables"]
    df["implantation.tissue.source"] = df["implantation.tissue.source"].where(
        df["implantation.tissue.source"].isin(allowed),
        ""
    )
    return df

def map_model(raw_data: str, map_file: str) -> pd.DataFrame:
    """
    Reads raw data from a TSV file, and drug data from a tsv file, updates drug information.

    Args:
        raw_data (str): Path to the input TSV file with raw data.
        drug_map (str): Path to the input drug mapping file.

    Returns:
        pd.DataFrame: The updated DataFrame.
    """
    # Read the raw data from the TSV file
    raw_df = pd.read_csv(raw_data, sep="\t", low_memory=False)
    # drug_df = pd.read_csv(drug_map, sep="\t")
    model_df = pd.read_csv(map_file)

    merge = raw_df.merge(model_df, on="ModelID", how="left")

    merge = merge.rename(columns={"Parent PMLB_xenograft": "modelID",
                                  "Parent PMLB_specimen": "specimenID",
                                  "Parent PMLB_case": "patientID",
                                  "ModelID": "modelID.original"})
    merge = merge.drop(columns=["stem_PPID","alias.id"])

    merge["model.type"] =  "parental"

    merge["modelID.parental"] = merge["modelID"]

    return merge

def update_drug_data(raw_data: str, drug_map: str) -> pd.DataFrame:
    """
    Reads raw data from a TSV file, and drug data from a tsv file, updates drug information.

    Args:
        raw_data (str): Path to the input TSV file with raw data.
        drug_map (str): Path to the input drug mapping file.

    Returns:
        pd.DataFrame: The updated DataFrame.
    """
    # Read the raw data from the TSV file
    raw_df = pd.read_csv(raw_data, sep="\t", low_memory=False)
    raw_df = raw_df.rename(columns={"Drug": "drug"})
    raw_df = raw_df.astype({'drug': 'string'})
    # drug_df = pd.read_csv(drug_map, sep="\t")
    drug_df = pd.read_excel(drug_map, header=0)
    #drug_df = drug_df.drop(columns={'comment_Dave','comment_Mitchell'})
    drug_df = drug_df.astype({'drug': 'string'})

    update_df = drug_df.merge(raw_df, on='drug', how='right')
    update_df = update_df.rename(columns={'treatment.standardized': 'drug.id'})
    update_df['drug.alternativename'] = update_df['drug'].where(update_df['drug'] != update_df['drug.id'], "")
    update_df['drug.original'] = update_df['drug']
    update_df['drug'] = update_df['drug.id']
    update_df['dose.original'] = update_df['dose']
    update_df = update_df[update_df['keep'] == "yes"]
    update_df = update_df.drop(columns=['keep'])

    # keep the filtered out rows in a sperate dataframe
    # 2. Anti-join: rows in raw_df that are NOT in drug_df
    drug_df = drug_df[drug_df["keep"] == "yes"]
    raw_unmatched = raw_df.merge(drug_df, on=['drug'], how='left', indicator=True)
    raw_unmatched = raw_unmatched[raw_unmatched['_merge'] == 'left_only'].drop(columns=['_merge'])

    return update_df, raw_unmatched

def update_batch_data(raw_file: str) -> pd.DataFrame:
    """
    Processes the raw data by adding a batch column and assigning mouse IDs.
    
    Args:
        raw_file (str): Path to the input TSV file with raw data.
    
    Returns:
        pd.DataFrame: The updated DataFrame.
    """
    # Read raw data
    df = pd.read_csv(raw_file, sep="\t", low_memory=False)
    
    # Create a 'batch' column by concatenating 'ModelID' and 'drug_id'
    df["drug.id"] = df["drug"]
    df["file_number"] = "f" + (df.groupby("file_name").ngroup() + 1).astype(str) # need file name in batch infor to distinguis
    df["batch"] = df["file_number"] + "." + df["SampleID"].str.split("_").str[0]+ ".TPX." + df["drug.id"].astype(str)
    df["model.id"] = df["file_number"] + "." + df["SampleID"].str.split("_").str[0]+ ".TPX." + df["drug.id"].astype(str) + "." + df["SampleID"].str.split("_").str[1].astype(str)
    df = df.rename(columns={"Drug": "drug"})
    
    return df


def save_updated_tsv(df: pd.DataFrame, updated_file: str) -> None:
    """
    Save the updated DataFrame to a CSV file using a comma as the separator.

    Args:
        df (pd.DataFrame): The DataFrame to be saved.
        updated_file (str): The path to the output CSV file.
    """
    df.to_csv(updated_file, sep="\t", index=False)
    print(f"Updated data saved to: {updated_file}")

