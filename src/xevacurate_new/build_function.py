import pandas as pd
import re
import os

def build_model(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a model DataFrame from the updated dataframe.

    This function process dataframe, extracts unique 'model.id' values,
    assigns default values for tissue and tissue_name, and computes a 
    patient.id by splitting the model.id.

    Args:
        df (pd.DataFrame): dataframe with updated data.

    Returns:
        pd.DataFrame: A DataFrame with columns 'model.id', 'tissue', 
        'tissue_name', and 'patient.id'. Can be extracted from PDX-MI
    """
    # Extract relevant drug columns
    columns = ['model.id', 'drug', 'patientID', 'sex', 'oncotree_code',
                'primary_tumor_site', 'cancer_type_detailed', 'specimenID',
                'primary_recurrence_metastasis','TNM_stage', 'start_date',
                'specimen_collection_site', 'modelID', 'sample_class', 'passage',
                'host_strain_name', 'host_strain_nomenclature',	'engraftment_site',
                'engraftment_type', 'engraftment_material', 'engraftment_material_status',
                'pubmed_ID','model.type', 'modelID.parental']

    # Remove duplicated rows
    model_df = df[columns].drop_duplicates().reset_index(drop=True)

    model_df.columns = model_df.columns.str.replace("_", ".", regex=False)

    model_df = model_df.rename(columns={"patientID": "patient.id", 
                                        "start.date" : "engraftment.date",
                                        "pubmed.ID": "pubmed.id",
                                        "modelID.parental": "parental.model.id"})

    model_df =  model_df.drop_duplicates()
    
    return model_df

def build_drug(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a drug DataFrame from the updated dataframe.

    Args:
        df (pd.DataFrame): dataframe with updated data.

    Returns:
        pd.DataFrame: A DataFrame containing processed drug information.
    """
    # Extract relevant drug columns
    # columns = [
    #     'drug.id', 'drug.alternativename', 'dose', 'dose.unit', 'frequency', 
    #     'treatment.day', 'is.ADC', 'is.PDC', 'is.single', 'is.control', 'paired.treatment'
    # ]
    columns = [
        'drug.id', 'drugname.standardized', 'drug.alternativename', 'target', 'method', 'dose', 
        'dose.unit', 'treatment.type', 'drug.original', 'dose.original', 'frequency', 'is.single']
    drug_df = df[columns].copy()

    # Define a function to reduce multiple values
    def merge_values(series):
        unique_vals = series.dropna().unique()
        if len(unique_vals) == 1:
            return unique_vals[0]
        else:
            return ';'.join(sorted(set(map(str, unique_vals))))

    # Remove duplicates and sort by 'drug.id'
    drug_df = drug_df.sort_values(by="drug.id").drop_duplicates().reset_index(drop=True)

    # Group by 'drug.id' and aggregate all other columns
    drug_df_merged = drug_df.groupby("drug.id", dropna=False).agg(merge_values).reset_index()
    
    return drug_df_merged

def build_experiment(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a experiment DataFrame from the updated dataframe.

    Args:
        df (pd.DataFrame): dataframe with updated data.

    Returns:
        pd.DataFrame: A DataFrame containing processed drug information.
    """
    # Extract relevant drug columns
    columns = [
        'model.id', 'drug', 'drug.id', 'dose', 'time', 'volume',
        'modelID', 'model.type', 'start_date', "start_treatment_date", "end_date",
        'file_name'
    ]

    # Extract relevant columns for experiments
    exp_df = df[columns].copy()

    exp_df.columns = exp_df.columns.str.replace("_", ".", regex=False)

#    exp_df =  exp_df.drop_duplicates()
    
    return exp_df

def add_control_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Used in genrating expDesign.csv
    Adds a 'control' column to the DataFrame based on the treatment value from the control batch.
    
    The control batch is identified by taking the prefix (the part before the period in 'batch.name')
    and appending '.H2O'. The corresponding treatment value from the control batch is then assigned
    as the control for the row.
    
    Args:
        df (pd.DataFrame): DataFrame containing at least the columns 'batch.name' and 'treatment'.
        
    Returns:
        pd.DataFrame: DataFrame with the added 'control' column.
    """
    # Extract prefix from 'batch.name'
    df['prefix'] = df['batch.name'].str.split('.').str[0]
    
    # Create a control key by appending '.H2O'
    df['control_key'] = df['prefix'] + '.H2O'
    
    # Create a mapping from 'batch.name' to 'treatment'
    control_mapping = df.set_index('batch.name')['treatment'].to_dict()
    
    # Generate 'control' column using the control_mapping dictionary
    df['control'] = df['control_key'].map(control_mapping)
    
    # Drop rows where 'batch.name' ends with '.H2O' (these are control batches themselves)
    df = df[~df['batch.name'].str.endswith('.H2O')]
    
    # Drop unnecessary columns
    df = df.drop(columns=['prefix', 'control_key'])
    
    return df


def build_expDesign(df: pd.DataFrame) -> pd.DataFrame:
    """
    Builds the experiment design DataFrame from the updated data file.
    
    Args:
        df (pd.DataFrame): dataframe with updated data.
    
    Returns:
        pd.DataFrame: Processed DataFrame ready for saving.
    """
    # Extract relevant drug columns
    columns = ['file_name','batch', 'model.id']

    df = df[columns].copy()

    df = (df.drop_duplicates()
                    .sort_values(by=['file_name', 'batch', 'model.id'])
                    .reset_index(drop=True))
    
    #print(df)

    # Aggregate 'model.id' for each 'batch'
    expDesign_df = (df.groupby(['file_name','batch'])['model.id']
                    .apply(lambda x: ','.join(map(str, x)))
                    .reset_index())
    
    control_df = expDesign_df[expDesign_df['batch'].str.contains("control")]
    control_df = control_df.rename(columns={"model.id": "control", "batch":"batch.control"})
    exp_df = expDesign_df[~expDesign_df['batch'].str.contains("control")]
    exp_df = exp_df.rename(columns={"model.id": "treatment"})

    #print(expDesign_df)

    expDesign_df = exp_df.merge(control_df, on="file_name", how="left")
    col = ['batch', 'treatment', 'control']
    expDesign_df = expDesign_df[col]
    expDesign_df = expDesign_df.rename(columns={"batch":"batch.name"})

    # # Rename columns for clarity
    # expDesign_df.rename(columns={'batch': 'batch.name', 'model.id': 'treatment'}, inplace=True)
    
    # # Add control column
    # expDesign_df = add_control_column(expDesign_df)
    
    return expDesign_df
    #return expDesign_df


def build_modToBiobaseMap_pdata(df: pd.DataFrame, omics_data: str) -> pd.DataFrame:
    """
    Build a modToBiobaseMap DataFrame by loading updated dataframe and omics data directory, 
    processing them, and merging on 'biobase.id'.
    
    Args:
        df (pd.DataFrame): dataframe with updated data.
        omics_data_file (str): Directory containing omics data.
        
    Returns:
        pd.DataFrame: Merged DataFrame containing 'model.id', 'biobase.id', and 'mDataType'.
    """
    ## 1. Process updated data file
    # Extract relevant drug columns
    columns = ['model.id', 'modelID']

    # Extract relevant columns for experiments
    df = df[columns].copy()

    # Process dataframe to extract 'model.id' and 'biobase.id'
    df = (df.sort_values(by='modelID')
            .drop_duplicates()
            .reset_index(drop=True)
            .rename(columns={'modelID': 'biobase.id'}))
    
    # Add "S" prefix if 'biobase.id' starts with a digit
    # potentially specific for TNBC dataset
    # probably should reformat RNAseq matrix sample ID, and skip this step in the future
    df['biobase.id'] = df['biobase.id'].str.replace("PHLC00", "PHLC").str.replace("PHLC0", "PHLC").str.replace(".TPX", "")
    #print(df)
    ## 2. Process Omics Data
    # List all files in the directory
    omics_df = pd.DataFrame(columns=['biobase.id', 'mDataType']) # for modToBiobaseMap
    pdata_dict = {} # for pdata

    files = os.listdir(omics_data)

    for file in files:
        file_path = os.path.join(omics_data, file)
        pdata_df = pd.DataFrame()

        # Process RNA-related file
        if 'rna' in file.lower():
            type = 'RNASeq'
            try:
                # Load RNAseq data
                RNAseq_df = pd.read_csv(file_path, sep="\t")
                
                # Prepare the DataFrame for merging
                rnaseq_model_df = pd.DataFrame({'mDataType': type}, index=RNAseq_df.columns[1:])
                rnaseq_model_df = (rnaseq_model_df.reset_index()
                                .rename(columns={'index': 'biobase.id'})
                                .sort_values(by='biobase.id'))
                
                # Concatenate the two DataFrames
                omics_df = pd.concat([omics_df, rnaseq_model_df], axis=0, ignore_index=True) 

                # Extract pdata info
                pdata_df['biobase.id'] = rnaseq_model_df['biobase.id']

            except Exception as e:
                print(f"Failed to process RNAseq file '{file}': {e}")

        elif any(keyword in file.lower() for keyword in ['wes', 'exome', 'snp', 'mutation', 'indel', 'snv']):
            type = 'mutation'
            try:
                # Load mutation data
                mutation_df = pd.read_csv(file_path, sep='\t')
                
                # Prepare the DataFrame for merging
                mutation_model_df = pd.DataFrame({'mDataType': type}, index=mutation_df.columns[1:])
                mutation_model_df = (mutation_model_df.reset_index()
                                .rename(columns={'index': 'biobase.id'})
                                .sort_values(by='biobase.id'))
                
                # Concatenate the two DataFrames
                omics_df = pd.concat([omics_df, mutation_model_df], axis=0, ignore_index=True)

                # Extract pdata info
                pdata_df['biobase.id'] = mutation_model_df['biobase.id']

            except Exception as e:
                print(f"Failed to process mutation file '{file}': {e}")

        elif any(keyword in file.lower() for keyword in ['cnv', 'cna', 'copynumber', 'copy_number', 'seg']):
            type = 'CNV'
            try:
                # Load mutation data
                cnv_df = pd.read_csv(file_path, sep='\t')
                
                # Prepare the DataFrame for merging
                cnv_model_df = pd.DataFrame({'mDataType': type}, index=cnv_df.columns[1:])
                cnv_model_df = (cnv_model_df.reset_index()
                                .rename(columns={'index': 'biobase.id'})
                                .sort_values(by='biobase.id'))
                
                # Concatenate the two DataFrames
                omics_df = pd.concat([omics_df, cnv_model_df], axis=0, ignore_index=True)

                # Extract pdata info
                pdata_df['biobase.id'] = cnv_model_df['biobase.id']

            except Exception as e:
                print(f"Failed to process CNV file '{file}': {e}")

        elif any(keyword in file.lower() for keyword in ['atac', 'accessibility']):
            type = 'ATACseq'
            try:
                # Load mutation data
                atac_df = pd.read_csv(file_path, sep='\t')
                
                # Prepare the DataFrame for merging
                atac_model_df = pd.DataFrame({'mDataType': type}, index=atac_df.columns[1:])
                atac_model_df = (atac_model_df.reset_index()
                                .rename(columns={'index': 'biobase.id'})
                                .sort_values(by='biobase.id'))
                
                # Concatenate the two DataFrames
                omics_df = pd.concat([omics_df, atac_model_df], axis=0, ignore_index=True)

                # Extract pdata info
                pdata_df['biobase.id'] = atac_model_df['biobase.id']

            except Exception as e:
                print(f"Failed to process ATACseq file '{file}': {e}")

        else:
            print(f"""Input '{file}' cannot be categorized into:
                    - RNASeq ("rna")
                    - Mutation
                    - CNV
                    - ATACseq
                    based on name.""")
            
        pdata_df['patient.id'] = pdata_df['biobase.id']
        pdata_df['tissue'] = 'Lung' # specific to TNBC
        pdata_df['tissue_name'] = 'Lung Cancer' # specific to TNBC

        pdata_dict[type] = pdata_df
    
    ## 3. Merge the two DataFrames on 'biobase.id'
    merged_df = df.merge(omics_df, how="inner", on="biobase.id")
    #types_omics = merged_df['mDataType'].unique()
    
    return merged_df, pdata_dict


def build_pdata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a pdata based on loading updated dataframe and omics data.
    
    Args:
        df (pd.DataFrame): dataframe with updated data.
        rna_data_file (str): Path to the RNAseq data file (TSV format).
        
    Returns:
        pd.DataFrame: Merged DataFrame containing 'biobase.id', 'patient.id', 'tissue', and 'tissue.name'.
    """
    # Extract all types of omics data
    omics_list = sorted(df['mDataType'].unique())

    # Populate the other columns
    # specific to TNBC
    df['tissue'] = 'BRCA'
    df['tissue_name'] = 'Breast Cancer'
    
    # Compute patient.id by splitting 'model.id':
    #   - Split on period ('.') then split the first part on '_P\d+' to remove passage info.
    df['patient.id'] = df['model.id'].apply(
        lambda x: re.split(r'_P\d+', x.split('.')[0])[0]
    )

    return df


def save_csv(df: pd.DataFrame, updated_file: str) -> None:
    """
    Save the updated DataFrame to a CSV file using a comma as the separator.

    Args:
        df (pd.DataFrame): The DataFrame to be saved.
        updated_file (str): The path to the output CSV file.
    """
    df.to_csv(updated_file, sep=",", index=False)
    print(f"Updated data saved to: {updated_file}")
