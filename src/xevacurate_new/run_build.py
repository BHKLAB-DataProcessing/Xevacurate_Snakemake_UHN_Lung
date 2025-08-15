import os
import pandas as pd
from xevacurate_new.build_function import build_model, build_drug, build_experiment, build_expDesign, build_modToBiobaseMap_pdata, save_csv
from xevacurate_new.config import OMICS_DIR, RESULTS_DIR

# Ensure the output directory exists
OUTPUT_DIR = os.path.join(RESULTS_DIR, "build/")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# STEP3: generate_inputs.py Generate Input files for XevaSet Creation 
# input: update_model_drug_batch_data.tsv
updated_file = os.path.join(RESULTS_DIR, "update/update_model_drug_batch_data.tsv")
df = pd.read_csv(updated_file, sep="\t", low_memory=False)

# Step3.1: generate model.csv
model_data = build_model(df)
model_file = os.path.join(OUTPUT_DIR, "model.csv")
save_csv(model_data, model_file)

# Step3.2: generate drug.csv
drug_data = build_drug(df)
drug_file = os.path.join(OUTPUT_DIR, "drug.csv")
save_csv(drug_data, drug_file)

# Step3.3: generate experiment.csv
experiment_data = build_experiment(df)
experiment_file = os.path.join(OUTPUT_DIR, "experiment.csv")
save_csv(experiment_data, experiment_file)

# Step3.4:generate expDesign.csv
expDesign_data = build_expDesign(df)
expDesign_file = os.path.join(OUTPUT_DIR, "expDesign.csv")
save_csv(expDesign_data, expDesign_file)

# Step3.5: generate modToBiobaseMap.csv and {WES|RNAseq|CNV|ATACseq}_pdata.csv
modToBiobaseMap_data, pdata_dict = build_modToBiobaseMap_pdata(df, OMICS_DIR)

modToBiobaseMap_file = os.path.join(OUTPUT_DIR, "modToBiobaseMap.csv")
save_csv(modToBiobaseMap_data, modToBiobaseMap_file)

for type, pdata in pdata_dict.items():
    pdata_file = os.path.join(OUTPUT_DIR, f"{type}_pdata.csv")
    save_csv(pdata, pdata_file)
