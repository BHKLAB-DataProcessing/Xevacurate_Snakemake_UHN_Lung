import pandas as pd

with open("drug_model.txt", "r") as f:
    model = [line.strip() for line in f]

print(model)

df = pd.read_csv("Tsao_Lung_WES_mutation_2025_old.tsv", sep="\t", index_col=0)

# Keep only columns that exist in both the model and the dataframe
valid_cols = [col for col in model if col in df.columns]
df = df[valid_cols]

df.to_csv("Tsao_Lung_WES_mutation_2025.tsv", sep="\t")
