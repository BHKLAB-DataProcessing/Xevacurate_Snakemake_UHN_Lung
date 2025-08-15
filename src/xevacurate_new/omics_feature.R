# Load the biomaRt library
library(biomaRt)

# Connect to Ensembl v98
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 98)

# Your Ensembl gene IDs (with or without versions)
ensg_ids <- c("ENSG00000223972.5", "ENSG00000139618.13")

# Remove version suffix (if present)
ensg_ids_clean <- sub("\\..*", "", ensg_ids)

# Map to gene symbols and other identifiers
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "gene_biotype", "description"),
  filters = "ensembl_gene_id",
  values = ensg_ids_clean,
  mart = ensembl
)

# View results
print(gene_info)
