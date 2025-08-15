#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Biobase)
  library(Xeva)
  library(dplyr)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--build-dir",  type = "character", default = "results/build",
              help = "Directory containing model.csv, drug.csv, experiment.csv, expDesign.csv, modToBiobaseMap.csv [default: %default]"),
  make_option("--omics-dir",  type = "character", default = "data/input/omics",
              help = "Directory containing omics matrices (e.g., data_expression_mRNA.txt, data_CNA.tsv, Tsao_Lung_WES_mutation_2025.tsv) [default: %default]"),
  make_option("--out-dir",    type = "character", default = "results/xevaset",
              help = "Output directory for the XevaSet RDS and helper files [default: %default]"),
  make_option("--name",       type = "character", default = "Tsao_Lung_2022",
              help = "XevaSet name [default: %default]"),
  make_option("--annotate-drugs", action = "store_true", default = TRUE,
              help = "Try to annotate drugs via AnnotationGx/PubChem [default: %default]"),
  make_option("--skip-annotate-drugs", action = "store_true", default = FALSE,
              help = "Force skip PubChem drug annotation (offline runs)"),
  make_option("--rna-file",   type = "character", default = "data_expression_mRNA.txt",
              help = "Filename (inside --omics-dir) for RNA expression matrix [default: %default]"),
  make_option("--wes-file",   type = "character", default = "Tsao_Lung_WES_mutation_2025.tsv",
              help = "Filename (inside --omics-dir) for mutation matrix [default: %default]"),
  make_option("--cnv-file",   type = "character", default = "data_CNA.tsv",
              help = "Filename (inside --omics-dir) for CNV matrix [default: %default]"),
  make_option("--rds-name",   type = "character", default = "UHN_Tsao_Lung_DrugResponse_2022_v1.rds",
              help = "Output RDS filename [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

BUILD_DIR <- normalizePath(opt$`build-dir`, mustWork = TRUE)
OMICS_DIR <- normalizePath(opt$`omics-dir`, mustWork = TRUE)
OUT_DIR   <- opt$`out-dir`
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.csv(path, stringsAsFactors = FALSE, ...)
}

read_table_matrix <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("Missing file: ", path)
  as.matrix(read.table(path,
                       header = TRUE, row.names = 1, sep = sep,
                       quote = "", comment.char = "", check.names = FALSE,
                       as.is = TRUE))
}

require_pkg <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

# BiomaRt mapping (shared by RNA / mutation / CNV)
feature_map_from_symbols <- function(symbols) {
  require_pkg("biomaRt")
  biomaRt <- asNamespace("biomaRt")
  ensembl <- biomaRt$useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
  gm <- biomaRt$getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "gene_biotype", "description"),
    filters    = "hgnc_symbol",
    values     = symbols,
    mart       = ensembl
  )
  gm_ordered <- gm[match(symbols, gm$hgnc_symbol), , drop = FALSE]
  rownames(gm_ordered) <- symbols
  gm_ordered
}

ensure_matrix_pheno_match <- function(mat, pheno_adf) {
  m_ids <- colnames(mat)
  p_ids <- rownames(pData(pheno_adf))
  if (!setequal(m_ids, p_ids)) {
    cat("Only in matrix:\n"); print(setdiff(m_ids, p_ids))
    cat("Only in pheno :\n"); print(setdiff(p_ids, m_ids))
    stop("Sample sets differ between matrix columns and phenoData rows.")
  }
  pheno_adf[m_ids, , drop = FALSE]
}

make_eset <- function(mat, pdata_csv, title, abstract) {
  pdata_df <- read_csv_safe(pdata_csv, row.names = 1)
  pdata_adf <- AnnotatedDataFrame(pdata_df)

  # feature annotation via biomaRt on HGNC symbols (rownames of mat)
  hugo <- rownames(mat)
  fdata_df <- feature_map_from_symbols(hugo)
  rownames(fdata_df) <- hugo
  fdata_adf <- AnnotatedDataFrame(fdata_df)

  # order pheno to matrix columns
  pdata_adf <- ensure_matrix_pheno_match(mat, pdata_adf)

  # minimal protocol data
  proto <- AnnotatedDataFrame(data.frame(row.names = colnames(mat)))

  expData <- new("MIAME",
                 name  = "Tsao Lab",
                 lab   = "Tsao Lab",
                 title = title,
                 abstract = abstract,
                 other = list(notes = "genome-build hg19; genome-version hg19"))

  ExpressionSet(assayData = mat,
                phenoData = pdata_adf,
                featureData = fdata_adf,
                protocolData = proto,
                experimentData = expData)
}

# PubChem drug annotation (optional / best-effort)
annotate_drugs <- function(drug_df) {
  if (opt$`skip-annotate-drugs`) return(drug_df)
  ok <- suppressWarnings(requireNamespace("AnnotationGx", quietly = TRUE))
  if (!ok) {
    message("AnnotationGx not available; skipping drug annotation.")
    return(drug_df)
  }
  AGx <- asNamespace("AnnotationGx")
  singles <- drug_df$drugname.standardized[drug_df$is.single == "yes"]
  singles <- singles[!is.na(singles)]
  if (!length(singles)) return(drug_df)

  cid_map <- try(AGx$mapCompound2CID(singles, first = TRUE), silent = TRUE)
  if (inherits(cid_map, "try-error")) {
    message("mapCompound2CID failed; skipping drug annotation.")
    return(drug_df)
  }
  cid_map <- cid_map[!is.na(cid_map$cids), , drop = FALSE]
  if (!nrow(cid_map)) return(drug_df)
  names(cid_map)[names(cid_map) == "cids"] <- "CID"

  props <- try(AGx$mapCID2Properties(
    ids = cid_map$CID,
    properties = c("MolecularFormula", "MolecularWeight", "IUPACName", "InChI", "InChIKey", "CanonicalSMILES")
  ), silent = TRUE)
  if (inherits(props, "try-error")) {
    message("mapCID2Properties failed; returning base drug table.")
    return(drug_df)
  }

  compound_info <- dplyr::left_join(cid_map, props, by = "CID")
  dplyr::left_join(drug_df, compound_info, by = c("drugname.standardized" = "name"))
}

# ---------- IO: build CSVs ----------
model_path <- file.path(BUILD_DIR, "model.csv")
drug_path  <- file.path(BUILD_DIR, "drug.csv")
exp_path   <- file.path(BUILD_DIR, "experiment.csv")
expd_path  <- file.path(BUILD_DIR, "expDesign.csv")
map_path   <- file.path(BUILD_DIR, "modToBiobaseMap.csv")

model       <- read_csv_safe(model_path)
drug        <- read_csv_safe(drug_path)
experiment  <- read_csv_safe(exp_path)
exp_df      <- read_csv_safe(expd_path)
modToBiobaseMap <- read_csv_safe(map_path)

# ---------- drug annotation (optional) ----------
drug_annotated <- tryCatch(
  annotate_drugs(drug),
  error = function(e) { message("Drug annotation error: ", conditionMessage(e)); drug }
)
drug_annotated_csv <- file.path(BUILD_DIR, "drug_annotated.csv")
write.csv(drug_annotated, file = drug_annotated_csv, row.names = TRUE)

# ---------- expDesign (list form) ----------
expDesign <- lapply(seq_len(nrow(exp_df)), function(i) {
  list(
    batch.name = exp_df$batch.name[i],
    treatment  = strsplit(exp_df$treatment[i], ",")[[1]],
    control    = strsplit(exp_df$control[i], ",")[[1]]
  )
})

# ---------- RNA ----------
rna_mat <- read_table_matrix(file.path(OMICS_DIR, opt$`rna-file`), sep = "\t")
rna_pdata_csv <- file.path(BUILD_DIR, "RNASeq_pdata.csv")
RNASeq <- make_eset(
  mat = rna_mat,
  pdata_csv = rna_pdata_csv,
  title = "Tsao Lab Lung Cancer RNA Expression",
  abstract = "HGNC symbols as gene IDs"
)
# write feature data for reference
rna_fdata_csv <- file.path(BUILD_DIR, "RNAseq_fdata.csv")
write.csv(Biobase::fData(RNASeq), file = rna_fdata_csv, row.names = TRUE)

# ---------- Mutation ----------
wes_mat <- read_table_matrix(file.path(OMICS_DIR, opt$`wes-file`), sep = "\t")
wes_pdata_csv <- file.path(BUILD_DIR, "mutation_pdata.csv")
mutation <- make_eset(
  mat = wes_mat,
  pdata_csv = wes_pdata_csv,
  title = "Tsao Lab Lung Cancer WES Mutation",
  abstract = "HGNC symbols as gene IDs"
)
mutation_fdata_csv <- file.path(BUILD_DIR, "mutation_fdata.csv")
write.csv(Biobase::fData(mutation), file = mutation_fdata_csv, row.names = TRUE)

# ---------- CNV ----------
cnv_mat <- read_table_matrix(file.path(OMICS_DIR, opt$`cnv-file`), sep = "\t")
cnv_pdata_csv <- file.path(BUILD_DIR, "CNV_pdata.csv")
cnv <- make_eset(
  mat = cnv_mat,
  pdata_csv = cnv_pdata_csv,
  title = "Tsao Lab Lung Cancer CNV",
  abstract = "HGNC symbols as gene IDs"
)
cnv_fdata_csv <- file.path(BUILD_DIR, "CNV_fdata.csv")
write.csv(Biobase::fData(cnv), file = cnv_fdata_csv, row.names = TRUE)

# ---------- experiment cleaning ----------
experiment <- as.data.frame(experiment, stringsAsFactors = FALSE)
# prefer drug.id when present/non-empty, otherwise fall back to 'drug'
if ("drug.id" %in% names(experiment)) {
  has_id <- !is.na(experiment$drug.id) & nzchar(experiment$drug.id)
  if (!"drug" %in% names(experiment)) experiment$drug <- NA_character_
  experiment$drug[has_id] <- experiment$drug.id[has_id]
  experiment$drug <- as.character(experiment$drug)
  experiment$drug.id <- NULL
} else if ("drug" %in% names(experiment)) {
  experiment$drug <- as.character(experiment$drug)
} else {
  stop("experiment.csv must have a 'drug' or 'drug.id' column.")
}

# ---------- create XevaSet ----------
xeva.set <- createXevaSet(
  name = opt$name,
  model = model,
  drug = drug_annotated,
  experiment = experiment,
  expDesign = expDesign,
  molecularProfiles = list(RNASeq = RNASeq, mutation = mutation, cnv = cnv),
  modToBiobaseMap = modToBiobaseMap
)

# ---------- set responses (best-effort) ----------
# helper that never stops the pipeline; logs warnings/errors and returns x unchanged on error
safe_setResponse <- function(x, measure) {
  withCallingHandlers(
    tryCatch(
      setResponse(x, res.measure = measure),
      error = function(e) {
        message("setResponse('", measure, "') error: ", conditionMessage(e))
        x
      }
    ),
    warning = function(w) {
      message("setResponse('", measure, "') warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

# explicit sequential calls
xeva.set <- safe_setResponse(xeva.set, "slope")    # no error
xeva.set <- safe_setResponse(xeva.set, "AUC")      # no error
xeva.set <- safe_setResponse(xeva.set, "angle")    # no error
xeva.set <- safe_setResponse(xeva.set, "abc")      # no error
xeva.set <- safe_setResponse(xeva.set, "mRECIST")  # >=50 warnings
xeva.set <- safe_setResponse(xeva.set, "TGI")      # error: replacement has length zero
xeva.set <- safe_setResponse(xeva.set, "lmm")      # error: NA/NaN/Inf in foreign function call (arg 1)

# ---------- save ----------
out_rds <- file.path(OUT_DIR, opt$`rds-name`)
saveRDS(xeva.set, file = out_rds)
cat("Wrote drug annotation data:\n",
    "  ", drug_annotated_csv, "\n")
cat("Wrote feature data:\n",
    "  ", rna_fdata_csv, "\n",
    "  ", mutation_fdata_csv, "\n",
    "  ", cnv_fdata_csv, "\n", sep = "")
cat("Wrote XevaSet to:", out_rds, "\n")
