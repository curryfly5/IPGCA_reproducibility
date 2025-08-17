library(tidyverse)
library(coloc)
library(data.table)
library(future)
library(furrr)
library(dplyr)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "Single gene colocalization analysis using GWAS and eQTL data")

# Required arguments
parser$add_argument("--query", required = TRUE,
                   help = "Gene ID to analyze")
parser$add_argument("--gwas", required = TRUE,
                   help = "Path to GWAS summary statistics file")
parser$add_argument("--eqtl", required = TRUE,
                   help = "Path to eQTL summary statistics file")
parser$add_argument("--tis", required = TRUE,
                   help = "Tissue name")
parser$add_argument("--trt", required = TRUE,
                   help = "Trait name")
parser$add_argument("--outdir", required = TRUE,
                   help = "Output directory")

# Optional arguments for eQTL column names
parser$add_argument("--eqtl_pval", default = "pval_nominal",
                   help = "Column name for eQTL p-values (default: pval_nominal)")
parser$add_argument("--eqtl_beta", default = "slope",
                   help = "Column name for eQTL effect sizes (default: slope)")
parser$add_argument("--eqtl_se", default = "slope_se",
                   help = "Column name for eQTL standard errors (default: slope_se)")
parser$add_argument("--eqtl_gene", default = "phenotype_id",
                   help = "Column name for gene IDs (default: phenotype_id)")
parser$add_argument("--eqtl_rsid", default = "rsid",
                   help = "Column name for SNP IDs (default: rsid)")
parser$add_argument("--eqtl_n", default = "ma_samples",
                   help = "Column name for sample sizes (default: ma_samples)")
parser$add_argument("--eqtl_maf", default = "af",
                   help = "Column name for allele frequencies (default: af)")

# Optional arguments for GWAS column names
parser$add_argument("--gwas_pval", default = "pvalue",
                   help = "Column name for GWAS p-values (default: pvalue)")
parser$add_argument("--gwas_beta", default = "effect_size",
                   help = "Column name for GWAS effect sizes (default: effect_size)")
parser$add_argument("--gwas_se", default = "standard_error",
                   help = "Column name for GWAS standard errors (default: standard_error)")
parser$add_argument("--gwas_chr", default = "chromosome",
                   help = "Column name for chromosome (default: chromosome)")
parser$add_argument("--gwas_pos", default = "position",
                   help = "Column name for position (default: position)")
parser$add_argument("--gwas_ea", default = "effect_allele",
                   help = "Column name for effect allele (default: effect_allele)")
parser$add_argument("--gwas_nea", default = "non_effect_allele",
                   help = "Column name for non-effect allele (default: non_effect_allele)")
parser$add_argument("--gwas_n", default = "sample_size",
                   help = "Column name for sample size (default: sample_size)")
parser$add_argument("--gwas_maf", default = "frequency",
                   help = "Column name for allele frequency (default: frequency)")
parser$add_argument("--gwas_varid", default = "variant_id",
                   help = "Column name for variant ID (default: variant_id)")

# Analysis parameters
parser$add_argument("--pp4_threshold", type = "double", default = 0.01,
                   help = "PP4 threshold for reporting SNPs (default: 0.01)")
parser$add_argument("--verbose", action = "store_true",
                   help = "Enable verbose output")

# Parse arguments
args <- parser$parse_args()

# Validate input files
if (!file.exists(args$gwas)) {
  stop("GWAS file not found: ", args$gwas)
}

if (!file.exists(args$eqtl)) {
  stop("eQTL file not found: ", args$eqtl)
}

# Create output directory
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

# Print parameters if verbose
if (args$verbose) {
  cat("=== Single Gene Colocalization Analysis ===\n")
  cat("Query gene:", args$query, "\n")
  cat("GWAS file:", args$gwas, "\n")
  cat("eQTL file:", args$eqtl, "\n")
  cat("Tissue:", args$tis, "\n")
  cat("Trait:", args$trt, "\n")
  cat("Output directory:", args$outdir, "\n")
  cat("PP4 threshold:", args$pp4_threshold, "\n")
  cat("===========================================\n")
}

# Function to read data safely
read_data_safely <- function(file_path, data_type, verbose = FALSE) {
  tryCatch({
    if (verbose) cat("Reading", data_type, "data from:", file_path, "\n")
    data <- fread(file_path)
    if (verbose) cat("Successfully read", nrow(data), "rows for", data_type, "\n")
    return(data)
  }, error = function(e) {
    stop("Error reading ", data_type, " file: ", e$message)
  })
}

# Function to check required columns
check_columns <- function(data, required_cols, data_type) {
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", data_type, " data: ", 
         paste(missing_cols, collapse = ", "), 
         "\nAvailable columns: ", paste(colnames(data), collapse = ", "))
  }
}

# Read GWAS and eQTL data
gwas <- read_data_safely(args$gwas, "GWAS", args$verbose)
eqtl <- read_data_safely(args$eqtl, "eQTL", args$verbose)

# Check required columns
gwas_required_cols <- c(args$gwas_chr, args$gwas_pos, args$gwas_ea, args$gwas_nea, 
                       args$gwas_pval, args$gwas_beta, args$gwas_se, args$gwas_n, 
                       args$gwas_maf, args$gwas_varid)
check_columns(gwas, gwas_required_cols, "GWAS")

eqtl_required_cols <- c(args$eqtl_rsid, args$eqtl_gene, args$eqtl_pval, args$eqtl_beta, 
                       args$eqtl_se, args$eqtl_n, args$eqtl_maf)
check_columns(eqtl, eqtl_required_cols, "eQTL")

# Standardize column names for processing
gwas_std <- gwas %>%
  rename(
    chromosome = !!args$gwas_chr,
    position = !!args$gwas_pos,
    effect_allele = !!args$gwas_ea,
    non_effect_allele = !!args$gwas_nea,
    pvalue = !!args$gwas_pval,
    effect_size = !!args$gwas_beta,
    standard_error = !!args$gwas_se,
    sample_size = !!args$gwas_n,
    frequency = !!args$gwas_maf,
    variant_id = !!args$gwas_varid
  ) %>%
  mutate(
    chromosome = as.numeric(gsub("chr", "", chromosome)),
    position = as.numeric(position),
    effect_size = as.numeric(effect_size),
    standard_error = as.numeric(standard_error),
    pvalue = as.numeric(pvalue),
    sample_size = as.numeric(sample_size),
    frequency = as.numeric(frequency)
  )

eqtl_std <- eqtl %>%
  rename(
    rsid = !!args$eqtl_rsid,
    phenotype_id = !!args$eqtl_gene,
    pval_nominal = !!args$eqtl_pval,
    slope = !!args$eqtl_beta,
    slope_se = !!args$eqtl_se,
    ma_samples = !!args$eqtl_n,
    af = !!args$eqtl_maf
  ) %>%
  mutate(
    slope = as.numeric(slope),
    slope_se = as.numeric(slope_se),
    pval_nominal = as.numeric(pval_nominal),
    ma_samples = as.numeric(ma_samples),
    af = as.numeric(af)
  )

# Check if query gene exists in eQTL data
if (!args$query %in% eqtl_std$phenotype_id) {
  stop("Query gene '", args$query, "' not found in eQTL data. Available genes: ", 
       paste(head(unique(eqtl_std$phenotype_id), 10), collapse = ", "))
}

if (args$verbose) {
  cat("Processing allele alignment...\n")
}

# Extract position information from eQTL rsid
pos_eqtl <- str_split(eqtl_std$rsid, "_", simplify = TRUE) %>%
  as.data.frame() %>%
  unique()

if (ncol(pos_eqtl) >= 4) {
  colnames(pos_eqtl) <- c("chromosome", "position", "ref", "alt")
  pos_eqtl <- pos_eqtl %>%
    mutate(
      chromosome = as.numeric(chromosome),
      position = as.numeric(position)
    )
  
  # Merge GWAS with position information
  gwas_merged <- merge(gwas_std, pos_eqtl, by = c("chromosome", "position"), all.x = TRUE)
  
  # Handle allele inconsistencies
  if ("ref" %in% colnames(gwas_merged)) {
    not_cons <- apply(gwas_merged, 1, function(x) {
      !is.na(x[["ref"]]) && !is.na(x[["effect_allele"]]) && x[["effect_allele"]] != x[["ref"]]
    })
    
    gwas_coloc <- gwas_merged
    
    if (sum(not_cons, na.rm = TRUE) > 0) {
      if (args$verbose) {
        cat("Flipping", sum(not_cons, na.rm = TRUE), "allele-inconsistent variants\n")
      }
      
      # Flip effect sizes for inconsistent variants
      inconsistent_idx <- which(not_cons)
      gwas_coloc$effect_size[inconsistent_idx] <- -gwas_coloc$effect_size[inconsistent_idx]
      
      # Swap effect and non-effect alleles
      temp_allele <- gwas_coloc$effect_allele[inconsistent_idx]
      gwas_coloc$effect_allele[inconsistent_idx] <- gwas_coloc$non_effect_allele[inconsistent_idx]
      gwas_coloc$non_effect_allele[inconsistent_idx] <- temp_allele
    }
  } else {
    gwas_coloc <- gwas_merged
  }
} else {
  if (args$verbose) {
    cat("Warning: Cannot parse allele information from eQTL rsid format\n")
  }
  gwas_coloc <- gwas_std
}

# Create SNP IDs for matching
snpid <- paste0(gwas_coloc$chromosome, "_", gwas_coloc$position, "_", 
                gwas_coloc$effect_allele, "_", gwas_coloc$non_effect_allele)
gwas_coloc$rsid <- snpid

# Filter eQTL data for the query gene
gene <- eqtl_std[eqtl_std$phenotype_id == args$query, ]

if (nrow(gene) == 0) {
  stop("No eQTL data found for gene: ", args$query)
}

if (args$verbose) {
  cat("Found", nrow(gene), "variants for gene", args$query, "\n")
}

# Prepare eQTL data for coloc
gene.coloc <- list(
  pvalues = gene$pval_nominal,
  beta = gene$slope,
  varbeta = gene$slope_se^2,
  snp = gene$rsid,
  type = "quant",
  N = unique(gene$ma_samples)[1],  # Use first unique value
  MAF = gene$af
)

# Find overlapping SNPs between GWAS and eQTL
gwas_coloc_sig <- gwas_coloc[gwas_coloc$rsid %in% gene$rsid, ]

if (nrow(gwas_coloc_sig) == 0) {
  stop("No overlapping SNPs found between GWAS and eQTL data for gene: ", args$query)
}

if (args$verbose) {
  cat("Found", nrow(gwas_coloc_sig), "overlapping SNPs\n")
}

# Prepare GWAS data for coloc
snpid2 <- paste0(gwas_coloc_sig$chromosome, "_", gwas_coloc_sig$position, "_", 
                 gwas_coloc_sig$effect_allele, "_", gwas_coloc_sig$non_effect_allele)

gwas.coloc <- list(
  pvalues = gwas_coloc_sig$pvalue,
  beta = gwas_coloc_sig$effect_size,
  varbeta = gwas_coloc_sig$standard_error^2,
  snp = snpid2,
  type = "quant",
  N = unique(gwas_coloc_sig$sample_size)[1],  # Use first unique value
  MAF = gwas_coloc_sig$frequency
)

# Run colocalization analysis
if (args$verbose) {
  cat("Running colocalization analysis...\n")
}

coloc_result <- coloc.abf(gene.coloc, gwas.coloc)

if (args$verbose) {
  cat("Colocalization summary:\n")
  cat("PP0 (no association):", round(coloc_result$summary[2], 4), "\n")
  cat("PP1 (eQTL only):", round(coloc_result$summary[3], 4), "\n")
  cat("PP2 (GWAS only):", round(coloc_result$summary[4], 4), "\n")
  cat("PP3 (both, different causal variants):", round(coloc_result$summary[5], 4), "\n")
  cat("PP4 (both, same causal variant):", round(coloc_result$summary[6], 4), "\n")
}

# Extract detailed results
coloc_result$results$rsid <- coloc_result$results$snp

# Merge with original p-values
gwasp <- gwas_coloc %>%
  select(pvalue, rsid) %>%
  rename(pval_gwas = pvalue)

# 修复：使用标准化后的列名
eqtlp <- eqtl_std %>%
  select(pval_nominal, rsid) %>%
  rename(pval_eqtl = pval_nominal)

# Merge p-values with colocalization results
coloc_detailed <- coloc_result$results %>%
  left_join(eqtlp, by = "rsid") %>%
  left_join(gwasp, by = "rsid") %>%
  mutate(
    Gene = args$query,
    Tissue = args$tis,
    Trait = args$trt
  ) %>%
  arrange(desc(SNP.PP.H4))

# Filter results based on PP4 threshold
significant_snps <- coloc_detailed %>%
  filter(SNP.PP.H4 > args$pp4_threshold)

if (args$verbose) {
  cat("Found", nrow(significant_snps), "SNPs with PP4 >", args$pp4_threshold, "\n")
}

# Write results
output_file <- file.path(args$outdir, paste0("ALLSNP.", args$query, ".", args$tis, ".coloc.res"))
fwrite(coloc_detailed, output_file, sep = "\t")

# Write summary results
summary_file <- file.path(args$outdir, paste0("SUMMARY.", args$query, ".", args$tis, ".coloc.res"))
summary_result <- data.frame(
  Tissue = args$tis,
  Gene = args$query,
  Trait = args$trt,
  nsnps = coloc_result$summary[1],
  PP0 = coloc_result$summary[2],
  PP1 = coloc_result$summary[3],
  PP2 = coloc_result$summary[4],
  PP3 = coloc_result$summary[5],
  PP4 = coloc_result$summary[6],
  n_significant_snps = nrow(significant_snps),
  top_snp = ifelse(nrow(significant_snps) > 0, significant_snps$rsid[1], NA),
  top_snp_pp4 = ifelse(nrow(significant_snps) > 0, significant_snps$SNP.PP.H4[1], NA)
)
fwrite(summary_result, summary_file, sep = "\t")

if (args$verbose) {
  cat("Results written to:\n")
  cat("- Detailed results:", output_file, "\n")
  cat("- Summary results:", summary_file, "\n")
  cat("Analysis completed successfully!\n")
}