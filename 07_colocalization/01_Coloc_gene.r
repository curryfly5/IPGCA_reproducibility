library(tidyverse)
library(coloc)
library(data.table)
library(future)
library(furrr)
library(dplyr)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "Colocalization analysis using GWAS and eQTL data")

# Required arguments
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
parser$add_argument("--workers", type = "integer", default = 8,
                   help = "Number of parallel workers (default: 8)")
parser$add_argument("--min_overlap", type = "integer", default = 10,
                   help = "Minimum number of overlapping SNPs (default: 10)")
parser$add_argument("--chunk_size", type = "integer", default = 100,
                   help = "Number of genes to process in each chunk (default: 100)")
parser$add_argument("--verbose", action = "store_true",
                   help = "Enable verbose output")

# Parse arguments
args <- parser$parse_args()

# Optimize worker count
args$workers <- min(args$workers, parallel::detectCores() - 1, 8)

# Print parameters if verbose
if (args$verbose) {
  cat("=== Colocalization Analysis Parameters ===\n")
  cat("GWAS file:", args$gwas, "\n")
  cat("eQTL file:", args$eqtl, "\n")
  cat("Tissue:", args$tis, "\n")
  cat("Trait:", args$trt, "\n")
  cat("Output directory:", args$outdir, "\n")
  cat("Workers:", args$workers, "\n")
  cat("Minimum overlap:", args$min_overlap, "\n")
  cat("Chunk size:", args$chunk_size, "\n")
  cat("==========================================\n")
}

# Create output directory
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

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

# Function to standardize column names
standardize_gwas_cols <- function(gwas_data, args) {
  gwas_data %>%
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
    ) %>%
    filter(!is.na(chromosome), !is.na(position), !is.na(pvalue), 
           !is.na(effect_size), !is.na(standard_error),
           pvalue > 0, standard_error > 0, sample_size > 0, 
           frequency > 0, frequency < 1)
}

standardize_eqtl_cols <- function(eqtl_data, args) {
  eqtl_data %>%
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
    ) %>%
    filter(!is.na(pval_nominal), !is.na(slope), !is.na(slope_se), 
           !is.na(ma_samples), !is.na(af),
           pval_nominal > 0, slope_se > 0, ma_samples > 0, 
           af > 0, af < 1)
}

# =============================================================================
# MAIN DATA LOADING AND PREPROCESSING (DONE ONCE)
# =============================================================================

if (args$verbose) cat("Loading and preprocessing data...\n")

# Read data ONCE
gwas <- read_data_safely(args$gwas, "GWAS", args$verbose)
eqtl <- read_data_safely(args$eqtl, "eQTL", args$verbose)

# Check required columns
gwas_required_cols <- c(args$gwas_chr, args$gwas_pos, args$gwas_ea, args$gwas_nea, 
                       args$gwas_pval, args$gwas_beta, args$gwas_se, args$gwas_n, args$gwas_maf)
check_columns(gwas, gwas_required_cols, "GWAS")

eqtl_required_cols <- c(args$eqtl_rsid, args$eqtl_gene, args$eqtl_pval, args$eqtl_beta, 
                       args$eqtl_se, args$eqtl_n, args$eqtl_maf)
check_columns(eqtl, eqtl_required_cols, "eQTL")

# Standardize columns ONCE
gwas <- standardize_gwas_cols(gwas, args)
eqtl <- standardize_eqtl_cols(eqtl, args)

if (args$verbose) {
  cat("After filtering - GWAS:", nrow(gwas), "variants\n")
  cat("After filtering - eQTL:", nrow(eqtl), "variants\n")
}

# Process allele alignment ONCE
if (args$verbose) cat("Processing allele alignment...\n")

tryCatch({
  pos_eqtl <- str_split(eqtl$rsid, "_", simplify = TRUE) %>%
    as.data.frame() %>%
    unique()
  
  if (ncol(pos_eqtl) >= 4) {
    colnames(pos_eqtl) <- c("chromosome", "position", "ref", "alt")
    pos_eqtl <- pos_eqtl %>%
      mutate(
        chromosome = as.numeric(chromosome),
        position = as.numeric(position)
      ) %>%
      filter(!is.na(chromosome), !is.na(position))
    
    # Convert to data.table for efficient merging
    gwas <- as.data.table(gwas)
    pos_eqtl <- as.data.table(pos_eqtl)
    
    # Ensure matching data types
    gwas[, chromosome := as.numeric(chromosome)]
    gwas[, position := as.numeric(position)]
    pos_eqtl[, chromosome := as.numeric(chromosome)]
    pos_eqtl[, position := as.numeric(position)]
    
    # Merge and align alleles
    gwas <- merge(gwas, pos_eqtl, by = c("chromosome", "position"), all.x = TRUE)
    
    # Handle allele flipping
    if ("ref" %in% colnames(gwas)) {
      inconsistent_idx <- which(!is.na(gwas$ref) & !is.na(gwas$effect_allele) & 
                               gwas$effect_allele != gwas$ref)
      
      if (length(inconsistent_idx) > 0) {
        if (args$verbose) cat("Flipping", length(inconsistent_idx), "inconsistent alleles\n")
        
        # Flip effect sizes and swap alleles
        gwas$effect_size[inconsistent_idx] <- -gwas$effect_size[inconsistent_idx]
        temp_allele <- gwas$effect_allele[inconsistent_idx]
        gwas$effect_allele[inconsistent_idx] <- gwas$non_effect_allele[inconsistent_idx]
        gwas$non_effect_allele[inconsistent_idx] <- temp_allele
      }
    }
  }
}, error = function(e) {
  if (args$verbose) cat("Warning: Could not process allele alignment:", e$message, "\n")
  gwas <- as.data.table(gwas)
})

# Create SNP matching IDs ONCE
if (args$verbose) cat("Creating SNP lookup tables...\n")

gwas_snp_ids <- paste0(gwas$chromosome, "_", gwas$position, "_", 
                      gwas$effect_allele, "_", gwas$non_effect_allele)

# Create efficient lookup structures
gwas_lookup <- gwas %>%
  mutate(snp_id = gwas_snp_ids) %>%
  select(snp_id, pvalue, effect_size, standard_error, sample_size, frequency) %>%
  as.data.table()

# Set key for fast lookup
setkey(gwas_lookup, snp_id)

# Split eQTL data by gene ONCE
if (args$verbose) cat("Organizing eQTL data by gene...\n")
eqtl <- as.data.table(eqtl)
eqtl_by_gene <- split(eqtl, eqtl$phenotype_id)

# Clean up large intermediate objects
rm(gwas, eqtl, pos_eqtl, gwas_snp_ids)
gc()

unique_genes <- names(eqtl_by_gene)
if (args$verbose) cat("Processing", length(unique_genes), "genes...\n")

# =============================================================================
# OPTIMIZED COLOCALIZATION FUNCTION (NO FILE I/O)
# =============================================================================

run_coloc_optimized <- function(gene_ids, gwas_lookup, eqtl_by_gene, args) {
  results <- list()
  
  for (gene_id in gene_ids) {
    tryCatch({
      # Get gene-specific eQTL data from pre-loaded split
      if (!gene_id %in% names(eqtl_by_gene)) {
        results[[gene_id]] <- data.frame(
          Tissue = args$tis, Gene = gene_id, Trait = args$trt,
          nsnps = 0, PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA
        )
        next
      }
      
      gene_data <- eqtl_by_gene[[gene_id]]
      
      if (nrow(gene_data) == 0) {
        results[[gene_id]] <- data.frame(
          Tissue = args$tis, Gene = gene_id, Trait = args$trt,
          nsnps = 0, PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA
        )
        next
      }
      
      # Find overlapping SNPs using fast lookup
      overlap_snps <- intersect(gene_data$rsid, gwas_lookup$snp_id)
      
      if (length(overlap_snps) < args$min_overlap) {
        results[[gene_id]] <- data.frame(
          Tissue = args$tis, Gene = gene_id, Trait = args$trt,
          nsnps = length(overlap_snps), PP0 = NA, PP1 = NA, 
          PP2 = NA, PP3 = NA, PP4 = NA
        )
        next
      }
      
      # Subset data for overlapping SNPs
      gene_subset <- gene_data[gene_data$rsid %in% overlap_snps, ]
      gwas_subset <- gwas_lookup[gwas_lookup$snp_id %in% overlap_snps, ]
      
      # Ensure we have matching data
      if (nrow(gene_subset) == 0 || nrow(gwas_subset) == 0) {
        results[[gene_id]] <- data.frame(
          Tissue = args$tis, Gene = gene_id, Trait = args$trt,
          nsnps = 0, PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA
        )
        next
      }
      
      # Prepare data for coloc
      gene_coloc <- list(
        pvalues = as.numeric(gene_subset$pval_nominal),
        beta = as.numeric(gene_subset$slope),
        varbeta = as.numeric(gene_subset$slope_se)^2,
        snp = as.character(gene_subset$rsid),
        type = "quant",
        N = as.numeric(median(gene_subset$ma_samples, na.rm = TRUE)),
        MAF = as.numeric(gene_subset$af)
      )
      
      gwas_coloc <- list(
        pvalues = as.numeric(gwas_subset$pvalue),
        beta = as.numeric(gwas_subset$effect_size),
        varbeta = as.numeric(gwas_subset$standard_error)^2,
        snp = as.character(gwas_subset$snp_id),
        type = "quant",
        N = as.numeric(median(gwas_subset$sample_size, na.rm = TRUE)),
        MAF = as.numeric(gwas_subset$frequency)
      )
      
      # Check for valid data
      if (any(is.na(gene_coloc$pvalues)) || any(is.na(gwas_coloc$pvalues)) ||
          any(gene_coloc$pvalues <= 0) || any(gwas_coloc$pvalues <= 0)) {
        results[[gene_id]] <- data.frame(
          Tissue = args$tis, Gene = gene_id, Trait = args$trt,
          nsnps = length(overlap_snps), PP0 = NA, PP1 = NA, 
          PP2 = NA, PP3 = NA, PP4 = NA
        )
        next
      }
      
      # Run colocalization
      coloc_result <- coloc.abf(gene_coloc, gwas_coloc)
      
      # Extract results
      results[[gene_id]] <- data.frame(
        Tissue = args$tis,
        Gene = gene_id,
        Trait = args$trt,
        nsnps = coloc_result$summary[1],
        PP0 = coloc_result$summary[2],
        PP1 = coloc_result$summary[3],
        PP2 = coloc_result$summary[4],
        PP3 = coloc_result$summary[5],
        PP4 = coloc_result$summary[6]
      )
      
    }, error = function(e) {
      results[[gene_id]] <- data.frame(
        Tissue = args$tis, Gene = gene_id, Trait = args$trt,
        nsnps = NA, PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA
      )
    })
  }
  
  return(do.call(rbind, results))
}

# =============================================================================
# PARALLEL PROCESSING
# =============================================================================

# Process genes in chunks
gene_chunks <- split(unique_genes, ceiling(seq_along(unique_genes) / args$chunk_size))

if (args$verbose) cat("Processing", length(gene_chunks), "chunks with", args$workers, "workers...\n")

# Set up parallel processing with reduced memory requirements
options(future.globals.maxSize = 4 * 1024^3)  # 4GB limit
plan(multisession, workers = args$workers)

# Process chunks in parallel - data is passed as arguments, not read from files
if (args$verbose) cat("Starting parallel processing...\n")

all_results <- future_map(gene_chunks, 
                         ~run_coloc_optimized(.x, gwas_lookup, eqtl_by_gene, args),
                         .options = furrr_options(seed = TRUE))

# Combine results
final_results <- rbindlist(all_results)

# Write output
output_file <- file.path(args$outdir, paste0(args$tis, ".coloc.res"))
fwrite(final_results, output_file, sep = "\t")

# Summary statistics
if (args$verbose) {
  cat("\n=== Analysis Summary ===\n")
  cat("Total genes analyzed:", nrow(final_results), "\n")
  cat("Genes with valid results:", sum(!is.na(final_results$nsnps)), "\n")
  cat("Strong colocalization signals (PP4 > 0.5):", 
      sum(final_results$PP4 > 0.5, na.rm = TRUE), "\n")
  cat("Medium colocalization signals (PP4 > 0.1):", 
      sum(final_results$PP4 > 0.1, na.rm = TRUE), "\n")
  cat("Results written to:", output_file, "\n")
}

# Clean up
plan(sequential)