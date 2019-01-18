### SMRI Array ACC SZ vs CTRL
### Cell type enrichment analysis
### Method: EWCE
### Prep cell data (single cell transcriptome): 
### Darmanis et al PNAS 2015 (human brain cell types; fresh tissue from surgery; multiple individuals)
### Alexis Norris
### Created: 2018-11-05
### Modified: 2019-01-05


# Input files --------------------------------------------------
#load("input/celldata/Darmanis/rpkmCounts_darmanisSingleCell.rda")
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)
library(EWCE)

### Parameters
analysis_name <- "darmanis_adult"
fdr_edgeR <- 0.1            # For edgeR DEGs
fdr_method <- "fdr"         # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")


# Load cell data -----------------------------------------------
### Using normalized RPKM data from Jaffe (LIBD)
### Includes: geneRpkm, geneMap, pd, exonRpkm, exonMap, jRpkm, jMap
load("input/celldata/Darmanis/rpkmCounts_darmanisSingleCell.rda")        

### Use exprsData for gene-level (exons and jxns also available)
### ExprsData rownames are ensembl IDS
exprs_all <- geneRpkm
pheno_all <- pd


# Filter samples and cell types --------------------------------
### Add ID column
pheno_all$id <- rownames(pheno_all)

### Summarize dataset
pheno_summary <- pheno_all %>%
  group_by(sub_tissue, Cell_type, AgeGroup, tissue_type) %>%
  summarize(
    n_cells = n(),
    n_subjects = length(unique(as.character(Age))),
    ages = toString(sort(unique(as.character(Age)))),
    age_min = min(Age),
    age_max = max(Age),
    age_mean = mean(Age),
    age_sd = sd(Age),
    age_median = median(Age)
  ) %>%
  # Add group name
  mutate(group = paste(Cell_type, Cell_type, AgeGroup, sep = "_")) 

### Filter
pheno <- pheno_all %>%
  # Adults only
  filter(AgeGroup == "postnatal") %>%
  # Celltypes
  filter(
    !(Cell_type %in% c(
      "Fetal_quiescent", "Fetal_replicating",
      # Remove mixed/hybrid cell types
      "Hybrid", 
      # Remove ifdata for <5 cells
      "OPC",                  # n = 2
      "Fetal_quiescent"       # n = 1
    ))) %>%
  mutate(group = paste(Cell_type, sub_tissue, AgeGroup, sep = "_")) 

### Note if sample group included in analysis
pheno_summary %>% 
  mutate(included_in_analysis = group %in% pheno$group) %>%
  select(group, included_in_analysis, everything()) %>%
  write_csv("input/celldata/Darmanis/Darmanis_summary.csv")

### Remove filtered samples from exprsData
exprs <- exprs_all[ , pheno$id]


# Filter features ----------------------------------------------
### Remove unexpressed
exprs <- exprs[rowSums(exprs) > 0, ]


# Add approved HGNC gene symbols -------------------------------
### If ensembl ID not in HGNC, make gene_symbol NA and remove

### Load HGNC anno (from anno_hgnc.R)
hgnc <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")$by_ensembl
names(hgnc) <- recode(
  names(hgnc),
  "symbol" = "gene_symbol"
)

### ensembl IDs to get gene symbols for
ensembl_ids <- rownames(exprs)

### Subset HGNC for ensembl_ids
df_hgnc <- data_frame("ensembl_gene_id" = ensembl_ids) %>%
  left_join(hgnc, by = "ensembl_gene_id") %>%
  filter(!is.na(gene_symbol))

### Make gene_symbols the exprsData rownames
### If ensembl ID not in HGNC, remove
stopifnot(!duplicated(df_hgnc$ensembl_gene_id))
stopifnot(!duplicated(df_hgnc$gene_symbol))
stopifnot(df_hgnc$ensembl_gene_id %in% rownames(exprs))
exprs <- exprs[which(rownames(exprs) %in% df_hgnc$ensembl_gene_id), ]
exprs <- exprs[match(df_hgnc$ensembl_gene_id, rownames(exprs)), ]
stopifnot(!duplicated(df_hgnc$gene_symbol))
rownames(exprs) <- df_hgnc$gene_symbol

### Save the new featureData
write_tsv(exprs, "input/celldata/Darmanis/Darmanis_featureData.csv")


# Assemble dataset ---------------------------------------------
stopifnot(colnames(exprs) == pheno$id)
celldata <- list(
  "exprs" = as.matrix(exprs),
  "pheno" = pheno
)

### Save
saveRDS(celldata, paste0("input/celldata/ewce_", analysis_name, ".rds"))


# Calculate specificity (all genes) ----------------------------
### Cell type column
celldata$pheno$group <- celldata$pheno$Cell_type

### Specifity Fxn (modified from EWCE package)
generate.celltype.data <- function(exprs, annotLevels) {
  require(parallel)
  no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  lapply(annotLevels, test <- function(x, exprs) {
    if (length(x) != dim(exprs)[2]) stop("Error: length of all annotation levels must equal the number of columns in exp matrix")
  }, exprs)
  exp2 <- suppressWarnings(apply(exprs, 2, function(x) {
    storage.mode(x) <- "double"; x
  }))
  ctd <- list()
  for (i in 1:length(annotLevels)) {
    ctd[[length(ctd) + 1]] <- list(annot = annotLevels[[i]])
  }
  aggregate.over.celltypes <- function(rowOfMeans, celltypes, func = "mean") {
    if (func == "mean") {
      exp_out <- as.matrix(data.frame(aggregate(
        rowOfMeans, 
        by = list(celltypes), 
        FUN = mean
      )))
    }
    else if (func == "median") {
      exp_out = as.matrix(data.frame(aggregate(
        rowOfMeans, 
        by = list(celltypes), 
        FUN = median
      )))
    }
    rownames(exp_out) <- exp_out[ , "Group.1"]
    exp_out <- exp_out[, 2]
    exp_out2 <- as.numeric(exp_out)
    names(exp_out2) <- names(exp_out)
    return(exp_out2)
  }
  calculate.meanexp.for.level <- function(ctd_oneLevel, expMatrix) {
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
      print(dim(expMatrix)[2])
      print(length(ctd_oneLevel$annot))
      if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) != 0) {
        stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
      }
      ctd_oneLevel$mean_exp <- expMatrix
    }
    else {
      mean_exp <- apply(expMatrix, 1, aggregate.over.celltypes, ctd_oneLevel$annot)
      ctd_oneLevel$mean_exp <- t(mean_exp)
    }
    return(ctd_oneLevel)
  }
  calculate.medianexp.for.level <- function(ctd_oneLevel, expMatrix) {
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
      print(dim(expMatrix)[2])
      print(length(ctd_oneLevel$annot))
      if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) != 0) {
        stop("There are an equal number of celltypes in expMatrix and ctd_oneLevel but the names do not match")
      }
      ctd_oneLevel$median_exp <- expMatrix
    }
    else {
      median_exp <- apply(expMatrix, 1, aggregate.over.celltypes, ctd_oneLevel$annot, func = "median")
      ctd_oneLevel$median_exp <- t(median_exp)
    }
    return(ctd_oneLevel)
  }
  calculate.specificity.for.level <- function(ctd_oneLevel) {
    normalised_meanExp <- t(t(ctd_oneLevel$mean_exp) * (1/colSums(ctd_oneLevel$mean_exp)))
    normalised_medianExp <- t(t(ctd_oneLevel$median_exp) * (1/colSums(ctd_oneLevel$mean_exp)))
    ctd_oneLevel$specificity <- normalised_meanExp/(apply(normalised_meanExp, 1, sum) + 1e-12)
    ctd_oneLevel$median_specificity <- normalised_medianExp/(apply(normalised_meanExp, 1, sum) + 1e-12)
    return(ctd_oneLevel)
  }
  ctd2 <- mclapply(ctd, calculate.meanexp.for.level, exp2)
  ctd2 <- mclapply(ctd2, calculate.medianexp.for.level, exp2)
  ctd3 <- mclapply(ctd2, calculate.specificity.for.level)
  stopCluster(cl)
  return(ctd3)  
}

### Get specificity for all genes 
set.seed(1116)
sct_data_all <- generate.celltype.data(
  exprs = celldata$exprs,
  annotLevels = list(l1 = celldata$pheno$group)
)

### Save
saveRDS(sct_data_all, paste0("input/celldata/ewce_", analysis_name, "_processed.rda"))


# Calculate specificity (only informative genes) ---------------
### Get specificity for informative genes only <- use this one in EWCE analysis!
  # Drop genes that don't differ across celltypes -
  # Modified for:
  # F.fdr_cutoff changeable (default 1e-05)
  # remove cell types expressed in less than 1 group of celltypes (min_celltypes)
  # remove lowly expressed genes (min 0.3 rpkm)
drop.uninformative.genes <- function (exprs, group, fdr_cutoff = 1e-05, e_min = 0.3) {
  group <- as.character(group)
  min_samples <- min(summary(factor(group)))
  n_all <- nrow(exprs)
  print(paste("Requiring min", e_min, "RPKM in at least", min_samples))
  exprs <- exprs[rowSums(exprs > e_min) >= min_samples, ]
  n_filt <- nrow(exprs)
  print(paste(n_all-n_filt, "of", n_all, "genes removed;", n_filt, "genes remain"))
  mod <- model.matrix(~group)
  fit <- limma::lmFit(exprs, mod)
  eb <- limma::eBayes(fit, robust = TRUE)
  fdr <- p.adjust(eb$F.p.value, method = "fdr")
  exprs <- exprs[fdr < fdr_cutoff, ]
  print(paste(nrow(exprs), "genes remain after F.FDR <", fdr_cutoff))
  return(exprs)
}
set.seed(1116)
write_lines(capture.output(
celldata_informative <- drop.uninformative.genes(
  exprs = celldata$exprs,
  group = celldata$pheno$group,
  fdr_cutoff = 1E-05, 
  e_min = 0.3
)), paste0("input/celldata/ewce_", analysis_name, "_informative_filter.txt")
)

### Generate specificity
set.seed(1116)
sct_data_informative <- generate.celltype.data(
  exprs = celldata_informative,
  annotLevels = list(l1 = celldata$pheno$group)
)

### Save
saveRDS(sct_data_informative, paste0("input/celldata/ewce_", analysis_name, "_processed_informative.rds"))


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
