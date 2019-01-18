### SMRI Array ACC SZ vs CTRL
### Cell type enrichment analysis
### Method: EWCE
### Prep cell data (single cell transcriptome): 
### Lake et al Science 2016 (human neuron subtypes; postmortem; 1 individual [female])
### Alexis Norris
### Created: 2018-11-05
### Modified: 2019-01-05


# Note ---------------------------------------------------------
### Cell types from multiple cortical regions
  # BA8 = FC   
  # BA10 = anterior PFC 
  # BA17 = V1 (primary visual cortex)
  # BA21 = mTC (middle temporal gyrus)
  # BA22 = sTC (superior temporal gyrus)
  # BA41 = aTC (auditory cortex)

# Input files --------------------------------------------------
#full_data <- read.csv("input/celldata/Lake/Lake-2016_Gene_TPM_filtered.csv", check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)
library(EWCE)

### Parameters
analysis_name <- "lake"


# Load cell data -----------------------------------------------
### Load data (TPM data from paper)
### One table contains both exprsData and phenoData
### Use check.names = FALSE to prevent R from renaming genes in colnames that have "-" --> "."
### e.g. ALDH1L1-AS2 --> ALDH1L1.AS2
full_data <- read.csv(
  "input/celldata/Lake/Lake-2016_Gene_TPM_filtered.csv", 
  check.names = FALSE,
  header = TRUE, stringsAsFactors = FALSE
)   

### Extract phenoData
pheno_cols <- c(
  "Published_Sample_Name", "Sample_Name",
  "SubGroup", "BA", "SCAP.T_ID", "Group"
)
pheno_all <- full_data %>%
  select(one_of(pheno_cols)) %>%
  mutate(id = Sample_Name) %>%
  as.data.frame()

### Extract exprsData
exprs_all <- full_data %>%
  select(-one_of(pheno_cols)) %>%
  as.matrix() 
rownames(exprs_all) <- pheno_all$id
exprs_all <- t(exprs_all)


# Filter samples and cell types --------------------------------
### Summarize dataset
pheno_summary <- pheno_all %>%
  group_by(Group, SubGroup, BA) %>%
  summarize(
    n_cells = n(),
    n_subjects = 1
  ) 

### Filter
### Remove if data for <5 cells
small_n_cells <- as.character(filter(pheno_summary, n_cells < 5)$Group)
pheno <- pheno_all %>%
  filter(!(Group %in% small_n_cells))

### Note if sample group included in analysis
pheno_summary %>% 
  mutate(included_in_analysis = Group %in% pheno$Group) %>%
  select(Group, included_in_analysis, everything()) %>%
  write_csv("input/celldata/Lake/Lake_dataset_summary.csv")

### Remove filtered samples from exprsData
exprs <- exprs_all[ , which(colnames(exprs_all) %in% pheno$id)]


# Filter features ----------------------------------------------
### Remove unexpressed
exprs <- exprs[rowSums(exprs) > 0, ]


# Update outdated HGNC gene symbols ----------------------------
### Check that gene_symbol is the current approved HGNC symbol
### if outdated (alias), update
### *using gene symbols* (not ensembl IDs!)

### Load HGNC anno (from anno_hgnc.R)
hgnc_keys <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")
hgnc_symbol <- hgnc_keys$by_symbol
hgnc_prev <- hgnc_keys$previous_to_symbol_noDuplicates
hgnc_alias <- hgnc_keys$alias_to_symbol_noDuplicates
names(hgnc_symbol) <- recode(names(hgnc_symbol), "symbol" = "approved_symbol")
names(hgnc_prev) <- recode(names(hgnc_prev), "symbol" = "approved_symbol")
names(hgnc_alias) <- recode(names(hgnc_alias), "symbol" = "approved_symbol")

### Current gene symbols
df <- data_frame(
  "gene_symbol" = rownames(exprs),
  # plus column to hold original gene symbol
  "gene_symbol_original" = rownames(exprs)
  )

### Join by gene symbols
# Outcome 1. Ok (is current approved symbol)
ok <- df %>%
  inner_join(hgnc_symbol, by = c("gene_symbol" = "approved_symbol")) %>%
  # Fix classes, for binding later
  mutate_if(is.factor, as.character) %>%
  mutate(gene_symbol_corrected_reason = "ok")
stopifnot(ok$mult_matches == "NO")

# Outcome 2. Updated gene symbol (was previous)
updated_prev <- df %>%
  # symbol isn't ok (approved symbol)
  filter(!(gene_symbol %in% ok$gene_symbol_original)) %>%
  # try to match with previous symbol
  inner_join(hgnc_prev, by = c("gene_symbol" = "prev_symbol")) %>%
  # update gene symbol
  mutate(gene_symbol = approved_symbol) %>%
  # Note that it was changed
  mutate(gene_symbol_corrected_reason = "previous symbol") %>%
  # Add full hgnc anno
  select(-approved_symbol) %>%
  left_join(select(hgnc_symbol, -mult_matches), by = c("gene_symbol" = "approved_symbol"))
# Check if there are multiple genes that previous symbol maps to
# If there are, then note it
if ("YES" %in% updated_prev$mult_matches)  {
  updated_prev <-updated_prev %>%
    filter(mult_matches == "YES") %>%
    mutate(
      gene_symbol_corrected_reason = "previous symbol with multiple approved gene symbols"
    )
}

# Outcome 2. Updated gene symbol (was alias)
updated_alias <- df %>%
  # symbol isn't ok (approved symbol)
  filter(!(gene_symbol %in% c(ok$gene_symbol_original, updated_prev$gene_symbol_original))) %>%
  # try to match with alias
  inner_join(hgnc_alias, by = c("gene_symbol" = "alias_symbol")) %>%
  # update gene symbol
  mutate(gene_symbol = approved_symbol) %>%
  # Note that it was changed
  mutate(gene_symbol_corrected_reason = "alias") %>%
  # Add full hgnc anno
  dplyr::select(-approved_symbol) %>%
  left_join(dplyr::select(hgnc_symbol, -mult_matches), by = c("gene_symbol" = "approved_symbol"))
# Check if there are multiple genes that alias maps to
# If there are, then note it
if ("YES" %in% updated_alias$mult_matches)  {
  updated_alias <-updated_alias %>%
    filter(mult_matches == "YES") %>%
    mutate(
      gene_symbol_corrected_reason = "alias with multiple approved gene symbols"
    )
}

### Outcome 4. old (ensembl gene id not in HGNC )
not_in_hgnc <- df %>%
  filter(!(gene_symbol %in% c(ok$gene_symbol_original, updated_prev$gene_symbol_original, updated_alias$gene_symbol_original))) %>%
  # If there are, then note it and make gene symbol NA
  mutate(
    gene_symbol = NA,
    gene_symbol_corrected_reason = "not in HGNC"
  )

### Combine outcomes
df_hgnc <- bind_rows(
  ok,
  updated_prev,
  updated_alias,
  not_in_hgnc
) %>% 
  select(
    gene_symbol, gene_symbol_original, 
    gene_symbol_corrected_reason, everything()
  ) 
stopifnot(nrow(df_hgnc) == nrow(df))

### Remove genes with ambigous or no approved gene symbol
reasons <- c(
  "previous symbol with multiple approved gene symbols",
  "alias with multiple approved gene symbols",
  "not in HGNC"
)
df_hgnc %>%
  filter(gene_symbol_corrected_reason %in% reasons) %>%
  mutate(approved_symbols = gene_symbol) %>%
  select(gene_symbol_original, approved_symbols, gene_symbol_corrected_reason) %>%
  write_csv("input/celldata/Lake/Lake_genes_removed.csv")
df_hgnc <- df_hgnc %>%
  filter(!(gene_symbol_corrected_reason %in% reasons))

### Make gene_symbols the exprsData rownames
### If ensembl ID not in HGNC, remove
stopifnot(!duplicated(df_hgnc$gene_symbol))
exprs <- exprs[which(rownames(exprs) %in% df_hgnc$gene_symbol_original), ]
exprs <- exprs[match(df_hgnc$gene_symbol_original, rownames(exprs)), ]
stopifnot(rownames(exprs) == df_hgnc$gene_symbol_original)
stopifnot(!duplicated(df_hgnc$gene_symbol))
rownames(exprs) <- df_hgnc$gene_symbol

### Save the new featureData
write_tsv(exprs,"input/celldata/Lake/Lake_featureData.csv")


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
celldata$pheno$group <- celldata$pheno$SubGroup
celldata$pheno$groupBA <- celldata$pheno$Group

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
  annotLevels = list(
    l1 = celldata$pheno$group,
    l2 = celldata$pheno$groupBA
  )
)

### Save
saveRDS(sct_data_all, paste0("input/celldata/ewce_", analysis_name, "_processed.rda"))


# Calculate specificity (only celltype-informative genes) ------
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
  print(paste(nrow(exprs), "genes remain after FDR <", fdr_cutoff))
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
  annotLevels = list(
    l1 = celldata$pheno$group,
    l2 = celldata$pheno$groupBA
  )
)

### Save
saveRDS(sct_data_informative, paste0("input/celldata/ewce_", analysis_name, "_processed_informative.rds"))


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
