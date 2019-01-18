### SMRI Array ACC SZ vs CTRL
### Filter low abundant features
### Alexis Norris
### Created: 2018-11-09
### Modified: 2019-01-02


# Input files --------------------------------------------------
### DGE list filtered for samples in final cohort
### And features not expressed in any sample were removed
### And repeatMasker overlap was added for region-level features
#dgeList <- readRDS("input/dge/ACC_PFC_HPC.rds")


# Setup --------------------------------------------------------
### Packages
library(edgeR)          # DGEs of exprs, pheno, and feature/anno data
library(tidyverse)      # plotting; wrangling; select fxn gets masked by GenomicRanges

### Parameters
analysis_name <- "exprsData_filter"
datasets <- c("ACC", "PFC", "HPC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_")


# Load DGEs ----------------------------------------------------
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC.rds"); beepr::beep()

### Want all dges
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# Normalize ----------------------------------------------------
### TMM normalization
dgeList <- lapply(dgeList, calcNormFactors); beepr::beep()


# Filter Features ----------------------------------------------
dgeList_filt <- list()

### Capture filtering stats
logList <- list()

### Filter -- "[]" means not used
  # 1. lowly expressed features (since they distort model)
    # We consider a gene to be expressed at a reasonable level in a sample if 
    # it has at least 10 counts in at least one diagnosis group
  # 2. chrY
  # [3. no HGNC approved gene_symbol for ensembl ID  (includes ensembl IDs that are no longer)]
  # 4. region far from gene
  # 5. region with high % repetitive overlap (added above)
  # [6. non-polyA (since ACC is mRNA library prep)]
min_counts <- 5              # counts
#min_cpm <- 0.3              # alternative to counts: cpm
max_distToGene <- 5000       # bp
max_rmskOverlap <- 0.50      # %
for (i in names(dgeList)) {
  print(i)
  
  # Subset
  brain_region <- gsub("_.*", "", i)
  feature_level <- gsub(".*_", "", i)
  y <- dgeList[[i]]
  
  # Log
  logList[[i]]$n_features_preFilter <- nrow(y$counts)
  
  # List of genes to remove
  removeList <- list()

  # 1. Low expression
    # Calculate minimum number of samples in a group
  n_CTRL <- length(y$samples$Dx[y$samples$Dx == "CTRL"])
  n_SZ <- length(y$samples$Dx[y$samples$Dx == "SZ"])
  min_sample <- min(n_CTRL, n_SZ)
    # Calculate libraries & read depth
  mean_libsize <- mean(y$samples$lib.size)/1E06
  min_libsize <- min(y$samples$lib.size)/1E06
  max_libsize <- max(y$samples$lib.size)/1E06
  mean_TotalLibSize <- mean(y$samples$TotalLibSize)/1E06
  min_TotalLibSize <- min(y$samples$TotalLibSize)/1E06
  max_TotalLibSize <- max(y$samples$TotalLibSize)/1E06
    # Alternative: Use counts & filter based on library size <------ did NOT use this; used preset cpm above
  libsize <- mean_libsize; lib_size_used <- "mean_libsize"
  min_cpm <- min_counts/libsize   
    # Subset
  keep <- rowSums(cpm(y) > min_cpm) >= min_sample 
  removeList <- append(removeList, list(y$genes$feature_id[!keep]))
      # Log
  logList[[i]]$min_counts <- min_counts
  logList[[i]]$lib_size_used <- lib_size_used
  logList[[i]]$mean_libsize <- mean_libsize
  logList[[i]]$min_libsize <- min_libsize
  logList[[i]]$max_libsize <- max_libsize
  logList[[i]]$mean_TotalLibSize <- mean_TotalLibSize
  logList[[i]]$min_TotalLibSize <- min_TotalLibSize
  logList[[i]]$max_TotalLibSize <- max_TotalLibSize
  logList[[i]]$n_CTRL <- n_CTRL 
  logList[[i]]$samples_CTRL <- paste(sort(as.character(y$samples$Sample[y$samples$Dx == "CTRL"])), collapse = ",")
  logList[[i]]$n_SZ <- n_SZ 
  logList[[i]]$samples_SZ <- paste(sort(as.character(y$samples$Sample[y$samples$Dx == "SZ"])), collapse = ",")
  logList[[i]]$min_sample <- min_sample
  logList[[i]]$min_cpm <- min_cpm
  logList[[i]]$n_features_lowExprs <- length(keep[!keep])

  # 2. chrY --> run post-edgeR
  keep <- y$genes$chr != "chrY"
  logList[[i]]$n_features_chrY <- length(keep[!keep])
  removeList <- append(removeList, list(y$genes$feature_id[!keep]))

  # 3. no HGNC approved gene_symbol --> run post-edgeR
  #if (!is.null(y$genes$gene_symbol)) {
  #  keep <- !is.na(y$genes$gene_symbol)
  #  logList[[i]]$n_features_noGeneSymbol <- length(keep[!keep])
  #  removeList <- append(removeList, list(y$genes$feature_id[!keep]))
  #}
  
  # Specific to region-level features --> run post-edgeR
  if (feature_level == "regions") {
    # 4. region far from gene
    stopifnot(y$genes$region_distToGene != "")
    stopifnot(!is.na(y$genes$region_distToGene))
    logList[[i]]$max_distToGene_regions <- max_distToGene
    keep <- y$genes$region_distToGene <= max_distToGene 
    logList[[i]]$n_features_noGeneNearby <- length(keep[!keep])
    removeList <- append(removeList, list(y$genes$feature_id[!keep]))

    # 5. region with high % repetitive overlap --> run post-edgeR
    stopifnot(y$genes$rmsk_overlap != "")
    stopifnot(!is.na(y$genes$rmsk_overlap))
    logList[[i]]$max_rmskOverlap_regions <- max_rmskOverlap
    keep <- y$genes$rmsk_overlap < max_rmskOverlap
    logList[[i]]$n_features_rmsk <- length(keep[!keep])
    removeList <- append(removeList, list(y$genes$feature_id[!keep]))
  }
  
  # Remove non-coding, since ACC data is from mRNA library --> don't do this; it's an assumption
  #if (brain_region == "ACC") {
  #  keep <- y$genes$gene_type %in% c("protein_coding", "processed_transcript")
  #  logList[[i]]$n_features_notmRNA <- length(keep[!keep])
  #  removeList <- removeList <- append(removeList, list(y$genes$feature_id[!keep]))
  #}

  # Remove all non-keep
  # Do this at end, to get actual num of total features removed for each reason in logList
  remove_features <- unique(as.character(unlist(removeList)))
  y <- y[!(y$genes$feature_id %in% remove_features), , keep.lib.sizes = FALSE]
  stopifnot(rownames(y$genes) == rownames(y$counts))
  
  # Add final num features (survived filtering; using for edgeR) and gene symbols
  logList[[i]]$n_features_postFilter <- nrow(y$counts)
  
  # Summarize genes for each level
  logList[[i]]$n_no_gene_symbol <- length(is.na(y$genes$gene_symbol))
  logList[[i]]$n_unique_genes <- length(unique(y$genes$gene_symbol[!is.na(y$genes$gene_symbol)]))
  logList[[i]]$n_feature_maps_multiple_hgnc_symbols <- length(!is.na(y$genes$mult_matches)[!is.na(y$genes$mult_matches) != "NO"])
  
  # Return
  dgeList_filt[[i]] <- y
}

### Write out filtering stats & summary
logList %>%
  bind_rows(.id = "dataset") %>%
  separate(dataset, into = c("brain_region", "feature_level"), remove = FALSE) %>%
  mutate(percent_remaining = n_features_postFilter/n_features_preFilter) %>%
  write_csv("logs/exprsData/exprsData_filterFeatures_stats.csv")


# Normalize ----------------------------------------------------
### Re-normalize filtered data
dgeList_filt <- lapply(dgeList_filt, calcNormFactors); beepr::beep()


# Save ---------------------------------------------------------
saveRDS(dgeList_filt, "input/dge/ACC_PFC_HPC_filtered.rds")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
