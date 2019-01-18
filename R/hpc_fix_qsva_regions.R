### SMRI Array ACC SZ vs CTRL
### HPC - fix qSVA regions
### Should be polyA not riboZero
### Alexis Norris
### Created: 2018-11-09
### Modified: 2018-12-29


# Notes --------------------------------------------------------
### DGEs created in Stanley_Array_ACC_PFC_HPC -- HPC has wrong degrad regions in y$samples
### Should be polyA!!! (fixed below)
### Original data source: /Volumes/Hopkins2018/Projects/Human/Stanley/ArrayCollection/ACC_PFC_HPC/20180605/input


# Input files --------------------------------------------------
#dgeList <- list(
# "HPC_genes" = readRDS("input/dge/dge_HPC_genes_ribozero.rds"),
# "HPC_exons" = readRDS("input/dge/dge_HPC_exons_ribozero.rds"),
# "HPC_jxns" = readRDS("input/dge/dge_HPC_jxns_ribozero.rds"),
# "HPC_regions" = readRDS("input/dge/dge_HPC_regions_ribozero.rds")
#)


# Setup --------------------------------------------------------
### Packages
library(edgeR)     # DGEs of exprs, pheno, and feature/anno data
library(tidyverse) # plotting; wrangling

### Parameters
analysis_name <- "HPC_fix-qsva_ribozero-to-polyA"
datasets <- c("HPC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_")


# Load DGEs ----------------------------------------------------
### feature levels: genes, exons, jxns, regions
### Original data source: /Volumes/Hopkins2018/Projects/Human/Stanley/ArrayCollection/ACC_PFC_HPC/20180605/input
dgeList <- list(
  "HPC_genes" = readRDS("input/dge/dge_HPC_genes_ribozero.rds"),
  "HPC_exons" = readRDS("input/dge/dge_HPC_exons_ribozero.rds"),
  "HPC_jxns" = readRDS("input/dge/dge_HPC_jxns_ribozero.rds"),
  "HPC_regions" = readRDS("input/dge/dge_HPC_regions_ribozero.rds")
)

### Only want HPC
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# Fix HPC degCov to polyA regions ------------------------------
### Only need to use one feature level's phenoData to get degcov
### Then will add to all feature levels' DGEs

### PhenoData
pd <- dgeList$HPC_genes$samples

### Variables
ids <- "FileName"     # col in y$samples (pd) that is filename used in degStats txt file
lib_type <- "polyA"   # vs ribozero
read_length <- 100

### List degStats files for samples
deg_files <- c(paste0(
  "input/degStats/", 
  pd[[ids]], 
  paste0(".degradeStats_", lib_type, ".txt") 
))
names(deg_files) <- pd[[ids]]

### Check that all samples have a degradeStats file
stopifnot(file.exists(deg_files)) 

### Create degradation coverage matrix
deg_cov <- sapply(
  deg_files, function (x) {
    # Load degStats file
    read.delim(
      # generate matrix, adjusting by read length
      pipe(paste("cut -f10", x)), 
      as.is = TRUE
    )$sum/read_length
  }
)

### Get the region positions by grabbing the first one
deg_regions <- read.delim(deg_files[1])

### Make region position the rowname
rownames(deg_cov) <- paste0(
  deg_regions$X.chrom, ":", 
  deg_regions$start, "-", deg_regions$end
)

### Check order
stopifnot(colnames(deg_cov) == pd[[ids]])

### Normalize by libSize 
### which was calculated on cluster using Jaffe script + chrsizes txt file from UCSC
median_libsize <- median(pd$TotalLibSize, na.rm = TRUE) 
deg_cov_adj <- deg_cov/matrix(
  rep(
    # Match sample order
    pd[match(colnames(deg_cov), pd[[ids]]), 
       # Divide by median Libsize -- to reduce #s
       ]$TotalLibSize/median_libsize
  ), 
  nc = ncol(deg_cov), nr = nrow(deg_cov), byrow = TRUE
)

### Add degcov to phenoDatas in all feature levels
deg_cov_adj_t <- as.data.frame(t(deg_cov_adj))
names(deg_cov_adj_t) <- paste0("deg", names(deg_cov_adj_t))
dgeList <- lapply(dgeList, function (y) {
  # Subset
  pheno_ribo <- y$samples
  
  # Check
  stopifnot(pheno_ribo$FileName == rownames(deg_cov_adj_t))
  
  # Replace ribozero with polyA
  pheno_polyA <- cbind.data.frame(
    pheno_ribo[ , -grep("^degchr", names(pheno_ribo))], 
    deg_cov_adj_t
  )
  
  # Add back to dgeList
  y$samples <- pheno_polyA
  
  # Return
  y
})


# Save ---------------------------------------------------------
saveRDS(dgeList$HPC_genes, "input/dge/dge_HPC_genes_polyA.rds")
saveRDS(dgeList$HPC_exons, "input/dge/dge_HPC_exons_polyA.rds")
saveRDS(dgeList$HPC_jxns, "input/dge/dge_HPC_jxns_polyA.rds")
saveRDS(dgeList$HPC_regions, "input/dge/dge_HPC_regions_polyA.rds")

# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
