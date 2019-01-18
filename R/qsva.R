### SMRI Array ACC SZ vs CTRL
### qSVA
### Alexis Norris
### Created: 2018-12-05
### Modified: 2019-01-02


# Input files --------------------------------------------------
#dgeList <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved.rds")


# Setup --------------------------------------------------------
### Packages
library(sva)       # run qSVA
library(edgeR)     # DGEs
library(tidyverse) # wrangling

### Parameters
analysis_name <- "qsva"
datasets <- c("ACC", "PFC", "HPC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_")


# Load DGEs ----------------------------------------------------
### datasets for brain regions: ACC, PFC, HPC
### feature levels: genes, exons, jxns, regions
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved.rds")

### Want all DEGs
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# qSVA ---------------------------------------------------------
### Check that each feature level has the correct deg regions
# ACC & HPC: polyA (n = 1000 regions)
polyA <- 1000
# PFC: ribozero (n = 515 regions)
ribozero <- 515
# Check
write_lines(capture.output(
  for (i in datasets) {
    # Number of degradation regions there should be for dataset
    n_qsva <- ifelse(
      i %in% c("ACC", "HPC"), polyA, ifelse(
        i == "PFC", ribozero, "error"
      )
    )
    
    # Check for each feature level
    for (j in c(paste(i, feature_levels, sep = "_"))) {
      
      # Number of deg regions in the dataset's phenodata
      y <- dgeList[[j]]
      n_pheno <- length(
        y$samples[ , grep("^degchr", names(y$samples))]
      )
      
      # Check
      stopifnot(n_pheno == n_qsva)
      
      # Print out to file
      print(paste0(
        j, ": n degradation regions = ", n_pheno
      ))
    }
  }), "logs/qsva/qsva_num_degradationRegions.txt"
)

### Only have to run once for dataset, and use that for all feature levels
### Write log out of # qSVs ID'd
write_lines(capture.output(
for (i in datasets) {
  # Subset
  pheno <- dgeList[[paste0(i, "_genes")]]$samples
  
  # PCA to find all qSVs
  deg_cov <- t(pheno[, grep("degchr", names(pheno))])
  deg_pca <- prcomp(t(log2(deg_cov + 1)))
  write_lines(capture.output(
    summary(deg_pca)), 
    paste0("output/qsva/qsva_pca_", i, "_summary.txt")
  )
  
  # Calculate qSVs
  design <- model.matrix(~Age + Sex + Dx, pheno)
  k <- num.sv(log2(deg_cov + 1), design)
  qsv <- data.frame(deg_pca$x[ , 1:k])
  names(qsv) <- gsub("PC", "qSV", names(qsv))
  
  # Print out num qSVs
  print(paste0(
    i, ": n qSVs = ", k
  ))
  
  # Add qSV cols to DGEs' phenoData (for all 4 feature levels)
  for (j in c(paste(i, feature_levels, sep = "_"))) {
    # Subset
    y <- dgeList[[j]]
    
    # Check
    stopifnot(rownames(y$samples) == rownames(qsv))
    
    # Add qSVs to phenoData
    for (k_num in 1:k) {
      y$samples[[paste0("qSV", k_num)]] <- qsv[[paste0("qSV", k_num)]]
    }
    
    # Save
    dgeList[[j]]$samples <- y$samples
    print(j)
  }
}
), "logs/qsva/qsva_num_qSVs.txt")


# Save ---------------------------------------------------------
saveRDS(dgeList, "input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds"); beepr::beep()


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
