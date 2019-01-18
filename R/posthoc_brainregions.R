### Stanley ACC SZ vs CTRL
### edgeR: posthoc brain regions
### models include all features
### which are then collapsed
### Alexis Norris
### Created: 2018-12-05
### Modified: 2019-01-16


# Input files --------------------------------------------------
#dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")


# Setup --------------------------------------------------------
### Packages
library(edgeR)          # for DGE object
library(tidyverse)      # wrangling

### Parameters
analysis_name <- "posthoc_brainregions"


# Load DGEs ----------------------------------------------------
### datasets for brain regions: ACC, PFC, HPC
### feature levels: genes, exons, jxns, regions
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")

### Only PFC and HPC
dgeList <- dgeList_all[-grep("ACC", names(dgeList_all))]


# Remove outlier by qSVA ---------------------------------------
### HPC's S15 is qSVA outlier
for (i in names(dgeList)[grep("HPC", names(dgeList))]) {
  # Subset
  y <- dgeList[[i]]
  
  # Remove outlier
  y <- y[ , which(y$samples$Sample != "S15"), keep.lib.sizes = TRUE]
  
  # Return
  dgeList[[i]] <- y
}


# Run ----------------------------------------------------------
edgeRList <- list()
for (i in names(dgeList)) {
  print(i)
  
  # Subset for dge
  y <- dgeList[[i]]
  
  # TMM normalization
  y <- calcNormFactors(y)
  
  # Model for ACC
  #design <- model.matrix(~Age + Sex + qSV1 + qSV2 + Dx, y$samples)
  
  # Use model based on # qSVs identified
  brain_region <- gsub("_.*", "", i)
  if (brain_region == "PFC") {
    design <- model.matrix(~Age + Sex + Dx + qSV1 + qSV2 + qSV3 + qSV4, y$samples)
  }
  if (brain_region == "HPC") {
    design <- model.matrix(~Age + Sex + Dx + qSV1 + qSV2 + qSV3 + qSV4 + qSV5, y$samples)
  }
  
  # Fit
  y <- estimateDisp(y, design, verbose = TRUE, robust = TRUE) 
  fit <- glmFit(y, design, robust = TRUE)
  
  # Run
  edgeRList[[i]] <- data.frame(topTags(glmLRT(fit, coef = "DxSZ"), n = Inf))
}

### Save full results
saveRDS(edgeRList, "output/edgeR/posthoc/edgeR_PFC_HPC_resList.rds")


# Collapse genes -----------------------------------------------
### Using same method used as ACC
### Collapse genes, take most significant feature for each
### Separate up vs down, in case one transcript is up and another down
for (i in c("PFC", "HPC")) {
  res_all <- edgeRList[grep(i, names(edgeRList))] %>%
    bind_rows(.id = "feature_level") 
  resList <-  list(
    "down" = filter(res_all, logFC < 0),
    "up" = filter(res_all, logFC >= 0)
  )
  lapply(resList, function (df) {
    df %>%
      filter(!is.na(gene_symbol)) %>%
      arrange(PValue) %>%
      .[!duplicated(.$gene_symbol), ]
  }) %>%
    bind_rows(.id = "direction") %>%
    mutate(FDR = p.adjust(.$PValue, method = "fdr")) %>%
    write_tsv(paste0("output/edger/posthoc/edgeR_", i, "_SZvsCTRL_res.tsv")) # contains dups, since allowed both up and down
}


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
