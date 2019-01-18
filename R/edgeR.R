### Stanley ACC SZ vs CTRL
### edgeR (ACC)
### Features filtered in exprsData_filter.R
### Alexis Norris
### Created: 2018-12-17
### Modified: 2019-01-02


# Input files --------------------------------------------------
#dgeList <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")


# Setup --------------------------------------------------------
### Packages
library(edgeR)          # for DGE object
library(tidyverse)      # wrangling

### Parameters
analysis_name <- "edgeR"
datasets <- c("ACC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_")
fdr_cutoff <- 0.1


# Load DGEs ----------------------------------------------------
### datasets for brain regions: ACC, PFC, HPC
### feature levels: genes, exons, jxns, regions
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")

### Only want ACC
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# Run for each feature level -----------------------------------
### Hold results in list
glmList <- list()         # DGEGLMs
fitList <- list()         # fits
ResList_Full <- list()    # full results
resList <- list()         # results table
modCheck <- list()        # model/fit stats

### Loop
for (i in names(dgeList)) {
  feature_level <- gsub(".*_", "", i)
  # Subset dge
  y <- dgeList[[i]]
  
  # TMM normalization
  y <- calcNormFactors(y)
  
  # Model
  design <- model.matrix(~Age + Sex + qSV1 + qSV2 + Dx, y$samples)
  
  # Run
  y <- estimateDisp(y, design, verbose = TRUE, robust = TRUE) 
  glmList[[i]] <- y
  fit <- glmFit(y, design, robust = TRUE)
  fitList[[i]] <- fit
  ResList_Full[[i]] <- glmLRT(fit, coef = "DxSZ")
  res <- data.frame(topTags(ResList_Full[[i]], n = Inf))
  
  # BCV (Biological Coefficient of Variation) test
  # Dispersion means biological coeffient of variation (BCV) squared 
  # e.g. if expression typically differs by 20% its BCV is 0.2
  # and its dispersion is 0.04. 
  # edgeR estimates dispersion using quantile-adjusted conditional max likelyhood method (qCML)
  # Common dispersion = common dispersion value for all features
  # Tagwise = feature-specific dispersions, using empirical Bayes strategy to squeeze original feature-wise dispersions towards global, abundance-dependent trend 
  # Typical BCV values; BCV = squareroot(common dispersion) 
  # Unrelated samples (e.g. humans) = 40%
  # Related samples (e.g. mouse studies) = 10%
  modCheck[[i]]$common_dispersion <- glmList[[i]]$common.dispersion
  modCheck[[i]]$BCV <- sqrt(modCheck[[i]]$common_dispersion)
  pdf(paste0("graphics/edgeR/bcv_plot_", i, ".pdf"))
  plotBCV(
    glmList[[i]], 
    xlab = "Average log CPM", 
    ylab = "Biological coefficient of variation",
    pch = 16, cex = 0.2, 
    col.common = "red", 
    col.trend = "blue", 
    col.tagwise = "black"
  ); dev.off
  
  # Save
  saveRDS(y, paste0("output/edgeR/edgeR_ACC_SZvsCTRL_DGEGLM_", feature_level, ".rds"))
  
  # Return
  resList[[i]] <- res
}

### Save
saveRDS(glmList, "output/edgeR/edgeR_ACC_SZvsCTRL_dgeGLMs.rds"); beepr::beep()
saveRDS(fitList, "output/edgeR/edgeR_ACC_SZvsCTRL_dgeFits.rds"); beepr::beep()
saveRDS(ResList_Full, "output/edgeR/edgeR_ACC_SZvsCTRL_resList_FULL.rds"); beepr::beep()
saveRDS(resList, "output/edgeR/edgeR_ACC_SZvsCTRL_resList.rds"); beepr::beep()

### Write out BCV vaues
modCheck %>%
  bind_rows(.id = "dataset") %>%
  separate(dataset, into = c("brain_region", "feature_level"), remove = FALSE) %>%
  write_csv("logs/edgeR/edgeR_ACC_SZvsCTRL_BCV_values.csv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))