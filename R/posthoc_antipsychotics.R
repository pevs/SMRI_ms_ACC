### Stanley ACC SZ vs CTRL
### edgeR: posthoc
### ACC DEGs for antipsychotics
### Alexis Norris
### Created: 2018-12-05
### Modified: 2019-01-16


# Input files --------------------------------------------------
#dgeList <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")
#acc <- read_tsv("output/edgeR/edgeR_ACC_SZvsCTRL_res.tsv")


# Setup --------------------------------------------------------
### Packages
library(edgeR)          # for DGE object
library(tidyverse)      # wrangling

### Parameters
analysis_name <- "posthoc_antipsychotics"


# Load DGEs ----------------------------------------------------
### datasets for brain regions: ACC, PFC, HPC
### feature levels: genes, exons, jxns, regions
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered_outliersRemoved_qSVA.rds")

### ACC only
dgeList <- dgeList_all[grep("ACC", names(dgeList_all))]


# Load ACC features --------------------------------------------
### Best feature used for each gene
acc <- read_tsv("output/edgeR/edgeR_ACC_SZvsCTRL_res.tsv")
acc_down <- filter(acc, logFC < 0)
acc_up <- filter(acc, logFC >= 0)


# Lifetime Antipsychotics --------------------------------------
edgeRList <- list()
for (i in names(dgeList)) {
  print(i)
  
  # Subset DGE
  y <- dgeList[[i]]
  
  # Check that NAs are 0s
  stopifnot(!is.na(y$samples$LifetimeAntipsychotics))
  
  # Use per kg instead of mg to make logFC numbers more readable
  y$samples$LifetimeAntipsychotics_kg <- y$samples$LifetimeAntipsychotics/1E06
  
  # Model
  design <- model.matrix(~Age + Sex + qSV1 + qSV2 + Dx + LifetimeAntipsychotics_kg, y$samples)
  
  # Fit
  y <- estimateDisp(y, design, verbose = TRUE, robust = TRUE) 
  fit <- glmFit(y, design, robust = TRUE)
  
  # Run
  edgeRList[[i]] <- data.frame(topTags(glmLRT(fit), n = Inf))
}

### Save full results
saveRDS(edgeRList, "output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_resList.rds")


### Extract posthoc results:
### Subset for features used in ACC SZvsCTRL
res_all <- edgeRList %>%
  bind_rows(.id = "feature_level") %>%
  filter(!is.na(gene_symbol))
resList <-  list(
  "down" = filter(res_all, feature_id %in% acc_down$feature_id),
  "up" = filter(res_all, feature_id %in% acc_up$feature_id)
)
lapply(resList, function (df) {
  df %>%
    arrange(PValue) %>%
    .[!duplicated(.$gene_symbol), ]
}) %>%
  bind_rows(.id = "direction") %>%
  mutate(FDR = p.adjust(.$PValue, method = "fdr")) %>%
  write_tsv("output/edger/posthoc/edgeR_ACC_AntipsychoticsKG_res.tsv") # contains dups, since allowed both up and down


# Lifetime Antipsychotics (SZ only) ----------------------------
edgeRList <- list()
for (i in names(dgeList)) {
  print(i)
  
  # Subset DGE
  y <- dgeList[[i]]

  # Subset samples
  y <- y[ , which(y$samples$Dx == "SZ"), keep.lib.sizes = TRUE]
  
  # Remove any unexpressed features
  print(table(rowSums(cpm(y)) > 0))
  y <- y[rowSums(cpm(y)) > 0, , keep.lib.sizes = TRUE]
  
  # Check that NAs are 0s
  stopifnot(!is.na(y$samples$LifetimeAntipsychotics))
  
  # Use per kg instead of mg to make logFC numbers more readable
  y$samples$LifetimeAntipsychotics_kg <- y$samples$LifetimeAntipsychotics/1E06
  
  # Model
  design <- model.matrix(~Age + Sex + qSV1 + qSV2 + LifetimeAntipsychotics_kg, y$samples)

  # Fit
  y <- estimateDisp(y, design, verbose = TRUE, robust = TRUE) 
  fit <- glmFit(y, design, robust = TRUE)
  
  # Run
  edgeRList[[i]] <- data.frame(topTags(glmLRT(fit), n = Inf))
}

### Save full results
saveRDS(edgeRList, "output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_SZonly_resList.rds")


### Extract posthoc results:
### Subset for features used in ACC SZvsCTRL
res_all <- edgeRList %>%
  bind_rows(.id = "feature_level") %>%
  filter(!is.na(gene_symbol))
resList <-  list(
  "down" = filter(res_all, feature_id %in% acc_down$feature_id),
  "up" = filter(res_all, feature_id %in% acc_up$feature_id)
)
lapply(resList, function (df) {
  df %>%
    arrange(PValue) %>%
    .[!duplicated(.$gene_symbol), ]
}) %>%
  bind_rows(.id = "direction") %>%
  mutate(FDR = p.adjust(.$PValue, method = "fdr")) %>%
  write_tsv("output/edger/posthoc/edgeR_ACC_AntipsychoticsKG_SZonly_res.tsv") # contains dups, since allowed both up and down


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
