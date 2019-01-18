### Stanley ACC SZ vs CTRL
### Cell type enrichment analysis
### Method: EWCE
### Cell data (single cell transcriptome): 
###   # Darmanis et al PNAS 2015 (human brain cell types; fresh tissue from surgery; multiple individuals)
###   # Lake et al Science 2016 (human neuron subtypes; postmortem; 1 individual [female])
### Alexis Norris
### Created: 2018-11-29
### Modified: 2019-01-07


# Note ---------------------------------------------------------
### Darmanis: cell types from anterior temporal lobe (BA15)
### Lake: cell types from multiple cortical regions
# BA8 = FC   
# BA10 = anterior PFC 
# BA17 = V1 (primary visual cortex)
# BA21 = mTC (middle temporal gyrus)
# BA22 = sTC (superior temporal gyrus)
# BA41 = aTC (auditory cortex)

# Input files --------------------------------------------------
#geneList_all <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_genelist.rds")
#celldataList[[i]] <- readRDS(paste0("input/celldata/ewce_celldata_", datasets[[i]], "_processed_informative.rds"))


# Setup --------------------------------------------------
### Packages
library(tidyverse)
library(EWCE)

### Parameters
analysis_name <- "ewce_celltypes"
reps <- 20000               # num bootstraps used for ewce


# Load celldata ------------------------------------------------
### From celldata_darmanis.R and celldata_Lake.R
cell_files <- list.files(path = "input/celldata/", pattern = "processed_informative")
celldataList <- list()
for (i in cell_files) {
  celldataList[[i]] <- readRDS(paste0("input/celldata/", i))
}

### Cleanup names
datasets <- gsub("darmanis_adult", "darmanis", gsub("ewce_", "", gsub("_processed_informative", "", gsub(".rds", "", cell_files))))
names(celldataList) <- datasets

# Load genelist ------------------------------------------------
### Load Genelist (from genelist.R)
geneList <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_genelist.rds")


# Run enrichment -----------------------------------------------
### Stattests for celltypes
### From vignette/help page (new version ofEWCE pkg)

### Log
overlapStats <- list()
resList <- list()
resTableList <- list()

### Run
for (i in names(celldataList)) for (j in names(geneList)) {
  # Genelist
  genes_all <- geneList[[j]]
  
  # Cell data
  sct_data <- celldataList[[i]]
  annotLevel <- 1 # cell type is level 1 for both datasets
  cell_genes <- rownames(sct_data[[annotLevel]]$mean_exp)
  
  # Genelist: remove genes that are not in cell data
  genes <- genes_all[which(genes_all %in% cell_genes)]
  
  # Summary of genes kept/removed for analysis
  overlapStats[[paste(j, i, sep = "_")]] <- data_frame(
    "Genelist" = j,
    "Celldata" = i,
    "n_genes_total" = length(genes_all),
    "n_genes_in_celldata" = length(genes),
    "genes_in_celldata" = toString(genes) 
  )
  seed <- 1116
  res <- bootstrap.enrichment.test(
    sct_data = sct_data,
    hits = genes,
    bg = cell_genes,
    reps = reps,
    annotLevel = annotLevel,
    sctSpecies = "human",
    genelistSpecies = "human"
  )
  
  # Fix classes (for bind_rows later)
  res$results$CellType <- as.character(res$results$CellType)
  
  # Save full results object
  resList[[paste(j, i, sep = "_")]] <- res
  
  # Save results table separately
  resTableList[[paste(j, i, sep = "_")]] <- res$results
}

### Save full results 
saveRDS(resList, paste0("output/ewce/ewce_res_", reps, "reps.rds"))

### Save results table (collapsed for all tests)
resTableList %>%
  bind_rows(.id = "analysis") %>%
  separate(analysis, into = c("genelist", "celldata"), sep = "_", remove = FALSE) %>%
  # Calculate FDR
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  # Remove duplicate column
  select(everything(), p, FDR) %>%
  write_csv(paste0("output/ewce/ewce_res_", reps, "reps.csv"))

### Write out stats on genes removed 
overlapStats %>%
  bind_rows(.id = "Celldata") %>%
  write_csv("output/ewce/geneList_geneOverlap_stats.csv")


# Gene ~ cell type specificity ---------------------------------
### Extract
specList <- lapply(celldataList, function (celldata) {
  df <- as.data.frame(celldata[[1]]$specificity)
  df$gene_symbol <- rownames(df)
  df
})

### Save
saveRDS(specList, "output/ewce/ewce_gene_specificityList.rds")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))