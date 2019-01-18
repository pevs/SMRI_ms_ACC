### SMRI Array ACC SZ vs CTRL
### GSEA with KEGG
### Alexis Norris
### Created: 2018-05-25
### Modified: 2019-01-13


# Sources ------------------------------------------------------
### clusterProfiler: http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html


# Input files --------------------------------------------------
### ACC significant DEGs
#sig <- read_tsv(paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_edgeR, "FDR.tsv"))


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling; plotting
library(patchwork)         # combining plots
library(clusterProfiler)   # GSEA

### Parameters
analysis_name <- "gsea"
fdr_edgeR <- 0.05
fdr_gsea <- 0.05
nPerm <- 10000


# Prep ranked genelist of ACC DEGs -----------------------------
### ACC significant DEGs
sig <- read_tsv(paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_edgeR, "FDR.tsv"))

### Use entrez IDs as gene identifiers
stopifnot(!duplicated(sig$entrez_id))

### Ranking metric: using logFC*log10(PValue)
sig$rank <- sig$logFC * log10(sig$PValue)

### Handle duplicate rank caused by one signal having two mappings/annotations
### here it is "chr17:79670387-79670605" = both MRPL12 (e994318) and SLC25A10 (e994338)
### Can't have duplicate ranks for GSEA --> so add small value to one
if (length(duplicated(sig$rank) > 0)) {
  # Sort by ascending gene symbol, add 0.0001 to it
  dups <- sig %>%
    filter(rank %in% sig$rank[duplicated(sig$rank)]) %>%
    arrange(gene_symbol)
  dups_first <- dups[!duplicated(dups$rank), ]
  sig$rank <- ifelse(
    sig$feature_id %in% dups_first$feature_id,
    sig$rank + 0.0001, 
    sig$rank
  )
}

### Entrez IDs with ranking metric, sorted by decreasing ranking metric
rankList_entrez <- sig$rank
names(rankList_entrez) <- sig$entrez_id
rankList_entrez <- rankList_entrez[order(rankList_entrez, decreasing = TRUE)]
stopifnot(!duplicated(rankList_entrez))

### Make ranked geneLists
geneList_ranked <- list(
  "all" = rankList_entrez,
  "down" = rankList_entrez[rankList_entrez > 0],
  "up" = rankList_entrez[rankList_entrez < 0]
)
saveRDS(geneList_ranked, "output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_ranked_geneList.rds")


# Geneset enrichment analysis (GSEA) ---------------------------
### Run for each ranked DEG geneList 
resList <- list()
for (i in c("down", "up")) {
  set.seed(1116)
  resList[[i]] <- gseKEGG(
    geneList = geneList_ranked[[i]],
    organism = "hsa",
    keyType = "ncbi-geneid",
    nPerm = nPerm,
    # defaults
    exponent = 1,
    pvalueCutoff = fdr_gsea, pAdjustMethod = "BH",
    minGSSize = 10, maxGSSize = 500, 
    use_internal_data = FALSE,
    seed = TRUE, 
    by = "fgsea",
    verbose = TRUE
  )
}

### Save
saveRDS(resList, "output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_KEGG.rds")

### Extract tables
res_gsea <- lapply(resList, as.data.frame) %>%
  bind_rows(.id = "DEGs")

### Add gene_symbols from entrez IDs
for (i in 1:nrow(res_gsea)) {
  hits <- as.vector(str_split(
    gsub("\\/", ",", res_gsea$core_enrichment[i]), 
    ",", simplify = TRUE
  ))
  res_gsea$geneID_symbols[i] <- paste(sig$gene_symbol[sig$entrez_id %in% hits], collapse = ";")
}

### Save tables
write_csv(res_gsea, "output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_KEGG.csv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))