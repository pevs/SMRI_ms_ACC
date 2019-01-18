### Stanley ACC SZ vs CTRL
### edgeR (ACC) combine
### Create genelist
### Summarize edgeR results
### Alexis Norris
### Created: 2018-11-05
### Modified: 2019-01-06


# Input files --------------------------------------------------
### Load edgeR results
#resList <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_resList.rds")

### qSVA regions (from anno_qsva.R)
#degrad <- read_tsv("input/qsva/qSVA_deg_regions_anno.tsv")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling

### Parameters
analysis_name <- "genelists"
fdr_cutoff <- 0.05


# Load ACC edgeR results ---------------------------------------
### Load edgeR results
resList <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_resList.rds")

### Unlist
res_all <- resList %>%
  bind_rows(.id = "data_set") %>%
  separate(data_set, into = c("brain_region", "feature_level"), remove = FALSE) %>%
  select(-gene_strand) %>% # Same as strand
  # Reorder cols
  select(
    gene_symbol, feature_level, feature_id, 
    logFC, LR, PValue, FDR, logCPM,
    pos_hg19, chr, start, end, strand,
    gene_id, gene_type, 
    exon_id, exon_number,
    transcript_id,
    jxn_class, jxn_fusion, jxn_exonStart_novel, jxn_exonEnd_novel,
    jxn_leftSeq, jxn_rightSeq, jxn_intron_length,
    region_class, region_distToGene, 
    rmsk_overlap, rmsk_repName, rmsk_repClass, rmsk_repFamily, 
    rmsk_pos_hg19, rmsk_strand,
    gene_symbol_uncorrected, approved_symbol, 
    everything()
  )


# Collapse genes -----------------------------------------------
### Collapse genes, take most significant feature for each
### Separate up vs down, in case one transcript is up and another down
resList <-  list(
  "down" = filter(res_all, logFC < 0),
  "up" = filter(res_all, logFC >= 0)
)
res <- lapply(resList, function (df) {
  df %>%
    filter(!is.na(gene_symbol)) %>%
    arrange(PValue) %>%
    .[!duplicated(.$gene_symbol), ]
}) %>%
  bind_rows(.id = "direction") %>%
  mutate(FDR = p.adjust(.$PValue, method = "fdr"))
sig <- res %>%
  filter(FDR < fdr_cutoff)

### Save as tables
write_tsv(res, "output/edger/edgeR_ACC_SZvsCTRL_res.tsv") # contains dups, since allowed both up and down
write_tsv(sig, paste0("output/edger/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR.tsv"))


# Genelist -----------------------------------------------------
### Combine into list -- just genes
geneList <- list(
  "all" = sig$gene_symbol,
  "down" = filter(sig, direction == "down")$gene_symbol,
  "up" = filter(sig, direction == "up")$gene_symbol
)


# Add qSVA degradation regions to Genelist ---------------------
### Add genes in qSVA regions
### (see anno_qsva_degregions.R)
degrad <- read_tsv("input/qsva/qSVA_degradation_regions_anno.tsv") %>%
  filter(
    method == "polyA",
    distToGene == 0,        # must overlap the gene
    !is.na(gene_symbol)     # must have HGNC approved symbol
  )
geneList$degrad <- degrad$gene_symbol


# Save genelist ------------------------------------------------
### Clean up
geneList <- lapply(geneList, function (x) {
  sort(unique(as.character(x)))
})

### Save
saveRDS(geneList, "output/edgeR/edgeR_ACC_SZvsCTRL_genelist.rds")

### Export
write_lines(capture.output(
  sapply(geneList, length),
  sapply(geneList, toString)
), "output/edger/edgeR_ACC_SZvsCTRL_genelist_summary.txt")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
