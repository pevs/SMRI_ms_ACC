### Stanley ACC SZ vs CTRL
### Results overlap
### Alexis Norris
### Created: 2018-11-05
### Modified: 2019-01-17


# Input files --------------------------------------------------
### ACC significant DEGs
### Start with one that already has cell type specificity data (ewce_celltypes.R)
#sig <- read_tsv(paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR.tsv"))

### Load posthoc test results (from posthoc_ACC_covs.R, posthoc_brainregions.R)
### Don't use read_tsv because it messes up formatting
#posthoc <- read.delim("output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_SZonly_res.tsv", stringsAsFactors = FALSE)
#posthocList <- list(
#  "PFC" = read.delim("output/edgeR/posthoc/edgeR_PFC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE),
#  "HPC" = read.delim("output/edgeR/posthoc/edgeR_HPC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE)
#)

### TaqMan validation results
#valid <- read_csv("input/TaqManValidation_072915_IssyNewkirk.txt")

### qSVA degradation-associated regions' genes (from anno_qsva.R)
#degrad <- read_tsv("input/qsva/qSVA_degradation_regions_anno.tsv")

### Other genesets
#zhao <- read_tsv("input/genesets/Zhao_ACC/Zhao_ACC_TableS4_gene_symbols.tsv")
#darby <- read_tsv("input/genesets/Darby_HPC/Darby_HPC_Table2_gene_symbols.tsv")
#df <- read_tsv("input/genesets/PGC_Ripke_2014/PGC_108loci_Ripke_2014_TableS3_gene_symbols.tsv")  # from the paper's TableS3
#pgcList <- read_tsv("input/genesets/PGC_scz2.regions/scz2_genes.rds") # from the PGC website (SZ1 matches TableS3 results)
#df <- read_tsv("input/genesets/PheGenI_Association_gene_symbols.tsv")

### Significant GSEA genesets
#oxphos <- read_tsv("input/genesets/Pathways/KEGG_hsa00190.txt") 


# Setup --------------------------------------------------------
### Packages
library(tidyverse)          # wrangling; plotting
library(janitor)            # data cleaning

### Parameters
analysis_name <- "overlaps"
fdr_cutoff <- 0.05


# Load ACC significant results ---------------------------------
sig <- read_tsv(paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR.tsv")) %>%
  # Remove empty columns
  # (transcript_id, mirbase, snornabase, horde_id, imgt, intermediate_filament_db)
  remove_empty(which = "cols")

# Add EWCE cell type specificity -------------------------------
### Load cell type specificty (from ewce_celltypes.R)
specList <- readRDS("output/ewce/ewce_gene_specificityList.rds")

### Add cell type specificty to DEGs table
sig <- sig %>%
  left_join(specList$darmanis, by = "gene_symbol") %>%
  left_join(specList$lake, by = "gene_symbol")


# Create list for all genesets ---------------------------------
### For all genesets, symbols converted to current HGNC approved symbol
### Create list to hold overlap
genesets <- list()


# Add posthoc antipsychotics results ---------------------------
### Get edgeR results for DE *features* ACC SZvsCTRL

### Load edgeR results (from posthoc_ACC_covs.R)
### Don't use read_tsv because it messes up formatting
df <- read.delim("output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_res.tsv", stringsAsFactors = FALSE)

### Subset
df <- df %>%
  # ACC SZvsCTRL DE features
  filter(feature_id %in% sig$feature_id) %>%
  
  # Just edgeR result columns
  select(gene_symbol, feature_id, logCPM, logFC, LR, PValue, FDR) 

### Add significant genes to list of genelists
genesets$ACC_Antipsychotics_05FDR_all <- filter(df, FDR < fdr_cutoff)
genesets$ACC_Antipsychotics_05FDR_down <- filter(df, FDR < fdr_cutoff, logFC < 0)
genesets$ACC_Antipsychotics_05FDR_up <- filter(df, FDR < fdr_cutoff, logFC > 0)

### Check
stopifnot(!duplicated(df$feature_id))

### Add to DEG table
sig <- sig %>%
  left_join(df, by = c("feature_id", "gene_symbol"), suffix = c("", "_AntipsychoticsKG"))


# Add posthoc brain regions results ----------------------------
### Load edgeR results (from posthoc_brainregions.R)
### Don't use read_tsv because it messes up formatting
posthocList <- list(
  "PFC" = read.delim("output/edgeR/posthoc/edgeR_PFC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE),
  "HPC" = read.delim("output/edgeR/posthoc/edgeR_HPC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE)
)

### Add posthoc results to DEG table
for (i in names(posthocList)) {
  ### Subset
  df <- posthocList[[i]] %>%
    select(direction, gene_symbol, feature_level, feature_id, logCPM, logFC, LR, PValue, FDR) 
  
  ### Add down and up separately (separate cols) - alternatively could just keep best p-value, but it's possible that one gene has alternative transcripts
  down <- df %>%
    filter(direction == "down") %>%
    select(-direction)
  up <- df %>%
    filter(direction == "up") %>%
    select(-direction)
  
  ### Add significant genes to list of genelists, plus directional
  genesets[[paste0(i, "_SZ_05FDR_all")]] <- filter(df, FDR < fdr_cutoff)
  genesets[[paste0(i, "_SZ_05FDR_down")]] <- filter(down, FDR < fdr_cutoff)
  genesets[[paste0(i, "_SZ_05FDR_up")]] <- filter(up, FDR < fdr_cutoff)

  ### Add to DEG table
    # Best down-regulated feature
  stopifnot(!duplicated(down$gene_symbol))
  sig <- sig %>%
    left_join(down, by = "gene_symbol", suffix = c("", paste0("_", i, "_down")))
    # Best up-regulated feature
  stopifnot(!duplicated(up$gene_symbol))
  sig <- sig %>%
    left_join(up, by = "gene_symbol", suffix = c("", paste0("_", i, "_up")))
  
  ### Check
  stopifnot(!duplicated(sig$feature_id))
}


# TaqMan qPCR validation results -------------------------------
### Load
df <- read_tsv("input/TaqManValidation_072915_IssyNewkirk.txt") %>%
  # Remove BP results
  select(gene_symbol, TaqMan_probe_ID, TaqMan_SZ_FC, TaqMan_SZ_p)

### Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 

### Add significant genes to list of genelists
genesets$TaqMan_Validation_p05 <- filter(df, TaqMan_SZ_p < 0.05)


# Overlap with other genesets ----------------------------------
### qSVA degradation-associated regions' genes (for polyA data; from anno_qsva.R)
  # Load
df <- read_tsv("input/qsva/qSVA_degradation_regions_anno.tsv") %>%
  filter(
    !is.na(gene_symbol),
    distToGene == 0,
    method == "polyA"
  ) %>%
  select(gene_symbol, region) %>%
  # collapse multiple degrad regions for the same gene
  group_by(gene_symbol) %>%
  summarise(
    qSVA_polyA_n = n(),
    qSVA_polyA_pos_hg19 = paste(region, collapse = "; ")
  ) %>%
  ungroup()
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genesets
genesets$qSVA_polyA_degradation_regions <- df

### Zhao, et al Mol Psych 2015 - TableS4 (from anno_Zhao_ACC.R); FDR < 0.25 for SZ or BP
  # Load
df <- read_tsv("input/genesets/Zhao_ACC/Zhao_ACC_TableS4_gene_symbols.tsv") %>%
  mutate(
    Zhao_ACC_SZ_p = SCZ_p,
    Zhao_ACC_SZ_q = SCZ_q
  ) %>%
  select(gene_symbol, Zhao_ACC_SZ_p, Zhao_ACC_SZ_q)
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genesets
  # Cutoff for TableS4 is 0.25FDR, so do subsets for 0.05 (FDR I use), 0.25 (theirs) and 0.10 (middle ground)
genesets$Zhao_ACC_SZ_05q <- df %>%
  filter(Zhao_ACC_SZ_q < 0.05)
genesets$Zhao_ACC_SZ_10q <- df %>%
  filter(Zhao_ACC_SZ_q < 0.10)
genesets$Zhao_ACC_SZ_20q <- df %>%
  filter(Zhao_ACC_SZ_q < 0.20)

### Darby, et al Transl Psych 2016 DEGs - Table2 (from anno_Darby_HPC.R); FDR < 0.1 for SZ or BP
  # Load
df <- read_tsv("input/genesets/Darby_HPC/Darby_HPC_Table2_gene_symbols.tsv") %>%
  mutate(
    Darby_HPC_SZ_logFC = SZ_logFC, 
    Darby_HPC_SZ_lfcSE = SZ_lfcSE, 
    Darby_HPC_SZ_Walk = SZ_Walk, 
    Darby_HPC_SZ_MTCP = SZ_MTCP 
  ) %>%
  select(gene_symbol, Darby_HPC_SZ_logFC, Darby_HPC_SZ_lfcSE, Darby_HPC_SZ_Walk, Darby_HPC_SZ_MTCP)
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genelists
  # Cutoff for TableS4 is 0.10FDR, so do subsets for 0.05 (FDR I use), 0.10 (theirs)
genesets$Darby_HPC_SZ_05FDR_all <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.05)
genesets$Darby_HPC_SZ_10FDR_all <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.10)
genesets$Darby_HPC_SZ_05FDR_down <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.05, Darby_HPC_SZ_logFC < 0)
genesets$Darby_HPC_SZ_10FDR_down <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.10, Darby_HPC_SZ_logFC < 0)
genesets$Darby_HPC_SZ_05FDR_up <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.05, Darby_HPC_SZ_logFC > 0)
genesets$Darby_HPC_SZ_10FDR_up <- df %>%
  filter(Darby_HPC_SZ_MTCP < 0.10, Darby_HPC_SZ_logFC > 0)


### PGC SZ risk loci - from original paper's TableS3 (from anno_PGC_108loci_Ripke2014.R)
df <- read_tsv("input/genesets/PGC_Ripke_2014/PGC_108loci_Ripke_2014_TableS3_gene_symbols.tsv") %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(
    PGC_Ripke2014_locus_pos_hg19 = locus_pos_hg19,
    PGC_Ripke2014_locus_rank = locus_rank,
    PGC_Ripke2014_locus_pvalue = locus_pvalue,
    PGC_Ripke2014_locus_genes_fromTableS3 = locus_genes
  ) %>%
  select(gene_symbol, PGC_Ripke2014_locus_pos_hg19, PGC_Ripke2014_locus_rank, PGC_Ripke2014_locus_pvalue, PGC_Ripke2014_locus_genes_fromTableS3)
# Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
# Add to genesets
genesets$PGC_Ripke2014_108Loci <- df


### PGC SZ risk loci - both original 108 and newer 128 datasets (from anno_PGC_SZ2_loci.R)
  # Load both 108 and 128 loci datasets
pgcList <- readRDS("input/genesets/PGC_scz2.regions/scz2_genes.rds")
  # 108 loci
df <- pgcList$scz1 %>%
  mutate(PGC_SZ2_108loci = loci_pos_hg19) %>%
  filter(!is.na(gene_symbol)) %>%
  select(gene_symbol, PGC_SZ2_108loci)
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genesets
genesets$PGC_SZ2_108Loci <- df
  # 128 loci
df <- pgcList$scz2 %>%
  mutate(PGC_SZ2_128loci = loci_pos_hg19) %>%
  filter(!is.na(gene_symbol)) %>%
  select(gene_symbol, PGC_SZ2_128loci)
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genesets
genesets$PGC_SZ2_128Loci <- df

### PheGenI - GWAS studies of SZ
### from https://www.ncbi.nlm.nih.gov/gap/phegeni (downloaded 2019-01-14)
  # Load
df <- read_tsv("input/genesets/PheGenI_Association_gene_symbols.tsv") %>%
  select(gene_symbol, PheGenI_n, PheGenI_evidence)
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol")
  
  # Add to genesets
genesets$PheGenI_SZ_GWAS <- df

### Oxidative phosphorylation KEGG pathway (hsa00190)
### Downloaded from BROAD GSEA 01-19-2019
  # Load
df <- read_tsv("input/genesets/Pathways/KEGG_hsa00190.txt") %>%
  mutate(OxPhos_kegg_hsa00190 = "YES")
  # Add to DEG table
stopifnot(!duplicated(df$gene_symbol))
sig <- sig %>%
  left_join(df, by = "gene_symbol") 
  # Add to genesets
genesets$KEGG_OxPhos_hsa00190 <- df


# Summarise DEG overlaps ---------------------------------------
### Separate down- and up-regulated
sig_down <- sort(filter(sig, direction == "down")$gene_symbol)
sig_up <- sort(filter(sig, direction == "up")$gene_symbol)

### Add ACC DEGs
genesets$ACC_SZ_05FDR_all <- sig
genesets$ACC_SZ_05FDR_down <- filter(sig, direction == "down")
genesets$ACC_SZ_05FDR_up <- filter(sig, direction == "up")

### Dataframe to hold info
overlaps <- data_frame(
  "geneset_name" = rep(NA, length(names(genesets))),
  "geneset_n" = rep(0, length(names(genesets))),
  "geneset_genes" = rep(NA, length(names(genesets))),
  "overlap_all_n" = rep(0, length(names(genesets))),
  "overlap_down_n" = rep(0, length(names(genesets))),
  "overlap_down_genes" = rep(NA, length(names(genesets))),
  "overlap_up_n" = rep(0, length(names(genesets))),
  "overlap_up_genes" = rep(NA, length(names(genesets)))
)

### Add overlap of nums and genes to one df
for (i in 1:length(names(genesets))) {
  ### Print geneset info
  overlaps$geneset_name[i] <- names(genesets)[i]
  genes <- sort(unique(genesets[[i]]$gene_symbol))
  
  ### Overlap with ACC DEGs (SZvsCTRL)
  if (length(genes) != 0) {
    overlaps$geneset_n[i] <- length(genes)
    overlaps$geneset_genes[i] <- toString(sort(genes))
    
    ### Down-regulated DEGs
    down <- genes[which(genes %in% sig_down)]
    if (length(down) != 0) {
      overlaps$overlap_down_n[i] <- length(down)
      overlaps$overlap_down_genes[i] <- toString(sort(down))
    }
    
    ### Up-regulated DEGs
    up <- genes[which(genes %in% sig_up)]
    if (length(up) != 0) {
      overlaps$overlap_up_n[i] <- length(up)
      overlaps$overlap_up_genes[i] <- toString(sort(up))
    }
    
    overlaps$overlap_all_n[i] <- length(c(down, up))
  }
}

### Save
write_tsv(overlaps, paste0("output/edger/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR_overlaps_Summary.tsv"))

### Add overlaps as a column in DEG table
genes_setsList <- list()
for (gene in unique(sig$gene_symbol)) {
  hitsList <- list()
  # Remove the ACC DEG results first
  for(j in sort(names(genesets)[-grep("ACC_SZ_05FDR", names(genesets))])) {
    hitsList[[j]] <- ifelse(gene %in% genesets[[j]]$gene_symbol, j, NA)
  } 
  hits <- unlist(hitsList)[!is.na(unlist(hitsList))]
  genes_setsList[[gene]] <- data_frame(
    "geneset_overlaps_n" = length(hits),
    "geneset_overlaps" = paste(hits, collapse = "; ")
  )
  
}
genes_sets <- genes_setsList %>% bind_rows(.id = "gene_symbol")
genes_sets$geneset_overlaps[genes_sets$geneset_overlaps == ""] <- NA
sig <- sig %>%
  left_join(genes_sets, by = "gene_symbol")


# Save ---------------------------------------------------------
### Clean up
stopifnot(!duplicated(sig$feature_id))
stopifnot(!is.na(sig$feature_id))
sig_clean <- sig %>%
  # Reorder; remove unnecessary columns
  select(
    direction, gene_symbol, feature_level, feature_id, pos_hg19, 
    logCPM, logFC, LR, PValue, FDR,
    geneset_overlaps_n, geneset_overlaps,
    one_of(
      paste0(c("logCPM", "logFC", "LR", "PValue", "FDR"), "_AntipsychoticsKG"),
      paste0(c("logCPM", "logFC", "LR", "PValue", "FDR", "feature_level", "feature_id"), "_PFC_down"),
      paste0(c("logCPM", "logFC", "LR", "PValue", "FDR", "feature_level", "feature_id"), "_PFC_up"),
      paste0(c("logCPM", "logFC", "LR", "PValue", "FDR", "feature_level", "feature_id"), "_HPC_down"),
      paste0(c("logCPM", "logFC", "LR", "PValue", "FDR", "feature_level", "feature_id"), "_HPC_up"),
      paste0("TaqMan_", c("probe_ID", "SZ_FC", "SZ_p")),
      paste0("qSVA_polyA_", c("n", "pos_hg19")),
      paste0("Zhao_ACC_SZ_", c("p", "q")),
      paste0("Darby_HPC_SZ_", c("logFC", "lfcSE", "Walk", "MTCP")),
      paste0("PGC_Ripke2014_locus_", c("pos_hg19", "rank", "pvalue", "genes_fromTableS3")),
      paste0("PGC_SZ2_", c("108loci", "128loci")),
      paste0("PheGenI_", c("n", "evidence"))
    ),
    OxPhos_kegg_hsa00190,
    Neurons, Oligodendrocytes, Astrocytes, Endothelial, Microglia, 
    one_of(
      paste0("In", 1:8),
      paste0("Ex", 1:8)
    ), 
    everything()
  ) %>%
  ### Unnecessary columns
  select(-one_of(c(
    "mult_matches", "gene_symbol_corrected_reason", "gene_symbol_uncorrected",
    "data_set", "brain_region"
    )))


### Save
stopifnot(!duplicated(sig_clean$feature_id))
stopifnot(!is.na(sig_clean$feature_id))
write_tsv(sig_clean, paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR_overlaps.tsv"))


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
