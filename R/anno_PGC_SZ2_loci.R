### Stanley ACC SZ vs CTRL
### Get genes in PGC2 loci
### and update symbols
### Alexis Norris
### Created: 2019-01-02
### Modified: 2019-01-16


# Input files --------------------------------------------------
### Source of PGC data: https://www.med.unc.edu/pgc/results-and-downloads/downloads
#pgc_pos <- list(
#  # 108 loci (SCZ1)
#  "SCZ1" = read_tsv("input/anno/scz2.regions/scz2.anneal.108.txt"),
#  # 128 loci (SCZ2)
#  "SCZ2" = pgc_pos2 <- read_tsv("input/anno/scz2.regions/scz2.rep.128.txt")
#)   

### geneMap 
#geneMap <- readRDS("input/genesets/geneMap_GRCh37.p13_gencode.v19.rds")

### HGNC approved gene symbol check (from anno_hgnc.R)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling; warning: GenomicRanges masks dplyr::select
library(GenomicRanges)     # annotate DEGs with overlapping repeatMasker & all HGNC genes

### Parameters
analysis_name <- "anno_PGC_SZ2_loci"


# Load PGC positions -------------------------------------------
### Load
### Source of PGC data: https://www.med.unc.edu/pgc/results-and-downloads/downloads
pgc_pos <- list(
  # 108 loci (SCZ1)
  "scz1" = read_tsv("input/genesets/PGC_scz2.regions/scz2.anneal.108.txt"),
  # 128 loci (SCZ2)
  "scz2" = pgc_pos2 <- read_tsv("input/genesets/PGC_scz2.regions/scz2.rep.128.txt")
) 

### Make colnames uniform
pgc_pos$scz2$chr <- NULL
names(pgc_pos$scz1) <- recode(
  names(pgc_pos$scz1),
  "hg19chrc" = "chr", "anneal1" = "start", "anneal2" = "end", "spananneal" = "loci_size_bp",
  "anneal.rank" = "loci_rank", "bestsnp" = "snp_best", "pmin" = "p_min"
)
names(pgc_pos$scz2) <- recode(
  names(pgc_pos$scz2),
  "hg19chrc" = "chr", "six1" = "start", "six2" = "end", "spananneal" = "loci_size_bp",
  "region.rank" = "loci_rank", "snpid" = "snp", "bp" = "snp_pos" 
)

### Add pos
pgc_pos <- lapply(pgc_pos, function (df) {
  df$pos_hg19 <- paste0(df$chr, ":", df$start, "-", df$end)
  df
})

### Save original version
saveRDS(pgc_pos, "input/genesets/PGC_scz2.regions/scz2_positions_hg19.rds")


# Get genes ----------------------------------------------------
### Get all genes in these loci

### Load geneMap 
geneMap <- readRDS("input/anno/geneMap_GRCh37.p13_gencode.v19.rds")

### Get genes
pgc_genes <- list()
for (i in names(pgc_pos)) {
  # Subset
  df <- pgc_pos[[i]]
  
  # Convert to GRanges object
  loci <- makeGRangesFromDataFrame(
    df, 
    keep.extra.columns = TRUE,
    ignore.strand = TRUE
  )
  
  # Get all genes
  hits <- findOverlaps(loci, geneMap, ignore.strand = TRUE, select = "all")
  
  ### Get overlaps
  hits <- findOverlaps(                     
    geneMap, loci,
    ignore.strand = TRUE,         
    type = "any", select = "all"
  )
  overlaps <- pintersect(
    geneMap[queryHits(hits)], 
    loci[subjectHits(hits)]
  )
  
   # Combine geneMap anno and loci data 
  mcols(overlaps) <- c(
    mcols(overlaps),
    mcols(loci[subjectHits(hits)])
  )
  overlaps$loci_pos_hg19 <- paste0(
    seqnames(loci[subjectHits(hits)]), ":",
    start(loci[subjectHits(hits)]), "-", end(loci[subjectHits(hits)])
  )

  # Clean up
  overlaps_df <- overlaps %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(-one_of(c("hit", "pos_hg19"))) 
  names(overlaps_df) <- recode(
    names(overlaps_df),
    "seqnames" = "chr"
  )
  # Return
  pgc_genes[[i]] <- overlaps_df
}


# Update outdated HGNC gene symbols ----------------------------
### Check that gene_symbol is the current approved HGNC symbol
### if outdated (alias), update
### If ensembl ID not in HGNC, remove gene_symbol (make NA)

### Load HGNC anno (from anno_hgnc.R)
hgnc <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")$by_ensembl
names(hgnc) <- recode(
  names(hgnc),
  "symbol" = "approved_symbol"
)

### Get current approved symbols
pgc_genes_updated <- list()
annoStats <- list()
for (i in names(pgc_genes)) {
  ### Subset
  original <- pgc_genes[[i]]
  
  ### Conform naming
  names(original) <- recode(names(original), "hgnc_symbol" = "gene_symbol")

  ### Update gene symbols
  df <- original %>%
    mutate(
      # column to match on; need to remove version suffix
      by_col = gsub("\\..*", "", gene_id),
      # write out original
      gene_symbol_original = gene_symbol
    ) %>%
    # Fix classes, for binding later
    mutate_if(is.factor, as.character)
  
  ### Check that all have ensembl ID
  stopifnot(!is.na(df$gene_id))

  ### Join by ensembl IDs
  df_hgnc <- df %>%
    left_join(hgnc, by = c("by_col" = "ensembl_gene_id")) %>%
    # Fix classes, for binding later
    mutate_if(is.factor, as.character)
  
  ### Outcome 1. Ok (current approved symbol)
  ok <- df_hgnc %>%
    # symbols match
    filter(gene_symbol == approved_symbol) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "ok")
  
  ### Outcome 2. Updated gene symbol (previous/alias)
  updated <- df_hgnc %>%
    # symbols don't match
    filter(gene_symbol != approved_symbol) %>%
    # Replace with approved
    mutate(gene_symbol = approved_symbol)  %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "updated")
  
  ### Outcome 3. old (not in HGNC)
  not_in_hgnc <- df_hgnc %>%
    # ensembl ID is NOT in HGNC
    filter(is.na(approved_symbol)) %>%
    # Replace with NA
    mutate(gene_symbol = NA) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "not in HGNC")
  
  ### Combine outcomes
  new_anno <- bind_rows(
    ok,
    updated,
    not_in_hgnc
  ) 
  
  ### Outcome 4. Gene/feature doesn't have ensembl ID
  if (is.null(is.na(df$gene_id))) {
    df %>%
      filter(is.na(by_col)) %>%
      # Replace with NA
      mutate(gene_symbol = NA) %>%
      # Note that it was changed
      mutate(gene_symbol_corrected_reason = "no ensembl ID") %>% 
      bind_rows(df_anno)
  }
  
  new_anno <- new_anno %>%
    dplyr::select(
      gene_symbol, 
      gene_symbol_original, gene_symbol_corrected_reason, 
      mult_matches, everything()
    ) 
  
  ### Check
  new_anno[new_anno == ""] <- NA
  new_anno[new_anno == "NA"] <- NA
  stopifnot(nrow(new_anno) == nrow(original))
  
  ### Summarize conversions
  new_anno$gene_symbol_corrected_reason <- fct_relevel(
    new_anno$gene_symbol_corrected_reason,
    "ok", "updated", "not in HGNC"
  )
  annoStats[[i]] <- new_anno %>%
    group_by(gene_symbol_corrected_reason) %>%
    summarise(
      n = n(),
      original = toString(unique(gene_symbol_original[!is.na(gene_symbol_original)])),
      updated_symbols = toString(unique(gene_symbol[!is.na(gene_symbol)]))
    ) 
  
  ### Return
  pgc_genes_updated[[i]] <- new_anno 
}


# Collapse genes -----------------------------------------------
### If one gene hit for multiple regions/SNPs
### (aka if multiple SNPs-evidence for a gene) -- similar to code used for anno_pheGenI.R
pgc_genes_updated <- lapply(pgc_genes_updated, function (df_full) {
  ### Get number of duplicates (if any)
  df_notNA <- df_full %>%
    filter(!is.na(gene_symbol))
  gene_dups <- df_notNA$gene_symbol[duplicated(df_notNA$gene_symbol)]
  
  ### If no duplicates
  if (length(gene_dups) == 0) {
    return(df_full)
  }
    
  ### If there ARE duplicates
  if (length(gene_dups) > 0) {
  uniques_char <- df_notNA %>%
    filter(gene_symbol %in% gene_dups) %>%
    group_by(gene_symbol) %>%
    summarise_if(is.character, function(x) paste(unique(as.character(x[!is.na(x)])), collapse = "; "))
  uniques_num <- df_notNA %>%
    filter(gene_symbol %in% gene_dups) %>%
    group_by(gene_symbol) %>%
    summarise_if(is.numeric, sum) 
  uniques <- inner_join(uniques_num, uniques_char, by = "gene_symbol")
  df_notNA_collapsed <- df_notNA %>%
    filter(!(gene_symbol %in% gene_dups)) %>%
    bind_rows(uniques)
  df_no_dups <- df_full %>%
    filter(is.na(gene_symbol)) %>%
    bind_rows(df_notNA_collapsed)
  stopifnot(!duplicated(df_no_dups$gene_symbol[!is.na(df_no_dups$gene_symbol)]))
  df_no_dups
  }
})

for (i in names(pgc_genes_updated)) {
  x <- pgc_genes_updated[[i]]$gene_symbol
  stopifnot(!duplicated(x[!is.na(x)]))
}

# Save ---------------------------------------------------------
saveRDS(pgc_genes_updated, "input/genesets/PGC_scz2.regions/scz2_genes.rds")
annoStats %>%
  bind_rows(.id = "genelist") %>%
  write_tsv("input/genesets/PGC_scz2.regions/gene_symbols_updateStats.tsv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
