### Stanley ACC SZ vs CTRL
### PheGenI database for SZ
### Alexis Norris
### Created: 2019-01-14
### Modified: 2019-01-16


# Input files --------------------------------------------------
### PheGenI Schizophrenia (https://www.ncbi.nlm.nih.gov/gap/phegeni; downloaded 2019-01-14)
#pgi <- read_tsv("input/genesets/PheGenI_Association.tab") %>%
  
### HGNC approved gene symbol check (from anno_hgnc.R)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling

### Parameters
analysis_name <- "anno_phegeni"


# Load PheGenI SZ table ----------------------------------------
### PheGenI Schizophrenia (https://www.ncbi.nlm.nih.gov/gap/phegeni; downloaded 2019-01-14)
### Column names were cleaned up manually
pgi <- read_tsv("input/genesets/PheGenI_Association.tab") %>%
  mutate(
    SNP_rsID = paste0("rs", as.character(SNP_rsID)),
    PubMed = paste0("PMID:", as.character(PubMed))
  )

### Separate the SNPs with two mapped genes, to separate rows for each gene
pgi <- bind_rows(
  mutate(pgi, gene_symbol = gene1_symbol),
  mutate(pgi, gene_symbol = gene2_symbol)
) %>%
  distinct() 

### Collapse cases where  multiple genes for 1 SNP 
### Add PheGenI info (SNP id and PubMed ID)
pgi$PheGenI <- ifelse(
  pgi$gene1_symbol == pgi$gene2_symbol,
  paste0(pgi$SNP_rsID, ", ", pgi$gene1_symbol, ", ", pgi$PubMed),
  paste0(pgi$SNP_rsID, ", ", pgi$gene1_symbol, " and ", pgi$gene2_symbol, ", ", pgi$PubMed)
)

### Collapse multiple SNPs-evidence for a gene
pgi <- pgi %>%
  group_by(gene_symbol) %>%
  summarise(
    PheGenI_n = n(),
    PheGenI_evidence = paste(PheGenI, collapse = "; ")
  ) 
  


# Update outdated HGNC gene symbols ----------------------------
### Check that gene_symbol is the current approved HGNC symbol
### if outdated (alias), update
### *using gene symbols* (not ensembl IDs!)

### Load HGNC anno (from anno_hgnc.R)
hgnc_keys <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")
hgnc_symbol <- hgnc_keys$by_symbol
hgnc_prev <- hgnc_keys$previous_to_symbol_noDuplicates
hgnc_alias <- hgnc_keys$alias_to_symbol_noDuplicates
names(hgnc_symbol) <- recode(names(hgnc_symbol), "symbol" = "approved_symbol")
names(hgnc_prev) <- recode(names(hgnc_prev), "symbol" = "approved_symbol")
names(hgnc_alias) <- recode(names(hgnc_alias), "symbol" = "approved_symbol")


### Get current approved symbols
original <- pgi

### Update gene symbols
df <- original %>%
  # column to hold original gene symbol
  mutate(gene_symbol_original = gene_symbol) %>%
  # Fix classes, for binding later
  mutate_if(is.factor, as.character)


### Join by gene symbols
### Outcome 1. Ok (is current approved symbol)
ok <- df %>%
  inner_join(hgnc_symbol, by = c("gene_symbol" = "approved_symbol")) %>%
  # Fix classes, for binding later
  mutate_if(is.factor, as.character) %>%
  mutate(gene_symbol_corrected_reason = "ok")
#stopifnot(ok$mult_matches == "NO") # ok, it's actually one gene symbol that had 2 hgnc IDs (one has been withdrawn)

### Outcome 2. Updated gene symbol (was previous)
updated_prev <- df %>%
  # symbol isn't ok (approved symbol)
  filter(!(gene_symbol %in% ok$gene_symbol_original)) %>%
  # try to match with previous symbol
  inner_join(hgnc_prev, by = c("gene_symbol" = "prev_symbol")) %>%
  # update gene symbol
  mutate(gene_symbol = approved_symbol) %>%
  # Note that it was changed
  mutate(gene_symbol_corrected_reason = "previous symbol") %>%
  # Add full hgnc anno
  dplyr::select(-approved_symbol) %>%
  left_join(dplyr::select(hgnc_symbol, -mult_matches), by = c("gene_symbol" = "approved_symbol"))
# Check if there are multiple genes that previous symbol maps to
# If there are, then note it
if ("YES" %in% updated_prev$mult_matches)  {
  updated_prev <-updated_prev %>%
    filter(mult_matches == "YES") %>%
    mutate(
      gene_symbol_corrected_reason = "previous symbol with multiple approved gene symbols"
    )
}

### Outcome 2. Updated gene symbol (was alias)
updated_alias <- df %>%
  # symbol isn't ok (approved symbol)
  filter(!(gene_symbol %in% c(ok$gene_symbol_original, updated_prev$gene_symbol_original))) %>%
  # try to match with alias
  inner_join(hgnc_alias, by = c("gene_symbol" = "alias_symbol")) %>%
  # update gene symbol
  mutate(gene_symbol = approved_symbol) %>%
  # Note that it was changed
  mutate(gene_symbol_corrected_reason = "alias") %>%
  # Add full hgnc anno
  dplyr::select(-approved_symbol) %>%
  left_join(dplyr::select(hgnc_symbol, -mult_matches), by = c("gene_symbol" = "approved_symbol"))
# Check if there are multiple genes that alias maps to
# If there are, then note it
if ("YES" %in% updated_alias$mult_matches)  {
  updated_alias <-updated_alias %>%
    filter(mult_matches == "YES") %>%
    mutate(
      gene_symbol_corrected_reason = "alias with multiple approved gene symbols"
    )
}

### Outcome 4. old (ensembl gene id not in HGNC )
not_in_hgnc <- df %>%
  filter(!(gene_symbol %in% c(ok$gene_symbol_original, updated_prev$gene_symbol_original, updated_alias$gene_symbol_original))) %>%
  # If there are, then note it and make gene symbol NA
  mutate(
    gene_symbol = NA,
    gene_symbol_corrected_reason = "not in HGNC"
  )

### Combine outcomes
df_hgnc <- bind_rows(
  ok,
  updated_prev,
  updated_alias,
  not_in_hgnc
) %>% 
  dplyr::select(
    gene_symbol, 
    gene_symbol_original, gene_symbol_corrected_reason, 
    mult_matches, everything()
  ) 
stopifnot(nrow(df_hgnc) == nrow(df))

### Remove genes with ambigous or no approved gene symbol
reasons <- c(
  "previous symbol with multiple approved gene symbols",
  "alias with multiple approved gene symbols",
  "not in HGNC"
)
df_hgnc %>%
  filter(gene_symbol_corrected_reason %in% reasons) %>%
  mutate(approved_symbols = gene_symbol) %>%
  select(gene_symbol_original, approved_symbols, gene_symbol_corrected_reason) %>%
  write_csv("input/genesets/PheGenI_Association_genes_removed.csv")
df_hgnc <- df_hgnc %>%
  filter(!(gene_symbol_corrected_reason %in% reasons))


### Collapse multiple SNPs-evidence for a gene
if (length(df_hgnc$gene_symbol) != length(unique((df_hgnc$gene_symbol)))) {
  gene_dups <- df_hgnc$gene_symbol[duplicated(df_hgnc$gene_symbol)]
  uniques_char <- df_hgnc %>%
    filter(gene_symbol %in% gene_dups) %>%
    group_by(gene_symbol) %>%
    summarise_if(is.character, function(x) paste(unique(as.character(x[!is.na(x)])), collapse = "; "))
  uniques_num <- df_hgnc %>%
    filter(gene_symbol %in% gene_dups) %>%
    group_by(gene_symbol) %>%
    summarise_if(is.numeric, sum) 
  uniques <- inner_join(uniques_num, uniques_char, by = "gene_symbol")
  df_hgnc <- df_hgnc %>%
    filter(!(gene_symbol %in% gene_dups)) %>%
    bind_rows(uniques)
}

### Check
stopifnot(!duplicated(df_hgnc$gene_symbol))


# Save ---------------------------------------------------------
### Write out
write_tsv(df_hgnc, "input/genesets/PheGenI_Association_gene_symbols.tsv")
df_hgnc %>%
  group_by(gene_symbol_corrected_reason) %>%
  summarise(
    n = n(),
    original = toString(unique(gene_symbol_original[!is.na(gene_symbol_original)])),
    updated_symbols = toString(unique(gene_symbol[!is.na(gene_symbol)]))
  ) %>%
  write_tsv("input/genesets/PheGenI_Association_gene_symbols_updateStats.tsv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
