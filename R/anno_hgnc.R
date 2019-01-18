### Stanley ACC SZ vs CTRL
### Get gene annotation
### Current HGNC symbols + their aliases
### Alexis Norris
### Created: 2018-11-05
### Modified: 2019-01-03




### [1/16/19 - TO DO] -- for multiple matches, sometimes it's NOT duplicate gene symbols, but instead duplicate HGNC IDs where one has been withdrawn
### I need to check that and remove those before calling something "mult_matches"



# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling; warning: dplyr::select masked by GenomicRanges
library(DBI)               # to query org.Hs.eg.db
library(org.Hs.eg.db)      # annotate DEGs with gene symbols

### Parameters
analysis_name <- "anno_hgnc"


# Load files ---------------------------------------------------
### Downloaded 2018-12-31
#download.file(
#  "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt",
#  destfile = "input/anno/hgnc_complete_set_ftp.tsv"
#)

### Load
hgnc <- read.delim("input/anno/hgnc_complete_set_ftp.tsv", stringsAsFactors = FALSE) # don't use readr::read_tsv

### Add column to indicate date accessed/downloaded and method
hgnc$download_url <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
hgnc$download_date <- "2018_12_31"

### Convert all missing data to NA
hgnc[hgnc == ""] <- NA
hgnc[hgnc == "NA"] <- NA


# By HGNC ID ---------------------------------------------------
### For conversion/matching using HGNC ID
by_id <- hgnc
stopifnot(!duplicated(by_id$hgnc_id))


# By ensembl ID ------------------------------------------------
### For conversion/matching using ensembl ID

### Some have more than one ID
table(duplicated(hgnc$ensembl_gene_id))

### Some have no ensembl ID
table(is.na(hgnc$ensembl_gene_id))

### So need to collapse 
hgnc_ensembl <- hgnc %>%
  filter(!is.na(ensembl_gene_id)) %>%
  # Column indicating multiple matches
  mutate(mult_matches = "NO") %>%
  # Convert all to character, to prevent conversion errors when collapse duplicates (e.g. ENTREZ was numerical)
  mutate_all(as.character)  
duplicates <- hgnc_ensembl$ensembl_gene_id[duplicated(hgnc_ensembl$ensembl_gene_id)]
by_ensembl <- hgnc_ensembl %>%
  filter(ensembl_gene_id %in% duplicates) %>%
  # collapse duplicates for each ensembl ID
  group_by(ensembl_gene_id) %>%
  mutate_all(.funs = function(x) paste0(unique(as.character(x)), collapse = ";")) %>%
  ungroup() %>%
  distinct() %>%
  # Multiple matches
  mutate(mult_matches = "YES") %>%
  # Add back non-duplicated data
  bind_rows(filter(hgnc_ensembl, !(ensembl_gene_id %in% duplicates))) %>%
  as.data.frame(stringsAsFactors = FALSE)
stopifnot(!duplicated(by_ensembl$ensembl_gene_id))
stopifnot(by_ensembl$ensembl_gene_id %in% hgnc_ensembl$ensembl_gene_id)


# By approved symbol -------------------------------------------
### For conversion/matching using HGNC approved symbol

### Some have more than one ID
table(duplicated(hgnc$symbol))

### All have gene symbol (no NAs)
table(is.na(hgnc$symbol))

### So need to collapse 
hgnc_symbol <- hgnc %>%
  filter(!is.na(symbol)) %>%                 
  # Column indicated in multiple mappings
  mutate(mult_matches = "NO") %>%
  # Convert all to character, to prevent conversion errors when collapse duplicates (e.g. ENTREZ was numerical)
  mutate_all(as.character)  
duplicates <- hgnc_symbol$symbol[duplicated(hgnc_symbol$symbol)]
by_symbol <- hgnc_symbol %>%
  filter(symbol %in% duplicates) %>%
  # collapse duplicates for each symbol
  group_by(symbol) %>%
  mutate_all(.funs = function(x) paste0(unique(as.character(x[!is.na(x)])), collapse = ";")) %>%
  ungroup() %>%
  distinct() %>%
  # Multiple mappings
  mutate(mult_matches = "YES") %>%
  # Add back non-duplicated data
  bind_rows(filter(hgnc_symbol, !(symbol %in% duplicates))) %>%
  as.data.frame(stringsAsFactors = FALSE)
stopifnot(!duplicated(by_symbol$symbol))
stopifnot(by_symbol$symbol %in% hgnc_symbol$symbol)


# By previous symbol -------------------------------------------
### For conversion/matching using HGNC alias_symbol
### Expand, to split aliases to one row each

### Alias ~ Symbol
prev_to_sym <- by_symbol %>%
  # split aliases
  separate_rows(prev_symbol, sep = "\\|") %>%   
  # remove empty
  filter(!is.na(prev_symbol)) %>%
  # remove if previous symbol is an approved symbol for another gene (~585)
  filter(!(prev_symbol %in% hgnc$symbol)) %>%
  as.data.frame()

### Previos ~ Currently Approved Symbol, removing duplicate matches
### Indicate if previous symbol matches to multiple symbols and which symbols they map to
dups <- prev_to_sym$prev_symbol[duplicated(prev_to_sym$prev_symbol)]
prev_to_sym_no_dups <- prev_to_sym %>%
  # only conversion columns
  dplyr::select(prev_symbol, symbol) %>%
  # collapse duplicate approved symbols for each previous
  group_by(prev_symbol) %>%
  mutate(symbol = paste0(symbol, collapse = "|")) %>%
  ungroup() %>%
  distinct() %>%
  # Column indicated in multiple approved symbols map to previous symbol
  mutate(mult_matches = ifelse(
    .$prev_symbol %in% dups, "YES", "NO"
  )) %>%
  as.data.frame()
stopifnot(!duplicated(prev_to_sym_no_dups$prev_symbol))


### Previous ~ ID
prev_to_id <- by_id %>%
  # only conversion columns
  dplyr::select(prev_symbol, hgnc_id) %>%
  # split previous symbols
  separate_rows(prev_symbol, sep = "\\|") %>%   
  # remove empty
  filter(!is.na(prev_symbol)) %>%
  # remove if previous is an approved symbol for another gene (~585)
  filter(!(prev_symbol %in% hgnc$symbol)) %>%
  as.data.frame()

### Previous ~ ID, removing duplicate matches
### Indicate if alias matches to multiple IDs and which IDs they map to
dups <- prev_to_id$prev_symbol[duplicated(prev_to_id$prev_symbol)]
prev_to_id_no_dups <- prev_to_id %>%
  # only conversion columns
  dplyr::select(prev_symbol, hgnc_id) %>%
  # collapse duplicate approved symbols for each previous
  group_by(prev_symbol) %>%
  mutate(hgnc_id = paste0(hgnc_id, collapse = "|")) %>%
  ungroup() %>%
  distinct() %>%
  # Column indicated in multiple approved symbols map to previous
  mutate(mult_matches = ifelse(
    .$prev_symbol %in% dups, "YES", "NO"
  )) %>%
  as.data.frame()
stopifnot(!duplicated(prev_to_id_no_dups$prev_symbol))


# By alias symbol ----------------------------------------------
### For conversion/matching using HGNC alias_symbol
### Expand, to split aliases to one row each

### Alias ~ Symbol
alias_to_sym <- by_symbol %>%
  # split aliases
  separate_rows(alias_symbol, sep = "\\|") %>%   
  # remove empty
  filter(!is.na(alias_symbol)) %>%
  # remove if alias is an approved or previous symbol for another gene (~585)
  filter(!(alias_symbol %in% hgnc$prev_symbol)) %>%
  filter(!(alias_symbol %in% hgnc$symbol)) %>%
  as.data.frame()

### Alias ~ Symbol, removing duplicate matches
### Indicate if alias matches to multiple symbols and which symbols they map to
dups <- alias_to_sym$alias_symbol[duplicated(alias_to_sym$alias_symbol)]
alias_to_sym_no_dups <- alias_to_sym %>%
  # only conversion columns
  dplyr::select(alias_symbol, symbol) %>%
  # collapse duplicate approved symbols for each alias
  group_by(alias_symbol) %>%
  mutate(symbol = paste0(symbol, collapse = "|")) %>%
  ungroup() %>%
  distinct() %>%
  # Column indicated in multiple approved symbols map to alias
  mutate(mult_matches = ifelse(
    .$alias_symbol %in% dups, "YES", "NO"
  )) %>%
  as.data.frame()
stopifnot(!duplicated(alias_to_sym_no_dups$alias_symbol))


### Alias ~ ID
alias_to_id <- by_id %>%
  # only conversion columns
  dplyr::select(alias_symbol, hgnc_id) %>%
  # split aliases
  separate_rows(alias_symbol, sep = "\\|") %>%   
  # remove empty
  filter(!is.na(alias_symbol)) %>%
  # remove if alias is an approved or previous symbol for another gene (~585)
  filter(!(alias_symbol %in% hgnc$prev_symbol)) %>%
  filter(!(alias_symbol %in% hgnc$symbol)) %>%
  as.data.frame()

### Alias ~ ID, removing duplicate matches
### Indicate if alias matches to multiple IDs and which IDs they map to
dups <- alias_to_id$alias_symbol[duplicated(alias_to_id$alias_symbol)]
alias_to_id_no_dups <- alias_to_id %>%
  # only conversion columns
  dplyr::select(alias_symbol, hgnc_id) %>%
  # collapse duplicate approved symbols for each alias
  group_by(alias_symbol) %>%
  mutate(hgnc_id = paste0(hgnc_id, collapse = "|")) %>%
  ungroup() %>%
  distinct() %>%
  # Column indicated in multiple approved symbols map to alias
  mutate(mult_matches = ifelse(
    .$alias_symbol %in% dups, "YES", "NO"
  )) %>%
  as.data.frame()
stopifnot(!duplicated(alias_to_id_no_dups$alias_symbol))


# Save ---------------------------------------------------------
hgncList <- list(
  "by_id" = by_id,
  "by_ensembl" = by_ensembl,
  "by_symbol" = by_symbol,
  
  "previous_to_id_duplicateIDs" = prev_to_id,
  "previous_to_id_noDuplicates" = prev_to_id_no_dups,
  "previous_to_symbol_duplicateSymbols" = prev_to_sym,
  "previous_to_symbol_noDuplicates" = prev_to_sym_no_dups,
  
  "alias_to_id_duplicateIDs" = alias_to_id,
  "alias_to_id_noDuplicates" = alias_to_id_no_dups,
  "alias_to_symbol_duplicateSymbols" = alias_to_sym,
  "alias_to_symbol_noDuplicates" = alias_to_sym_no_dups
)
hgncList <- lapply(hgncList, function (df) {
  df[df == "NA"] <- NA
  df[df == ""] <- NA
  df
})
saveRDS(hgncList, "input/anno/hgnc_complete_set_ftp_key.rds")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
