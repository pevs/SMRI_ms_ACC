### Stanley ACC SZ vs CTRL
### Get genes in PGC 108 loci
### From Ripke, et al. Nature 2014
### Get gene symbols from chrom pos 
### Rather than by updating gene symbols
### Alexis Norris
### Created: 2019-01-02
### Modified: 2019-01-17


# Input files --------------------------------------------------
### Source of data: Table S3 in Ripke, et al. Nature 2014
#pgc <- read_tsv("input/genesets/PGC_Ripke_2014/PGC_108loci_Ripke_2014_TableS3.txt")   

### geneMap 
#geneMap <- readRDS("input/genesets/geneMap_GRCh37.p13_gencode.v19.rds")

### HGNC approved gene symbol check (from anno_hgnc.R)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling; warning: GenomicRanges masks dplyr::select
library(janitor)           # wrangling
library(GenomicRanges)     # annotate DEGs with overlapping repeatMasker & all HGNC genes

### Parameters
analysis_name <- "anno_PGC_Ripke2014_loci"


# Load PGC positions -------------------------------------------
### Load
### Ripke, et al. Nature 2014's Table S3
df <- read_tsv("input/genesets/PGC_Ripke_2014/PGC_108loci_Ripke_2014_TableS3.txt") %>%
  remove_empty() %>%
  
  # Extract positions
  separate(locus_pos_hg19, into = c("chr", "start_end"), sep = ":", remove = FALSE, convert = TRUE) %>%
  separate(start_end, into = c("start", "end"), sep = "-", remove = TRUE, convert = TRUE) %>%
  
  # Clean up genes column 
  # Remove * in the genes
  # Change separation from space to ";'
  mutate(protein_coding_genes = gsub(" ", ";", gsub(" \\*", "", protein_coding_genes)))

# Clarify column names
names(df) <- recode(
  names(df), 
  protein_coding_genes = "locus_genes",
  SZ = "SZ_GWAS_prev",
  NHGRI_GWAS_catalog = "SZ_GWAS_prev_study"
) 

### Clean up gene column more
df$locus_genes <- gsub(";\\(micro-RNA\\)", "", df$locus_genes)
df$locus_genes <- gsub("Locus;too;broad", "(locus too broad)", df$locus_genes)


# Get genes ----------------------------------------------------
### Get all genes in these loci

### Load geneMap 
geneMap <- readRDS("input/anno/geneMap_GRCh37.p13_gencode.v19.rds")

### Convert to GRanges object
loci <- makeGRangesFromDataFrame(
  df, 
  keep.extra.columns = TRUE,
  ignore.strand = TRUE
)

### Get all genes
hits <- findOverlaps(loci, geneMap, ignore.strand = TRUE, select = "all")

### Get overlaps
hits <- findOverlaps(                     
  geneMap, loci,
  ignore.strand = TRUE,         
  type = "any", select = "all"
)
overlapsGR <- pintersect(
  geneMap[queryHits(hits)], 
  loci[subjectHits(hits)]
)

### Combine geneMap anno and loci data 
mcols(overlapsGR) <- c(
  mcols(overlapsGR),
  mcols(loci[subjectHits(hits)])
)
overlapsGR$gene_pos_hg19 <- paste0(
  seqnames(loci[subjectHits(hits)]), ":",
  start(loci[subjectHits(hits)]), "-", end(loci[subjectHits(hits)])
)

### Extract table from GR object & Clean up
### Only keep ensembl ID, since will add gene anno from HGNC next
### Also not keeping chr/start/end, since the pos could be for gene OR locus, so use explicit "gene_pos_hg19" and "locus_pos_hg19" instead
keep_cols <- names(df)[names(df) %in% names(mcols(overlapsGR))]
overlaps <- overlapsGR %>%
  as_data_frame() %>%
  dplyr::select(gene_id, gene_pos_hg19, keep_cols)


# Update outdated HGNC gene symbols ----------------------------
### Check that gene_symbol is the current approved HGNC symbol
### if outdated (alias), update
### If ensembl ID not in HGNC, remove gene_symbol (make NA)

### Load HGNC anno (from anno_hgnc.R)
hgnc <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")$by_ensembl
names(hgnc) <- recode(
  names(hgnc),
  "symbol" = "gene_symbol"
)

### Get gene symbols using ensembl
### Remove version suffix of ensembl gene id (to match hgnc)
overlaps_hgnc <- overlaps %>%
  mutate(ensembl_gene_id = gsub("\\..*", "", gene_id)) %>%
  left_join(hgnc, by = "ensembl_gene_id")

### Check
overlaps_hgnc[overlaps_hgnc == ""] <- NA
overlaps_hgnc[overlaps_hgnc == "NA"] <- NA
stopifnot(nrow(overlaps_hgnc) == nrow(overlaps))


# Collapse genes -----------------------------------------------
### If one gene hit for multiple regions/SNPs
### (aka if multiple SNPs-evidence for a gene) -- similar to code used for anno_pheGenI.R

### Get number of duplicates (if any)
df_notNA <- overlaps_hgnc %>%
  filter(!is.na(gene_symbol))
gene_dups <- df_notNA$gene_symbol[duplicated(df_notNA$gene_symbol)]

### If there are duplicates
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
  df_no_dups <- pgc_genes_updated %>%
    filter(is.na(gene_symbol)) %>%
    bind_rows(df_notNA_collapsed)
  stopifnot(!duplicated(df_no_dups$gene_symbol[!is.na(df_no_dups$gene_symbol)]))
  overlaps_hgnc <- df_no_dups
}

### Check
x <- overlaps_hgnc$gene_symbol
stopifnot(!duplicated(x[!is.na(x)]))


# Save ---------------------------------------------------------
write_tsv(overlaps_hgnc, "input/genesets/PGC_Ripke_2014/PGC_108loci_Ripke_2014_TableS3_gene_symbols.tsv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
