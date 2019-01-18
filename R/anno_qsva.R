### SMRI Array ACC SZ vs CTRL
### Annotation of qSVA degradation-associated regions 
### Get updated gene symbols
### Alexis Norris
### Created: 2018-11-07
### Modified : 2019-01-05


# Notes --------------------------------------------------------
# this includes both polyA and riboZero regions
# some gene annotation is very far from gene, 
# recommended to filter for <1kb on distToGene column


# Input files --------------------------------------------------
### From qSVA PNAS 2017 paper
#polyA_regions <- read.delim("input/qsva/qsva_polyA_regions.txt", header = FALSE)[,1:3]
#ribo_regions <- read.delim("input/qsva/qsva_polyA_regions.txt", header = FALSE)[,1:3]

### HGNC approved gene symbol check (from anno_hgnc.R)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")


# Setup --------------------------------------------------------
### Packages
library(tidyverse)     # wrangling
library(GenomicRanges) # annotation
library(derfinder)     # annotation

### Parameters
analysis_name <- "anno_qsva"


# Load degradation regions -------------------------------------
### Load regions (data from qSVA paper)
polyA_regions <- read.delim("input/qsva/qsva_polyA_regions.txt", header = FALSE)[ , 1:3]
ribo_regions <- read.delim("input/qsva/qsva_ribozero_regions.txt", header = FALSE)[ , 1:3]

### Combine 
polyA_regions$method <- "polyA"
ribo_regions$method <- "ribo"
deg_regions <- bind_rows(polyA_regions, ribo_regions)
names(deg_regions) <- c("chr", "start", "end", "method")
deg_regions$region <- paste0(deg_regions$chr, ":", deg_regions$start, "-", deg_regions$end)
rownames(deg_regions) <- deg_regions$region


# Get genes ----------------------------------------------------
### Get nearest gene names for degradation-associated regions
### Convert to GRanges object
regs <- makeGRangesFromDataFrame(deg_regions, keep.extra.columns = TRUE)

### Load geneMap 
geneMap <- readRDS("input/anno/geneMap_GRCh37.p13_gencode.v19.rds")

### Get nearest gene
hits <- distanceToNearest(regs, geneMap, ignore.strand = TRUE)
regs$gene_id <- geneMap$gene_id[subjectHits(hits)]
regs$gene_symbol <- geneMap$gene_name[subjectHits(hits)]
regs$gene_strand <- strand(geneMap)[subjectHits(hits)]

### Add column with distance to nearest gene
regs$distToGene <- mcols(hits)$distance

### Add annoClass using GenomicState 
gs <- readRDS("input/anno/genomicState_GRCh37.p13_gencode.v19.rds")
countTable <- annotateRegions(regs, gs$codingGenome)$countTable
mcols(regs) <- cbind(mcols(regs), countTable)
regs$annoClass <- NA
regs$annoClass[regs$X3UTR > 0] <- "3UTR"
regs$annoClass[regs$X5UTR > 0] <- "5UTR"
regs$annoClass[regs$promoter > 0] <- "promoter"
regs$annoClass[regs$exon > 0 & regs$intron == 0 & regs$intergenic == 0] <- "strictExonic"
regs$annoClass[regs$exon == 0 & regs$intron > 0 & regs$intergenic == 0] <- "strictIntronic"
regs$annoClass[regs$exon == 0 & regs$intron == 0 & regs$intergenic > 0] <- "strictIntergenic"
regs$annoClass[regs$exon > 0 & regs$intron > 0 & regs$intergenic == 0] <- "exonIntron"
regs$annoClass[regs$exon > 0 & regs$intergenic > 0] <- "extendUTR"
write_lines(capture.output(summary(factor(regs$annoClass))), "logs/qsva/anno_qsva_classes.txt")

## Extract annotation of degrad-assoc regions
regs <- as_data_frame(regs)


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
hgnc[hgnc == "NA"] <- NA
hgnc[hgnc == ""] <- NA

### Setup
original <- regs
df <- original %>%
  mutate(
    # id column (ensembl ID)
    id = region,
    # column to match on; need to remove version suffix
    by_col = gsub("\\..*", "", gene_id),
    # Column indicating original/uncorrected gene symbol
    gene_symbol_uncorrected = gene_symbol
  ) %>%
  # Fix classes, for binding later
  mutate_if(is.factor, as.character)

### Don't have ensembl ID (many jxns)
no_ensembl <- df %>%
  filter(is.na(by_col)) %>%
  # Replace with NA
  mutate(gene_symbol = NA) %>%
  # Note that it was changed
  mutate(gene_symbol_corrected_reason = "no ensembl ID")

### Join by ensembl IDs
df_hgnc <- df %>%
  left_join(hgnc, by = c("by_col" = "ensembl_gene_id")) %>%
  # Fix classes, for binding later
  mutate_if(is.factor, as.character)

### Outcome 1. Ok (is current approved symbol)
ok <- df_hgnc %>%
  # symbols match
  filter(gene_symbol == approved_symbol) %>%
  # no change
  mutate(gene_symbol_corrected_reason = "OK")

### Outcome 2. Updated gene symbol (previous/alias)
updated <- df_hgnc %>%
  # symbols don't match
  filter(gene_symbol != approved_symbol) %>%
  # Replace with approved
  mutate(gene_symbol = approved_symbol) %>%
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
  new_anno <- df %>%
    filter(is.na(by_col)) %>%
    # Replace with NA
    mutate(gene_symbol = NA) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "no ensembl ID") %>% 
    bind_rows(new_anno)
}

new_anno <- new_anno %>%
  # clarify colnames
  mutate(
    chr = seqnames,
    pos_hg19 = paste0(seqnames, ":", start, "-", end)
  ) %>%
  # remove extra columns
  dplyr::select(
    region, gene_symbol, mult_matches, pos_hg19, chr, start, end, 
    -one_of(c("id", "by_col", "strand", "seqnames")), # created for HGNC matching, NA, or replaced
    everything()
  )

### Check
new_anno[new_anno == ""] <- NA
new_anno[new_anno == "NA"] <- NA
stopifnot(nrow(new_anno) == nrow(original))

### Summarize
write_lines(
  capture.output(
    new_anno %>% 
      group_by(method) %>% 
      summarise(
        n_regions = n(), 
        n_genes = length(unique(gene_symbol)),
        n_region_maps_multiple_genes = length(!is.na(mult_matches)[!is.na(mult_matches) != "NO"])
    )),
  "logs/qsva/anno_qsva_summary.txt"
)


# Save ---------------------------------------------------------
write_tsv(new_anno, "input/qsva/qSVA_degradation_regions_anno.tsv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
