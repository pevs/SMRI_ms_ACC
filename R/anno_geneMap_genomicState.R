### Stanley ACC SZ vs CTRL
### Get gene annotation
### geneMap and genomicState
### (ensembl, hg19)
### Alexis Norris
### Created: 2018-11-05
### Modified: 2018-12-30


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling
library(rtracklayer)       # get and load anno
library(GenomicFeatures)   # get geneMap
library(derfinder)         # get GenomicState
library(biomaRt)           # add gene IDs

### Parameters
analysis_name <- "anno_geneMap_genomicState_GRCh37.p13_gencode.v19"


# Load files ---------------------------------------------------
### Chromosome sizes
download.file(
  url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
  destfile = "input/anno/hg19.chrom.sizes.txt"
)
chrom_sizes <- read.delim(
  "input/anno/hg19.chrom.sizes.txt", 
  header = FALSE,
  col.names = c("chr", "size")
)

### Gencode annotation
### Note:
  # gtf format is 1-based start: http://www.ensembl.org/info/website/upload/gff.html
  # bed format is 0-based start: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# Download
# Prev url: "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
download.file(
  url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
  destfile = "input/anno/gencode.v19.annotation.gtf.gz"
)
# Load
map <- rtracklayer::import(
  "input/anno/gencode.v19.annotation.gtf.gz"
)

### Add chrom.sizes info to genome map
seqlengths(map) <- chrom_sizes[match(seqlevels(map), chrom_sizes$chr), ]$size

### Add genome version to map
genome(map) <- "hg19"

### Bed correction
map$start <-map$start + 1 # bed correction

### Save entire map
### "TranscriptomeMap" since it includes genes, transcripts, exons, introns, ...
saveRDS(map, "input/anno/transcriptomeMap_GRCh37.p13_gencode.v19.rds") 


# geneMap ------------------------------------------------------
### Subset transcriptomeMap for gene-level annotations
geneMap <- map[map$type == "gene"]

### Remove columns without data for gene-level annotation
mcols(geneMap) <- mcols(geneMap)[sapply(mcols(geneMap), function (x) {
  length(x[is.na(x)]) != length(x)
})]

### Remove transcript-specific columns
mcols(geneMap) <- mcols(geneMap)[-grep("transcript_", names(mcols(geneMap)))]

### Add additional gene annotation (HGNC, entrez, description)
### Code from run_preprocess_counts_genes_exons.R
  # Get anno from ensembl using biomaRt
mart <- useMart(        # VERSION 75, hg19
  "ENSEMBL_MART_ENSEMBL",               
  dataset = "hsapiens_gene_ensembl",
  host = "feb2014.archive.ensembl.org"
)
attributes <- c(
  "hgnc_id", "hgnc_symbol", 
  "entrezgene",
  "description", 
  "chromosome_name", "band"
)
mart_anno <- getBM(
  attributes = c("ensembl_gene_id", attributes), 
  values = geneMap$gene_id, 
  mart = mart
)
  # Clean up
table(duplicated(mart_anno$ensembl_gene_id))
anno <- mart_anno %>%
  # combine chrband
  mutate(chrBand = paste(chromosome_name, band)) %>%
  dplyr::select(-one_of(c("chromosome_name", "band"))) %>%
  # Collapse duplicates
  group_by(ensembl_gene_id) %>%
  mutate(
    hgnc_id = paste0(unique(hgnc_id), collapse = "|"),
    hgnc_symbol = paste0(unique(hgnc_symbol), collapse = "|"),
    entrezgene = paste0(unique(entrezgene), collapse = "|"),
    description = paste0(unique(description), collapse = "|"),
    chrBand = paste0(unique(chrBand), collapse = "|")
  ) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame()
  # account for geneMap's version control (.# suffix) that biomaRt anno doesn't have
geneMap$abbrev <- gsub("\\..*", "", geneMap$gene_id)
  # Add any gene id's without biomaRt anno
anno <- anno %>%
  bind_rows(data_frame(
    "ensembl_gene_id" = as.character(
      geneMap$abbrev[!(geneMap$abbrev %in% anno$ensembl_gene_id)]
    )
  ))
table(geneMap$abbrev %in% anno$ensembl_gene_id)
  # Add to geneMap's metadata
stopifnot(!duplicated(anno$ensembl_gene_id))
stopifnot(!duplicated(geneMap$gene_id))
stopifnot(geneMap$abbrev %in% anno$ensembl_gene_id)
for (i in attributes) {
  mcols(geneMap)[[i]] <- anno[[i]][match(geneMap$abbrev, anno$ensembl_gene_id)]
}
  # Remove temporary column for matching
geneMap$abbrev <- NULL

### Save full geneMap
saveRDS(geneMap, "input/anno/geneMap_GRCh37.p13_gencode.v19.rds") 


### Coding genes only = geneMap_coding
### Remove gencode features that don't have a gene symbol
# summary(factor(geneMap$gene_type)) # print all gene types
geneMap_coding <- geneMap[geneMap$gene_type == "protein_coding"]
saveRDS(geneMap_coding, "input/anno/geneMap_GRCh37.p13_gencode.v19_coding.rds") 


# genomicState -------------------------------------------------
### Using approach described here: https://support.bioconductor.org/p/67680/#67766
txdb <- GenomicFeatures::makeTxDbFromGRanges(map)
genomicState <- derfinder::makeGenomicState(
  txdb, 
  species = "Homo_sapiens", 
  chrs = paste0("chr", c(1:22, "M", "X", "Y")), 
  keytype = "transcript_id"
)

### Save
saveRDS(genomicState, "input/anno/genomicState_GRCh37.p13_gencode.v19.rds") 
  
  
# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
