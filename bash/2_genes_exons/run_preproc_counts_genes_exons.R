### Preprocessing featureCounts results: genes & exons
### Alexis Norris
### Created: 2018-02-28 (Andrew Jaffe)
### Modified: 2018-03-02
### For HPC!

# IMPORTANT: wd should be "hisat2_hg19"
setwd("../hisat2_hg19")
Sys.time()

# Packages ----------------------------------
#library(tidyverse) # install failed on cluster
library(magrittr)
library(readr)


# Load lieber_functions_ajr.R ----------------
# jaffelab pkg requires R >= 3.3.0, so needed to download that 
#library("devtools")
#install_github("LieberInstitute/jaffelab") 
#library("jaffelab")
### Issues loading the library on cluster, because updated R won't configure because of zlib not being detected as updated (it should be in path though!)
### So load the functions here instead
ss <- function (x, pattern, slot = 1, ...) {
    sapply(strsplit(x = x, split = pattern, ...), "[", slot)
}

# Load phenoData -----------------------------
pd <- read_tsv("featureCounts_phenoData.txt") # must have "totalMapped" column
pd$full_path <- paste0("bams/", pd$FileName, ".bam")

### Make sure no duplicates
pd <- pd[!duplicated(pd), ]  

# LibSize ------------------------------------
### Use totalLibSize calculated for qSVA
print(head(pd))

# Packages -----------------------------------
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)

# Gene Counts --------------------------------
### Load gene count files 
geneFn <- paste0("counts/genes/", 
                 pd$FileName, 
                 ".counts")
names(geneFn) <- pd$FileName
table(file.exists(geneFn))

### Load gene annotation 
geneMap <- read.delim(geneFn[1], skip = 1, as.is = TRUE)[ ,1:6]
rownames(geneMap) <- geneMap$Geneid
geneMap$Chr <- paste0("chr", ss(geneMap$Chr, ";"))
geneMap$Start <- as.numeric(ss(geneMap$Start, ";"))
tmp <- strsplit(geneMap$End, ";")
geneMap$End <- as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand <- ss(geneMap$Strand, ";")
geneMap$Geneid <- NULL

### Annotate with biomaRt  
ensembl <- useMart(
  "ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
	dataset = "hsapiens_gene_ensembl",
	host = "feb2014.archive.ensembl.org"
)
sym <- getBM(
  attributes = c("ensembl_gene_id",
                 "hgnc_symbol",
                 "entrezgene"), 
	values = rownames(geneMap), 
  mart = ensembl
)
geneMap$Symbol <- sym$hgnc_symbol[match(rownames(geneMap), 
                                        sym$ensembl_gene_id)]
geneMap$EntrezID <- sym$entrezgene[match(rownames(geneMap), 
                                         sym$ensembl_gene_id)]
	
### Get counts 
geneCountList <- mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[ ,1]
}, mc.cores = 12)
geneCounts <- do.call("cbind", geneCountList)
rownames(geneCounts) <- rownames(geneMap)
head(geneCounts)

### Put in order
geneCounts <- geneCounts[ ,pd$FileName] 

### Check
head(pd$totalMapped)
dim(pd)
dim(geneCounts)
head(geneCounts)

### Calculate RPKM 
bg <- matrix(
  rep(pd$totalMapped), 
  nc = nrow(pd), 
  nr = nrow(geneCounts),	
  byrow = TRUE)
widG <- matrix(
  rep(geneMap$Length), 
  nr = nrow(geneCounts), 
	nc = nrow(pd),	
  byrow = FALSE)
geneRpkm <- geneCounts / (widG/1000) / (bg/1e6)

### number of reads assigned 
geneStatList <- lapply(
  paste0(geneFn, ".summary"), 
	read.delim,
  row.names = 1)
geneStats <- do.call("cbind", geneStatList)
colnames(geneStats) <- pd$FileName 
pd$totalAssignedGene <- as.numeric(geneStats[1, ] / colSums(geneStats))


# Exon counts ---------------------------------
### Filenames
exonFn <- paste0("counts/exons/", pd$FileName, ".counts")
names(exonFn) = pd$FileName
all(file.exists(exonFn))

### Load exon annotation 
exonMap <- read.delim(exonFn[1], skip = 1, as.is = TRUE)[ ,1:6]
exonMap$Symbol <- sym$hgnc_symbol[match(exonMap$Geneid, 
                                        sym$ensembl_gene_id)]
exonMap$EntrezID <- sym$entrezgene[match(exonMap$Geneid, 
                                         sym$ensembl_gene_id)]
rownames(exonMap) <- paste0("e", rownames(exonMap))
exonMap$Chr <- paste0("chr", exonMap$Chr)

### Load counts
exonCountList <- mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[ ,1]
}, mc.cores = 12)
exonCounts <- do.call("cbind", exonCountList)
rownames(exonCounts) <- rownames(exonMap)

### Put in order
exonCounts <- exonCounts[ , pd$FileName] 

### Remove duplicates
eMap <- GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex <- which(!duplicated(eMap))
exonCounts <- exonCounts[keepIndex, ]
exonMap <- exonMap[keepIndex, ]

### Calculate RPKM
bgE <- matrix(
  rep(pd$totalMapped), 
  nc = nrow(pd), 
	nr = nrow(exonCounts),	
  byrow = TRUE)
widE <- matrix(
  rep(exonMap$Length), 
  nr = nrow(exonCounts), 
	nc = nrow(pd),	
  byrow = FALSE)
exonRpkm = exonCounts / (widE/1000) / (bgE/1e6)

 
# Save counts ----------------------------
save(pd, 
     geneCounts, geneMap, 
     exonCounts, exonMap, 
     compress = TRUE,
     file = "rawCounts_n20_genes_exons_hpc.rda")
save(pd, 
     geneRpkm,	geneMap, 
     exonRpkm, exonMap, 
     compress = TRUE,
     file = "rpkmCounts_n20_genes_exons_hpc.rda")

# Write out for coverage ------------------
write_tsv(pd, "samples_with_bams_hpc.txt")

# Session info ---------------------------
library(sessioninfo)
sessioninfo::session_info()
Sys.time()
proc.time()

