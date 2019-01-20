#!/usr/bin/env r

### Annotate derfinder regions
### Alexis Norris
### Created: 2018-03-12
### Modified: 2018-03-16 March
### Run in bash script with "module add R/3.2.3"


# Setup --------------------------------------------------------
### Packages
library(derfinder)
library(GenomicRanges)
library(magrittr)
library(dplyr)

### Log
Sys.time()


# Setup --------------------------------------------------------
# Load regions -------------------------------------------------
regionMat <- readRDS("regionMat.rds")


# Get regions --------------------------------------------------
### Extract regions from GRanges object 
regions <- unlist(GRangesList(lapply(regionMat, "[[", "regions")))
names(regions) <- paste0(
  seqnames(regions), ":", 
  start(regions), "-", 
  end(regions)
) 


# GeneMap ------------------------------------------------------
### Check regions for overlapping genes
geneMap <- readRDS("GeneMap_Hsapiens.gencode.v19.GRCh37.rds")
geneMapGR <- makeGRangesFromDataFrame(
  geneMap, 
  keep.extra.columns = TRUE
)
dA <- distanceToNearest(regions, geneMapGR)
regions$feature_name <- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_symbol<- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_type <- geneMapGR$gene_type[subjectHits(dA)]
regions$gene_id <- geneMapGR$gene_id[subjectHits(dA)]
regions$distToGene <- mcols(dA)$distance


# GenomicState -------------------------------------------------
### Count number of region anno types
gs <- readRDS("GenomicState_Hsapiens.gencode.v19.GRCh37.p13.rds")
ensemblAnno <-annotateRegions(regions, gs$fullGenome)
mcols(regions) <- cbind(mcols(regions), ensemblAnno$countTable)
regions$anno_class <- NA
regions$anno_class[regions$exon > 0 & 
                     regions$intron == 0 & 
                     regions$intergenic == 0] <- "strictExonic"
regions$anno_class[regions$exon == 0 & 
                     regions$intron > 0 & 
                     regions$intergenic == 0] <- "strictIntronic"
regions$anno_class[regions$exon == 0 & 
                     regions$intron == 0 & 
                     regions$intergenic > 0] <- "strictIntergenic"
regions$anno_class[regions$exon > 0 & 
                     regions$intron > 0 & 
                     regions$intergenic == 0] <- "exonIntron"
regions$anno_class[regions$exon > 0 & 
                     regions$intergenic > 0] <- "extendUTR"
summary(factor(regions$anno_class))


# Extract anno -------------------------------------------------
### Extract featureData (annotation) from granges object
features <- as.data.frame(regions) %>%
  mutate(chr = seqnames) %>%
  mutate(
    pos_hg19 = paste0(chr, ":", start, "-", end),
    feature_id = paste0(chr, ":", start, "-", end)
  ) %>% 
  select(
    feature_id, feature_name, 
    pos_hg19, chr, start, end, strand, 
    gene_symbol, gene_type, gene_id, 
    distToGene, anno_class
  )

### feature_id = Region (chr_:_-_) = pos_hg19 = rownames
features <- as.data.frame(features)
rownames(features) <- features$feature_id


# Save ---------------------------------------------------------
### Full annotation
saveRDS(features, "features_regions.rds")

###  bed
rtracklayer::export(regions, "derfinder__ERs.bed")


# Log ----------------------------------------------------------
Sys.time()
