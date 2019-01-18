#!/usr/bin/env r

### Annotate derfinder regions on cluster
### Alexis Norris
### Created: 2018-03-12
### Modified: 16 March 2018 (PFC -> ACC)

### Libraries
library(derfinder)
library(GenomicRanges)
library(magrittr)
library(dplyr)

### Dir paths
REFS <- "/mnt/data/RNASeq/methods/references/hg19/"
brain_region <- "ACC"

### Load regions
Sys.time()
regionMat <- readRDS("regionMat_ACC_75bp.rds")
Sys.time()

### Extract regions from GRanges object 
regions <- unlist(GRangesList(lapply(regionMat, "[[", "regions")))
names(regions) <- paste0(seqnames(regions), ":", 
                         start(regions), "-", 
                         end(regions)) 

### GeneMap: Check regions for overlapping genes
geneMap <- readRDS(paste0(REFS, "GeneMap_Hsapiens.gencode.v19.GRCh37.rds"))
geneMapGR <- GenomicRanges::makeGRangesFromDataFrame(
  geneMap, 
  keep.extra.columns = TRUE
)
rm(geneMap)

Sys.time()
dA <- GenomicRanges::distanceToNearest(regions, geneMapGR)
Sys.time()
regions$feature_name <- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_symbol<- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_type <- geneMapGR$gene_type[subjectHits(dA)]
regions$gene_id <- geneMapGR$gene_id[subjectHits(dA)]
regions$distToGene <- mcols(dA)$distance

### GenomicState: Count number of region anno types
gs <- readRDS(paste0(REFS, "GenomicState_Hsapiens.gencode.v19.GRCh37.p13.rds"))
ensemblAnno <- derfinder::annotateRegions(regions, gs$fullGenome)
colSums(ensemblAnno$countTable)
mcols(regions) <- cbind(mcols(regions), ensemblAnno$countTable)
regions$anno_class <- NA
regions$anno_class[regions$exon > 0 & 
                     regions$intron == 0 & 
                     regions$intergenic == 0] = "strictExonic"
regions$anno_class[regions$exon == 0 & 
                     regions$intron > 0 & 
                     regions$intergenic == 0] = "strictIntronic"
regions$anno_class[regions$exon == 0 & 
                     regions$intron == 0 & 
                     regions$intergenic > 0] = "strictIntergenic"
regions$anno_class[regions$exon > 0 & 
                     regions$intron > 0 & 
                     regions$intergenic == 0] = "exonIntron"
regions$anno_class[regions$exon > 0 & 
                     regions$intergenic > 0] = "extendUTR"
summary(factor(regions$anno_class))

### Save bed
rtracklayer::export(regions, 
                    paste0("derfinder_", brain_region, "_ERs.bed"))

### Extract featureData (annotation) from granges object
class(regions)
head(regions)
head(as.data.frame(regions))

features <- as.data.frame(regions) %>%
  mutate(chr = seqnames) %>%
  mutate(pos_hg19 = paste0(chr, ":", start, "-", end),
         feature_id = paste0(chr, ":", start, "-", end)) %>% 
  select(feature_id, feature_name, 
         pos_hg19, chr, start, end, strand, 
         gene_symbol, gene_type, gene_id, 
         distToGene, anno_class)

### Summarise strand
summary(factor(features$strand))

### feature_id = Region (chr_:_-_) = pos_hg19 = rownames
features <- as.data.frame(features)
rownames(features) <- features$feature_id

### Save
Sys.time()
saveRDS(features, paste0("features_", brain_region, "_regions.rds")) 

### Done
Sys.time()




