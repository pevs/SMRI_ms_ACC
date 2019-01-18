#!/usr/bin/env r

### Alexis Norris
### created: 1 August 2016
### modified: 12 March 2018 *for HPC, which needs more threads*

print(Sys.time())
library(derfinder)

### Generate full coverage
# Load chromosome information used in assembly
chrom.sizes <- read.delim("/mnt/data/RNASeq/methods/references/hg19/hg19.chrom.sizes", header = FALSE)
print(head(chrom.sizes))

# Remove unmapped contigs
chrom.sizes <- chrom.sizes[-grep("_", chrom.sizes$V1), ]

# HPC (100bp) ------------------------------------------
# Locate bigwig files
# .bai files must also be in the directory and have identical names
files_HPC <- rawFiles("/mnt/data/RNASeq/data/Stanley/all_hisat_hg19/derfinder/bigwigs/HPC", samplepatt = "bw", fileterm = NULL)
names(files_HPC) <- gsub(".bw", "", names(files_HPC))

# Load full coverage from bigwig files
fullCov_HPC <- fullCoverage(files = files_HPC, chrs = as.character(chrom.sizes$V1), chrlens = chrom.sizes$V2, verbose = TRUE, mc.cores.load = 10)
saveRDS(fullCov_HPC, "fullCov_HPC.rds")

### Filter coverage to reduce computational time of `regionMatrix`
# Using mean coverage of 5 reads (mean across all samples), for a region to be included as an expressed region  
filteredCov_HPC <- lapply(fullCov_HPC, filterData, returnMean = TRUE, filter = "mean", cutoff = 5, verbose = TRUE)
saveRDS(filteredCov_HPC, "filteredCov_HPC.rds")

### Calculate expressed region counts  
# L = read length (here, 100bp)
regionMat100_HPC <- regionMatrix(
	filteredCov_HPC, 
	L = 100, 
	runFilter = FALSE, 
	returnBP = FALSE, # don't need this, and might speed it up?!
	verbose = TRUE)

### Save
print("Saving regionMat..."); Sys.time()
saveRDS(regionMat100_HPC, "regionMat_HPC_100bp.rds")
print("Saving regionMat...DONE"); Sys.time()

### Clean-up
rm(fullCov_HPC, 
   filteredCov_HPC,
   regionMat100_HPC)

# Annotate regions -------------------------------
### Libraries
library(derfinder)
library(GenomicRanges)
library(magrittr)
library(dplyr)

### Dir paths # had typo here initially, so then ran separately as run_regionAnno_hpc.r on 3018-03-14
REFS <- "/mnt/data/RNASeq/methods/references/hg19/"
brain_region <- "HPC"

### Load regions (already loaded, since combined with script generating the regionMat!)
#regionMat <- readRDS("regionMat_HPC_100bp.rds")

### Extract regions from GRanges object 
print("Extracting regions from regionMatrix"); Sys.time()
regions <- unlist(GRangesList(lapply(regionMat, "[[", "regions")))
names(regions) <- paste0(seqnames(regions), ":", 
                         start(regions), "-", 
                         end(regions)) 

### GenomicState: Count number of region anno types (pre-filtered)
print("Loading GenomicState_Hsapiens.gencode.v19.GRCh37.p13.rds"); Sys.time()
gs <- readRDS(paste0(REFS, "GenomicState_Hsapiens.gencode.v19.GRCh37.p13.rds"))
print("derfinder::annotateRegions()"); Sys.time()
ensemblAnno <- derfinder::annotateRegions(regions, gs$fullGenome)
print(countTable <- colSums(ensemblAnno$countTable))

### GeneMap: Check regions for overlapping genes
print("Loading GeneMap_Hsapiens.gencode.v19.GRCh37.rds..."); Sys.time()
geneMap <- readRDS(paste0(REFS, "GeneMap_Hsapiens.gencode.v19.GRCh37.rds"))
geneMapGR <- GenomicRanges::makeGRangesFromDataFrame(
  geneMap, 
  keep.extra.columns = TRUE
); rm(geneMap)

print("adding nearest gene..."); Sys.time()
dA <- GenomicRanges::distanceToNearest(regions, geneMapGR)
regions$feature_name <- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_symbol<- geneMapGR$gene_symbol[subjectHits(dA)]
regions$gene_type <- geneMapGR$gene_type[subjectHits(dA)]
regions$gene_id <- geneMapGR$gene_id[subjectHits(dA)]
regions$distToGene <- mcols(dA)$distance
mcols(regions) <- cbind(mcols(regions), ensemblAnno$countTable)
print("adding nearest gene...DONE"); Sys.time()

### Add additional annotation
print("Adding anno class"); Sys.time()
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
print("Exporting regions as bed file."); Sys.time()
rtracklayer::export(regions, 
                    paste0("derfinder_", brain_region, "_ERs.bed"))

### Extract featureData (annotation) from granges object
features <- data.frame(regions) %>%
  mutate(chr = seqnames) %>%
  mutate(pos_hg19 = paste0(chr, ":", start, "-", end),
         feature_id = paste0(chr, ":", start, "-", end)) %>% 
  select(feature_id, feature_name, 
         pos_hg19, chr, start, end, strand, 
         gene_symbol, gene_type, gene_id, 
         distToGene, anno_class)
summary(factor(features$strand))

### feature_id = Region (chr_:_-_) = pos_hg19 = rownames
features <- as.data.frame(features)
rownames(features) <- features$feature_id

### Save
print("Saving annotation..."); Sys.time()
saveRDS(features, paste0("features_", brain_region, "_regions.rds")) 

### Done
print("Annotation complete."); Sys.time()


### Session info
print("Script finished.")
library(sessioninfo)
sessioninfo::session_info()
Sys.time()
proc.time()

