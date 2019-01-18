#!/usr/bin/env r

### Alexis Norris
### created: 1 August 2016
### modified: 3 March 2018 *for HPC, which needs more threads*

print(Sys.time())
library(derfinder)

### Generate full coverage
# Load chromosome information used in assembly
chrom.sizes <- read.delim("/mnt/data/RNASeq/methods/references/hg19/hg19.chrom.sizes", header = FALSE)
print(head(chrom.sizes))

# Remove unmapped contigs
chrom.sizes <- chrom.sizes[-grep("_", chrom.sizes$V1), ]


# ACC (75bp) ------------------------------------------
# Locate bigwig files
# .bai files must also be in the directory and have identical names
#files_ACC <- rawFiles("/mnt/data/RNASeq/data/Stanley/all_hisat_hg19/derfinder/bigwigs/ACC", samplepatt = "bw", fileterm = NULL)
#names(files_ACC) <- gsub(".bw", "", names(files_ACC))

# Load full coverage from bigwig files
#fullCov_ACC <- fullCoverage(files = files_ACC, chrs = as.character(chrom.sizes$V1), chrlens = chrom.sizes$V2, verbose = TRUE, mc.cores.load = 10)
#saveRDS(fullCov_ACC, "fullCov_ACC.rds")

### Filter coverage to reduce computational time of `regionMatrix`
# Using mean coverage of 5 reads (mean across all samples), for a region to be included as an expressed region  
#filteredCov_ACC <- lapply(fullCov_ACC, filterData, returnMean = TRUE, filter = "mean", cutoff = 5, verbose = TRUE)
#saveRDS(filteredCov_ACC, "filteredCov_ACC.rds")

### Calculate expressed region counts  
# L = read length (here, 75bp)
#regionMat75_ACC <- regionMatrix(filteredCov_ACC, L = 75, runFilter = FALSE, verbose = TRUE)
#saveRDS(regionMat75_ACC, "regionMat_ACC_75bp.rds")
#rm(fullCov_ACC)
#rm(filteredCov_ACC)
#rm(regionMat75_ACC)

# HPC (100bp) ------------------------------------------
# Locate bigwig files
# .bai files must also be in the directory and have identical names
#files_HPC <- rawFiles("/mnt/data/RNASeq/data/Stanley/all_hisat_hg19/derfinder/bigwigs/HPC", samplepatt = "bw", fileterm = NULL)
#names(files_HPC) <- gsub(".bw", "", names(files_HPC))

# Load full coverage from bigwig files
#fullCov_HPC <- fullCoverage(files = files_HPC, chrs = as.character(chrom.sizes$V1), chrlens = chrom.sizes$V2, verbose = TRUE, mc.cores.load = 10)
#saveRDS(fullCov_HPC, "fullCov_HPC.rds")

### Filter coverage to reduce computational time of `regionMatrix`
# Using mean coverage of 5 reads (mean across all samples), for a region to be included as an expressed region  
#filteredCov_HPC <- lapply(fullCov_HPC, filterData, returnMean = TRUE, filter = "mean", cutoff = 5, verbose = TRUE)
#saveRDS(filteredCov_HPC, "filteredCov_HPC.rds")

# read in previously saved
filteredCov_HPC <- readRDS("filteredCov_HPC.rds")

# HPC previously freezed here, at chr7, last message was:
#2018-02-28 05:13:47 getRegionCoverage: processing chr7

### Calculate expressed region counts  
# L = read length (here, 100bp)
regionMat100_HPC <- regionMatrix(
	filteredCov_HPC, 
	L = 100, 
	runFilter = FALSE, 
	returnBP = FALSE, # don't need this, and might speed it up?!
	verbose = TRUE)

saveRDS(regionMat100_HPC, "regionMat_HPC_100bp.rds")

#rm(fullCov_HPC)
#rm(filteredCov_HPC)
#rm(regionMat100_HPC)

# PEC (100bp) ------------------------------------------
# Locate bigwig files
# .bai files must also be in the directory and have identical names
#files_PFC <- rawFiles("PFC", samplepatt = "bw", fileterm = NULL)
#names(files_PFC) <- gsub(".bw", "", names(files_PFC))

# Load full coverage from bigwig files
#fullCov_PFC <- fullCoverage(files = files_PFC, chrs = as.character(chrom.sizes$V1), chrlens = chrom.sizes$V2, verbose = TRUE, mc.cores.load = 10)
#saveRDS(fullCov_PFC, file = "fullCov_PFC.rds")

### Filter coverage to reduce computational time of `regionMatrix`
# Using mean coverage of 5 reads (mean across all samples), for a region to be included as an expressed region  
#filteredCov_PFC <- lapply(fullCov_PFC, filterData, returnMean = TRUE, filter = "mean", cutoff = 5, verbose = TRUE)
#saveRDS(filteredCov_PFC, file = "filteredCov_PFC.rds")

### Calculate expressed region counts  
# L = read length (here, 100bp)
#regionMat100_PFC <- regionMatrix(filteredCov_PFC, L = 100, runFilter = FALSE, verbose = TRUE)
#saveRDS(regionMat100_PFC, file = "regionMat_PFC_100bp.rds")
#rm(fullCov_PFC)
#rm(filteredCov_PFC)
#rm(regionMat100_PFC)

library(sessioninfo)
sessioninfo::session_info()
Sys.time()
proc.time()

