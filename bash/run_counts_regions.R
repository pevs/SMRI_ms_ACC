#!/usr/bin/env r

## Alexis Norris
## Created: 2016-08-01
## Modified: 2018-02-21
### Run in bash script with "module add R/3.2.3"


# Setup --------------------------------------------------------
### Packages
library(derfinder)

### Log
Sys.time()


# Load chromosome sizes ----------------------------------------
### Load chromosome information used in assembly
chrom.sizes <- read.delim("hg19.chrom.sizes", header = FALSE)

### Remove unmapped contigs
chrom.sizes <- chrom.sizes[-grep("_", chrom.sizes$V1), ]


# Load coverage data -------------------------------------------
### Locate bigwig files
### .bai files must also be in the directory and have identical names
bw_files <- rawFiles(".", samplepatt = "bw", fileterm = NULL)
names(bw_files) <- gsub(".bw", "", names(bw_files))

### Load full coverage from bigwig files
fullCov <- fullCoverage(
  files = bw_files, 
  chrs = as.character(chrom.sizes$V1), 
  chrlens = chrom.sizes$V2, 
  verbose = TRUE, 
  mc.cores.load = 10
)


# Filter -------------------------------------------------------
### Filter coverage to reduce computational time of `regionMatrix`
### Using mean coverage of 5 reads (mean across all samples), for a region to be included as an expressed region  
filteredCov <- lapply(
  fullCov, 
  filterData, 
  returnMean = TRUE, 
  filter = "mean", 
  cutoff = 5, 
  verbose = TRUE
)


# Get counts ---------------------------------------------------
### Read length for dataset (ACC = 75bp, PFC = 100bp, HPC = 100bp)
read_length <- 100

### Calculate expressed region counts  
regionMat <- regionMatrix(
  filteredCov, 
  L = read_length, 
  runFilter = FALSE, 
  verbose = TRUE
)


# Save ---------------------------------------------------------
### Full coverage
saveRDS(fullCov_ACC, "fullCov.rds")

### Filtered coverage
saveRDS(filteredCov, "filteredCov.rds")

### Final exprsData (feature [region] counts)
saveRDS(regionMat, file = "regionMat.rds")


# Log ----------------------------------------------------------
Sys.time()
