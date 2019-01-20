### Stanley ACC SZ vs CTRL
### Import exprsData, featureData, and phenoData 
### Combine into DGEs
### Junctions
### Alexis Norris
### Created: 2016-08-27
### Modified: 2018-06-01



#### NOT RUN!!! PATHS NEED UPDATED, but this the current ones to find the data on external HD
#### Need to add phenoData too "phenoList"
### Specific for ACC???


# Input files --------------------------------------------------
### exprsData (counts) 
#regionMat <- readRDS(paste0("input/counts/regionMat_", brain_region, "_", read_length, "bp.rds"))

### featureData
### Load from cluster
#features <- readRDS(paste0("input/refs/features_", brain_region, "_regions.rds"))

### phenoData
#phenoList <-


# Setup --------------------------------------------------------
### Packages
library(GenomicRanges)

### Parameters
brain_region <- "ACC"
read_length <- 75
analysis_name <- paste0("exprsData_import_regions_", brain_region)


# [NEED TO ADD] Load phenoData -----------------------------------------------
### Load
#phenoList <- 

### Subset
pheno <- phenoList[[brain_region]]
rownames(pheno) <- pheno$FileName_counts 


# Load counts --------------------------------------------------
### exprsData
regionMat <- readRDS(paste0("input/counts/regionMat_", brain_region, "_", read_length, "bp.rds"))

### Extract regions from GRanges object 
regions <- unlist(GRangesList(lapply(regionMat, "[[", "regions")))
names(regions) <- paste0(
  seqnames(regions), ":", 
  start(regions), "-", 
  end(regions)
) 

### Extract counts from GRanges object
counts <- do.call("rbind", lapply(regionMat, "[[", "coverageMatrix"))

### Put sample columns in order
### First check that FileName used for rownames(phenoData) was correct
### For ACC, use ACC_key
if (brain_region == "ACC") {
  colnames(counts) <- plyr::revalue(colnames(counts), ACC_key) 
}
colnames(counts) <- plyr::revalue(colnames(counts), ACC_key) 
stopifnot(colnames(counts) %in% rownames(pheno))
stopifnot(rownames(pheno) %in% colnames(counts))
counts <- counts[ , rownames(pheno)]
rownames(counts) <- names(regions)


# featureData --------------------------------------------------
### Load from cluster
features <- readRDS(paste0("input/refs/features_", brain_region, "_regions.rds"))
stopifnot(rownames(counts) %in% rownames(features))
stopifnot(rownames(features) %in% rownames(counts))


# Create dge ---------------------------------------------------
### Match feature order: exprsData rownames = featureData rownames
counts <- counts[which(rownames(counts) %in% rownames(features)), ]
features <- features[match(rownames(counts), rownames(features)), ]
stopifnot(rownames(counts) == rownames(features))

### Match sample order: exprsData colnames = phenoData rownames
pheno <- pheno[which(rownames(pheno) %in% colnames(counts)), ]
counts <- counts[ , match(rownames(pheno), colnames(counts))]
stopifnot(colnames(counts) == rownames(pheno))

### Set sample_id to rowname
rownames(pheno) <- pheno$sample_id
colnames(counts) <- rownames(pheno)

### Change Dx column name from group
pheno$Dx <- pheno$group
pheno$group <- NULL

### Combine into dge
dge <- edgeR::DGEList(
  counts = counts,
  genes = features,
  samples = pheno,
  group = pheno$Dx
)


# Save ---------------------------------------------------------
saveRDS(dge, "input/dge_", brain_region, "_regions.rds") 


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
