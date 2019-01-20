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
### from `featureCounts`: raw counts; contains:
#   jCounts - exprsData (matrix)      
#   jMap    - featureData  
#   pd      - phenoData (only filenames) 
#load("input/counts/rawCounts_n20_jxns_acc.rda") 

### featureData
### GeneMap 
#geneMap <- readRDS("input/refs/GeneMap_Hsapiens.gencode.v19.GRCh37.rds")

### phenoData
#phenoList <-


# Setup --------------------------------------------------------
### Packages
library(GenomicRanges)

### Parameters
brain_region <- "ACC"
analysis_name <- paste0("exprsData_import_jxns_", brain_region)


# [NEED TO ADD] Load phenoData -----------------------------------------------
### Load
#phenoList <- 

### Subset
pheno <- phenoList[[brain_region]]
rownames(pheno) <- pheno$FileName_counts 

# Load counts --------------------------------------------------
### exprsData (counts)
load("input/counts/rawCounts_n20_jxns_acc.rda") 

### Only expressed
counts <- jCounts[which(rowSums(jCounts) > 0), ] 

### Check that FileName used for rownames(phenoData) was correct
### For ACC, use ACC_key
if (brain_region == "ACC") {
  colnames(counts) <- plyr::revalue(colnames(counts), ACC_key) 
}
stopifnot(colnames(counts) %in% rownames(pheno))
stopifnot(rownames(pheno) %in% colnames(counts))


# featureData --------------------------------------------------
### Extract 
features <- as.data.frame(jMap) 

### names --> feature_id & pos_hg19
features$feature_id <- names(jMap)

### Remove strand info
features$pos_hg19 <- gsub("\\(\\*\\)", "", features$feature_id) 


### Get gene_symbol using gene_id (ENSG) 
geneMap <- readRDS("input/refs/GeneMap_Hsapiens.gencode.v19.GRCh37.rds")
features$gene_id <- features$newGeneID
table(!is.na(features$gene_id) %in% geneMap$gene_id)
table(duplicated(geneMap$gene_id))
features <- geneMap %>%
  mutate(gene_strand = strand) %>%
  select(gene_id, gene_symbol, gene_type, gene_strand) %>%
  plyr::join(
    features,
    by = "gene_id",
    "right", "all"
  ) %>% 
  mutate(
    feature_name = gene_symbol,
    chr = seqnames
  ) %>%
  select(-one_of(c("newGeneID", "newGeneSymbol", "seqnames")))

### Rownames = feature_id (must be df, not tibble!)
features <- as.data.frame(features)
rownames(features) <- features$feature_id

### Save
saveRDS(features, paste0("input/refs/jxnMap_", brain_region, ".rds"))


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
  group = pheno$Dx # or else it makes it col1 (sample_ids); but if already have "group" col, then you end up with columns "group" AND "group.1"
)


# Save ---------------------------------------------------------
saveRDS(dge, paste0("input/dge_", brain_region, "_jxns.rds")) 


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
