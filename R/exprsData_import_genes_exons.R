### Stanley ACC SZ vs CTRL
### Import exprsData, featureData, and phenoData 
### Combine into DGEs
### Genes & exons
### Alexis Norris
### Created: 2016-08-27
### Modified: 2018-06-01



# Input files --------------------------------------------------
### exprsData (counts) 
### from `featureCounts`: raw counts; contains:
#   geneCounts - exprsData (matrix)      
#   geneMap    - featureData (empty, use my own one in refs/)    
#   exonCounts - exprsData (matrix)   
#   exonMap    - featureData (only e#s, merge with my own one in refs/)    
#   pd         - phenoData (only filenames) 
#load("input/counts/rawCounts_n20_genes_exons_*.rda")
#load("input/counts/rawCounts_n20_genes_exons_*.rda")

### featureData
### TranscriptomeMap(for exons)
#load("~/References/TranscriptomeMaps/Hsapiens.gencode.v19.GRCh37.p13.rda")
### GeneMap (for genes) 
#features <- readRDS("input/refs/GeneMap_Hsapiens.gencode.v19.GRCh37.rds")

### phenoData
#phenoList <- readRDS("input/pheno/phenoData_full_list.rds")


# Setup --------------------------------------------------------
### Packages
library(GenomicRanges)

### Parameters
### All brain regions run in the same script (unlike jxns and regions)
analysis_name <- "exprsData_import_genes_exons"


# Load phenoData -----------------------------------------------
### Load
phenoList <- readRDS("input/pheno/phenoData_full_list.rds")

### Subset
pheno <- phenoList[[brain_region]]
rownames(pheno) <- pheno$FileName_counts 


# Prepare exon featureData -------------------------------------
### get exonMap from TranscriptomeMap 
load("~/References/TranscriptomeMaps/Hsapiens.gencode.v19.GRCh37.p13.rda")
txMap <- data.frame(
  "type" = map$type,
  "feature_id" = map$exon_id,
  "exon_id" = map$exon_id,
  "exon_number" = map$exon_number, # makes duplicates
  "gene_symbol" = map$gene_name, 
  "feature_name" = map$gene_name, 
  "gene_id" = map$gene_id, 
  "gene_type" = map$gene_type,
  "chr" = seqnames(map), 
  "start" = start(map), 
  "end" = end(map), 
  "strand" = strand(map),
  "pos_hg19" = paste0(seqnames(map), ":", 
                      start(map), "-", end(map))
)
summary(factor(txMap$type))
exonMap <- txMap[txMap$type == "exon", ]
exonMap$type <- NULL
exonMap <- exonMap[!duplicated(exonMap), ]
saveRDS(exonMap, "input/refs/ExonMap_Hsapiens.gencode.v19.GRCh37.rds")

### Add with old exonMap, to get e# (feature_id in counts) 
features_counts$Chr <- str_replace_all(features_counts$Chr, "chrchr", "chr")
features_counts$pos_hg19 <- paste0(
  features_counts$Chr, ":", 
  features_counts$Start, "-", 
  features_counts$End
)
features_counts$e_id <- rownames(features_counts)
features <- features_counts %>%
  mutate(gene_id = Geneid) %>%
  select(gene_id, e_id, pos_hg19) %>%
  plyr::join(
    exonMap, 
    by = c("pos_hg19", "gene_id"), 
    type = "left",   # want everything in counts data
    match = "first"  # still 512 dups in exonMap
  ) %>% as.data.frame()

### Make e_id (e#) rowname for now, since that's what counts rownames are
features$feature_id <- features$e_id # not exon_id, because then dups exist!
rownames(features) <- features$feature_id

### Save
saveRDS(features, "input/refs/ExonMap_Hsapiens.gencode.v19.GRCh37_with_eID_.rds")


# Prepare gene featureData -------------------------------------
### Load geneMap
features <- readRDS("input/refs/GeneMap_Hsapiens.gencode.v19.GRCh37.rds")

### Make gene_id (ENSG) --> rownames & feature_id
features$feature_id <- features$gene_id
rownames(features) <- features$feature_id
features$feature_name <- features$gene_symbol

### Add pos_hg19 column (match Regions)
features$pos_hg19 <- paste0(
  features$chr, ":", 
  features$start, "-", 
  features$end
)
saveRDS(features, "input/refs/GeneMap_Hsapiens.gencode.v19.GRCh37_features.rds")


# Load counts --------------------------------------------------
### ACC
load("input/counts/rawCounts_n20_genes_exons_acc.rda") 
geneCounts_acc <- geneCounts
exonCounts_acc <- exonCounts
colnames(geneCounts_acc) <- plyr::revalue(colnames(geneCounts_acc), ACC_key) # fix
colnames(exonCounts_acc) <- plyr::revalue(colnames(exonCounts_acc), ACC_key) # fix

### PFC
load("input/counts/rawCounts_n20_genes_exons_pfc.rda") 
geneCounts_pfc <- geneCounts
exonCounts_pfc <- exonCounts

### HPC
load("input/counts/rawCounts_n20_genes_exons_hpc.rda") 
geneCounts_hpc <- geneCounts
exonCounts_hpc <- exonCounts


# Combine all --------------------------------------------------
### Check
stopifnot(rownames(geneCounts_acc) == rownames(geneCounts_pfc)) 
stopifnot(rownames(geneCounts_acc) == rownames(geneCounts_hpc))
stopifnot(rownames(exonCounts_acc) == rownames(exonCounts_pfc)) 
stopifnot(rownames(exonCounts_acc) == rownames(exonCounts_hpc))

### cbind
geneCounts <- cbind(geneCounts_acc, geneCounts_pfc, geneCounts_hpc)
exonCounts <- cbind(exonCounts_acc, exonCounts_pfc, exonCounts_hpc)
save(exonCounts, geneCounts, file = "input/counts/rawCounts_n20_genes_exons.rda", compress = TRUE)


# Check sample IDs ---------------------------------------------
stopifnot(colnames(geneCounts) == colnames(exonCounts))
stopifnot(colnames(geneCounts) %in% pheno$sample_id)
stopifnot(pheno$sample_id %in% colnames(geneCounts))
rownames(pheno) <- pheno$sample_id


# Create dge (genes) -------------------------------------------
### Match feature order: exprsData rownames = featureData rownames
geneCounts <- geneCounts[which(rownames(geneCounts) %in% rownames(geneMap)), ]
geneMap <- geneMap[match(rownames(geneCounts), rownames(geneMap)), ]
stopifnot(rownames(geneCounts) == rownames(geneMap))

### Match sample order: exprsData colnames = phenoData rownames
pheno <- pheno[which(rownames(pheno) %in% colnames(geneCounts)), ]
geneCounts <- geneCounts[ , match(rownames(pheno), colnames(geneCounts))]
stopifnot(colnames(geneCounts) == rownames(pheno))

### Change Dx column name from group
pheno$Dx <- pheno$group
pheno$group <- NULL

### Combine into dge
dge <- edgeR::DGEList(
  counts = geneCounts,
  genes = geneMap,
  samples = pheno,
  group = pheno$Dx # or else it makes it col1 (sample_ids); but if already have "group" col, then you end up with columns "group" AND "group.1"
)


# Create dge (exons) -------------------------------------------
### Match feature order: exprsData rownames = featureData rownames
exonCounts <- exonCounts[which(rownames(exonCounts) %in% rownames(exonMap)), ]
exonMap <- exonMap[match(rownames(exonCounts), rownames(exonMap)), ]
table(rownames(exonCounts) == rownames(exonMap))

### Match sample order: exprsData colnames = phenoData rownames
pheno <- pheno[which(rownames(pheno) %in% colnames(exonCounts)), ]
exonCounts <- exonCounts[ , match(rownames(pheno), colnames(exonCounts))]
table(colnames(exonCounts) == rownames(pheno))

### Combine into dge
dge <- edgeR::DGEList(
  counts = exonCounts,
  genes = exonMap,
  samples = pheno,
  group = pheno$Dx
)


# Save ---------------------------------------------------------
saveRDS(dge, "input/dge_all_genes.rds")
saveRDS(dge, "input/dge_all_exons.rds")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
