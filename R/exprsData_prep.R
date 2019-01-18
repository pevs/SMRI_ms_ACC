### SMRI Array ACC SZ vs CTRL
### Combine DGEs (exprsData) for all feature levels and brain regions
### Remove samples not in final cohort
### Remove features not expressed in any sample in the final cohort
### Alexis Norris
### Created: 2018-11-09
### Modified: 2019-01-03


# Notes --------------------------------------------------------
### Libsizes
# TotalReads = total # reads sequenced
# LibSize = total # reads aligned
# lib.size = edgeR calculation of reads (*or relative reads?*) in dge (so less when filter out low abundance genes/features)

### Not all samples in ACC also have data in PFC and/or HPC; 
### and there are samples with PFC and/or HPC data that do not have PFC data
### Restricting to PFC/HPC samples in ACC; but not restricting ACC to only samples with PFC and HPC!


# Input files --------------------------------------------------
### PhenoData cleaned up and Pevs Codes added
#pheno_final <- read_tsv("input/phenoData/phenoData_array_ACC.txt") %>%

### DGEs created in Stanley_Array_ACC_PFC_HPC 
### ! Use polyA version of HPC DGEs -- the ribozero was erroneously used originally
### hpc_fix_qsva_regions.R fixed that
#dgeList <- list(
# "ACC_genes" = readRDS("input/dge/dge_ACC_genes_polyA.rds"),
# "ACC_exons" = readRDS("input/dge/dge_ACC_exons_polyA.rds"),
# "ACC_jxns" = readRDS("input/dge/dge_ACC_jxns_polyA.rds"),
# "ACC_regions" = readRDS("input/dge/dge_ACC_regions_polyA.rds"),
# "PFC_genes" = readRDS("input/dge/dge_PFC_genes_ribozero.rds"),
# "PFC_exons" = readRDS("input/dge/dge_PFC_exons_ribozero.rds"),
# "PFC_jxns" = readRDS("input/dge/dge_PFC_jxns_ribozero.rds"),
# "PFC_regions" = readRDS("input/dge/dge_PFC_regions_ribozero.rds"),
# "HPC_genes" = readRDS("input/dge/dge_HPC_genes_polyA.rds"),
# "HPC_exons" = readRDS("input/dge/dge_HPC_exons_polyA.rds"),
# "HPC_jxns" = readRDS("input/dge/dge_HPC_jxns_polyA.rds"),
# "HPC_regions" = readRDS("input/dge/dge_HPC_regions_polyA.rds")
#)

### RepeatMasker anno (for region-level features)
#rmsk <- readRDS("input/anno/rmsk_hg19.rds")

### Load HGNC anno (from anno_hgnc.R)
#hgncList <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")

### Cleaned up Array phenoData for ACC final cohort
#pheno_final <- read_tsv("input/pheno/phenoData_array_ACC.txt")

# Setup --------------------------------------------------------
### Packages
library(edgeR)     # DGEs of exprs, pheno, and feature/anno data
library(tidyverse)      # plotting; wrangling; select fxn gets masked by GenomicRanges
library(GenomicRanges)  # rmsk anno

### Parameters
analysis_name <- "exprsData_prep"
datasets <- c("ACC", "PFC", "HPC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_")


# Load all dges ------------------------------------------------
### brain regions: ACC, PFC, HPC
### feature levels: genes, exons, jxns, regions
### Original data source: /Volumes/Hopkins2018/Projects/Human/Stanley/ArrayCollection/ACC_PFC_HPC/20180605/input
dgeList_all <- list(
  "ACC_genes" = readRDS("input/dge/dge_ACC_genes_polyA.rds"),
  "ACC_exons" = readRDS("input/dge/dge_ACC_exons_polyA.rds"),
  "ACC_jxns" = readRDS("input/dge/dge_ACC_jxns_polyA.rds"),
  "ACC_regions" = readRDS("input/dge/dge_ACC_regions_polyA.rds"),
  "PFC_genes" = readRDS("input/dge/dge_PFC_genes_ribozero.rds"),
  "PFC_exons" = readRDS("input/dge/dge_PFC_exons_ribozero.rds"),
  "PFC_jxns" = readRDS("input/dge/dge_PFC_jxns_ribozero.rds"),
  "PFC_regions" = readRDS("input/dge/dge_PFC_regions_ribozero.rds"),
  "HPC_genes" = readRDS("input/dge/dge_HPC_genes_polyA.rds"),
  "HPC_exons" = readRDS("input/dge/dge_HPC_exons_polyA.rds"),
  "HPC_jxns" = readRDS("input/dge/dge_HPC_jxns_polyA.rds"),
  "HPC_regions" = readRDS("input/dge/dge_HPC_regions_polyA.rds")
)

### Want all dges
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# Load final ACC phenoData -------------------------------------
### Cols to remove from generic phenoData ("pheno_final" below), but to keep in each dataset's phenoData
### (These are dataset-specific values that we want to keep for each dataset)
specific_cols <- c(
  "FileName",
  "lib.size", "norm.factors",    # edgeR calculations
  "RIN", 
  "Batch", "ReplicateBatch", 
  "rnaIsolation_batch", "Flowcell", "Instrument",
  "TotalLibSize", "AlignStats", "AlignStats_hisat2",
  "rRNA_rate", "Reads_chrM",
  "gene_assignment_rate", "exon_assignment_rate"
)

### Only includes samples used for final ACC analysis in paper
### Filtered samples (e.g., B3, B5, only 1yr SZ duration, and non-EA) not in this pheno
pheno_final <- read_tsv("input/pheno/phenoData_array_ACC.txt") %>%
  # Remove ACC-specific cols for sample filenames -- just use actual FileName (dge$samples' rownames)
  # (these are ACC-specific since this phenoData originated from ACC phenoData)
  # Change ACC FileName col to match format of PFC & HPC --> "ACC_FileName"
  mutate("ACC_FileName" = FilenameFull) %>%
  dplyr::select(-one_of(c(
    # Redundant with A_Code
    "A_Num",
    # ACC-specific
    "FilenameFull", "FilenameShort", "FilenameAbbrev", "TotalReads", "rRNA",
    # Dataset specific, but want to keep for dges (see above)
    specific_cols                                      
  ))) %>%
  # rearrange
  dplyr::select(Sample, ACC_FileName, PFC_FileName, HPC_FileName, everything())


# Filter Samples -----------------------------------------------
### 1. Remove samples not in final ACC cohort
### 2. Update sample IDs to those coded ones used in paper (add PevsCode), by merging pheno_final/sub with y$samples
### 3. Remove unexpressed features

### Run for each dataset
### use for loop; lapply hits memory limit
for (i in names(dgeList)) {
  # Subset dge
  y <- dgeList[[i]]
  
  # Subset pheno for samples with data for region
  # make filename the one for that region, to join later with y$samples
  pheno_sub <- bind_rows(
    filter(pheno_final, ACC_FileName %in% rownames(y$samples)) %>% mutate(FileName = ACC_FileName),
    filter(pheno_final, PFC_FileName %in% rownames(y$samples)) %>% mutate(FileName = PFC_FileName),
    filter(pheno_final, HPC_FileName %in% rownames(y$samples)) %>% mutate(FileName = HPC_FileName)
  )
  
  # Remove samples that aren't in final sample cohort used in paper, reorder to match pheno_final/sub
  y <- y[ , pheno_sub$FileName, keep.lib.sizes = TRUE]
  
  # Check
  stopifnot(rownames(y$samples) == pheno_sub$FileName)
  
  # Add to existing phenoData (e.g., add PevsCode) --> use left join so that row order in y doesn't change!
  y$samples <- cbind.data.frame(
    y$samples[ , which(names(y$samples) %in% specific_cols)],
    dplyr::select(pheno_sub, -FileName),
    y$samples[ , grep("degchr", names(y$samples))] # add to end
  )
  
  # Change rownames to annonymized "Sample"
  stopifnot(rownames(y$samples) == colnames(y$counts)) # check first
  rownames(y$samples) <- y$samples$Sample
  colnames(y$counts) <- rownames(y$samples)
  
  # Remove unexpressed features
  y <- y[rowSums(y$counts) > 0, , keep.lib.sizes = FALSE]
  
  # Return
  dgeList[[i]] <- y
}
saveRDS(dgeList, "temp/dge_ACC_PFC_HPC_temp2.rds"); beepr::beep()



# Get repeatMasker overlap for regions --------------------------
### Load region featureData
### These are low confidence annotation because of mapping impaired by repetitive nature
features <- lapply(
  # Regions only
  dgeList[grep("_regions", names(dgeList))], function(y) y$genes
) %>%
  # Combine all datasets
  bind_rows() %>%
  # Remove duplicates
  distinct() %>%
  # Only position columns and feature_id
  dplyr::select(chr, start, end, strand, feature_id) %>%
  # Convert to GR
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
# Make feature_id the name/rowname of GR object
stopifnot(!duplicated(features$feature_id))

### Load repeatMasker anno
rmsk <- readRDS("input/anno/rmsk_hg19.rds")
# Add prefix to rmsk IDs
names(mcols(rmsk)) <- paste0("rmsk_", names(mcols(rmsk)))

### Get overlaps
hits <- findOverlaps(                     
  features, rmsk,
  ignore.strand = TRUE,         # ACC is unstranded
  type = "any", select = "all"
)
overlaps <- pintersect(
  features[queryHits(hits)], 
  rmsk[subjectHits(hits)]
)
# Add percent overlap
overlaps$rmsk_overlap <- width(overlaps)/width(features[queryHits(hits)])
# Add rmsk anno 
mcols(overlaps) <- c(
  mcols(overlaps),
  mcols(rmsk[subjectHits(hits)])
)
overlaps$rmsk_pos_hg19 <- paste0(
  seqnames(rmsk[subjectHits(hits)]), ":",
  start(rmsk[subjectHits(hits)]), "-", end(rmsk[subjectHits(hits)])
)
overlaps$rmsk_strand <- as.character(strand(rmsk[subjectHits(hits)]))

### Convert to df
overlaps_df <- overlaps %>%
  as_data_frame() %>%
  # Only columns that will be added to DGEs' featureData
  dplyr::select(
    feature_id,
    rmsk_overlap, 
    rmsk_repName, rmsk_repClass, rmsk_repFamily, 
    rmsk_pos_hg19, rmsk_strand                       
  ) %>%
  # sort by decreasing % overlap
  arrange(-rmsk_overlap)
# Save full
write_tsv(overlaps_df, "logs/exprs/ACC_PFC_HPC_rmskOverlap.tsv")
# For featureData, only keep the first (most overlapping) repeat for each feature_id
overlaps_df <- overlaps_df[!duplicated(overlaps_df$feature_id), ]
# Add back region features that have no overlap with rmsk
overlaps_df <- overlaps_df %>%
  bind_rows(data_frame(
    "feature_id" = features$feature_id[!(features$feature_id %in% overlaps_df$feature_id)],
    "rmsk_overlap" = 0
  )) %>%
  as.data.frame()
# Fix classes
overlaps_df$rmsk_repName <- as.character(overlaps_df$rmsk_repName)
overlaps_df$rmsk_repClass <- as.character(overlaps_df$rmsk_repClass)
overlaps_df$rmsk_repFamily <- as.character(overlaps_df$rmsk_repFamily)

### Add rmsk overlap anno to featureData
stopifnot(!is.na(overlaps_df$feature_id))
stopifnot(!duplicated(overlaps_df$feature_id))
stopifnot(features$feature_id %in% overlaps_df$feature_id)
for (dge in paste0(datasets, "_regions")) {
  print(dge)
  # Subset DGE
  y <- dgeList[[dge]]
  
  # Combine
  new_anno <- left_join(
    y$genes, overlaps_df, 
    by = "feature_id"
  )
  
  # Add rownames
  stopifnot(new_anno$feature_id == rownames(y$counts))
  rownames(new_anno) <- new_anno$feature_id
  
  # Return
  dgeList[[dge]]$genes <- new_anno
}
saveRDS(dgeList, "temp/dge_ACC_PFC_HPC_temp3.rds"); beepr::beep()


# Update outdated HGNC gene symbols ----------------------------
### Check that gene_symbol is the current approved HGNC symbol
### if outdated (alias), update
### If ensembl ID not in HGNC, remove gene_symbol (make NA)

### Load HGNC anno (from anno_hgnc.R)
hgnc <- readRDS("input/anno/hgnc_complete_set_ftp_key.rds")$by_ensembl
names(hgnc) <- recode(
  names(hgnc),
  "symbol" = "approved_symbol"
)

### Update featureData
### Note some novel jxn/read-through/fusion features at jxn-level
for (dge in names(dgeList)) {
  # Subset featureData
  print(dge)
  y <- dgeList[[dge]]
  ### Account for novel jxn/read-through/fusion features at jxn-level
  if (dge %in% paste0(datasets, "_jxns")) {
    # Get featureData
    df <- y$genes %>%
      mutate(
        # id column (ensembl ID)
        id = feature_id,
        # column to match on; need to remove version suffix --> will remove fusion partner!
        by_col = gsub("\\..*", "", gene_id),
        # Column indicating original/uncorrected gene symbol
        gene_symbol_uncorrected = gene_symbol
      ) %>%
      # Fix classes, for binding later
      mutate_if(is.factor, as.character)
    
    # Find fusions
    two_ensembl <- df$gene_id[grep("-ENS", df$gene_id)]
    
    # Just fusions
    # Will add this to no_ensembl df later/below
    fusions <- df %>%
      # Has 2 ensembl ids (is a fusion/novel jxn)
      filter(gene_id %in% two_ensembl) %>%
      # Note that it was changed
      mutate(gene_symbol_corrected_reason = "novel fusion (2 ensembl IDs)")
    # Write out fusions (since not annotating; they wouldn't be in HGNC)
    write_tsv(fusions, paste0("logs/exprsData/fusion_jxns_", dge, ".tsv"))
    
    # Remove fusions from df
    df <- df %>%
      filter(!(gene_id %in% fusions$gene_id))
  } else {
    df <- y$genes %>%
      mutate(
        # id column (ensembl ID)
        id = feature_id,
        # column to match on; need to remove version suffix
        by_col = gsub("\\..*", "", gene_id),
        # Column indicating original/uncorrected gene symbol
        gene_symbol_uncorrected = gene_symbol
      ) %>%
      # Fix classes, for binding later
      mutate_if(is.factor, as.character)
  }
  
  ### Don't have ensembl ID (many jxns)
  no_ensembl <- df %>%
    filter(is.na(by_col)) %>%
    # Replace with NA
    mutate(gene_symbol = NA) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "no ensembl ID")

  ### Join by ensembl IDs
  df_hgnc <- df %>%
    filter(!is.na(by_col)) %>%
    left_join(hgnc, by = c("by_col" = "ensembl_gene_id")) %>%
    # Fix classes, for binding later
    mutate_if(is.factor, as.character)
  
  ### Outcome 1. Ok (is current approved symbol)
  ok <- df_hgnc %>%
    filter(!is.na(by_col)) %>%
    # ensembl ID is in HGNC
    # and symbols match
    filter(gene_symbol == approved_symbol) %>%
    # no change
    mutate(gene_symbol_corrected_reason = "OK; uncorrected")
  
  ### Outcome 2. Updated gene symbol (was outdated/alias)
  updated <- df_hgnc %>%
    filter(!is.na(by_col)) %>%
    # ensembl ID is in HGNC
    # but symbols don't match
    filter(gene_symbol != approved_symbol) %>%
    # Replace with approved
    mutate(gene_symbol = approved_symbol) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "symbol outdated")
  
  ### Outcome 3. old (ensembl gene id not in HGNC )
  old_ensembl <- df_hgnc %>%
    filter(!is.na(by_col)) %>%
    # ensembl ID is NOT in HGNC
    filter(is.na(approved_symbol)) %>%
    # Replace with NA
    mutate(gene_symbol = NA) %>%
    # Note that it was changed
    mutate(gene_symbol_corrected_reason = "ensembl ID not in HGNC")
  
  ### Combine outcomes
  df_anno <- bind_rows(
    ok,
    updated,
    old_ensembl,
    no_ensembl
  ) 
  # Account for novel jxn/read-through/fusion features at jxn-level (from above)
  if (dge %in% paste0(datasets, "_jxns")) {
    df_anno <- bind_rows(df_anno, fusions)
  }
  stopifnot(y$genes$feature_id %in% df_anno$id)
  stopifnot(!duplicated(df_anno$id))

  ### Reorder to match y$genes
  df_anno <- df_anno[match(y$genes$feature_id, df_anno$id), ]
  
  ### Add rownames
  stopifnot(df_anno$id == rownames(y$genes))
  rownames(df_anno) <- df_anno$id
  
  ### Remove temp columns just for matching
  df_anno <- df_anno %>%
    dplyr::select(-one_of(c("id", "by_col")))
  
  ### Return
  stopifnot(rownames(df_anno) == rownames(y$counts))
  dgeList[[dge]]$genes <- df_anno
}
saveRDS(dgeList, "temp/dge_ACC_PFC_HPC_temp4.rds"); beepr::beep()


# Clean up FeatureData -----------------------------------------
### Fix classes
for (dge in names(dgeList)) {
  # Subset featureData
  df <- dgeList[[dge]]$genes
  
  # Convert exon feature level's exon_number from factor to integer
  if (!is.null(df$exon_number)) {
    df$exon_number <- as.integer(df$exon_number)
  }
  
  # Convert factors --> characters (affects all feature levels)
  i <- sapply(df, is.factor) 
  df[i] <- lapply(df[i], as.character)
  
  # Convert logical to YES/NO (affects exons & jxns)
  j <- sapply(df, is.logical)
  df[j] <- lapply(df[j], function (x) {
    as.character(ifelse(x == FALSE, "NO", ifelse(x == TRUE, "YES", NA)))
  })
  
  # Remove columns
  remove_cols <- c(
    # same as gene_symbol (all feature levels)
    "feature_name",
    # Erroneously converted because of wrong class (jxns)
    "startExon", "endExon", "numTx",
    # Duplicates/unnnecessary
    "ensemblGeneID",             # Same as gene_id (jxns)
    "ensemblStrand",             # Same as strand (jxns)
    "inEnsembl",                 # if TRUE, there is gene_id (jxns)
    "e_id"                       # Same as feature_id (exons)
  )
  for (k in names(df)[names(df) %in% remove_cols]) {
    df[k] <- NULL
  }
  
  # Add back rownames
  stopifnot(df$feature_id == rownames(dgeList[[dge]]$counts))
  rownames(df) <- df$feature_id
  
  # Return
  dgeList[[dge]]$genes <- df
}

### Specific to jxn feature level
for (dge in paste0(datasets, "_jxns")) {
  # Subset featureData
  df <- dgeList[[dge]]$genes
  
  # Convert ensemblTx AsIs (list) --> character
  # Required for flattening resList into df table
  df$ensemblTx <- as.character(df$ensemblTx)
  df$ensemblTx[df$ensemblTx == "character(0)"] <- "NA"
  
  # Strand: replace unstranded info of "*" with gene_strand
  df$strand <- as.character(df$gene_strand)
  
  # Feature ID: remove strand suffix of "(*)"
  df$feature_id <- gsub("\\(\\*\\)", "", df$feature_id)
  
  # Align with nomenclature used in other feature levels
  names(df) <- recode(
    names(df),
    "width" = "jxn_intron_length",
    "ensemblTx" = "transcript_id",
    "isFusion" = "jxn_fusion", 
    "leftSeq" = "jxn_leftSeq",
    "rightSeq" = "jxn_rightSeq"
  )
  
  # Clarify
  df$jxn_exonStart_novel <- ifelse(
    df$inEnsemblStart == "YES", "NO", ifelse(
      df$inEnsemblStart == "NO", "YES", NA
    )
  ); df$inEnsemblStart <- NULL
  df$jxn_exonEnd_novel <- ifelse(
    df$inEnsemblEnd == "YES", "NO", ifelse(
      df$inEnsemblEnd == "NO", "YES", NA
    )
  ); df$inEnsemblEnd <- NULL
  
  # Return
  stopifnot(rownames(df) == rownames(dgeList[[dge]]$counts))
  dgeList[[dge]]$genes <- df
}

### Specific to region feature level
for (dge in paste0(datasets, "_regions")) {
  # Subset featureData
  df <- dgeList[[dge]]$genes
  
  # Strand: change from "*" to NA
  df$strand <- NA
  
  # Align with nomenclature used in other feature levels
  names(df) <- recode(
    names(df),
    "anno_class" = "region_class",
    "distToGene" = "region_distToGene"
  )
  
  # Return
  stopifnot(rownames(df) == rownames(dgeList[[dge]]$counts))
  dgeList[[dge]]$genes <- df
}
saveRDS(dgeList, "temp/dge_ACC_PFC_HPC_temp5.rds"); beepr::beep()


# Fix structure ------------------------------------------------
dgeList <- lapply(dgeList, function (y) {
  DGEList(
    counts = y$counts,
    samples = y$samples[ , -grep("^lib.size|^norm.factors", names(y$samples))],
    genes = y$genes
  )
})
saveRDS(dgeList, "temp/dge_ACC_PFC_HPC_temp6.rds"); beepr::beep()


# Normalize ----------------------------------------------------
### TMM normalization
for (i in names(dgeList)) {
  dgeList[[i]] <- calcNormFactors(dgeList[[i]])
}; beepr::beep()


# Save ---------------------------------------------------------
saveRDS(dgeList, "input/dge/ACC_PFC_HPC.rds"); beepr::beep()


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
