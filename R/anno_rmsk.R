### Stanley ACC SZ vs CTRL
### Get repeatMasker annotation for hg19
### Alexis Norris
### Created: 2018-11-05
### Modified: 2018-12-30


# Setup --------------------------------------------------------
### Packages
library(rtracklayer)       # get repeatMasker from ucsc
library(GenomicRanges)     # convert to GR for anno

### Parameters
analysis_name <- "anno_rmsk"


# Get anno -----------------------------------------------------
### Connect to ucsc
getucsc <- browserSession("UCSC")
genome(getucsc) <- "hg19"

### Get available track names and their table names
tracks <- trackNames(ucscTableQuery(getucsc))        # "rmsk"
tableNames(ucscTableQuery(getucsc, track = "rmsk"))  # "rmsk"

### Get repeatMasker 
rmsk_table <- getTable(ucscTableQuery(
  getucsc, 
  track = "rmsk", 
  table = "rmsk"
))

### Put rmsk in GR object
rmskGR <- makeGRangesFromDataFrame(
  rmsk_table, 
  keep.extra.columns = TRUE, 
  seqnames.field = "genoName",
  start.field = "genoStart", 
  end.field = "genoEnd"
)


# Save ---------------------------------------------------------
saveRDS(rmskGR, "input/anno/rmsk_hg19.rds")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
