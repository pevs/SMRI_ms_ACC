#!/usr/bin/env r

### Preprocessing featureCounts counts
### Junctions
### Alexis Norris
### Created: 2018-02-28 (Andrew Jaffe)
### Modified: 2018-03-09
### Run in bash script with "module add R/3.2.3"


# Setup --------------------------------------------------------
### Packages
library(plyr)
library(magrittr)
library(readr)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(GenomicRanges)
library(sessioninfo)


# Function -----------------------------------------------------
### From LieberInstitute/jaffelab
junctionCount <- function (
  junctionFiles, sampleNames = names(junctionFiles), 
  output = "Count", minOverhang = 0, strandSpecific = FALSE, 
  illuminaStranded = FALSE, minCount = 1, maxCores = 12
) 
  {
  stopifnot(length(junctionFiles) == length(sampleNames))
  names(junctionFiles) <- sampleNames
  message(paste(Sys.time(), "reading in data"))
  theData <- mclapply(junctionFiles, function(x) {
    y <- read.delim(
      x, skip = 1, header = FALSE, sep = "\t",
      col.names = c("chr", "start", "end", "strand", "count"), 
      colClasses = c("character", "integer", "integer", "character", "integer")
    )

    y <- y[y$count >= minCount, ]
    weird <- which(y$strand == "?")
    if (length(weird) > 0) y <- y[-weird, ]
    
    gr <- GRanges(
      y$chr,
      IRanges(y$start, y$end), 
      strand = y$strand, 
      count = y$count
    )
    
	return(gr)
  }, mc.cores = maxCores)
  
message(paste(Sys.time(), "creating master table of junctions"))
  grList <- GRangesList(theData)
  
  if (illuminaStranded & strandSpecific) {
    grList <- GRangesList(mclapply(grList, function(x) {
      strand(x) = ifelse(strand(x) == "+", "-", "+")
      return(x)
    }, mc.cores = maxCores))
  }
  
  fullGR <- unlist(grList)
  
  if (!strandSpecific)strand(fullGR) <- "*"
  
  fullGR <- fullGR[!duplicated(fullGR)]
  fullGR <- sort(fullGR)
  fullGR$count <- NULL
  message(paste(Sys.time(), "there are", length(fullGR), "total junctions"))
  message(paste(Sys.time(), "populating count matrix"))
  jNames <- paste0(
    as.character(seqnames(fullGR)), ":", 
    start(fullGR), "-", end(fullGR), 
    "(", as.character(strand(fullGR)), ")"
  )
  options(warn = -1)
  mList <- mclapply(
    grList, match, fullGR, 
    ignore.strand = !strandSpecific, 
    mc.cores = maxCores
  )
  options(warn = 0)
  countList <- mList
  M <- length(jNames)
  message(paste(Sys.time(), "filling in the count matrix"))
  for (i in seq(along = grList)) {
    if (i%%25 == 0) 
      cat(".")
    cc <- rep(0, M)
    cc[mList[[i]]] <- theData[[i]]$count
    countList[[i]] <- Rle(cc)
  }
  
  countDF <- DataFrame(countList, row.names = jNames, check.names = FALSE)
  names(fullGR) <- jNames
  out <- list(countDF = countDF, anno = fullGR)
  
  return(out)
}


# Load feature annotation for genes, exons ---------------------
### Gene and exons maps (gencode annotation)
### (from featureCounts that I used gencode maps to add gene info to)
exonMap <- readRDS("ExonMap_Hsapiens.gencode.v19.GRCh37_with_eID.rds")
geneMap <- readRDS("GeneMap_Hsapiens.gencode.v19.GRCh37.rds")

### Junctions (ensembl annotation)
jxnMap <- readRDS("ensembl_hg19_v75_junction_annotation.rds")


# Load phenoData -----------------------------------------------
pd <- read_tsv("featureCounts_phenoData.txt")


# Get jxn counts ----------------------------------------------- 
### From primary alignments
junctionFiles <- paste0(pd$FileName, ".jxns.counts")
names(junctionFiles) <- pd$FileName
stopifnot(file.exists(junctionFiles)) 
juncCounts <- junctionCount(junctionFiles, pd$FileName, maxCores = 12)


# Annotate junctions -------------------------------------------
### Get anno (genomic coordinates) of sample junctions
anno <- juncCounts$anno

### Clean up their chromosome names
seqlevels(anno, force = TRUE) <- paste0("chr", c(1:22, "X", "Y", "M"))

## add additional annotation
anno$inEnsembl <- countOverlaps(anno, jxnMap, type = "equal") > 0
anno$inEnsemblStart <- countOverlaps(anno, jxnMap, type = "start") > 0
anno$inEnsemblEnd <- countOverlaps(anno, jxnMap, type = "end") > 0

oo <- findOverlaps(anno, jxnMap, type = "equal")
anno$ensemblGeneID <- NA
anno$ensemblGeneID[queryHits(oo)] <- as.character(jxnMap$ensemblID[subjectHits(oo)])
anno$ensemblSymbol <- NA
anno$ensemblSymbol[queryHits(oo)] <- jxnMap$symbol[subjectHits(oo)]
anno$ensemblStrand <- NA
anno$ensemblStrand[queryHits(oo)] <- as.character(strand(jxnMap)[subjectHits(oo)])
anno$ensemblTx <- CharacterList(vector("list", length(anno)))
anno$ensemblTx[queryHits(oo)] <- jxnMap$tx[subjectHits(oo)]
anno$numTx <- elementLengths(anno$ensemblTx)

### Clean up
anno$ensemblSymbol <- geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

### Junction class (jaffe calls "jxn_class" "code" instead)
anno$jxn_class <- ifelse(
  anno$inEnsembl, "InEns", 
  ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
         ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")
  )
)

### Exons and junctions
exonGR <- GRanges(exonMap$chr,	IRanges(exonMap$start, exonMap$end))
anno$startExon <- match(
  paste0(seqnames(anno), ":", start(anno) - 1), 
	paste0(seqnames(exonGR), ":", end(exonGR))
)
anno$endExon <- match(
  paste0(seqnames(anno), ":", end(anno) + 1),
	paste0(seqnames(exonGR), ":", start(exonGR))
)

### Get genes
g <- data.frame(
	leftGene = exonMap$gene_id[anno$startExon],
	rightGene = exonMap$gene_id[anno$endExon],
	leftGeneSym = exonMap$gene_symbol[anno$startExon],
	rightGeneSym = exonMap$gene_symbol[anno$endExon],
	stringsAsFactors = FALSE
)
g$newGene <- NA
g$newGeneSym <- NA
g$newGene[which(g$leftGene == g$rightGene)] <- g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene == g$rightGene)] <- g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene != g$rightGene)] <- paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene != g$rightGene)] <- paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] <- g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] <- g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] <- g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] <- g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym == ""] <- NA
g$newGeneSym[g$newGeneSym == "-"] <- NA
anno$newGeneID <- g$newGene
anno$newGeneSymbol <- g$newGeneSym
anno$isFusion <- grepl("-", anno$newGeneID)

### Extract
jMap <- anno
jCounts <- juncCounts$countDF
jCounts <- as.data.frame(jCounts[names(jMap), pd$FileName])

### Reduce
mappedPer80M <- pd$totalMapped/80e6 # normalize to 80M libsize
countsM <- DataFrame(mapply(function(x,d) x/d, jCounts, mappedPer80M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm <- as.data.frame(countsM)

### Get sequence of acceptor/donor sites
left = right = jMap
end(left) <- start(left)+1
start(right) <- end(right)-1
# Using library(BSgenome.Hsapiens.UCSC.hg19)
jMap$leftSeq <- getSeq(Hsapiens, left)
jMap$rightSeq <- getSeq(Hsapiens, right)

### save counts
save(
  jMap, jCounts, pd,
  compress = TRUE,
  file = "rawCounts_n20_jxns_hpc.rda"
)
save(
  jMap, jRpkm, pd,
  compress = TRUE,
  file = "rpkmCounts_n20_jxns_hpc.rda"
)

# Session info -------------------------------------------------
library(sessioninfo)
Sys.time()
