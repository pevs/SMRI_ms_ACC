#!/bin/bash
### Run R script, using updated R version
### Alexis Norris
### 01 Aug 2016
### 03 Mar 2018

module add R/3.2.3


#srun R CMD BATCH ./run_regionMatrix.r # triggers "cannot create R_TempDIR" error

### Derfinder
#env TMPDIR=/mnt/data/RNASeq/_scratch R CMD BATCH ./run_regionMatrix.r 

### Andrew Jaffe's featureCounts (stringtie alternative)
env TMPDIR=/mnt/data/RNASeq/_scratch R CMD BATCH run_preproc_counts_genes_exons.R


