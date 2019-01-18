#!/bin/bash

### Run R script for derfinder
### Alexis Norris
### 1 August 2016
### 10 March 2018

module add R/3.2.3

#srun R CMD BATCH ./run_regionMatrix.r
env TMPDIR=/mnt/data/RNASeq/_scratch R CMD BATCH --max-mem-size=190G ./run_regionAnno_pfc.r #./run_regionMatrix_pfc.r 

