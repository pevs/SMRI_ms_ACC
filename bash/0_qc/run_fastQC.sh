#!/bin/bash
### Alexis Norris
### 19 June 2017
### FastQC
### Check if re-transferred files are corrupt
### If so, need to get from SNCID.org or Sarven
### If files are OK, then unzip
### Set to use 12 cores

echo " FastQC on re-transferred ArrayCollection Hippocampus RNA-Seq files (from Sarven) and a subset of OFC RNA-Seq files from NPC that also were corrupt, starting at:"
date

echo "Directory path is:"
pwd

#find . -name "*.fastq.gz" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 


### Additional samples [2017-07-06]
find . -name "S201_NoIndex_L002_R1_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 
find . -name "S201_NoIndex_L002_R2_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 
find . -name "S_359_NoIndex_L008_R1_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 
find . -name "S_359_NoIndex_L008_R2_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 
find . -name "S192_NoIndex_L008_R1_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 
find . -name "S192_NoIndex_L008_R2_001*" | xargs -n 1 /mnt/data/RNASeq/methods/FastQC/fastqc 


echo "FastQC complete at:"
date


