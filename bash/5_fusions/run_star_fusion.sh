#!/usr/bin/env bash

### Alexis Norris
### Created: 21 Dec 2017
### Modified: 09 Mar 2018
### Run star fusion with while read line (not simple-job-array)

#sbatch -n 20 --mem=100000 --nice=1000 -J sfuse3 run_star_fusion_whilereadline.sh

DIR=/mnt/data/RNASeq/data/Stanley/ArrayCollection_Hippocampus
cd $DIR/fastq
ls *_R1_001.fastq | sed 's/_R1_001.fastq//g' > $DIR/star_fusion/samples_star.txt
cd $DIR/star_fusion

### For while read line (not array) job
while read line
do 

REF=/mnt/data/RNASeq/methods/references/starfusionindex

echo "Starting star fusion for $line at:"
date

ID=$line
mkdir ${ID}_star_fusion
cd $DIR/star_fusion/${ID}_star_fusion

echo "Directory path is:"
pwd

STAR-Fusion --genome_lib_dir ${REF} --left_fq $DIR/fastq/${ID}_R1_001.fastq --right_fq $DIR/fastq/${ID}_R2_001.fastq --runThreadN 20 --limitBAMsortRAM 91532137230 --output_dir ${ID}_star_fusion

echo "Star fusion for $line is completed at:"
date

### WARNING: Need to rename file in the Sample directory before combining, since the output files do not have sample name!!!
cp star-fusion.fusion_candidates.final ../${ID}_fusion_candidates_final.txt
cd $DIR/star_fusion

done < "samples_star.txt"

