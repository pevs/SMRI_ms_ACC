#!/bin/bash

### Alexis Norris
### HISAT2 alignment
### 2017-06-25


### Log
echo "Script started at $(date), using 48 cores."
echo "TOOLS:"
echo "HISAT2 version 2.0.4"
echo "samtools version 1.3.1"
REF=hg19-hisat
echo "Using hg19 (UCSC) reference genome"

### Setup directories
mkdir bams
mkdir logs


### Get sample list
ls *_R1_001.fastq | sed s'/_R1_001.fastq//'g > samples.txt


### Run alignment
echo "HISAT2 alignment using hisat2-2.0.4/hisat2 --dta -x ${REF}"
while read line
do
    ID=$line
    FWD=$(echo "$ID" | sed 's/$/_R1_001.fastq/')
    REV=$(echo "$ID" | sed 's/$/_R2_001.fastq/')
    Sample=$(basename $ID)

	echo "$Sample alignment started at $(date)."
    echo "Using $FWD and $REV files."
    echo "Written out as $(Sample)."
	hisat2 -p 48 --dta -x $REF -1 fastq/$FWD -2 fastq/$REV -S ${Sample}.sam

    echo "$Sample sort and index started at $(date)."
    samtools sort -@ 48 -o ${Sample}.bam ${Sample}.sam
    samtools index ${Sample}.bam

done < "samples.txt"


### Completed
echo "Script completed at $(date)."
