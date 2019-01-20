#!/bin/bash

### Alexis Norris
### HISAT2 alignment
### Created: 2017-06-25
### Modified: 2018-01-04


### Log
echo "Script started at $(date), using 48 cores."
echo "TOOLS:"
echo "HISAT2 version 2.0.4"
echo "samtools version 1.3.1"
REF=hg19-hisat
CHRSIZE=hg19.chrom.sizes
echo "Using hg19 (UCSC) reference genome"


### Get sample list
ls *_R1_001.fastq | sed s'/_R1_001.fastq//'g > samples.txt


### Run alignment
echo "HISAT2 alignment using hisat2-2.0.4/hisat2 --dta -x ${REF}"
while read line
do
    FQ=$line
    ID=$(basename $FQ)
    FWD=$(echo "$FQ" | sed 's/$/_R1_001.fastq/')
    REV=$(echo "$FQ" | sed 's/$/_R2_001.fastq/')

    ### Align
	echo "$ID alignment started at $(date)."
    echo "Using $FWD and $REV files."
    echo "Written out as ${ID}."
	hisat2 -p 48 --dta -x $REF -1 $FWD -2 $REV -S ${Sample}.sam

    ### Sort & index
    echo "$ID sort and index started at $(date)."
    samtools sort -@ 48 -o ${Sample}.bam ${Sample}.sam
    samtools index ${Sample}.bam

done < "samples.txt"


### Completed
echo "Script completed at $(date)."
