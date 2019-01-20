#!/bin/bash

### Alexis Norris
### Estimate postmortem degradation using qSVA
### Created: 2017-03-30
### Modified: 2017-07-06


### Log
echo "Script started at $(date), using 48 cores."
echo "TOOLS:"
echo "Bedtools version 2.25.0"
echo "Samtools version 1.1"
echo "Bedtools version 2.25.0"
echo "kentUtils-302.1.0"
echo "bwtools downloaded Jul 12 2016"
CHRSIZE=hg19.chrom.sizes
echo "Using hg19 (UCSC) reference genome"


### Create list of samples
ls *.bam > samples.txt


### Degradation-associated regions bed file (from qSVA publication)
    # if polyA library (n = 1000 regions)
BED=polyA_degradation_regions_biorxivS1.bed
    # if RiboZero library (n = 515 regions)
#BED=ribozero_degradation_regions_biorxivS2.bed


### Get coverage for degradation regions via bigwig
while read line
do

    BAM=${line}
    ID=$(basename $line)
    OUT=${ID}.degStats.txt
    BW=${ID}.bw
    BG=${ID}.bedGraph
    echo "Starting for $line file, for $ID Sample."

    ### Convert bam to bigwig
    bedtools genomecov -ibam $BAM -bga -split > $BG
    bedGraphToBigWig $BG $CHRSIZE $BW

    ### Get degradation coverage for polyA or riboZero regions
    bwtool summary $BED $BW $OUT -header -fill=0 -with-sum

    ### Calculate library size (total number of mapped reads)
    samtools idxstats $BAM > ${ID}.libsize.txt

done < "samples.txt"


### Completed
echo "Script completed at $(date)."
