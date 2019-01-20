#!/bin/sh

### Quantify expression of genes, exons, and junctions
### Alexis Norris
### Created: 2018-02-28
### Modified: 2018-03-01


### Log
echo "Script started at $(date)."
echo "TOOLS:"
echo "featureCounts from subread linux version 1.6.0"
echo "samtools version 1.3.1"
echo "bedtools version 2.25.0"
echo "subread version 1.6.0"
GTF=gencode.v19.annotation.gtf
echo "Using $GTF annotation"


### Get sample list
ls *.bam | sed 's/\.[^.]*$//' > samples.txt

### Run with loop
while read line
do

  ### Sample files
  ID=$(basename $line)  
  BAM=${ID}.bam
  echo "Script started for $ID at $(date)."

  ### featureCounts: Genes
  echo "Started genes at $(date), using featureCounts -p -s 1."
  date
  featureCounts -p -s 1 -T 10 -a $GTF -o ${ID}.gene.counts $BAM
  
  ### featureCounts: Exons
  echo "Started exons at $(date)."
  featureCounts -p -s 1 -T 10 -O -f -a $GTF -o ${ID}.exon.counts $BAM

  ### Junctions from primary alignments (from cigar strings in bam)
  echo "Started junctions at $(date)."
  TEMP=${ID}.temp.bam
  BED=${ID}.jxns.bed
  COUNTS=${ID}.jxns.counts
  samtools view -bh -@ 10 -F 0x100 $BAM > $TEMP
  samtools index $TEMP
  regtools junctions extract -i 9 -o $BED $TEMP
  bed_to_juncs_tophat.py < $BED > $COUNTS

done < "samples.txt"


### Completed
echo "Script completed at $(date)."
