#!/usr/bin/env bash

### Alexis Norris
### Created: 10 July 2017
### Modified: 19 Jan 2018
### Estimate rRNA contamination in RNASeq data (Illumina, paired-end)
### Using refseq source for hg rRNA
### Currently using 8 threads (p=8)

echo "run_rRNA_refSeq.sh started at:"; date

### Sample list
ls fastq/*_R1_001.fastq | sed s'/_R1_001.fastq//'g > samples_rRNA.txt

### Make directory
mkdir rRNA_refSeq
#mkdir rRNA_hg38


### Reference = RefSeq rRNA sequences for Homo sapiens
REF=/mnt/data/RNASeq/methods/rRNA/refSeq.rRNA-hisat 
### Using RefSeq: sent all 44 sequences to file, then uploaded as Homo_sapiens_rRNA_refSeq.fa
#https://www.ncbi.nlm.nih.gov/nuccore/?term=rRNA+%22Homo+sapiens%22%5Bporgn%3A__txid9606%5D # Accessed 2018-Jan-19


### Reference = hg38 rRNA from ensembl
#REF=/mnt/data/RNASeq/methods/rRNA/hg38.rRNA-hisat 
### Using ENSEMBL and an awk trick I found on online forum
### Download ncRNA fasta from ensembl
#wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
	### Just rRNA
#cat Homo_sapiens.GRCh38.ncrna.fa | awk '/^>/ { p = ($0 ~ /gene_biotype:rRNA/)} p' > Homo_sapiens.GRCh38.ncrna_rRNA.fa


### Method
echo "Starting alignment to hg rRNA from RefSeq with hisat2-2.0.4 [options: --dta]; using 8 cores"

### Run alignment in loop
while read line
do

   FWD=$(echo "$line" | sed 's/$/_R1_001.fastq/')
   REV=$(echo "$line" | sed 's/$/_R2_001.fastq/')
   ID=$(basename $line)
   echo "file $line, writing out as $ID."
   
   /mnt/data/RNASeq/methods/packages/hisat2-2.0.4/hisat2 -p 8 --dta -x $REF -1 $FWD -2 $REV -S rRNA_refSeq/${ID}.sam

done < "samples_rRNA.txt"

echo "run_rRNA_refSeq.sh finished at:"; date

