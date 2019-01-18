#!/bin/bash

### Alexis Norris
### 25 June 2017
### HISAT2-Stringtie-Ballgown Workflow
### PART ONE
### Alignment and expression estimation
### Set to use 48 cores

echo "HISAT2-StringTie-Ballgown Workflow Step 2 [HISAT2 alignment and StringTie estimation] starting at:"
date

echo "Directory path is:"
pwd

echo "Samtools version: 1.3.1"
module add samtools/1.3.1

### Setup directories
mkdir hisat2_stringtie_hg19
mkdir hisat2_stringtie_hg19/bams
mkdir hisat2_stringtie_hg19/ballgown
mkdir hisat2_stringtie_hg19/stringtie_unmerged
mkdir hisat2_stringtie_hg19/logs

### Combine FastQC results
#echo "Tabulating FastQC results"
#cd FastQC
#module add python3/python3.3.6
#/mnt/data/RNASeq/methods/FastQC/tabulate_results.py
### Add header with information
#sed -i '1iFileName,Non_informative,Per_base_sequence_quality,Per_tile_sequence_quality,Per_sequence_quality_scores,Per_base_sequence_content,Per_sequence_GC_content,Per_base_N_content,Sequence_Length_Distribution,Sequence_Duplication_Levels,Overrepresented_Sequences,Adapter_Content,Kmer_Content' FastQC_results_table.csv
#echo "Finished tabulating FastQC results at:"
#date

### Get file sample list
cd fastq
ls *_R1_001.fastq | sed s'/_R1_001.fastq//'g > ../samples.txt
cd ..

### Alignment
echo "Starting alignment of paired-end reads."
while read line
do
        FWD=$(echo "$line" | sed 's/$/_R1_001.fastq/')
        REV=$(echo "$line" | sed 's/$/_R2_001.fastq/')
        Sample=$(basename $line)
	Ref=/mnt/data/RNASeq/methods/references/hg19/hg19-hisat
	#Ref=/mnt/data/RNASeq/methods/references/mm10/mm10-hisat
	echo "Alignment starting for $line, using files $FWD and $REV, and written out as $Sample."
	/mnt/data/RNASeq/methods/packages/hisat2-2.0.4/hisat2 -p 48 --dta -x $Ref -1 fastq/$FWD -2 fastq/$REV -S hisat2_stringtie_hg19/bams/${Sample}.sam
	echo "Alignment complete for $Sample at:"
	date
done < "samples.txt"

echo "Alignment to $Ref using hisat2-2.0.4 [options: --dta -x] is finished at:"
date


### StringTie assembly
echo "Starting StringTie assembly. StringTie version stringtie-1.2.4.Linux_x86_64, using the reference annotation:" 
	### Gencode for hg19
RefAnno=/mnt/data/RNASeq/methods/references/hg19/gencode.v19.annotation.gtf
	### Gencode for mm10
#RefAnno=/mnt/data/RNASeq/methods/references/mm10/gencode.vM10.annotation.gtf
echo ${RefAnno}

while read line
do
        Sample=$(basename $line)
        echo "Stringtie assembly starting for $line, written out as $Sample."
        
	samtools sort -@ 48 -o bams/${Sample}.bam bams/${Sample}.sam
        samtools index bams/${Sample}.bam bams/${Sample}.bai
	echo "View and sort complete for $Sample"
        
	/mnt/data/RNASeq/methods/packages/stringtie-1.2.4.Linux_x86_64/stringtie -p 48 -G ${RefAnno} -o stringtie_unmerged/${Sample}.gtf bams/${Sample}.bam
	echo "Stringtie assembly complete for $Sample at:"
	date
done < "samples.txt"

echo "StringTie estimation is complete. Script ended at:"
date



