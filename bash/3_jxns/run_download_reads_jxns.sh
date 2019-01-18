#!/bin/sh

### RNASeq: Download reads from bam
### Alexis Norris
### Created: 2018-02-28 (from Andrew Jaffe)
### Modified: 2018-03-09 

### JXNS only since that failed due to python synatx error
### bed files made, but not counts files for HPC & PFC!

### Must have already installed:
### regtools (requires cmake)
### subread (for featureCount fxn)

### Setup
echo "run_download_reads.sh starting at:"
date
DIR=/mnt/data/RNASeq/data/Stanley/ArrayCollection_Hippocampus/hisat2_hg19
PKGS=/mnt/data/RNASeq/methods/packages
GTF=/mnt/data/RNASeq/methods/references/hg19/gencode.v19.annotation.gtf
#CHRALIAS=references/chrAliases_GRCh37_to_hg19.csv #chrom names in bam vs gtf
SUBREAD=$PKGS/subread-1.6.0-Linux-x86_64/bin # install featureCounts (http://bioinf.wehi.edu.au/subread-package/)
REGTOOLS=$PKGS/regtools/build
#SCRIPTS=/mnt/data/RNASeq/methods/scripts
THREADS=10 # 1-32
echo "Directory path is: $DIR"
echo "GTF is: $GTF"
echo "Annotation is: $REF"
module add samtools/1.3.1
echo "samtools version 1.3.1"
module add bedtools/2.25.0
echo "bedtools version 2.25.0"
echo "subread version 1.6.0"
module add python2.7/2.7.9
echo "python2.7/2.7.9"


### More about ChrAliases
# Provide a chromosome name alias file to match chr names in
# annotation (gencode gtf chr format = "chr1")
# with those in the 
# reads (bam fasta ref file chr format = hg19, aka "chr1"). 
# This should be a twocolumn
# comma-delimited text file. Its first column should
# include chr names in the annotation and its second column
# should include chr names in the reads. Chr names are case
# sensitive. No column header should be included in the file

### This script is located in scripts directory in main directory, counts is also a 
### Create sub-directories for results
#mkdir $DIR/counts
cd $DIR/counts
#mkdir genes
#mkdir exons
#mkdir jxns


### Get file sample list (remove suffix) --> only files beginning with "A", so that not getting "2_" samples (partial samples)
### Then I manually removed the partners of "2_" samples, files A17.bam, A45.bam, A87.bam, A8.bam, and A91.bam --> full files are "Array....A##_combined.bam"
  #ls ${DIR}/bams/*.sam | 's/\.[^.]*$//' > ${DIR}/samples.txt # starting after hisat
#ls ${DIR}/bams/A*.bam | sed 's/\.[^.]*$//' > ${DIR}/samples.txt # starting with indexed/sorted bams (e.g., if already ran stringtie)

### Run with loop
while read line
do

  ### Sample files
  ID=$(basename $line)  
  BAM=$DIR/bams/${ID}.bam
  echo "Script starting for $line (BAM = $BAM) hisat output and writing out as ${ID}."
  date

  ### BAM must be indexed & sorted! (I already did this in stringtie)
  ### THIS IS NOT CODED CORRECTLY, just a template!
  #samtools view -b -S -o ${ID}.bam ${ID}.sam # bam > sam
	#samtools sort ../bams/${ID}.bam ${ID}.sorted # sorted --> doesn't actually add "sorted" suffix?!
	#samtools index ${ID}.sorted.bam # indexed
	#rm $SAM

	### Calculate library size
	### THIS IS NOT CODED CORRECTLY, just a template!
	#samtools idxstats ${ID}.sorted.bam > ${ID}.libsize.txt

  ### featureCounts: Genes (paired-end, unstranded)
  #echo "   genes..."
  
  #$SUBREAD/featureCounts -p -s 0 -T $THREADS -a $GTF -o genes/${ID}.counts $BAM 
    #-A $CHRALIAS # don't need?
  
  ### featureCounts: Exons
  #echo "      exons..."
  #date
  #$SUBREAD/featureCounts -p -s 0 -T $THREADS -O -f -a $GTF -o exons/${ID}.counts $BAM
    #-A $CHRALIAS # don't need?

  ### Junctions (directly from bam using cigar strings)
  echo "         junctions..."
  date
  OUTJXN=jxns/${ID}_jxns_primaryOnly_regtools.bed
  #TMPBAM=/mnt/data/RNASeq/_scratch/${ID}.bam
  #samtools view -bh -@ $THREADS -F 0x100 $BAM > $TMPBAM
  #samtools index $TMPBAM
  #$REGTOOLS/regtools junctions extract -i 9 -o $OUTJXN $TMPBAM 
  ### need to download regtools (https://github.com/griffithlab/regtools); requires cmake installed too (https://cmake.org/install/)
  OUTCOUNT=jxns/${ID}_jxns_primaryOnly_regtools.counts
  python $PKGS/bed_to_juncs_tophat.py < $OUTJXN > $OUTCOUNT ### Andrew uses: bed_to_junc_withCount (can't find it -- not in his github jaffelab R pkg; not in regtools...; might be synonymous with bed_to_juncs_tophat, it IS from: https://github.com/infphilo/tophat/blob/master/src/bed_to_juncs)
  ### need to get this fxn (is it in regtools?)
  
  ### Convert to bigwig (already did in stringtie)
  #echo "            convert to bw..."
  #date
  #CHRSIZE=/mnt/data/RNASeq/methods/references/hg19/hg19.chrom.sizes
  #BG=$DIR/bams/${ID}.bedGraph
  #BW=$DIR/bams/${ID}.bw
  #kentUtils=/mnt/data/RNASeq/methods/packages/kentUtils-302.1.0/bin
  
  #bedtools genomecov -ibam $BAM -bga -split > $BG
  #awk '$1 ~ /^chr/' $BG > ${BG}.tmp
  #$kentUtils/bedGraphToBigWig ${BG}.tmp $CHRSIZE $BW
  #rm $BG ${BG}.tmp

done < "$DIR/samples.txt"

echo "Script completed at:"
date

