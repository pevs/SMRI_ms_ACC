# SMRI_Array_ACC_SZ

***

## Description
Analysis of RNA-Seq from the anterior cingulate cortex (ACC) of patients with schizophrenia (SZ) and unaffected controls (CTRL), from the Stanley Medical Research Institute (SMRI) Array Collection, using qSVA methodology to adjust for postmortem RNA degradation and derfinder for flexible expression quantification. Results were compared with RNA-Seq data from prefrontal cortex (PFC) and hippocampus (HPC) of the same individuals.  

***

## Data

### Anterior cingulate cortex (ACC; BA24)   
[Zhao, et al. Transcriptome sequencing and genome-wide association analyses reveal lysosomal function and actin cytoskeleton remodeling in schizophrenia and bipolar disorder. Mol Psych 2015](http://www.nature.com/mp/journal/v20/n5/full/mp201482a.html). Data available from [SMRI](http://sncid.stanleyresearch.org/).

### Hippocampus (HPC)
[Darby, et al. Consistently altered expression of gene sets in postmortem brains of individuals with major psychiatric disorders. Transl Psych 2016](http://www.nature.com/tp/journal/v6/n9/full/tp2016173a.html). Data available from [SMRI](http://sncid.stanleyresearch.org/).

### Dorsolateral pre-frontal cortex (PFC; BA46)    
[PsychENCODE's BrainGVEx](https://www.synapse.org/#!Synapse:syn4590909). Data available from [PsychENCODE Consortium’s BrainGVEX project](https://www.synapse.org/#!Synapse:syn4590909).

### Non-RNA-Seq
Additional data (e.g., qRT-PCR, protein expression, and cell type counts) available for samples at [SMRI Database](sncid.org).

***

# Analysis workflow

## 1. Alignment  
- run_align_hisat2.sh

## 2. Quantify feature expression

### Genes, exons, and junctions (featureCounts)  
- run_counts_genes_exons_jxns.sh  
- run_anno_genes_exons.R   
- run_anno_jxns.R 

### Regions (derfinder)  
- run_counts_regions.sh  
- run_anno_regions.sh  

## 3. Preprocessing 

### Library size and degradation stats  
- import_libsize_degstats.R  

### Combine region datasets  
- import_genes_exons.R  
- import_jxns.R  
- import_regions.R  
- hpc_fix_qsva_regions.R  --> is this in exprsData_prep.R???

### Filtering samples and features
- exprsData_prep.R  
- exprsData_filter.R  

### Cohort descriptive stats (tableone)   
- cohort_stats.R 

## 4. Outlier analysis (PCA)   
- pca.R

## 5. Estimate postmortem degradation (qSVA)

### Get expression (coverage) for degradation-associated regions  
- run_qsva_degstats.sh

### Run qSVA
- qsva.R

### Plot
- plot_qsva.R

## 6. Differential expression analysis

### Run edgeR stattests  
- edgeR.R

### Run edgeR for posthoc tests

#### Antipsychotic medication in ACC dataset  
- posthoc_antipsychotics.R

#### Other brain regions (PFC and HPC)  
- posthoc_brainregions.R

### Plot  
- plot_edgeR.R

### Genelists of differentially expressed genes (DEGs)  
- genelists.R


## 7. Enrichment testing of DEGs

### Cell type enrichment analysis (EWCE)

#### Preprocessing single cell transcriptome datasets   
- celldata_Darmanis.R  
- celldata_Lake.R  

#### Run EWCE  
- ewce_celltypes.R  

#### Plot  
- plot_ewce_celltypes.R

### Geneset enrichment analysis (ClusterProfiler)    
- gsea_kegg.R  
- plot_gsea_kegg.R  


## 8. Overlap of DEGs with posthoc tests and other genesets

### Annotation references  
- anno_geneMap_genomicState.R  
- anno_hgnc.R  
- anno_rmsk.R  
- anno_qsva.R  


### Updating genesets' gene symbols

#### Previous analysis of RNA-Seq data  
- anno_Zhao_ACC.R  
- anno_Darby_HPC.R  

#### SZ-associated genes  
- anno_PGC_108loci_Ripke2014.R  
- anno_PheGenI.R  

### Add overlaps to DEG table  
- overlaps.R  

### Plot overlaps  
- _not finished_ plot_UpSet.R  


## 9. Correlation/validation with SNCID data  
- _not finished_ addl_tests.R  
