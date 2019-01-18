### SMRI Array ACC SZ vs CTRL
### Cohort: descriptive stats
### Alexis Norris
### Created: 2018-12-29
### Modified: 2019-01-10


# Notes --------------------------------------------------------
### Sources: 
  # Tutorial: https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html
  # Manual: https://cran.r-project.org/web/packages/tableone/tableone.pdf
  # Github: https://github.com/kaz-yos/tableone


# Input files --------------------------------------------------
### DGEs, for phenoData
#dgeList <- readRDS("input/dge/ACC_PFC_HPC.rds")


# Setup --------------------------------------------------------
### Packages
library(tableone)  # descriptive stats
library(tidyverse) # plotting; wrangling

### Parameters
analysis_name <- "cohort_stats"


# Load DGE -----------------------------------------------------
dgeList <- readRDS("input/dge/ACC_PFC_HPC.rds")


# Stattests used -----------------------------------------------

### Categorical variables
### If <5 counts in a group/category --> use exact (fisher.test); otherwise use approx (chisq.test)
# APPROX: A function used to perform the large sample approximation based tests. 
  # The default is chisq.test. 
  # This is not recommended when some of the cell have small counts like fewer than 5.
#testApprox = chisq.test

# EXACT: A function used to perform the exact tests. 
  # The default is fisher.test. 
  # If the cells have large numbers, it will fail because of memory limitation. 
  # In this situation, the large sample approximation based should suffice.
#testExact = fisher.test      



### Continuous variable
### If highly right/left skewed (>2 in a group) --> use nonnormal, otherwise use normal
# NORMALANOVA:
  # A function used to perform the normal assumption based tests. 
  # The default is oneway.test
  # This is equivalent of the t-test when there are only two groups.
#testNormal = oneway.test        

# NON-NORMAL(RANKED) ANOVA: 
  # A function used to perform the nonparametric tests. 
  # The default is kruskal.test (Kruskal-Wallis Rank Sum Test). 
  # This is equivalent of the wilcox.test (Man-Whitney U test) when there are only two groups.
#testNonNormal = kruskal.test  


# PhenoData ----------------------------------------------------
### Extract phenoData
pheno_full <- dgeList$ACC_genes$samples

### Subset for covariates of interest
### Removing sample ID
keep <- c(
  "Age", "Sex", 
  "BrainHemisphere", 
  "RIN", "Brain_pH", "PMI", 
  #"RefrigInterval",
  "Smoking", "Alcohol", "Drugs"
)
pheno <- pheno_full %>% 
  select(Dx, keep)


# Setup variables ---------------------------------------------
### Variables
vars <- names(pheno)[-grep("Dx", names(pheno))]

### Variable classes
### Get all non-numeric cols (these are factor/categorical comparisons)
vars_factors <- names(pheno)[!sapply(pheno, is.numeric)]


# Run tableOne -------------------------------------------------
## Create Table 1 stratified by Dx
tableOne <- CreateTableOne(
  vars = vars, 
  strata = c("Dx"), 
  data = pheno,
  factorVars = vars_factors,
  
  # If TRUE, NA is handled as a regular factor level rather than missing. 
  # NA is shown as the last factor level in the table. 
  # Only effective for categorical variables.
  includeNA = FALSE,              # default is FALSE
  
  ### Categorical variables
  testApprox = chisq.test,        # default
  testExact = fisher.test,        # use if groups < 5
  
  # Continuous variables
  testNormal = oneway.test,       # default
  testNonNormal = kruskal.test    # non-normal (ranked); use if abs(skew) > 2
)


# Decide stattests to use --------------------------------------
### Decide whether to use exact or non-normal tests for each variable
### https://rpubs.com/kaz_yos/tableone-demo-e
  # Continuous variable: if highly right/left skewed (>2 in a group) --> use nonnormal
summary(tableOne$ContTable)

  # Categorical: if small #s in some categories --> use exact
summary(tableOne$CatTable)


### Export 
### Source: https://github.com/kaz-yos/tableone/blob/master/R/print.TableOne.R (raw script)
tableOne_mat <- print(
#tableone:::print.TableOne(
  tableOne,
  
  # Number of digits to show
  catDigits = 1, 
  contDigits = 2, 
  pDigits = 3, 
  
  # Whether to show quotes; use quotes if copying/pasting into Excel
  quote = FALSE,       
  
  ## Common options
  missing = TRUE,      # Not implemented yet --> actually seems to work!
  explain = TRUE,      # Whether to show explanation in variable names
  printToggle = TRUE,  # Whether to print the result visibly
  test = TRUE,         # Whether to add p-values
  smd = TRUE,          # Whether to add standardized mean differences
  noSpaces = TRUE,     # Whether to remove spaces for alignments
  padColnames = FALSE, # Whether to pad column names for alignments
  varLabels = FALSE,   # Whether to show variable labels instead of names
  
  ### Categorical options
    # Use exact (fisher) if category sizes small (< 5)
    # if variable not explicitly listed here for exact, approx (chisq) will be used
    # HERE, only brain hemisphere has min 5
  exact = c("Sex", "Smoking", "Drugs", "Alcohol"),
  
  format = c("fp","f","p","pf")[1], # Format f_requency and/or p_ercent
  showAllLevels = FALSE, # Show all levels of a categorical variable
  
    # Which 2-level variables to show both levels in one row
  cramVars = c("Sex", "BrainHemisphere"),
  
    # Do not show " = second level" for two-level variables
  dropEqual = FALSE, 
  
  
  ### Continuous options
    # If highly right/left skewed, if abs(skew) > 2 (for either Dx group) --> use nonnormal ("ranked")
    # HERE, all abs(skew) < 2. Closest is refrigerator interval is 1.96 for ctrl but only 1.34 for sz
  nonnormal = NULL, 
  
    # Whether to show median
  minMax = FALSE  
) 

### Convert to dataframe
tableOne_df <- as.data.frame(tableOne_mat)
tableOne_df$covariate <- rownames(tableOne_df)
tableOne_df$covariate <- gsub("\\.\\.mean\\.\\.sd\\.\\.", "\nmean (SD)", tableOne_df$covariate)
tableOne_df$covariate <- gsub("\\.\\.\\.\\.", "", tableOne_df$covariate)
tableOne_df$covariate <- gsub("X\\.\\.\\.", "", tableOne_df$covariate)
  # Remove the number suffix added when categorical variables have overlapping groups (here alcohol and drug use scores)
tableOne_df$covariate <- gsub("\\.[0-9]$", "", tableOne_df$covariate)
tableOne_df$covariate <- gsub("\\.\\.\\.", "\n", tableOne_df$covariate)
tableOne_df$covariate <- gsub("\\.", "\\/", tableOne_df$covariate)
tableOne_df$covariate <- gsub("_", " ", tableOne_df$covariate)

### Specific to these covariates
tableOne_df$covariate <- gsub("BrainHemisphere", "Brain Hemisphere", tableOne_df$covariate)

### Clean up and write out  
tableOne_df %>%
  as_data_frame() %>%
  select(covariate, everything()) %>%
  # Add comment at bottom
  bind_rows(
    data_frame(
      ### If test empty
      # Continous = one.way
      # Categorical = absolute (with chisq test)
      "covariate" = "Note",
      "test" = "exact (fisher.test) used if categorical and <5 counts in one category for either Dx; nonnormal (Wilcox) non-parametric ranked test used if continuous and absolute skew > 2.00 for either Dx group (note: this is equivalent to wilcox.test (Man-Whitney U test) when there are only 2 groups); if empty, default tests were used: for continuous, a normal (one-way) parametric test, and for categorical, N absolute (chisq)"
    )
  ) %>%
  write_csv("output/phenoData/cohort_summary_ACC_byDx.csv")






# Do summary of SNCID data -------------------------------------s

### Subset for SNCID data of interest
keep <- c(
  # To test for correlation Dx ~ Sx, to relate to EWCE results
  "Celltype_",  "Beasley", "Cotter",
  
  # Expression of genes of interest, with respective housekeeping genes
  "PVALB", "GAD", "Calbindin", "vWF", "BDNF", "GFAP", "RELN",
  "GAPDH_Expression_112003Yolken_FrontalCortex",
  
  # Bahn (protein?)
  "Endothelin1", 
  
  # Other (Mischa)
  "JUNB" 
)
keep_cols <- names(pheno_full)[grep(paste(keep, collapse = "|"), names(pheno_full))]
  # Remove "SD" values
keep_cols <- keep_cols[-grep("cellspermm2SD", keep_cols)]
  # Get phenoData
pheno <- pheno_full %>%     
  select(Dx, keep_cols)

### Normalize
pheno <- pheno %>%     
  mutate(
    # Beasley cell counts for both area and section thickness
    "Astro_Beasley_BA9" = Celltype_Astro_Count_052010Beasley_BA9/(AreaAnalyzed_mm2_052010Beasley_BA9*SectionThickness_um_052010Beasley_BA9),
    "Oligo_Beasley_BA9" = Celltype_Oligo_Count_052010Beasley_BA9/(AreaAnalyzed_mm2_052010Beasley_BA9*SectionThickness_um_052010Beasley_BA9),
    "Neuron_Beasley_BA9" = Celltype_Neuron_Count_052010Beasley_BA9/(AreaAnalyzed_mm2_052010Beasley_BA9*SectionThickness_um_052010Beasley_BA9),
    
    # Yolken expression data
    "RELN_112003Yolken_FC" = RELN_Expression_112003Yolken_FrontalCortex/GAPDH_Expression_112003Yolken_FrontalCortex,
    "GAD1akaGAD67_112003Yolken_FC" = GAD1akaGAD67_Expression_112003Yolken_FrontalCortex/GAPDH_Expression_112003Yolken_FrontalCortex,
    "PVALB_112003Yolken_FC" = PVALB_Expression_112003Yolken_FrontalCortex/GAPDH_Expression_112003Yolken_FrontalCortex,
    "GFAP_112003Yolken_FC"= GFAP_Expression_112003Yolken_FrontalCortex/GAPDH_Expression_112003Yolken_FrontalCortex,
    "BDNF_112003Yolken_FC" = BDNF_Expression_112003Yolken_FrontalCortex/GAPDH_Expression_112003Yolken_FrontalCortex
  ) %>%
  select(-one_of(c(
    "AreaAnalyzed_mm2_052010Beasley_BA9", "SectionThickness_um_052010Beasley_BA9",
    "Celltype_Astro_Count_052010Beasley_BA9", 
    "Celltype_Oligo_Count_052010Beasley_BA9",
    "Celltype_Neuron_Count_052010Beasley_BA9",
    
    
    "Celltype_Astro_Count_052010Beasley_BA9",
    "AreaAnalyzed_mm2_052010Beasley_BA9",
    "SectionThickness_um_052010Beasley_BA9",
    "Celltype_Oligo_Count_052010Beasley_BA9",
    "Celltype_Neuron_Count_052010Beasley_BA9"
  )))

### Setup variables
### Variables
vars <- names(pheno)[-grep("Dx", names(pheno))]

### Variable classes
  # Fix factors

  # Get all non-numeric cols (these are factor/categorical comparisons)
vars_factors <- names(pheno)[!sapply(pheno, is.numeric)]


### Run tableOne
## Create Table 1 stratified by Dx
tableOne <- CreateTableOne(
  vars = vars, 
  strata = c("Dx"), 
  data = pheno,
  factorVars = vars_factors,
  
  # If TRUE, NA is handled as a regular factor level rather than missing. 
  # NA is shown as the last factor level in the table. 
  # Only effective for categorical variables.
  includeNA = FALSE,              # default is FALSE
  
  ### Categorical variables
  testApprox = chisq.test,        # default
  testExact = fisher.test,        # use if groups < 5
  
  # Continuous variables
  testNormal = oneway.test,       # default
  testNonNormal = kruskal.test    # non-normal (ranked); use if abs(skew) > 2
)


### Decide stattests to use 
  # Continuous variable: if highly right/left skewed (>2 in a group) --> use nonnormal
summary(tableOne$ContTable)
  # Categorical: if small #s in some categories --> use exact
summary(tableOne$CatTable)


### Export
tableOne_mat <- print(
  #tableone:::print.TableOne(
  tableOne,
  
  # Number of digits to show
  catDigits = 1, 
  contDigits = 2, 
  pDigits = 3, 
  
  # Whether to show quotes; use quotes if copying/pasting into Excel
  quote = FALSE,       
  
  ## Common options
  missing = TRUE,      # Not implemented yet --> actually seems to work!
  explain = TRUE,      # Whether to show explanation in variable names
  printToggle = TRUE,  # Whether to print the result visibly
  test = TRUE,         # Whether to add p-values
  smd = TRUE,          # Whether to add standardized mean differences
  noSpaces = TRUE,     # Whether to remove spaces for alignments
  padColnames = FALSE, # Whether to pad column names for alignments
  varLabels = FALSE,   # Whether to show variable labels instead of names
  
  
  ### Continuous options
  # If highly right/left skewed, if abs(skew) > 2 (for either Dx group) --> use nonnormal ("ranked")
  nonnormal = c(
    "Celltype_Astro_GFAP_meanOD_012007Cotter_FrontalDeepWhiteMatter",
    "vWF_082012Bahn_Serum",
    "GFAP_112003Yolken_FC",
    "GAD1akaGAD67_112003Yolken_FC", "PVALB_112003Yolken_FC"
  ), 
  
  # Whether to show median
  minMax = FALSE  
) 

### Save
write.csv(tableOne_mat, "output/phenoData/cohort_summary_sncidData.csv")


# Do summary of SZ-only covariates -----------------------------
### Subset
sz_cols <- c("Suicide", "Age_Onset", "LifetimeAntipsychotics", "Lithium")
phenoSZ <- pheno_full %>%
  filter(Dx == "SZ") %>%
  select(sz_cols) 

### Summarize continuous
sz_cont <-names(phenoSZ)[sapply(phenoSZ, is.numeric)]
df_cont <- phenoSZ %>%
  select(sz_cont) %>%
  gather(covariate, value) %>%
  group_by(covariate) %>%
  summarise_all(funs(mean, sd, min, max, median))

### Summarize categorical
sz_cat <-names(phenoSZ)[!sapply(phenoSZ, is.numeric)]
phenoSZ_cat <- phenoSZ %>%
  select(sz_cat) 
  # Summarize
df_cat <- data_frame(
  "covariate" = c("Suicide", "Lithium"),
  "n_No" = c(
    nrow(filter(phenoSZ_cat, Suicide == "NO")),
    nrow(filter(phenoSZ_cat, Lithium == "NO"))
  ),
  "n_Yes" = c(
    nrow(filter(phenoSZ_cat, Suicide == "YES")),
    nrow(filter(phenoSZ_cat, Lithium == "YES"))
  ),
  "n_NA" = c(
    nrow(filter(phenoSZ_cat, is.na(Suicide))),
    nrow(filter(phenoSZ_cat, is.na(Lithium)))
  )
)

### Combine and write out
bind_rows(df_cont, df_cat) %>%
  write_csv("output/phenoData/cohort_summary_ACC_SZspecific.csv")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
