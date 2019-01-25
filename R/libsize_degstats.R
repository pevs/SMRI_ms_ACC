### SMRI Array ACC SZ vs CTRL
### qSVA
### Alexis Norris
### Created: 
### Modified: 


# Input files --------------------------------------------------


# Setup --------------------------------------------------------
### Packages
library(data.table)
library(readr)
library(dplyr)
library(plyr)
library(ggplot)
library(magrittr)
library(reshape2)

### Parameters
analysis_name <- "import_libsize_degstats"
figcolors <- c("CTRL" = "black", "SCZ" = "blue2", "BPD" = "green")


# Library Sizes ------------------------------------------------
### From hisat2 hg19 bam files; calculated during degradeStats script.   
### Count read pairs separately, so the value is 2x the "Total Readpairs" 

### Load library sizes
lib_files <- paste0("input/libSize/", pd$FileName, ".libsize.txt")

### Match naming to files
names(lib_files) <- pd$FileName
stopifnot(file.exists(lib_files))

### Load library sizes
lib_size <- data.frame(
  "FileName" = names(lib_files),
  "Reads_chrM" = sapply(lib_files, function (x) {
    read.delim(
      pipe(paste("cut -f3", x)), 
      header = FALSE, as.is = TRUE
    )[25, 1]
  }),
  "TotalLibSize" = sapply(lib_files, function (x) {
    colSums(read.delim(
      pipe(paste("cut -f3", x)), 
      header = FALSE, as.is = TRUE
    ))
  })
)

### Add to phenoData
stopifnot(pd$FileName %in% lib_size$FileName)
pd <- plyr::join(pd, lib_size, by = "FileName", "left", "all") 


# Plot library size --------------------------------------------
### Prep
df <- pd %>%
  filter(group %in% c("CTRL", "SCZD", "BPAD"))
df$group <- factor(
  df$group, 
  levels = c("CTRL", "SCZD", "BPAD")
)

### Density plot
g <- df %>%
  ggplot(aes(x = TotalLibSize/2E06, color = group, linetype = cohort)) +
  geom_density(lwd = 1.2, alpha = 1/20, na.rm = TRUE) +
  scale_color_manual(values = figcolors) + 
  scale_fill_manual(values = figcolors) + 
  theme_bw() + theme(
    legend.justification = c(1,1), 
    legend.position = c(1,1), 
    legend.title = element_blank(),
    legend.background = element_blank(), 
    legend.spacing = unit(1, "line")
  ) +
  labs(x = "Total Read Pairs (Millions)", title = NULL)  +
  facet_wrap(~subgroup, scales = "free")
print(g)
ggsave("graphics/LibSize_Density.png", height = 4, width = 3, dpi = 300)

### Boxplot
g <- df %>%
  ggplot(aes(
    x = group, 
    y = TotalLibSize/2E06, 
    color = group
  )) +
  geom_boxplot(
    lpha = 1/20, outlier.color = NA, 
    show.legend = FALSE, na.rm = TRUE
  ) +
  geom_point(alpha = 1/5, show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = figcolors) + 
  scale_fill_manual(values = figcolors) + 
  theme_bw() + theme(axis.text.x = element_text(size = 7)) +
  labs(
    title = NULL, 
    x = NULL, 
    y = "Total Read Pairs (Millions)"
  )  +
  facet_grid(~subgroup, scales = "free_x")
print(g)
ggsave("graphics/LibSize_Boxplot.png", height = 3.5, width = 4, dpi = 300)


# Degradation coverage -----------------------------------------
### Import files, adjusting for read-length
### Create list to hold degcov
degCovList <- vector(mode = "list", length = length(phenoList))
names(degCovList) <- names(phenoList)

### Import files, normalizing for readLength
for(i in names(degFiles)) {
  adj <- studyKeyList[[i]]$readLength
  degCovList[[i]] <- sapply(degFiles[[i]], function (x) {
    read.delim(pipe(paste("cut -f10", x)), as.is = TRUE)$sum/adj
  })
}

### Add genomic coordinates as rownames
for (i in names(degCovList)) {
  if (studyKeyList[[i]]$regions == "polyA") x <- degRegionsList$polyA$Region
  if (studyKeyList[[i]]$regions == "ribozero") x <- degRegionsList$ribozero$Region
  rownames(degCovList[[i]]) <- x
}

### Save
saveRDS(degCovList, "output/qsva/deg_cov.rds")

### Normalize by libSize 
### Check and match order of samples to that in pd (if not already matching)
for (i in names(degCovList)) {
  stopifnotprint(colnames(degCovList[[i]]) %in% phenoList[[i]]$FileName)
  phenoList[[i]] <- phenoList[[i]][which(phenoList[[i]]$FileName %in% colnames(degCovList[[i]])), ]
  degCovList[[i]] <- degCovList[[i]][ , match(phenoList[[i]]$FileName, colnames(degCovList[[i]]))]
}

### Divide by median Libsize -- to reduce scale of values
degCovAdjList <- degCovList
for (i in names(degCovList)) {
  adj <- studyKeyList[[i]]$medianLibSize
  stopifnot(degCovList[[i]]) %in% phenoList[[i]]$FileName0
  degCovAdjList[[i]] <- degCovList[[i]]/matrix(
    rep(phenoList[[i]][match(
      colnames(degCovList[[i]]), 
      phenoList[[i]]$FileName), ]$TotalLibSize/adj), 
    nc = ncol(degCovList[[i]]), 
    nr = nrow(degCovList[[i]]), 
    byrow = TRUE
  ) 
}

### Check order
for (i in names(degCovAdjList)) {
  stopifnot(colnames(degCovAdjList[[i]]) == phenoList[[i]]$FileName)
}


### Save
saveRDS(degCovAdjList, "output/qsva/degradationMatrix.rds")


# Add degstats to pheno ----------------------------------------
### First need to transpose and add shared col with pd (FileName)
degMatList <- lapply(degCovAdjList, function(degCovAdj) {
  df <- as.data.frame(t(degCovAdj))
  df$FileName <- rownames(df)
  df
})

### Join with merge
degMatList_pd <- vector(mode = "list", length = length(names(phenoList)))
names(degMatList_pd) <- names(phenoList)
for(i in names(phenoList)) {
  stopifnot(phenoList[[i]]$FileName %in% rownames(degMatList[[i]]))
  degMatList_pd[[i]] <- plyr::join(
    phenoList[[i]], 
    degMatList[[i]], 
    by = "FileName"
  )
  rownames(degMatList_pd[[i]]) <- degMatList_pd[[i]]$sample_id
}

### Melt
degMatList_molten <- lapply(degMatList_pd, function (df) {
  molten <- melt.data.table(
    data.table::data.table(df), 
    id.vars = names(df)[grep("^chr", names(df), invert = TRUE)], 
    measure.vars = names(df)[grep("^chr", names(df))],
    variable.name = "Region",
    value.name = "adjCov"
  )
  
  ### Remove 0 expr values
  data.frame(molten[which(molten$adjCov > 0), ])
})

### Unlist
### Note: these are by region, so duplicate FileNames exist
degMat_df <- do.call("rbind", degMatList_molten)
degMat_df$cohort_model <- gsub("\\..*", "", rownames(degMat_df))
summary(factor(degMat_df$cohort_model))
rownames(degMat_df) <- NULL

### Add deg_ to degMat colnames inphenoList 
degMatList_pd <- lapply(degMatList_pd, function (df) {
  names(df)[grep("^chr", names(df))] <- paste0("deg", names(df)[grep("^chr", names(df))]) 
  df
})
saveRDS(
  degMatList_pd,
  "input/pheno/phenoData_full_list_degMat.rds"
)


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
