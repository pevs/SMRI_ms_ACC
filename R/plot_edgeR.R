### Stanley ACC SZ vs CTRL
### edgeR plots
### Alexis Norris
### Created: 2018-11-15
### Modified: 2019-01-16


# Notes --------------------------------------------------------
### Plot volcano and p-value histogram for edgeR results
### Comparison of feature sizes for feature levels


# Input files --------------------------------------------------
### DFull results 
#res <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_resList.rds") 

### Significant results
#sig <- read_tsv("output/edger/edgeR_ACC_SZvsCTRL_res_0.05FDR.tsv")

### Posthoc results -- full edgeR feature data
#antipsych <- readRDS("output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_resList.rds")
#pfc_hpc <- readRDS("output/edgeR/posthoc/edgeR_PFC_HPC_resList.rds")

### Posthoc results -- gene-collapsed
#read.delim(paste0("output/edgeR/posthoc/edgeR_+_res.tsv"), stringsAsFactors = FALSE)


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling
library(patchwork)         # combine plots

### Parameters
analysis_name <- "plot_edgeR"
fdr_cutoff <- 0.05


# Load ACC edgeR results ---------------------------------------
res <- readRDS("output/edgeR/edgeR_ACC_SZvsCTRL_resList.rds") %>%
  # Collapse levels
  bind_rows(.id = "feature_level") %>%
  select(feature_level, gene_symbol, feature_id, logFC, PValue, start, end)
res$feature_level <- gsub(".*_", "", res$feature_level)

# Plotting order
res$feature_level <- fct_relevel(
  res$feature_level,
  "genes", "exons", "jxns", "regions"
)

### Significant
sig <- read_tsv(paste0("output/edgeR/edgeR_ACC_SZvsCTRL_res_", fdr_cutoff, "FDR.tsv"))

### Split, to highlight direction of ACC DEGs
sig_down <- sig %>%
  filter(logFC < 0)
sig_up <- sig %>%
  filter(logFC > 0)


# Volcano plot -------------------------------------------------
### Color if DE feature
res$ACC_DEG <- ifelse(
  res$feature_id %in% sig_down$feature_id, "Down", ifelse(
    res$feature_id %in% sig_up$feature_id, "Up", "n.s."
  )
)

### Plot
volcano <- res %>%
  ggplot(aes(x = logFC, y = -log10(PValue))) +
  geom_point(data = filter(res, ACC_DEG == "n.s."), shape = "o", color = "grey50") + 
  geom_point(data = filter(res, ACC_DEG == "Down"), shape = "o", color = "blue2") + 
  geom_point(data = filter(res, ACC_DEG == "Up"), shape = "o", color = "firebrick2") +
  facet_grid(feature_level~.) +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.position = "none", 
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  ) +
  labs(
    x = expression(paste(log[2], " Fold Change")), 
    y = expression(paste(-log[10], " p-value")),
    title = "Volcano plot"
  ) +
  scale_y_continuous(expand = c(0,0,0.05,0))


# Distribution of p-values -------------------------------------
phisto <- res %>%
  ggplot(aes(x = PValue)) + 
  geom_histogram(binwidth = 0.05, na.rm = TRUE, color = "black", fill = "grey25") +
  facet_grid(feature_level~., scales = "free_y") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    strip.background.y = element_rect(fill = "white"),
    strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8)
  ) +
  labs(
    x = "p-value", 
    y = "Number of features",
    title = "Distribution of p-values"
  ) +
  scale_y_continuous(expand = c(0,0,0.1,0), breaks = scales::pretty_breaks(n = 2)) +
  scale_x_continuous(expand = c(0,0,0,0), breaks = c(0, 0.25, 0.5, 0.75, 1))


# Volcano + p-value histogram ----------------------------------
### With patchwork
plotList <- list()   # to hold all edgeR test results' plots
plotList$ACC <- volcano + phisto
print(plotList$ACC)
ggsave("graphics/edgeR/edgeR_volcano_pvalueHisto_ACC.png", width = 5, height = 5)
rm(volcano, phisto)


# Feature sizes ------------------------------------------------
### Size distribution of features --> especially regions! (do they tend to be smaller thatn other feature levels?)
### Remove jxns, since their coordinates are for exon-exon boundaries, not the length of expressed read/transcript
feature_sizes <- res %>%
  filter(feature_level != "jxns") %>%
  mutate(size = end - start) %>%
  group_by(feature_level) %>%
  summarise(min = min(size), max = max(size), mean = mean(size), sd = sd(size), median = median(size))
feature_sizes_sig <- sig %>%
  filter(feature_level != "jxns") %>%
  mutate(size = end - start) %>%
  group_by(feature_level) %>%
  summarise(min = min(size), max = max(size), mean = mean(size), sd = sd(size), median = median(size))
bind_rows(
  mutate(feature_sizes, genelist = "All features"),
  mutate(feature_sizes_sig, genelist = paste("FDR <", fdr_cutoff))
) %>%
  write_tsv("output/edger/feature_sizes_by_featureLevel.tsv")

### Plot
sz_cols <- c("genes" = "#984ea3", "exons" = "#377eb8", "regions" = "#4daf4a")
keep <- c("feature_level", "start", "end")
sizes <- bind_rows(
  mutate(select(res, one_of(keep)), features = "All"),
  mutate(select(sig, one_of(keep)), features = "Significant")
) %>%
  filter(feature_level != "jxns") %>%
  mutate(size = end - start) %>%
  ggplot(aes(
    x = reorder(feature_level, -size), y = log(size),
    color = feature_level, fill = feature_level, alpha = 0.9
  )) + 
  geom_violin() +
  facet_grid(~features) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "none",
    strip.background.x = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 10)
  ) +
  scale_color_manual(values = sz_cols) +
  scale_fill_manual(values = sz_cols) + 
  labs(
    x = NULL,
    y = expression(paste(log[10], " Size (bp)")),
    title = "Feature sizes by feature level"
  )
print(sizes)
ggsave("graphics/edger/feature_sizes_by_featureLevel.png", height = 5, width = 5)

### Size vs p-value
sizes_pval <- res %>%
  filter(feature_level != "jxns") %>%
  mutate(size = end - start) %>%
  ggplot(aes(
    x = log10(size), 
    y = -log10(PValue), 
    group = feature_level,
    color = feature_level, fill = feature_level
  )) + 
  geom_point(alpha = 0.7, shape = ".") + 
  geom_smooth() +
  geom_rug() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  scale_color_manual(values = sz_cols) +
  scale_fill_manual(values = sz_cols) + 
  labs(
    x = expression(paste(log[10], " size (bp)")),
    y = expression(paste(-log[10], " p-value")),
    title = "Feature Size ~ Significance"
  ) +
  scale_x_continuous(expand = c(0,0,0,0)) +
  scale_y_continuous(expand = c(0,0,0,0))
print(sizes_pval)
ggsave("graphics/edger/feature_sizes_vs_pvalue.png", height = 5, width = 5)


# Plot edgeR results for posthoc -------------------------------
### Load full results (not gene-collapsed)
antipsych <- readRDS("output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_resList.rds")
pfc_hpc <- readRDS("output/edgeR/posthoc/edgeR_PFC_HPC_resList.rds")
post_resList <- list(
  "ACC_Antipsychotics" = antipsych,
  "PFC" = pfc_hpc[grep("PFC", names(pfc_hpc))],
  "HPC" = pfc_hpc[grep("HPC", names(pfc_hpc))]
)
posthocTables <- list(
  "ACC_Antipsychotics" = read.delim("output/edgeR/posthoc/edgeR_ACC_AntipsychoticsKG_res.tsv", stringsAsFactors = FALSE),
  "PFC" = read.delim("output/edgeR/posthoc/edgeR_PFC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE),
  "HPC" = read.delim("output/edgeR/posthoc/edgeR_HPC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE)
)

### Plot for each
for (i in names(posthocList)) {
  print(i)
  
  ### Prep
  # Collapse levels
  post_res <- post_resList[[i]] %>%
    bind_rows(.id = "feature_level") %>%
    select(feature_level, gene_symbol, feature_id, logFC, PValue)
  post_res$feature_level <- gsub(".*_", "", post_res$feature_level)
  
  # Plotting order
  post_res$feature_level <- fct_relevel(
    post_res$feature_level,
    "genes", "exons", "jxns", "regions"
  )
  
  # Color if ACC DEG (FEATURE for antipsychotics)
  if (i == "ACC_Antipsychotics") {
    post_res$ACC_DEG <- ifelse(
      post_res$feature_id %in% sig_down$feature_id, "Down", ifelse(
        post_res$feature_id %in% sig_up$feature_id, "Up", "n.s."
      )
    )
  }
  if (i %in% c("PFC", "HPC")) { 
    ### Get posthoc features for ACC DEGs
    degs_down <- posthocTables[[i]] %>%
      filter(gene_symbol %in% sig_down$gene_symbol)
    degs_up <- posthocTables[[i]] %>%
      filter(gene_symbol %in% sig_up$gene_symbol)
    
    ### Get results for the best feature
    post_res$ACC_DEG <- ifelse(
      post_res$feature_id %in% degs_down$feature_id, "Down", ifelse(
        post_res$feature_id %in% degs_up$feature_id, "Up", "n.s."
      )
    )
    rm(degs_down, degs_up)
  }
  print(table(post_res$ACC_DEG))
  
  ### Volcano plot
  post_volcano <- post_res %>%
    ggplot(aes(x = logFC, y = -log10(PValue))) +
    geom_point(data = filter(post_res, ACC_DEG == "n.s."), shape = "o", color = "grey50") + 
    geom_point(data = filter(post_res, ACC_DEG == "Down"), shape = "o", color = "blue2") + 
    geom_point(data = filter(post_res, ACC_DEG == "Up"), shape = "o", color = "firebrick2") + 
    # Highlight discordant byb plotting them in front
    #geom_point(data = filter(post_res, ACC_DEG == "Down", logFC > 0), shape = "o", color = "blue2") + 
    #geom_point(data = filter(post_res, ACC_DEG == "Up", logFC < 0), shape = "o", color = "firebrick2") + 
    facet_grid(feature_level~.) +
    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      legend.position = "none", 
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8)
    ) +
    labs(
      x = expression(paste(log[2], " Fold Change")), 
      y = expression(paste(-log[10], " p-value")),
      title = "Volcano plot"
    ) +
    scale_y_continuous(expand = c(0,0,0.05,0))
  
  ### Distribution of p-values
  post_phisto <- post_res %>%
    ggplot(aes(x = PValue)) + 
    geom_histogram(binwidth = 0.05, na.rm = TRUE, color = "black", fill = "grey25") +
    facet_grid(feature_level~., scales = "free_y") +
    theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      strip.background.y = element_rect(fill = "white"),
      strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8)
    ) +
    labs(
      x = "p-value", 
      y = "Number of features",
      title = "Distribution of p-values"
    ) +
    scale_y_continuous(expand = c(0,0,0.1,0), breaks = scales::pretty_breaks(n = 2)) +
    scale_x_continuous(expand = c(0,0,0,0), breaks = c(0, 0.25, 0.5, 0.75, 1))
  
  
  ### Combine
  plotList[[i]] <- post_volcano + post_phisto
  print(plotList[[i]])
  ggsave(paste0("graphics/edgeR/edgeR_volcano_pvalueHisto_", i, ".png"), width = 5, height = 5)

  ### Clean-up
  rm(post_res, post_volcano, post_phisto)
}


# Combine all --------------------------------------------------
### Save
saveRDS(plotList, "graphics/edger/edgeR_volcano_pvalueHisto_all_objects.rds")

#fig_all <- readRDS("graphics/edger/edgeR_volcano_pvalueHisto_all_objects.rds")
### Add figure parts (A-D)
fig_a <- wrap_elements(plotList$ACC) + ggtitle("A") + theme(plot.title = element_text(face = "bold", size = 20))
fig_b <- wrap_elements(plotList$ACC_Antipsychotics) + ggtitle("B") + theme(plot.title = element_text(face = "bold", size = 20))
fig_c <- wrap_elements(plotList$PFC) + ggtitle("C") + theme(plot.title = element_text(face = "bold", size = 20))
fig_d <- wrap_elements(plotList$HPC) + ggtitle("D") + theme(plot.title = element_text(face = "bold", size = 20))

### Print combined
fig_all <- (fig_a | fig_b) / (fig_c | fig_d) 
print(fig_all)
ggsave("graphics/edgeR/edgeR_volcano_pvalueHisto_all.png", width = 10, height = 10)
ggsave("graphics/edgeR/edgeR_volcano_pvalueHisto_all.eps", width = 10, height = 10)


# Compare dataseets' gene-collapsed res ------------------------
### Should I collapse the genes further??? (e.g. remove the down AND up feature, only keeping lowest p?)

### Create fxn
### Source for overplotting reduction methods: https://stackoverflow.com/questions/7714677/scatterplot-with-too-many-points
plot_compare_edger <- function (
  # data options = c("ACC_SZvsCTRL", "ACC_AntipsychoticsKG","PFC_SZvsCTRL", "HPC_SZvsCTRL")
  data1 = "ACC_SZvsCTRL", 
  data2 = "ACC_AntipsychoticsKG", 
  dir1 = "output/edgeR/",
  dir2 = "output/edgeR/posthoc/",
  
  # edgeR stats to plot, options: c("logCPM", "logFC", "neg_log10_PValue", "logCPM",  "LR", "PValue", "FDR", "signed_LR")
  statvals = c("logCPM", "logFC", "neg_log10_PValue", "signed_LR"), 
  
  # Options to reduce overplotting: c("alpha", "hex", "density_raster", "density_contours)
  geom_reduce = "none",  
  shape = ".", alpha = 0.4, color = "black", fill = "black",
  
  # Add points
  geom_point = TRUE, 
  point_shape = ".", point_alpha = 0.7, point_color = "black", point_fill = "black",
  
  # Add line
  line_stat = "auto", # c("auto", "lm", "glm", "gam", "loess")
  line_se = TRUE, line_color = "blue4", line_fill = "blue1",
  
  # Add rug
  rug = FALSE,
  
  # Add margin (with ggExtra package) options: c("density", "histogram", "boxplot", "violin") 
  margin = NULL, 
  margin_size = 5,
  margin_axis = "both", # c("both", "x", "y")
  margin_color = "black", margin_fill = "white",
  
  #### Format options
  title_face = "bold", title_size = 10,
  ...
) {

  ### Load gene-collapsed edgeR results
  dataList <- list(
    read.delim(paste0(dir1, "edgeR_", data1, "_res.tsv"), stringsAsFactors = FALSE),
    read.delim(paste0(dir2, "edgeR_", data2, "_res.tsv"), stringsAsFactors = FALSE)
  ) %>%
    lapply(., function (df) {
      ### Do before select cols because need logFC
      df %>% 
        ### Log-transform
        mutate(
          neg_log10_PValue = -log10(PValue),
          signed_LR = ifelse(logFC < 0, -LR, LR)
        ) %>%
        ### Remove unnecessay columns
        select(gene_symbol, direction, one_of(statvals))
    })
  
  ### Join datasets
  both <- dataList[[1]] %>%
    inner_join(
      dataList[[2]],
      by = c("gene_symbol", "direction"), 
      suffix = c(paste0(".", data1), paste0(".", data2))
    )
  
  ### Plot
  pList <- lapply(statvals, function (val) {
    
    ### Get x and y values
    p <- both %>%
      ggplot(
        aes_string(
          x = paste0(val, ".", data1), 
          y = paste0(val, ".", data2)
        )
      )
    
    ### Reducing over-plotting
    # Option 1: use hexbins (requires hexbin pkg from CRAN)
    if (geom_reduce == "hex") {
      p <- p + 
        geom_hex() +
        scale_fill_viridis_c(direction = -1) #+ geom_point(shape = shape, col = "white") 
    }
    
    # Option 2:Density heatmap
    if (geom_reduce == "density_raster") {
      p <- p + 
        stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +       
        scale_fill_viridis_c(direction = -1) 
    }
    
    # Opt 3: Density contours
    if (geom_reduce == "density_contours") {
      p <- p +
        stat_density_2d(aes(fill = ..level..), geom = "polygon") +
        scale_fill_viridis_c(name = "density") 
    }
    
    ### Plot points
    if (geom_point) {
      p <- p +
        geom_point(alpha = point_alpha, shape = point_shape, color = point_color, fill = point_fill)
    }
    
    
    ### Add fit line
    if (!is.null(line_stat)) {
      p <- p +
        geom_smooth(method = line_stat, se = line_se, color = line_color, fill = line_fill)
    }
    
    ### Add rug
    if (rug) {
      p <- p + 
        geom_rug(alpha = alpha, color = color)
    }
    
    ### Formatting
    p <- p + theme_bw() +
      theme(legend.position = "bottom", legend.key.width = unit(5, "line")) +
      coord_cartesian() +
      scale_x_continuous(expand = c(0,0.05,0,0.05)) +
      scale_y_continuous(expand = c(0,0.05,0,0.05)) 
    
    ### Return
    p
  })
  
  ### Add margin with ggExtra
  ### ! can't be added to plotting step before this; MUST be separate
  if (!is.null(margin)) {
    pList <- lapply(pList, function (p) {
      ggExtra::ggMarginal(
        p,
        type = margin, margins = margin_axis,
        size = margin_size, col = margin_color, fill = margin_fill
      )
    })
  }
    
  ### Combine
  ### ggMarginal requires using wrap_plots instead of "+"
  if (length(statvals) == 1) {
    compare_plot <- pList[[1]]
  }
  if (length(statvals) == 2) {
    compare_plot <- wrap_plots(pList[[1]], pList[[2]])
  }
  if (length(statvals) == 3) {
    compare_plot <- wrap_plots(
      pList[[1]], 
      pList[[2]], 
      pList[[3]]
    ) 
  }
  if (length(statvals) == 4) {
    compare_plot <- wrap_plots(pList[[1]], pList[[2]]) /
        wrap_plots(pList[[3]], pList[[4]]) 
  }
  if (length(statvals) == 5) {
    compare_plot <- wrap_plots(pList[[1]], pList[[2]], pList[[3]]) /
      wrap_plots(pList[[4]], pList[[5]], plot_spacer()) 
  }
  if (length(statvals) == 6) {
    compare_plot <- wrap_plots(pList[[1]], pList[[2]], pList[[3]]) /
      wrap_plots(pList[[4]], pList[[5]], pList[[6]]) 
  }
  compare_plot <- wrap_elements(compare_plot) +
    ggtitle(paste(data1, "~", data2, "edgeR results for each's best feature for a gene")) +
    theme(plot.title = element_text(face = title_face, size = title_size, hjust = 0.5))
  
  ### Return
  return(compare_plot)
  
}

### Plot scatter with rug
for (i in c("ACC_AntipsychoticsKG","PFC_SZvsCTRL", "HPC_SZvsCTRL")) {
  data1 <- "ACC_SZvsCTRL"
  data2 <- i # "ACC_AntipsychoticsKG"
  compare_plot <- plot_compare_edger(
    # data options = c("ACC_SZvsCTRL", "ACC_AntipsychoticsKG","PFC_SZvsCTRL", "HPC_SZvsCTRL")
    data1 = data1, 
    data2 = data2, 
    dir1 = "output/edgeR/",
    dir2 = "output/edgeR/posthoc/",
    
    # edgeR stats to plot, options: c("logCPM", "logFC", "neg_log10_PValue", "logCPM",  "LR", "PValue", "FDR", "signed_LR")
    statvals = c("logCPM", "logFC", "signed_LR", "neg_log10_PValue"), 
    
    # method options: c("alpha", "hex", "density", "hex", "contours", "none")
    geom_reduce = "none", 
    alpha = 0.2, color = "black", fill = "black", shape = ".",
    
    # Points
    geom_point = TRUE,
    point_alpha = 0.2, point_color = "black", point_fill = "black", point_shape = ".",
    
    # Line
    line_stat = "lm", # c(NULL, "auto", "lm", "glm", "gam", "loess")
    line_se = TRUE, line_color = "blue2", line_fill = "grey20",
    
    # Rug/margin
    rug = TRUE,
    margin = NULL, # c(NULL, "density", "histogram", "boxplot", "violin") - from ggExtra
    margin_size = 5,
    margin_axis = "both", # c("both", "x", "y")
    
    #### Format options
    title_face = "bold",
    title_size = 10
  )
  png(filename = paste0("graphics/edgeR/edgeR_stats_", data1, "-", data2, "_point_rug.png"), width = 12, height = 12, units = "in", res = 600)
  print(compare_plot); dev.off()
  
  ### not working
  #ggsave(paste0("graphics/edgeR/edgeR_stats_", data1, "-", data2, "_ggsave.png"), width = 12, height = 12)
}

### Plot hex
for (i in c("ACC_AntipsychoticsKG","PFC_SZvsCTRL", "HPC_SZvsCTRL")) {
  data1 <- "ACC_SZvsCTRL"
  data2 <- i # "ACC_AntipsychoticsKG"
  compare_plot <- plot_compare_edger(
    # data options = c("ACC_SZvsCTRL", "ACC_AntipsychoticsKG","PFC_SZvsCTRL", "HPC_SZvsCTRL")
    data1 = data1, 
    data2 = data2, 
    dir1 = "output/edgeR/",
    dir2 = "output/edgeR/posthoc/",
    
    # edgeR stats to plot, options: c("logCPM", "logFC", "neg_log10_PValue", "logCPM",  "LR", "PValue", "FDR", "signed_LR")
    statvals = c("logCPM", "logFC", "signed_LR", "neg_log10_PValue"), 

    # method options: c("alpha", "hex", "density", "hex", "contours)
    geom_reduce = "hex", 
    alpha = 0.25, color = "black", fill = "black", shape = ".",
    
    # Points
    geom_point = TRUE,
    point_alpha = 0.2, point_color = "black", point_fill = "black", point_shape = ".",
    
    # Line
    line_stat = NULL, # c(NULL, "auto", "lm", "glm", "gam", "loess")
    line_se = TRUE, line_color = "blue4", line_fill = "blue1",
    
    # Rug/margin
    rug = FALSE,
    margin = "density", # c(NULL, "density", "histogram", "boxplot", "violin") - from ggExtra
    margin_size = 5,
    margin_axis = "both", # c("both", "x", "y")
    
    #### Format options
    title_face = "bold",
    title_size = 10
  )
  png(filename = paste0("graphics/edgeR/edgeR_stats_", data1, "-", data2, "_hex.png"), width = 12, height = 12, units = "in", res = 600)
  print(compare_plot); dev.off()
  
  ### not working
  #ggsave(paste0("graphics/edgeR/edgeR_stats_", data1, "-", data2, "_ggsave.png"), width = 12, height = 12)
}


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
