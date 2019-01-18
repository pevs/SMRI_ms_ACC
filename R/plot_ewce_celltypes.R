### Stanley ACC SZ vs CTRL
### Plot cell type enrichment analysis results
### Method: EWCE
### Cell data (single cell transcriptome): 
###   # Darmanis et al PNAS 2015 (human brain cell types; fresh tissue from surgery; multiple individuals)
###   # Lake et al Science 2016 (human neuron subtypes; postmortem; 1 individual [female])
### Alexis Norris
### Created: 2018-11-29
### Modified: 2019-01-14


# Input files --------------------------------------------
### EWCE enrichment results
#res <- read_csv(paste0("output/ewce/ewce_res_", reps, "reps.csv"))

### EWCE specificity data
#specList <- readRDS("output/ewce/ewce_res_gene_specificityList.rds")

### DEGs from edgeR
#sig <- read_tsv(paste0("output/edger/edgeR_ACC_SZvsCTRL_res_", fdr_edger, "FDR.tsv"))


# Setup --------------------------------------------------
### Packages
library(tidyverse)          # wrangling/plotting
library(patchwork)          # combine plots

### Parameters
analysis_name <- "plot_ewce_celltypes"
reps <- 20000               # num bootstraps used for ewce
fdr_ewce <- 0.05            # for EWCE analysis
fdr_edger <- 0.05           # for edgeR DE analysis


# Load EWCE results --------------------------------------
res <- read_csv(paste0("output/ewce/ewce_res_", reps, "reps.csv"))


# Plot enrichment ----------------------------------------
df <- res

### Note significant cell types
ast <- rep("", dim(df)[1])
ast[df$FDR < fdr_ewce] <- "*"
df$ast <- ast

### Remove negative SD values (since they are essentially 0?) !!!!!!!!!!!!
### per ewce vignette, no reasoning given
df$sd_from_mean[df$sd_from_mean < 0] <- 0

### Make room for * to denote sig cell type
upperLim <- max(df$sd_from_mean)
df$y_ast <- df$sd_from_mean*1.05
df$sd <- df$sd_from_mean

### Reorder genelists
df$genelist <- fct_relevel(
  df$genelist,
  "up", "down", "all", "degrad"
)

### Add column to cleanly separate for plot facets 
df$Class <- ifelse(
  df$celldata == "darmanis", "Cell type",
  ifelse(
    df$CellType %in% df$CellType[grep("^In", df$CellType)], "Inhibitory Neuron",
    ifelse(
      df$CellType %in% df$CellType[grep("^Ex", df$CellType)], "Excitatory Neuron",
      NA # nothing should be NA!
    )
  )
)

### Shorten cell type names
df$CellType <- recode(
  df$CellType,
  # Darmanis
  Astrocytes = "Astro",
  Endothelial = "Endo",
  Microglia = "Micro",
  Neurons = "Neuron",
  Oligodendrocytes = "Oligo"
)

### Reorder cell types
df$CellType <- fct_relevel(
  df$CellType,
  # Darmanis
  "Neuron", "Oligo",
  "Astro",  "Micro", "Endo",
  # Lake
  c(paste0("In", 1:8), paste0("Ex", 1:8))
)

### For colors
df$Color <- ifelse(
  df$CellType %in% df$CellType[grep("^In", df$CellType)], "In", 
  ifelse(
    df$CellType %in% df$CellType[grep("^Ex", df$CellType)], "Ex",
    as.character(df$CellType)
  )
)

### Cell type colors
cols <- c(
  # Darmanis
  "Neuron" = "black", 
  "Oligo" = "grey15", 
  "Astro" = "grey30", 
  "Micro"= "grey45",
  "Endo" = "grey60",
  # Lake             <--  Not done
  # color by cell type (neuron) class
  "In" = "blue3", 
  "Ex" = "firebrick3"
)

### Plot
cell_plot <- df %>%
  filter(genelist != "all") %>%
  ggplot(aes(x = CellType, label = ast)) + 
  geom_bar(aes(y = sd, color = Color, fill = Color), stat = "identity") + 
  geom_text(aes(y = y_ast), size = 10, color = "black") + 
  facet_grid(genelist~Class, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_y_continuous(limits = c(0, 1.4*upperLim), expand = c(0,0)) + 
  theme_bw() + theme(
    axis.text.x = element_text(size = 8.5, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold")
  ) +
  labs(
    x = "",   # "Cell type"
    y = "SDs from mean"
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

### Print 
print(cell_plot)
ggsave("graphics/ewce/ewce_celltype.png", height = 6, width = 10)
ggsave("graphics/ewce/ewce_celltype.eps", height = 6, width = 10)


# Plot gene specificity (heatmap) ------------------------
### Add cell type specificty to DEG table
# Load cell type specificty (from ewce_celltypes.R)
specList <- readRDS("output/ewce/ewce_gene_specificityList.rds")

# Load DEG table
sig <- read_tsv(paste0("output/edger/edgeR_ACC_SZvsCTRL_res_", fdr_edger, "FDR.tsv"))

# Add cell type specificty to DEGs table
sig <- sig %>%
  left_join(specList$darmanis, by = "gene_symbol") %>%
  left_join(specList$lake, by = "gene_symbol")

### Get celltypes enriched for SZ DEGs
cells <- res %>%
  filter(genelist != "degrad") %>%
  filter(FDR < fdr_ewce) %>%
  distinct(CellType) %>%
  as_vector() %>%
  as.character() 

### Require a minimum % specificity for one of the enriched cell types
spec_stats <- sig %>%
  select(gene_symbol, one_of(cells)) %>%
  gather(CellType, Specificity, -gene_symbol, na.rm = TRUE) %>%
  group_by(gene_symbol) %>%
  summarise(
    spec_min = min(Specificity),
    spec_max= max(Specificity)
  ) 
write_csv(spec_stats, "output/ewce/ewce_gene_specificity_DEGs_stats.csv")
spec_cutoff <- 1/3
spec_genes <- sig %>%
  filter(gene_symbol %in% filter(spec_stats, spec_max > spec_cutoff)$gene_symbol) 

### Use clustering to order genes in plot
gene_mat <- spec_genes %>%
  select(gene_symbol, one_of(cells)) %>%
  as.data.frame()
rownames(gene_mat) <- gene_mat$gene_symbol
gene_mat[is.na(gene_mat)] <- 0
gene_mat <- gene_mat %>%
  select(-gene_symbol)
spec_genes$GeneOrder <- hclust(dist(as.matrix(gene_mat)))$order
cell_order <- data_frame(
  "CellType" = colnames(gene_mat),
  "CellOrder" = hclust(dist(t(as.matrix(gene_mat))))$order
)

### Wrangle for plotting
spec_genes <- spec_genes %>%
  select(direction, gene_symbol, GeneOrder, one_of(cells)) %>%
  gather(CellType, Specificity, -direction, -gene_symbol, -GeneOrder, na.rm = FALSE) %>%
  left_join(cell_order, by = "CellType")

### Abbrev CellTypes for plot
spec_genes$CellType <- recode(
  spec_genes$CellType,
  Oligodendrocytes = "Oligo.",
  Neurons = "Neuron"
)

### Plot
spec_plot <- spec_genes %>%
  filter(direction != "all") %>%
  ggplot(aes(
    # Plot my clustering of genes
    #x = fct_reorder(gene_symbol, GeneOrder, .desc = FALSE),
    # Plot alphabetically
    x = gene_symbol,
    y = CellType,
    fill = Specificity
  )) + 
  geom_tile(color = "white") + 
  facet_grid(~direction, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 6),
    axis.text.y = element_text(face = "bold"),
    strip.background.x = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_blank()
  ) +
  labs(
    title = "Specificity of DEGs in Enriched Cell Types",
    subtitle = paste0("min specificity in an enriched celltype = ", round(spec_cutoff, 3)*100, "%")
  ) +
  #scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90") +
  scale_fill_gradient2(
    low = "white", 
    mid = "white",
    high = "firebrick", 
    midpoint = 0,
    space = "Lab",
    na.value = "grey95", 
    guide = "colourbar", 
    aesthetics = "fill"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))  
print(spec_plot)
ggsave("graphics/ewce/ewce_deg_celltype_specificity.png", height = 5, width = 20)
ggsave("graphics/ewce/ewce_deg_celltype_specificity.eps", height = 5, width = 20)


# Plot gene specificity (stacked) ------------------------
### Wrangle for plotting
spec_genes2 <- sig %>%
  select(direction, gene_symbol, logFC, PValue, FDR) %>%
  left_join(specList$darmanis, by = "gene_symbol") %>%
  left_join(specList$lake, by = "gene_symbol") %>%
  filter(gene_symbol %in% filter(spec_stats, spec_max > spec_cutoff)$gene_symbol)  %>%
  gather(CellType, Specificity, -direction, -gene_symbol, -logFC, -PValue, -FDR, na.rm = FALSE) %>%
  mutate(dataset = ifelse(
    CellType %in% c(paste0("In", 1:8), paste0("Ex", 1:8)), 
    "Lake", "Darmanis"
  ))

### Celltype ordering for plot
spec_genes2$CellType <- fct_relevel(
  spec_genes2$CellType, c(
    "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia", "Endothelial",
    "In1", "In2", "In3", "In4", "In5", "In6", "In7", "In8",
    "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8"
  )
)

### Abbrev CellTypes for plot
spec_genes2$CellType <- recode(
  spec_genes2$CellType,
  Neurons = "Neuron",
  Oligodendrocytes = "Oligo",
  Astrocytes = "Astro",
  Microglia = "Micro",
  Endothelial = "Endo"
)

### Colors
cols2 <- c(
  # Darmanis
  "Neuron" = "blue3", 
  "Oligo" = "firebrick3", 
  "Astro" = "grey30", 
  "Micro"= "grey45",
  "Endo" = "grey60",
  # Lake             <--  Not done
  # color by cell type (neuron) class
  "In1" = "blue3", 
  "In2" = "blue3", 
  "In3" = "blue3", 
  "In4" = "blue3", 
  "In5" = "blue3", 
  "In6" = "blue3", 
  "In7" = "blue3", 
  "In8" = "blue3", 
  "Ex1" = "firebrick3",
  "Ex2" = "firebrick3",
  "Ex3" = "firebrick3",
  "Ex4" = "firebrick3",
  "Ex5" = "firebrick3",
  "Ex6" = "firebrick3",
  "Ex7" = "firebrick3",
  "Ex8" = "firebrick3"
)

### Plot
spec_plot_stacked <- spec_genes2 %>%
  ggplot(aes(
    # Plot my clustering of genes
    #x = fct_reorder(gene_symbol, GeneOrder, .desc = FALSE),
    # Plot alphabetically
    x = fct_reorder(gene_symbol, Specificity, .desc = TRUE),
    y = Specificity,
    group = CellType, #order(CellType, decreasing = TRUE),
    color = CellType,
    fill = CellType
  )) + 
  geom_bar(stat = "identity", position = "stack") + 
  facet_grid(direction~dataset, scales = "free", space = "free") + 
  theme_bw() + theme(
    legend.position = "right",
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
    panel.spacing = unit(0, "lines"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Specificity of DEGs in Enriched Cell Types",
    subtitle = paste0("min specificity in an enriched celltype = ", round(spec_cutoff, 3)*100, "%")
  ) +
  scale_color_manual(values = cols2) +
  scale_fill_manual(values = cols2) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))  +
  coord_flip() +
  guides(col = guide_legend(ncol = 1))
print(spec_plot_stacked)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_stacked.png", height = 20, width = 10)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_stacked.eps", height = 20, width = 10)


# Plot gene specificity for just In6 ---------------------
In6_plot <- sig %>%
  select(gene_symbol, direction, In6) %>%
  filter(In6 > 0.1, !is.na(In6)) %>%
  ggplot(aes(
    x = fct_reorder(gene_symbol, -In6),
    y = In6,
    fill = In6
  )) + 
  geom_bar(stat = "identity") + 
  facet_grid(~direction, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    strip.background.x = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    panel.spacing.x = unit(0, "lines"), 
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene", y = "Specificity",   
    title = "Specificity of DEGs in In6"
  ) +
  scale_fill_gradient2(
    low = "white", 
    mid = "white",
    high = "firebrick", 
    midpoint = 0,
    space = "Lab",
    na.value = "grey95", 
    guide = "colourbar", 
    aesthetics = "fill"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))  
print(In6_plot)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_In6.png", height = 5, width = 12)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_In6.eps", height = 5, width = 12)


# Plot gene specificity for just In4 ---------------------
In4_plot <- sig %>%
  select(gene_symbol, direction, In4) %>%
  filter(In4 > 0.1, !is.na(In4)) %>%
  ggplot(aes(
    x = fct_reorder(gene_symbol, -In4),
    y = In4,
    fill = In4
  )) + 
  geom_bar(stat = "identity") + 
  facet_grid(~direction, scales = "free_x", space = "free_x") + 
  theme_bw() + theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    strip.background.x = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    panel.spacing.x = unit(0, "lines"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Gene", y = "Specificity",   
    title = "Specificity of DEGs in In4"
  ) +
  scale_fill_gradient2(
    low = "white", 
    mid = "white",
    high = "firebrick", 
    midpoint = 0,
    space = "Lab",
    na.value = "grey95", 
    guide = "colourbar", 
    aesthetics = "fill"
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))  
print(In4_plot)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_In4.png", height = 5, width = 12)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_In4.eps", height = 5, width = 12)


# Combine specificity plots ------------------------------
### Specificity heatmap for all cell types plus In4 and In6
specs_In6_In4 <- spec_plot * theme(
  legend.position = "bottom", 
  legend.title = element_text(size = 20, face = "bold", vjust = 1),
  axis.text.y = element_text(size = 14, face = "bold")
) + {
  In6_plot * theme(
    axis.title.y = element_text(vjust = -5, size = 20, face = "bold")
  ) + 
     In4_plot * theme(
       axis.title.y = element_blank(),
       axis.text.y = element_blank()
     ) +
    plot_layout(ncol = 2) & theme(
      legend.position = "none", 
      axis.text.x = element_text(size = 7)
    ) 
} + 
  plot_layout(nrow = 2) & theme(
    plot.title = element_blank(), 
    plot.subtitle = element_blank(), 
    axis.title.x = element_blank()
  )
print(specs_In6_In4)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_plus_In6_In4.png", height = 10, width = 20)
ggsave("graphics/ewce/ewce_deg_celltype_specificity_plus_In6_In4.eps", height = 10, width = 20)


# Write out packages and versions ------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
