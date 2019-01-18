### SMRI Array ACC SZ vs CTRL
### Plot GSEA results
### Alexis Norris
### Created: 2018-05-25
### Modified: 2019-01-12


# Sources ------------------------------------------------------
### clusterProfiler: http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html


# Input files --------------------------------------------------
### Genelist, ranked
### Ranking metric = logFC * log10(PValue)
#geneList_ranked <- readRDS("output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_geneList_ranked.rds")

### GSEA results
#kegg_gseaList <- readRDS("output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_KEGG.rds")

### Full ACC edgeR results
### Don't use read_tsv because it messes up formatting
#full_res <- read.delim("output/edgeR/edgeR_ACC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE)


# Setup --------------------------------------------------------
### Packages
library(tidyverse)         # wrangling; plotting
library(clusterProfiler)   # GSEA
library(org.Hs.eg.db)      # Convert Entrez IDs to gene symbols
library(enrichplot)        # for network plot
library(ggraph)            # for network plot
library(igraph)            # for network plot
library(pathview)          # KEGG pathway plot

### Parameters
analysis_name <- "plot_gsea"
#fdr_edgeR <- 0.05
nPerm <- 10000


# Function: network plot ---------------------------------------
### enrichplot::cnetplot.enrichResult
### Requires enrichplot, ggraph, and igraph
cnetplot_mod <- function (x, showCategory = 5, foldChange = NULL, layout = "kk", 
                          colorEdge = FALSE, circular = FALSE, 
                          node_label = TRUE,  # pathway?
                          geneLabelSize = 2, genePointSize = 2, 
                          pathwayLabelSize = 4,  ...) 
{
  require(enrichplot)
  require(ggraph)
  require(igraph)
  
  if (circular) {
    layout <- "linear"
    geom_edge <- ggraph::geom_edge_arc
  }
  else {
    geom_edge <- ggraph::geom_edge_link
  }
  
  geneSets <- enrichplot:::extract_geneSets(x, showCategory)
  
  # Added set colors
  #colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")[1:length(names(geneSets))]
  #names(colors) <- names(geneSets)
  
  g <- enrichplot:::list2graph(geneSets)
  foldChange <- enrichplot:::fc_readable(x, foldChange)
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~category), alpha = 0.8)
  }
  else {
    edge_layer <- geom_edge(alpha = 0.8, color = "darkgrey")
  }
  if (!is.null(foldChange)) {   # Plotting gene nodes
    # Use blue-red heatmap instead
    fc <- foldChange[V(g)$name[(n + 1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n + 1):length(V(g))] <- fc
    palette <- enrichplot:::fc_palette(fc)
    # For categories, added this -- doesn't work
    #V(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    # Plot
    p <- ggraph(g, layout = layout, circular = circular) + 
      edge_layer + 
      #geom_node_point(aes_(color = ~as.numeric(as.character(color)), size = ~size)) + 
      geom_node_point(aes_(color = ~as.numeric(as.character(color))), size = genePointSize) + # genes
      #scale_color_gradientn(name = "fold change", colors = palette, na.value = "#E5C494")
      scale_fill_gradient2(
        low = "blue", 
        mid = "white",
        high = "red", 
        midpoint = 0,
        space = "Lab",
        na.value = "grey95", 
        guide = "colourbar", 
        aesthetics = "fill"
      ) +
      labs(fill = "logFC")
  }
  else {
    V(g)$color <- "#B3B3B3"               # gray
    V(g)$color[1:n] <- "#E5C494"          # tan
    p <- ggraph(g, layout = layout, circular = circular) + 
      edge_layer + 
      #geom_node_point(aes_(color = ~I(color), size = ~size))
      geom_node_point(size = genePointSize)   # would like to add color here to match ~category above!!
  }
  p <- p + 
    scale_size(range = c(3, 10), breaks = unique(round(seq(min(size), max(size), length.out = 4)))) + 
    theme_void()
  if (node_label) {  ### Plotting pathway names
    p <- p + 
      geom_node_text(
        aes_(label = ~name), 
        repel = TRUE, 
        #fontface = "bold",
        size = pathwayLabelSize
      )
  }
  p <- p + theme(legend.position = "none")
  return(p)
}

# Load ACC DEG geneLists used ----------------------------------
### Genelist, ranked
### Ranking metric = logFC * log10(PValue)
geneList_ranked <- readRDS("output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_ranked_geneList.rds")


# Load ACC edgeR results ---------------------------------------
### Load full ACC edgeR results
### to plot logFC (rather than ranking metric of logFC*log10(PValue))
full_res <- read.delim("output/edgeR/edgeR_ACC_SZvsCTRL_res.tsv", stringsAsFactors = FALSE)

### Genes with logFC
genes_fc <- full_res$logFC
names(genes_fc) <- full_res$entrez_id


# Load GSEA results --------------------------------------------
### Load
resList <- readRDS("output/gsea/gsea_ACC_SZvsCTRL_0.05FDR_KEGG.rds")

### Only results for down-regulated, so subset
res <- resList$down


# Dot plot -----------------------------------------------------
dot <- dotplot(res) + 
  labs(x = "Gene Ratio", title = "GSEA (down DEGs)") +
  theme(plot.title = element_text(hjust = 0.5))
print(dot)
ggsave("graphics/gsea/KEGG_down_dotPlot.png", height = 4, width = 8)


# Enrichment map -----------------------------------------------
enrichmap <- emapplot(res) + 
  labs(title = "Enrichment Map") +
  theme(plot.title = element_text(hjust = 0.5))
print(enrichmap)
ggsave("graphics/gsea/KEGG_down_enrichmentMap.png", height = 3, width = 5)


# Network plot -------------------------------------------------
### cnetplot: consider potentially biological complexities in which gene is in multiple pathways
### enrichplot::cnetplot.enrichResult
### Requires enrichplot, ggraph, and igraph

### I created modified version
network <- cnetplot_mod(
  # Convert Entrez IDs to gene symbols
  DOSE::setReadable(res, OrgDb = "org.Hs.eg.db", keytype = "ENTREZID"),
  colorEdge = TRUE,
  categorySize = FALSE,
  node_label = TRUE,
  geneLabelSize = 3,     # I think this is plotted same as pathway...
  genePointSize = 4,
  pathwayLabelSize = 4,
  #circular = TRUE, 
  #categorySize = "pvalue",                # alternative: geneNum'
  foldChange = genes_fc #geneList_ranked$down
) + 
  labs(title = "GSEA Network Plot") +
  theme(plot.title = element_text(hjust = 0.5))
print(network)
ggsave("graphics/gsea/KEGG_down_networkPlot.png", height = 12, width = 12)
ggsave("graphics/gsea/KEGG_down_networkPlot.eps", height = 12, width = 12)


# Plot KEGG pathway (via pathview) -----------------------------
### Change directory because pathview plots to cwd
setwd("graphics/gsea")

### Plot
for (i in as.data.frame(res)$ID) {
  pathview(
    gene.data = genes_fc,
    gene.idtype = "entrez",
    pathway.id = i,
    species = "hsa",
    limit = list(gene = max(abs(genes_fc)), cpd = 1),
    kegg.dir = ".",
    low = list(gene = "blue", cpd = "blue"), 
    mid = list(gene = "white", cpd = "gray"), 
    high = list(gene = "red", cpd = "yellow"), 
    na.col = "gray" # "transparent"
  )
}

### Fix cwd
setwd("../..")


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))