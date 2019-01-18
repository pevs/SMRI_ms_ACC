### SMRI Array SZ vs CTRL
### Outlier analysis with PCA
### Alexis Norris
### Created: 2018-11-18
### Modified: 2019-01-02


# Input files --------------------------------------------------
#dgeList <- readRDS("input/dge/ACC_PFC_HPC_filtered.rds")


# Setup --------------------------------------------------------  
### Packages
library(tidyverse)    # wrangling/plotting
library(ggrepel)      # plotting
library(edgeR)        # exprsData (DGEs)
library(patchwork)    # Combine plots
### Note that dplyr fxns are replaced by other pkg's fxns, so must call them with dplyr::<fxn_name>

### Parameters
analysis_name <- "pca"
cols <- c("CTRL" = "black", "SZ" = "blue2")
datasets <- c("ACC", "PFC", "HPC")
feature_levels <- c("genes", "exons", "jxns", "regions")
dges <- paste(rep(datasets, length(feature_levels)), feature_levels, sep = "_") # names(dgeList)


# PCA plot fxn -------------------------------------------------
### Modified from ggbiplot pkg
### Modified 2017-07-31; 2018-10-12 (include ggrepel of labels); 2018-12-16 (control plotting size of x vs y)
plot_pca <- function (
  pcobj, choices = 1:2, 
  scale = 1, obs.scale = 1 - scale, var.scale = scale, 
  groups = NULL, ellipse = TRUE, ellipse.prob = 0.69, 
  labels = NULL, labels.size = 3, labels.only = FALSE,
  alpha = 1, var.axes = FALSE, circle = FALSE, circle.prob = 0.69, 
  varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
  out.prob = 0.95, out.text = FALSE, out.label = FALSE, out.size = 3, out.names = NULL, ...) {
  
  stopifnot(length(choices) == 2)
  
  ### Get PCA data, based on the pca function used
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(
      pcobj$x, 
      2, 
      1/(d * nobs.factor), 
      FUN = "*"
    )
    v <- pcobj$rotation
  }
  if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(
      pcobj$scores, 
      2, 
      1/(d * nobs.factor), 
      FUN = "*"
    )
    v <- pcobj$loadings
  }
  if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(
      pcobj$ind$coord, 
      2, 
      1/(d * nobs.factor), 
      FUN = "*"
    )
    v <- sweep(
      pcobj$var$coord, 
      2, 
      sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 1]), 
      FUN = "/"
    )
  }
  if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  
  ### Combine data
  #choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(
    u[ , choices], 
    2, 
    d[choices]^obs.scale, 
    FUN = "*"
  ))
  v <- sweep(
    v, 
    2, 
    d^var.scale, 
    FUN = "*"
  )
  df.v <- as.data.frame(v[ , choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  df.u <- df.u * nobs.factor
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  
  ### Axis labels
    # Add prefix "PC"
  if (obs.scale == 0) u.axis.labs <- paste0("PC", choices)
  if (obs.scale != 0)  u.axis.labs <- paste0("PC", choices)
    # Add % variance to PC axis label
  u.axis.labs <- paste(
    u.axis.labs, 
    sprintf(
      "(%0.1f%%)", 
      100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)
    )
  )
  if (!is.null(labels)) df.u$labels <- labels
  if (!is.null(groups)) df.u$Group <- groups
  
  ### 
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else df.v$varname <- rownames(v)
  
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust <- with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  g <- ggplot(df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) #+ coord_equal()
  
  if (var.axes) {
    if (circle) {
      theta <- c(
        seq(-pi, pi, length = 50), 
        seq(pi, -pi, length = 50)
      )
      circle <- data.frame(
        xvar = r * cos(theta), 
        yvar = r * sin(theta)
      )
      g <- g + geom_path(
        data = circle, 
        color = muted("white"), size = 1/2, alpha = 1/3,
        show.legend = FALSE
      )
    }
    g <- g + geom_segment(
      data = df.v, 
      aes(x = 0, y = 0, xend = xvar, yend = yvar), 
      arrow = arrow(length = unit(1/2, "picas")), 
      color = muted("black"),
      show.legend = FALSE
    )
  }
  if (!is.null(df.u$labels)) if (labels.only) {
    if (!is.null(df.u$Group)) {
      g <- g + geom_text(
        aes(label = labels, color = Group), 
        size = labels.size,
        show.legend = FALSE
      )
    }
    if (is.null(df.u$Group)) {
      g <- g + geom_text(
        aes(label = labels), 
        size = labels.size,
        show.legend = FALSE
      )
    }
  }
  if (!is.null(df.u$labels)) if (!labels.only) {
    if (!is.null(df.u$Group)) {
      g <- g + 
        geom_point(
          aes(color = Group), 
          alpha = alpha, 
          show.legend = TRUE
        ) +
        ggrepel::geom_text_repel(
          aes(label = labels, color = Group), 
          size = labels.size,
          show.legend = FALSE
        ) 
    }
    if (is.null(df.u$Group)) {
      g <- g + 
        geom_point( 
          alpha = alpha, 
          show.legend = TRUE
        ) +
        ggrepel::geom_text_repel(
          aes(label = labels), 
          size = labels.size,
          show.legend = FALSE
        ) 
    }
  }
  if (is.null(df.u$labels)) {
    if (!is.null(df.u$Group)) {
      g <- g + geom_point(
        aes(color = Group), 
        alpha = alpha, 
        show.legend = TRUE
      ) 
    }
    if (is.null(df.u$Group)) {
      g <- g + geom_point(
        alpha = alpha,
        show.legend = FALSE
      )
    }
  }
  
  ### Plot ellipse for groups
  if (!is.null(df.u$Group) && ellipse) {
    theta <- c(
      seq(-pi, pi, length = 50), 
      seq(pi, -pi, length = 50)
    )
    circle <- cbind(
      cos(theta), 
      sin(theta))
    ellipse <- plyr::ddply(df.u, "Group", function (x) {
      if (nrow(x) <= 2) return(NULL)
      sigma <- var(cbind(
        x$xvar, 
        x$yvar
      ))
      mu <- c(
        mean(x$xvar), 
        mean(x$yvar)
      )
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(
        sweep(
          circle %*% chol(sigma) * ed, 
          2, 
          mu, 
          FUN = "+"
        ), 
        Group = x$Group[1]
      )
    })
    names(ellipse)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(
        data = ellipse, 
        aes(color = Group, group = Group),
        show.legend = FALSE
    )
    ellipse.out <- plyr::ddply(df.u, "Group", function(x) {
      if (nrow(x) <= 2) return(NULL)
      sigma <- var(cbind(
        x$xvar, 
        x$yvar
      ))
      mu <- c(
        mean(x$xvar), 
        mean(x$yvar)
      )
      ed <- sqrt(qchisq(out.prob, df = 2))
      data.frame(
        sweep(
          circle %*% chol(sigma) * ed, 
          2, 
          mu, 
          FUN = "+"
        ), 
        Group = x$Group[1]
      )
    })
    names(ellipse.out)[1:2] <- c("xvar", "yvar")
    edge <- plyr::ddply(ellipse.out, "Group", function (x) {
      data.frame(
        Group = x$Group[1], 
        xmin = min(x$xvar), 
        xmax = max(x$xvar), 
        ymin = min(x$yvar), 
        ymax = max(x$yvar)
      )
    })
    
    ### Add text/labels for samples
    if (!is.null(out.names)) df.u$names <- out.names
    if (is.null(out.names)) df.u$names <- rownames(df.u)
    
    outs <- list()
    for (i in 1:length(levels(df.u$Group))) {
      grp <- df.u[df.u$Group == levels(df.u$Group)[i], ]
      outs <- append(outs, c(
        rownames(grp[grp$xvar < edge$xmin[i] | grp$xvar > edge$xmax[i], ]), 
        rownames(grp[grp$yvar < edge$ymin[i] | grp$yvar > edge$ymax[i], ]))
      )
    }
    df.out <- df.u[unique(unlist(outs)), ]
 
    if (out.text) {
      g <- g + geom_text(
        data = df.out, 
        aes(label = names, color = Group), 
        position = "jitter",
        size = out.size, fontface = "bold",
        hjust = "inward", vjust = "inward",
        show.legend = FALSE
      )
    }
    
    if (out.label) {
      g <- g + geom_label(
        data = df.out, 
        aes(label = names, color = Group), 
        size = out.size, 
        show.legend = FALSE
      )
    }
  }
  if (var.axes) {
    g <- g + geom_text(
      data = df.v, 
      aes(
        label = varname, 
        x = xvar, y = yvar, 
        angle = angle, hjust = hjust
      ), 
      color = "black", size = varname.size,
      show.legend = FALSE
    )
  }
  g <- g + theme_bw() + theme(
    axis.title = element_text(face = "bold"), 
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
  return(g)
}



# Load DGEs ----------------------------------------------------
dgeList_all <- readRDS("input/dge/ACC_PFC_HPC_filtered.rds")

### Want all dges
dgeList <- dgeList_all[which(names(dgeList_all) %in% dges)]


# Plot PCA by Dx -----------------------------------------
### Plot in loop
pca_plot <- list()
for (i in names(dgeList)) {
  # Subset
  y <- dgeList[[i]]
  
  # Run PCA
  pca_obj <- prcomp(t(cpm(y, log = TRUE))) 
  
  # Plot
  pca_plot[[i]] <- plot_pca(
    pca_obj,
    groups = y$samples$Dx
  ) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 8, face = "plain"),
      axis.text = element_text(size = 5.5),
      legend.position = "none"
    ) + 
    scale_color_manual(values = cols) 
}

### Combine into one figure
pca_plot_all <- pca_plot$ACC_genes + 
  pca_plot$PFC_genes +
  pca_plot$HPC_genes + plot_layout(ncol = 3) +
  pca_plot$ACC_exons + 
  pca_plot$PFC_exons +
  pca_plot$HPC_exons + plot_layout(ncol = 3) +
  pca_plot$ACC_jxns + 
  pca_plot$PFC_jxns +
  pca_plot$HPC_jxns + plot_layout(ncol = 3) +
  pca_plot$ACC_regions + 
  pca_plot$PFC_regions +
  pca_plot$HPC_regions + plot_layout(ncol = 3) 
print(pca_plot_all)
ggsave("graphics/pca/pca_Dx.png", height = 8, width = 6)
ggsave("graphics/pca/pca_Dx.eps", height = 8, width = 6)

### Write out key for rows/cols
write_lines(
  "Rows (top to bottom) are: genes, exons, junctions, regions \nColumns (left to right) are: ACC, PFC, and HPC.",
  "graphics/pca/pca_Dx_legend.txt"
)


# Plot PCA by Dx (label samples) -------------------------------
### Plot in loop
pca_plot_lab <- list()
for (i in names(dgeList)) {
  # Subset
  y <- dgeList[[i]]
  
  # Run PCA
  pca_obj <- prcomp(t(cpm(y, log = TRUE))) 
  
  # Plot
  pca_plot_lab[[i]] <- plot_pca(
    pca_obj,
    groups = y$samples$Dx, 
    labels = y$samples$Sample,
    labels.size = 1
  ) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 8, face = "plain"),
      axis.text = element_text(size = 5.5),
      legend.position = "none"
    ) + 
    scale_color_manual(values = cols)
}

### Combine into one figure
pca_plot_lab_all <- pca_plot_lab$ACC_genes + 
  pca_plot_lab$PFC_genes +
  pca_plot_lab$HPC_genes + plot_layout(ncol = 3) +
  pca_plot_lab$ACC_exons + 
  pca_plot_lab$PFC_exons +
  pca_plot_lab$HPC_exons + plot_layout(ncol = 3) +
  pca_plot_lab$ACC_jxns + 
  pca_plot_lab$PFC_jxns +
  pca_plot_lab$HPC_jxns + plot_layout(ncol = 3) +
  pca_plot_lab$ACC_regions + 
  pca_plot_lab$PFC_regions +
  pca_plot_lab$HPC_regions + plot_layout(ncol = 3) 
print(pca_plot_lab_all)
ggsave("graphics/pca/pca_dx_labelled.png", height = 8, width = 6)
ggsave("graphics/pca/pca_dx_labelled.eps", height = 8, width = 6)


# Remove outliers ----------------------------------------------
### List outliers
beepr::beep()
outliers <- list(
  "ACC" = NA,
  "PFC" = "S03",      
  "HPC" = "S20"
)
### Log
outliers %>%
  bind_rows(.id = "brain_region") %>%
  write_csv("logs/pca/outliers_removed.txt")

### To avoid memory limit error, do this:
#cd ~
#touch .Renviron
#open .Renviron
#R_MAX_VSIZE=200Gb # add to first line of .Renviron
### Then restart RStudio
### Note: This limit includes both physical and virtual memory; so setting _MAX_VSIZE=16Gb on a machine with 16Gb of physical memory may not prevent this error. You may have to play with this parameter, depending on the specs of your machine

### Remove samples from DGE
dgeList_filt <- list()
for (i in datasets) for (j in feature_levels) {
  y <- dgeList[[paste(i, j, sep = "_")]]
  print(paste(i, j, paste(outliers[[i]], collapse = "."), sep = "_"))
  print(nrow(y$samples))
  print(paste(rownames(y$samples), collapse = ","))
  if (is.na(outliers[[i]])) {
    dgeList_filt[[paste(i, j, sep = "_")]] <- y
  }
  if (!is.na(outliers[[i]])) {
    y <- y[ , !(rownames(y$samples) %in% outliers[[i]]), keep.lib.sizes = TRUE]
    dgeList_filt[[paste(i, j, sep = "_")]] <- y
  }
  print(nrow(y$samples))
  print(paste(rownames(y$samples), collapse = ","))
}


# Save ---------------------------------------------------------
saveRDS(dgeList_filt, "input/dge/ACC_PFC_HPC_filtered_outliersRemoved.rds"); beepr::beep()


# Write out packages and versions ------------------------------
writeLines(capture.output(sessioninfo::session_info()), paste0("logs/versions/sessionInfo_", analysis_name, ".txt"))
