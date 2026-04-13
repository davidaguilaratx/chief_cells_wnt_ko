library(ggplot2)
library(cowplot)
library(patchwork)
library(ggthemes)
library(ggpubr)
library(viridis)
library(scales)
library(Seurat)
library(scCustomize)
library(scuttle)
library(scattermore)

# ============================================================
# Panel-saving helper ----
# ============================================================

#' Save a composed (patchwork/cowplot) plot to disk.
save_panel <- function(
      panel,
      fig_dir,
      name,
      width = 16,
      height = 12,
      dpi = 320
) {
      ggsave(
            file.path(fig_dir, name),
            plot = panel,
            width = width,
            height = height,
            dpi = dpi
      )
      message("[save_panel] Saved: ", file.path(fig_dir, name))
}


# ============================================================
# Calculate MAD thresholds ----
# ============================================================

# calculate MAD thresholds for a given metric, grouped by sample
get_mad_thresholds <- function(
      metadata,
      metric_col,
      sample_col = "sample",
      apply_log1p = FALSE,
      back_transform = FALSE
) {
      vals <- metadata[[metric_col]]
      batches <- metadata[[sample_col]]

      if (apply_log1p) {
            vals <- log1p(vals)
      }

      thresh3 <- attr(isOutlier(vals, batch = batches, nmads = 3), "thresholds")
      thresh5 <- attr(isOutlier(vals, batch = batches, nmads = 5), "thresholds")

      res <- data.frame(
            sample = colnames(thresh3),
            mad3_lower = thresh3["lower", ],
            mad3_upper = thresh3["higher", ],
            mad5_lower = thresh5["lower", ],
            mad5_upper = thresh5["higher", ]
      )

      # Exponentiate the thresholds back to original scale if requested
      if (back_transform) {
            res[, -1] <- expm1(res[, -1])
      }

      return(res)
}

# ============================================================
# MAD plotting function ----
# ============================================================

# plot a metric with MAD thresholds and optional hard cutoffs, with log10 scaling if desired
plot_metric_with_mad <- function(
      metadata,
      metric_col,
      mad_thresholds,
      title,
      use_log10_scale = FALSE,
      cutoff_low = NULL,
      cutoff_high = NULL,
      apply_log1p = FALSE
) {
      plot_data <- metadata[, c("sample", metric_col)]

      # Handle scale transformation and axis labeling
      if (apply_log1p) {
            plot_data[[metric_col]] <- log1p(plot_data[[metric_col]])
            if (!is.null(cutoff_low)) {
                  cutoff_low <- log1p(cutoff_low)
            }
            if (!is.null(cutoff_high)) {
                  cutoff_high <- log1p(cutoff_high)
            }
            y_label <- paste0("log1p(", metric_col, ")")
      } else {
            y_label <- metric_col
      }

      p <- ggplot(
            plot_data,
            aes(x = sample, y = !!sym(metric_col), fill = sample)
      ) +
            geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) +
            theme_classic2() +
            theme(
                  legend.position = "none",
                  axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(title = title, x = "Sample", y = y_label)

      # Only use log10 scale if the data isn't log1p transformed
      if (use_log10_scale && !apply_log1p) {
            p <- p +
                  scale_y_log10(labels = scales::comma_format()) +
                  annotation_logticks(sides = "l")
      }

      # Add MAD lines (both upper and lower bounds)
      thresh_long <- tidyr::pivot_longer(
            mad_thresholds,
            -sample,
            names_to = "type",
            values_to = "y"
      )
      thresh_long$color <- ifelse(
            grepl("mad3", thresh_long$type),
            "cyan",
            "magenta"
      )

      p <- p +
            geom_segment(
                  data = thresh_long,
                  mapping = aes(
                        x = as.numeric(factor(sample)) - 0.4,
                        xend = as.numeric(factor(sample)) + 0.4,
                        y = y,
                        yend = y,
                        color = color
                  ),
                  inherit.aes = FALSE,
                  linewidth = 1.0,
                  linetype = "dashed"
            ) +
            scale_color_identity()

      # Add Hard Cutoffs (Low & High)
      if (!is.null(cutoff_low)) {
            p <- p +
                  geom_hline(
                        yintercept = cutoff_low,
                        color = "red",
                        linetype = "dashed",
                        linewidth = 1.0
                  )
      }
      if (!is.null(cutoff_high)) {
            p <- p +
                  geom_hline(
                        yintercept = cutoff_high,
                        color = "red",
                        linetype = "dashed",
                        linewidth = 1.0
                  )
      }

      return(p)
}

# ============================================================
# Gene Filtering Impact plot ----
# ============================================================

# plot gene filtering impact with eCDF and zoomed histogram
plot_gene_filtering_impact <- function(seurat_obj, cutoffs = c(3, 20)) {
      # Extract counts matrix
      counts <- LayerData(seurat_obj, assay = 'RNA', layer = 'counts')

      # @x contains the non-zero values. We convert them all
      # to 1 and sum the rows directly.
      counts@x <- rep(1, length(counts@x))
      cells_per_gene <- Matrix::rowSums(counts)

      df <- data.frame(gene = names(cells_per_gene), num_cells = cells_per_gene)

      # Calculate summary stats for the subtitles
      stats <- sapply(cutoffs, function(c) {
            kept <- sum(df$num_cells >= c)
            dropped <- sum(df$num_cells < c)
            sprintf(
                  "Cutoff ≥ %d: Keep %d genes | Drop %d genes",
                  c,
                  kept,
                  dropped
            )
      })

      # eCDF Plot
      p1 <- ggplot(df, aes(x = num_cells)) +
            stat_ecdf(geom = "step", linewidth = 1) +
            scale_x_log10(labels = scales::comma_format()) +
            annotation_logticks(sides = "b") +
            geom_vline(
                  xintercept = cutoffs,
                  color = c("blue", "red"),
                  linetype = "dashed",
                  linewidth = 1
            ) +
            theme_classic() +
            labs(
                  title = "Cumulative Distribution of Gene Detection",
                  subtitle = "Shows the proportion of the transcriptome dropped at any cutoff",
                  x = "Number of cells expressing the gene (Log10 scale)",
                  y = "Cumulative Fraction of Genes"
            ) +
            annotate(
                  "text",
                  x = cutoffs,
                  y = 0.1,
                  label = paste("Cutoff =", cutoffs),
                  color = c("blue", "red"),
                  angle = 90,
                  vjust = -1
            )

      # Zoomed Histogram
      p2 <- ggplot(df[df$num_cells < 100, ], aes(x = num_cells)) +
            geom_histogram(
                  binwidth = 1,
                  fill = "grey70",
                  color = "black",
                  linewidth = 0.2
            ) +
            geom_vline(
                  xintercept = cutoffs,
                  color = c("blue", "red"),
                  linetype = "dashed",
                  linewidth = 1
            ) +
            theme_classic() +
            labs(
                  title = "Zoomed Histogram: Rare Genes (Expressed in < 100 cells)",
                  subtitle = paste(stats, collapse = "\n"),
                  x = "Number of cells expressing the gene",
                  y = "Number of Genes"
            )

      return(p1 / p2) # Return stitched plot
}

# ============================================================
# QC plots ----
# ============================================================

#' Generate QC plots for a Seurat object
#'
#' `plot_qcs()` creates a standard set of single-cell QC plots from a Seurat
#' object and saves them as PNG files to a directory.
#'
#' The function uses `seurat_object@meta.data` and assumes the metadata contains
#' columns named `sample`, `condition`, `percent.mt`, `complexity`,
#' `nCount_RNA`, and `nFeature_RNA`.
#'
#' @param seurat_object A Seurat object containing QC metadata.
#' @param percent.mt.cutoff Numeric cutoff used for mitochondrial percentage
#'   annotations and color scaling. Default `8`.
#' @param gene_cutoff_low Numeric lower gene-count cutoff. Default `100`.
#' @param gene_cutoff_high Numeric upper gene-count cutoff. Default `6000`.
#' @param umi_cutoff_low Numeric lower UMI-count cutoff. Default `100`.
#' @param umi_cutoff_high Numeric upper UMI-count cutoff. Default `10000`.
#' @param complexity_cutoff Numeric novelty score cutoff used for the gene/UMI
#'   complexity density plot. Default `0.8`.
#' @param min_cells_per_gene Minimum cells per gene threshold. This parameter is
#'   accepted for compatibility but is not used within the function.
#' @param fig_dir Character path to the folder where PNG files will be written.
#'
#' @return Invisibly returns `NULL`. The primary effect is to save QC plot
#'   images to `fig_dir`.
#'
#' @details
#' The function generates:
#' - bar plot of cell counts per sample
#' - violin and density plots for mitochondrial percentage
#' - density plot of complexity scores
#' - violin and density plots for UMI counts
#' - violin and density plots for gene counts
#' - scatter plots of `nCount_RNA` versus `nFeature_RNA` colored by sample,
#'   condition, and mitochondrial percentage
#' - Seurat QC plots by condition for genes, UMIs, and mitochondrial ratio
#'
#' It relies on Seurat plotting helpers like `VlnPlot()`, `FeatureScatter()`,
#' and custom helper functions such as `QC_Plots_Genes()`, `QC_Plots_UMIs()`,
#' and `QC_Plots_Mito()`.
#'
#' @examples
#' \dontrun{
#' fig_dir <- "path/to/qc/figures"
#' plot_qcs(merged_seurat,
#'          percent.mt.cutoff = 8,
#'          gene_cutoff_low = 100,
#'          gene_cutoff_high = 7500,
#'          umi_cutoff_low = 100,
#'          umi_cutoff_high = 50000,
#'          complexity_cutoff = 0.8,
#'          fig_dir = fig_dir)
#' }
#'
#' @export
plot_qcs = function(
      seurat_object,
      percent.mt.cutoff = 8,
      gene_cutoff_low = 100,
      gene_cutoff_high = 7500,
      umi_cutoff_low = 100,
      umi_cutoff_high = 50000,
      complexity_cutoff = 0.8,
      min_cells_per_gene = 20,
      fig_dir = NULL
) {
      meta = seurat_object@meta.data

      ## plot number of cells per dataset ----
      cat('\nPlotting number of cells per sample')
      cond_counts <- table(meta$condition)
      subtitle <- paste0(
            "Total cells = ",
            nrow(meta),
            ", ",
            paste0(
                  names(cond_counts),
                  " = ",
                  as.integer(cond_counts),
                  collapse = ", "
            )
      )

      ggplot(meta, aes(x = sample, fill = sample)) +
            geom_bar(position = 'stack') +
            theme_classic2() +
            NoLegend() +
            labs(
                  title = "Number of cells per experiment",
                  subtitle = subtitle
            ) +
            theme(
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            scale_y_continuous(n.breaks = 8)
      ggsave(paste0(fig_dir, '/ncells.png'), dpi = 320, width = 10, height = 10)

      ## plot histograms of UMIs and genes ----
      cat('\nPlotting histogram of UMIs and genes per cell')

      p_hist_umi <- QC_Histogram(
            seurat_object = seurat_object,
            features = "nCount_RNA",
            low_cutoff = umi_cutoff_low,
            high_cutoff = umi_cutoff_high
      ) +
            theme_classic2() +
            scale_x_log10(labels = comma_format()) +
            annotation_logticks(sides = "b") +
            labs(
                  title = "Histogram of UMIs per cell",
                  x = "nCount_RNA (log10)",
                  y = "Number of Cells"
            ) +
            theme(
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12)
            )

      p_hist_gene <- QC_Histogram(
            seurat_object = seurat_object,
            features = "nFeature_RNA",
            low_cutoff = gene_cutoff_low,
            high_cutoff = gene_cutoff_high
      ) +
            theme_classic2() +
            scale_x_log10(labels = comma_format()) +
            annotation_logticks(sides = "b") +
            labs(
                  title = "Histogram of Genes per cell",
                  x = "nFeature_RNA (log10)",
                  y = "Number of Cells"
            ) +
            theme(
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12)
            )

      p_hist_mt <- QC_Histogram(
            seurat_object = seurat_object,
            features = "percent.mt",
            high_cutoff = percent.mt.cutoff
      ) +
            theme_classic2() +
            annotation_logticks(sides = "b") +
            labs(
                  title = "Histogram of Mitochondrial Genes per cell",
                  x = "percent.mt",
                  y = "Number of Cells"
            ) +
            theme(
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12)
            )

      # Combine side-by-side using patchwork
      p_hists <- p_hist_umi | p_hist_gene | p_hist_mt

      ggsave(
            file.path(fig_dir, 'histogram_ncounts_ngenes_percent_mt.png'),
            plot = p_hists,
            dpi = 320,
            width = 21,
            height = 10
      )

      rm(p_hist_umi, p_hist_gene, p_hist_mt, p_hists)
      gc()

      ## plot proportion of mitochondrial genes per sample ----
      cat('\nplotting percent.mt per sample')
      VlnPlot(
            seurat_object,
            features = "percent.mt",
            group.by = 'sample',
            pt.size = 0.05,
            alpha = 0.3
      ) +
            NoLegend() +
            labs(
                  title = 'Percent mitochondrial genes',
                  subtitle = paste0('cutoff above ', percent.mt.cutoff, '%'),
                  x = 'sample',
                  y = 'percent.mt'
            ) +
            scale_y_continuous(n.breaks = 20, limits = c(0, NA)) +
            theme(
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            geom_hline(
                  yintercept = percent.mt.cutoff,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            annotate(
                  'text',
                  x = 7.5,
                  y = 10.5,
                  label = 'cutoff',
                  size = 5,
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/percent_mt.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot novelty score (complexity = log10genes / log10UMI) ----
      cat('\nPlotting complexity per sample')
      ggplot(meta, aes(color = sample, x = complexity, fill = sample)) +
            geom_density(alpha = 0.1) +
            theme_classic2() +
            labs(
                  title = 'Overal complexity via genes detected per UMI (novelty score)',
                  x = 'log10(nGenes) / log10(nUMI)',
                  y = 'Cell density'
            ) +
            theme(
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            geom_vline(
                  xintercept = complexity_cutoff,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            annotate(
                  'text',
                  x = complexity_cutoff - 0.02,
                  y = 40,
                  label = 'cutoff',
                  size = 5,
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/complexity_novelty_scores.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot number of reads per cell as vln plots ----
      cat('\nPlotting UMIs per cell per sample as vln plots')
      VlnPlot(
            seurat_object,
            features = "nCount_RNA",
            group.by = 'sample',
            pt.size = 0.05,
            alpha = 0.3
      ) +
            theme_classic2() +
            NoLegend() +
            labs(
                  title = 'Number of transcript UMIs per cell',
                  subtitle = paste0(
                        'cutoffs below ',
                        umi_cutoff_low,
                        ' and above ',
                        umi_cutoff_high,
                        ' UMIs'
                  ),
                  x = 'sample',
                  y = 'nUMI'
            ) +
            scale_y_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "l") +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            geom_hline(
                  yintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/umi_counts.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot number of reads per cell as density plots ----
      cat('\nPlotting UMIs per cell per sample as density plots')
      ggplot(meta, aes(color = sample, x = nCount_RNA, fill = sample)) +
            geom_density(alpha = 0.1) +
            labs(
                  title = 'Number of transcript UMIs per cell',
                  subtitle = paste0(
                        'cutoffs below ',
                        umi_cutoff_low,
                        ' and above ',
                        umi_cutoff_high,
                        ' UMIs'
                  ),
                  x = 'nUMI',
                  y = 'Cell density'
            ) +
            theme_classic2() +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            scale_x_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "b") +
            geom_vline(
                  xintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/umi_counts_density.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot number of genes per cell as vln plots ----
      cat('\nPlotting genes per cell per sample as vln plots')
      VlnPlot(
            seurat_object,
            features = "nFeature_RNA",
            group.by = 'sample',
            pt.size = 0.05,
            alpha = 0.3
      ) +
            theme_classic2() +
            NoLegend() +
            labs(
                  title = 'Number of genes detected per cell',
                  subtitle = paste0(
                        'cutoffs below ',
                        gene_cutoff_low,
                        ' and above ',
                        gene_cutoff_high,
                        ' genes'
                  ),
                  x = 'sample',
                  y = 'nGenes'
            ) +
            scale_y_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "l") +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            geom_hline(
                  yintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/gene_counts.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot number of genes per cell as density plots ----
      cat('\nPlotting genes per cell per sample as density plots')
      ggplot(meta, aes(color = sample, x = nFeature_RNA, fill = sample)) +
            geom_density(alpha = 0.1) +
            labs(
                  title = 'Number of genes detected per cell',
                  subtitle = paste0(
                        'cutoffs below ',
                        gene_cutoff_low,
                        ' and above ',
                        gene_cutoff_high,
                        ' genes'
                  ),
                  x = 'nGene',
                  y = 'Cell density'
            ) +
            theme_classic2() +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18)
            ) +
            scale_x_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "b") +
            geom_vline(
                  xintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )
      ggsave(
            paste0(fig_dir, '/gene_counts_density.png'),
            dpi = 320,
            width = 10,
            height = 10
      )

      ## plot count depth (number of UMIs) vs number of genes, colored by sample name (left) and percent mitochondrial gene counts (right) ----
      p1 = FeatureScatter(
            seurat_object,
            feature1 = 'nCount_RNA',
            feature2 = 'nFeature_RNA',
            group.by = 'sample',
            pt.size = 0.5,
            shuffle = TRUE
      ) +
            DarkTheme() +
            labs(
                  title = '# of UMIs vs # of genes',
                  subtitle = 'colored by sample',
                  x = '# of UMI counts (log scale)',
                  y = '# of genes (log scale)',
                  color = 'sample'
            ) +
            scale_x_log10(n.breaks = 8, labels = comma_format()) +
            scale_y_log10(
                  n.breaks = 8,
                  limits = c(NA, 10000),
                  labels = comma_format()
            ) +
            annotation_logticks(sides = "bl", colour = 'white', size = 1) +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18),
                  plot.background = element_rect(color = 'black')
            ) +
            NoGrid() +
            stat_smooth(method = lm, linewidth = 0.5, colour = '#61e4e8') +
            geom_hline(
                  yintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )

      p2 = ggplot(
            meta,
            aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)
      ) +
            geom_scattermore(size = 1.5) +
            DarkTheme() +
            labs(
                  title = '# of UMIs vs # of genes',
                  subtitle = paste0(
                        'colored by percent mitochondrial (mt) counts >',
                        percent.mt.cutoff
                  ),
                  x = '# of UMI counts (log scale)',
                  y = '# of genes (log scale)',
                  color = '% mt counts'
            ) +
            scale_x_log10(n.breaks = 8, labels = comma_format()) +
            scale_y_log10(
                  n.breaks = 8,
                  limits = c(NA, 10000),
                  labels = comma_format()
            ) +
            annotation_logticks(sides = "bl", colour = 'white', size = 1) +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18),
                  plot.background = element_rect(color = 'black')
            ) +
            NoGrid() +
            stat_smooth(method = lm, linewidth = 0.5, colour = '#61e4e8') +
            geom_hline(
                  yintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            scale_color_gradientn(
                  colors = viridis_plasma_light_high,
                  limits = c(percent.mt.cutoff, NA),
                  na.value = '#292a2d'
            )

      cat(
            '\nPlotting UMIs vs genes colored by sample (left) and mt.percent (right)'
      )
      plot_grid(p1, p2)
      ggsave(
            paste0(fig_dir, '/ncounts_vs_ngenes_coloredby_percentmt.png'),
            dpi = 320,
            width = 20,
            height = 10
      )

      ## plot count depth (number of reads) vs number of genes, colored by exposure condition and mitochondrial gene percentage ----
      p1 = FeatureScatter(
            seurat_object,
            feature1 = 'nCount_RNA',
            feature2 = 'nFeature_RNA',
            group.by = 'condition',
            pt.size = 0.5,
            shuffle = TRUE
      ) +
            DarkTheme() +
            labs(
                  title = '# of UMIs vs # of genes',
                  subtitle = 'colored by condition',
                  x = '# of UMI counts (log scale)',
                  y = '# of genes (log scale)',
                  color = 'condition'
            ) +
            scale_x_log10(n.breaks = 8, labels = comma_format()) +
            scale_y_log10(
                  n.breaks = 8,
                  limits = c(NA, 10000),
                  labels = comma_format()
            ) +
            annotation_logticks(sides = "bl", colour = 'white', size = 1) +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18),
                  plot.background = element_rect(color = 'black')
            ) +
            NoGrid() +
            stat_smooth(method = lm, linewidth = 0.5, colour = '#61e4e8') +
            geom_hline(
                  yintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            )

      p2 = ggplot(
            meta,
            aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)
      ) +
            geom_scattermore(size = 1.5) +
            DarkTheme() +
            labs(
                  title = '# of UMIs vs # of genes',
                  subtitle = paste0(
                        'colored by percent mitochondrial (mt) counts >',
                        percent.mt.cutoff
                  ),
                  x = '# of UMI counts (log scale)',
                  y = '# of genes (log scale)',
                  color = '% mt counts'
            ) +
            scale_x_log10(n.breaks = 8, labels = comma_format()) +
            scale_y_log10(
                  n.breaks = 8,
                  limits = c(NA, 10000),
                  labels = comma_format()
            ) +
            annotation_logticks(sides = "bl", colour = 'white', size = 1) +
            theme(
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.title = element_text(size = 12),
                  axis.text = element_text(size = 12),
                  title = element_text(size = 18),
                  plot.background = element_rect(color = 'black')
            ) +
            NoGrid() +
            stat_smooth(method = lm, linewidth = 0.5, colour = '#61e4e8') +
            geom_hline(
                  yintercept = gene_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_hline(
                  yintercept = gene_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_low,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            geom_vline(
                  xintercept = umi_cutoff_high,
                  linetype = 'dashed',
                  color = 'red'
            ) +
            scale_color_gradientn(
                  colors = viridis_plasma_light_high,
                  limits = c(percent.mt.cutoff, NA),
                  na.value = '#292a2d'
            )

      cat(
            '\nPlotting UMIs vs genes colored by condition (left) and mt.percent (right)'
      )
      plot_grid(p1, p2)
      ggsave(
            paste0(
                  fig_dir,
                  '/ncounts_vs_ngenes_coloredby_percentmt_and_condition.png'
            ),
            dpi = 320,
            width = 20,
            height = 10
      )

      ## more vln plots for UMIs, genes, percent.mt, and more by exposure condition ----
      p1 <- QC_Plots_UMIs(
            seurat_object,
            group.by = 'condition',
            low_cutoff = umi_cutoff_low,
            high_cutoff = umi_cutoff_high,
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = nCount_RNA),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            ) +
            scale_y_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "l", size = .5)

      p2 <- QC_Plots_Genes(
            seurat_object,
            group.by = 'condition',
            low_cutoff = gene_cutoff_low,
            high_cutoff = gene_cutoff_high,
            y_axis_log = TRUE,
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = nFeature_RNA),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            ) +
            scale_y_log10(n.breaks = 10, labels = comma_format()) +
            annotation_logticks(sides = "l", size = .5)

      p3 <- QC_Plots_Mito(
            seurat_object,
            mito_name = 'percent.mt',
            group.by = 'condition',
            low_cutoff = percent.mt.cutoff,
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = percent.mt),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            )

      p4 <- QC_Plots_Complexity(
            seurat_object,
            feature = 'complexity',
            group.by = 'condition',
            low_cutoff = complexity_cutoff,
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = complexity),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            )

      p5 <- QC_Plots_Feature(
            seurat_object,
            feature = 'pct_counts_in_top_20_genes',
            group.by = 'condition',
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = pct_counts_in_top_20_genes),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            )

      p6 <- QC_Plots_Feature(
            seurat_object,
            feature = 'percent.hb',
            group.by = 'condition',
            pt.size = 0
      ) +
            geom_scattermore(
                  aes(x = ident, y = percent.hb),
                  position = position_jitter(width = 0.2),
                  pointsize = 1,
                  alpha = 0.3
            )

      cat('\nPlotting nUMIs, nGenes, and percent.mt per condition')
      combined_p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, align = 'hv')
      ggsave(
            paste0(fig_dir, '/vln_metrics_by_condition.png'),
            dpi = 320,
            width = 20,
            height = 10
      )
      rm(p1, p2, p3, p4, p5, p6, combined_p)
      gc()

      ## plot Malat1 diagnostics (if MALAT1 is in data) ----
      # https://doi.org/10.1101/2024.07.14.603469
      if (
            length(grep(
                  "^Malat1",
                  rownames(seurat_object),
                  ignore.case = TRUE,
                  value = TRUE
            )) >
                  0
      ) {
            cat('\nPlotting Malat1_Threshold diagnostics')

            # Ensure it's treated as a discrete factor/character
            meta$Malat1_Threshold <- as.character(meta$Malat1_Threshold)

            # 1. UMIs vs Genes colored by Malat1 threshold
            p_malat1_scatter <- ggplot(
                  meta,
                  aes(
                        x = nCount_RNA,
                        y = nFeature_RNA,
                        color = Malat1_Threshold
                  )
            ) +
                  geom_scattermore(size = 1.5) +
                  DarkTheme() +
                  scale_x_log10(n.breaks = 8, labels = comma_format()) +
                  scale_y_log10(
                        n.breaks = 8,
                        limits = c(NA, 10000),
                        labels = comma_format()
                  ) +
                  annotation_logticks(
                        sides = "bl",
                        colour = 'white',
                        size = 1
                  ) +
                  # Explicitly map TRUE to dark grey and FALSE to a bright cyan/red
                  scale_color_manual(
                        values = c("TRUE" = "#292a2d", "FALSE" = "cyan")
                  ) +
                  theme(
                        axis.text.x = element_text(
                              angle = 45,
                              vjust = 1,
                              hjust = 1
                        ),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        plot.subtitle = element_text(hjust = 0.5),
                        plot.background = element_rect(color = 'black')
                  ) +
                  NoGrid() +
                  geom_hline(
                        yintercept = gene_cutoff_low,
                        linetype = 'dashed',
                        color = 'red'
                  ) +
                  geom_hline(
                        yintercept = gene_cutoff_high,
                        linetype = 'dashed',
                        color = 'red'
                  ) +
                  geom_vline(
                        xintercept = umi_cutoff_low,
                        linetype = 'dashed',
                        color = 'red'
                  ) +
                  geom_vline(
                        xintercept = umi_cutoff_high,
                        linetype = 'dashed',
                        color = 'red'
                  ) +
                  labs(
                        title = '# of UMIs vs # of genes',
                        subtitle = 'Colored by Malat1_Threshold',
                        x = '# of UMI counts (log scale)',
                        y = '# of genes (log scale)',
                        color = 'Malat1 Threshold'
                  )

            # 2. Percent MT vs Malat1 separated by sample (using jitter since Y is categorical)
            p_mt_vs_malat1 <- ggplot(
                  meta,
                  aes(x = percent.mt, y = Malat1_Threshold, color = condition)
            ) +
                  geom_jitter(size = 0.5, alpha = 0.6, height = 0.2) +
                  facet_wrap(~sample) +
                  theme_classic2() +
                  DarkTheme() +
                  theme(
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        strip.text = element_text(face = "bold")
                  ) +
                  geom_vline(
                        xintercept = percent.mt.cutoff,
                        linetype = "dashed",
                        color = "red"
                  ) +
                  labs(
                        title = "Percent MT vs Malat1 Threshold",
                        x = "Percent Mitochondrial Counts",
                        y = "Malat1 Threshold",
                        color = "Condition"
                  )

            # Combine and save
            malat1_combined <- p_malat1_scatter | p_mt_vs_malat1

            ggsave(
                  file.path(fig_dir, 'malat1_diagnostics.png'),
                  plot = malat1_combined,
                  dpi = 320,
                  width = 20,
                  height = 8
            )

            rm(p_malat1_scatter, p_mt_vs_malat1, malat1_combined)
            gc()
      }

      ## plot metrics with MAD cutoffs per sample ----
      metrics <- list(
            percent.mt = list(
                  log10 = FALSE,
                  low = NULL,
                  high = percent.mt.cutoff,
                  title = 'Percent MT'
            ),
            nCount_RNA = list(
                  log10 = TRUE,
                  low = umi_cutoff_low,
                  high = umi_cutoff_high,
                  title = 'UMIs per cell'
            ),
            nFeature_RNA = list(
                  log10 = TRUE,
                  low = gene_cutoff_low,
                  high = gene_cutoff_high,
                  title = 'Genes per cell'
            ),
            pct_counts_in_top_20_genes = list(
                  log10 = FALSE,
                  low = NULL,
                  high = NULL,
                  title = 'Percent Counts in Top 20 Genes'
            ),
            pct_counts_in_top_50_genes = list(
                  log10 = FALSE,
                  low = NULL,
                  high = NULL,
                  title = 'Percent Counts in Top 50 Genes'
            )
      )

      for (m in names(metrics)) {
            cat(sprintf('\nPlotting %s with MAD variations', m))
            c_low <- metrics[[m]]$low
            c_high <- metrics[[m]]$high
            t <- metrics[[m]]$title
            l <- metrics[[m]]$log10

            # Plot 1: Raw Data + Raw MADs
            thresh_raw <- get_mad_thresholds(
                  meta,
                  m,
                  apply_log1p = FALSE,
                  back_transform = FALSE
            )
            p1 <- plot_metric_with_mad(
                  meta,
                  m,
                  thresh_raw,
                  paste(t, "\nRaw MADs"),
                  l,
                  c_low,
                  c_high,
                  apply_log1p = FALSE
            )

            # Plot 2: Raw Data + Back-transformed Log1p MADs (expm1)
            thresh_back <- get_mad_thresholds(
                  meta,
                  m,
                  apply_log1p = TRUE,
                  back_transform = TRUE
            )
            p2 <- plot_metric_with_mad(
                  meta,
                  m,
                  thresh_back,
                  paste(t, "\nExpm1 MADs (on Raw Data)"),
                  l,
                  c_low,
                  c_high,
                  apply_log1p = FALSE
            )

            # Stitch side-by-side using patchwork
            combined_plot <- p1 | p2

            ggsave(
                  file.path(fig_dir, sprintf('%s_MAD_comparisons.png', m)),
                  plot = combined_plot,
                  dpi = 320,
                  width = 16,
                  height = 8
            )

            rm(p1, p2, combined_plot)
            gc()
      }

      ## plot min cells per gene cutoff impact ----
      cat('\nPlotting gene level filtering impact (min_cells_per_gene)\n')

      # Compare baseline standard (3 cells) vs proposed cutoff (20 cells)
      # or any value dynamically provided by min_cells_per_gene
      cutoffs_to_test <- unique(c(3, min_cells_per_gene))

      p_genes <- plot_gene_filtering_impact(
            seurat_object,
            cutoffs = cutoffs_to_test
      )

      ggsave(
            file.path(fig_dir, 'gene_filtering_impact_distributions.png'),
            plot = p_genes,
            dpi = 320,
            width = 10,
            height = 12
      )

      rm(p_genes)
      gc()

      ## plot percent counts in top n genes ----
      cat('\nPlotting percent of counts in top genes per sample')

      top_metrics <- c(
            "pct_counts_in_top_20_genes",
            "pct_counts_in_top_50_genes",
            "pct_counts_in_top_100_genes",
            "pct_counts_in_top_200_genes",
            "pct_counts_in_top_500_genes"
      )

      top_metrics <- intersect(top_metrics, colnames(meta))

      if (length(top_metrics) > 0) {
            # Reshape data to long format
            meta_long <- tidyr::pivot_longer(
                  meta,
                  cols = all_of(top_metrics),
                  names_to = "metric",
                  values_to = "pct_counts"
            )

            # Clean up metric names for the legend
            meta_long$metric <- factor(meta_long$metric, levels = top_metrics)
            levels(meta_long$metric) <- gsub(
                  "pct_counts_in_top_|_genes",
                  "",
                  levels(meta_long$metric)
            )
            levels(meta_long$metric) <- paste("Top", levels(meta_long$metric))

            # Combined Violin Plot
            p_vln_combined <- ggplot(
                  meta_long,
                  aes(x = sample, y = pct_counts, fill = sample)
            ) +
                  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) +
                  facet_wrap(~metric, ncol = 3) +
                  theme_classic2() +
                  theme(
                        legend.position = "none",
                        axis.text.x = element_text(
                              angle = 45,
                              hjust = 1,
                              size = 10
                        ),
                        strip.text = element_text(face = "bold", size = 12)
                  ) +
                  labs(
                        title = "Percent Counts in Top Genes",
                        x = NULL,
                        y = "Percent of Total Counts"
                  )

            ggsave(
                  file.path(
                        fig_dir,
                        'pct_counts_in_top_genes_vln_combined.png'
                  ),
                  plot = p_vln_combined,
                  dpi = 320,
                  width = 18,
                  height = 10
            )

            # Overlaid Scatter Plot vs Total UMIs
            p_scatter_combined <- ggplot(
                  meta_long,
                  aes(x = nCount_RNA, y = pct_counts, color = metric)
            ) +
                  geom_scattermore(size = 1.5, alpha = 0.5) +
                  facet_wrap(~sample) + # Facet only by sample
                  scale_x_log10(labels = comma_format()) +
                  scale_color_viridis_d(option = "plasma", end = 0.9) + # Distinct colors for the bands
                  DarkTheme() +
                  theme(
                        legend.position = "bottom",
                        legend.title = element_text(face = "bold"),
                        axis.text.x = element_text(
                              angle = 45,
                              hjust = 1,
                              vjust = 1,
                              size = 10
                        ),
                        strip.text = element_text(face = "bold", size = 12),
                        plot.background = element_rect(color = 'black')
                  ) +
                  guides(
                        color = guide_legend(
                              override.aes = list(size = 4, alpha = 1)
                        )
                  ) +
                  labs(
                        title = "Library Complexity by Sample",
                        subtitle = "Saturation bands indicate the percentage of reads consumed by the top N genes",
                        x = "nCount_RNA (log10)",
                        y = "Percent of Total Counts",
                        color = "Threshold"
                  )

            ggsave(
                  file.path(
                        fig_dir,
                        'pct_counts_in_top_genes_scatter_combined.png'
                  ),
                  plot = p_scatter_combined,
                  dpi = 320,
                  width = 12,
                  height = 10
            )

            rm(meta_long, p_vln_combined, p_scatter_combined)
            gc()
      }
}

# ============================================================
# Outlier diagnostic plots ----
# ============================================================

#' Plot outlier-flagged cells before filtering
#'
#' Generates three diagnostic plots after outlier columns have been added to
#' the Seurat object metadata (e.g. via `scuttle::isOutlier`), but before any
#' cells are removed. Saves each plot as a PNG to `fig_dir`.
#'
#' @param seurat_object A Seurat object with columns `outlier`, `mt_outlier`,
#'   `nCount_RNA`, `nFeature_RNA`, `percent.mt`, and `sample` in `@meta.data`.
#' @param mt_hard_cutoff Numeric hard mitochondrial percentage cutoff drawn as
#'   an orange dashed line. Default `8`.
#' @param fig_dir Character path to the folder where PNG files will be written.
#'   If `NULL`, plots are printed to the active device instead.
#'
#' @return Invisibly returns a named list of the three ggplot objects:
#'   `scatter_count_gene`, `mt_scatter`, and `removal_summary`.
#'
#' @export
plot_outliers <- function(
      seurat_object,
      mt_hard_cutoff = 8,
      fig_dir = NULL
) {
      meta <- seurat_object@meta.data

      # ----------------------------------------------------------
      ## 1. nCount vs nFeature coloured by outlier status ----
      # ----------------------------------------------------------
      p_scatter <- ggplot(
            meta,
            aes(x = nCount_RNA, y = nFeature_RNA, color = flagged)
      ) +
            geom_scattermore(pointsize = 2, alpha = 0.5) +
            scale_color_manual(
                  values = c("TRUE" = "red", "FALSE" = "grey60"),
                  name = "Flagged"
            ) +
            scale_x_log10(labels = scales::comma) +
            scale_y_log10(labels = scales::comma) +
            annotation_logticks() +
            facet_wrap(~sample) +
            theme_classic2() +
            theme(
                  strip.background = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(
                  title = "nCount vs nFeature — outliers flagged",
                  x = "Total UMIs (log10)",
                  y = "Genes detected (log10)"
            )

      # ----------------------------------------------------------
      ## 2. percent.mt vs nCount coloured by mt_outlier ----
      # ----------------------------------------------------------
      p_mt <- ggplot(
            meta,
            aes(x = nCount_RNA, y = percent.mt, color = mt_outlier)
      ) +
            geom_scattermore(pointsize = 2, alpha = 0.5) +
            scale_color_manual(
                  values = c("TRUE" = "red", "FALSE" = "grey60"),
                  name = "MT outlier"
            ) +
            scale_x_log10(labels = scales::comma) +
            annotation_logticks(sides = "b") +
            geom_hline(
                  yintercept = mt_hard_cutoff,
                  linetype = "dashed",
                  color = "orange",
                  linewidth = 0.8
            ) +
            facet_wrap(~sample) +
            theme_classic2() +
            theme(
                  strip.background = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(
                  title = "MT% vs UMIs — MT outliers flagged",
                  subtitle = paste0(
                        "Orange line = hard cutoff (",
                        mt_hard_cutoff,
                        "%)"
                  ),
                  x = "Total UMIs (log10)",
                  y = "% mitochondrial reads"
            )

      # ----------------------------------------------------------
      # Save or print ----
      # ----------------------------------------------------------
      plots <- list(
            scatter_count_gene = p_scatter,
            mt_scatter = p_mt
      )

      if (!is.null(fig_dir)) {
            ggsave(
                  file.path(fig_dir, "outlier_count_gene_scatter.png"),
                  plot = p_scatter,
                  dpi = 320,
                  width = 14,
                  height = 10
            )
            ggsave(
                  file.path(fig_dir, "outlier_mt_scatter.png"),
                  plot = p_mt,
                  dpi = 320,
                  width = 14,
                  height = 10
            )

            message("[plot_outliers] Saved 2 plots to: ", fig_dir)
      } else {
            print(p_scatter)
            print(p_mt)
      }

      invisible(plots)
}


# ============================================================
# Pre/post filtering comparison plots ----
# ============================================================

plot_before_and_after <- function(
      seurat_obj,
      filtered_seurat_obj,
      percent.mt.cutoff = 8,
      fig_dir = NULL,
      smooth_subsample = 5000L
) {
      metadata <- seurat_obj@meta.data
      filtered_metadata <- filtered_seurat_obj@meta.data

      # -- compute-once summaries --
      # condition levels pulled from the object, not hardcoded
      condition_levels <- sort(unique(as.character(c(
            metadata$condition,
            filtered_metadata$condition
      ))))
      condition_counts <- function(df) {
            tbl <- table(factor(df$condition, levels = condition_levels))
            setNames(as.integer(tbl), names(tbl))
      }

      n_total <- nrow(metadata)
      n_total_f <- nrow(filtered_metadata)
      cond_before <- condition_counts(metadata)
      cond_after <- condition_counts(filtered_metadata)
      max_y <- ceiling(
            max(table(metadata$sample), table(filtered_metadata$sample)) / 1000
      ) *
            1000

      format_subtitle <- function(total, counts) {
            parts <- c(
                  sprintf("Total cells = %d", total),
                  sprintf("# %s cells = %d", names(counts), counts)
            )
            paste(parts, collapse = ", ")
      }

      # -- helper: cells-per-sample bar plot --
      make_bar <- function(df, title, subtitle) {
            ggplot(df, aes(x = sample, fill = sample)) +
                  geom_bar() +
                  theme_classic2() +
                  NoLegend() +
                  labs(title = title, subtitle = subtitle) +
                  theme(
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        plot.subtitle = element_text(hjust = 0.5),
                        axis.title = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        title = element_text(size = 18)
                  ) +
                  scale_y_continuous(n.breaks = 10, limits = c(0, max_y))
      }

      bar_before <- make_bar(
            metadata,
            "Number of cells per sample",
            format_subtitle(n_total, cond_before)
      )
      bar_after <- make_bar(
            filtered_metadata,
            "Number of cells per sample after filtering",
            format_subtitle(n_total_f, cond_after)
      )
      bar_combined <- plot_grid(bar_before, bar_after)

      # -- helper: UMIs-vs-genes scatter --
      make_scatter <- function(df, title) {
            smooth_df <- if (nrow(df) > smooth_subsample) {
                  df[sample.int(nrow(df), smooth_subsample), ]
            } else {
                  df
            }
            ggplot(
                  df,
                  aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)
            ) +
                  geom_scattermore(pointsize = 1.2) +
                  DarkTheme() +
                  labs(
                        title = title,
                        subtitle = paste0(
                              "colored by percent mitochondrial (mt) counts >",
                              percent.mt.cutoff
                        ),
                        x = "# of UMI counts (log scale)",
                        y = "# of genes (log scale)",
                        color = "% mt counts"
                  ) +
                  scale_x_log10(n.breaks = 8, labels = comma_format()) +
                  scale_y_log10(n.breaks = 8, labels = comma_format()) +
                  coord_cartesian(ylim = c(NA, 10000)) +
                  annotation_logticks(
                        sides = "bl",
                        colour = "white",
                        size = 1
                  ) +
                  theme(
                        axis.text.x = element_text(
                              angle = 45,
                              vjust = 1,
                              hjust = 1
                        ),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        plot.subtitle = element_text(hjust = 0.5),
                        axis.title = element_text(size = 12),
                        axis.text = element_text(size = 12),
                        title = element_text(size = 18),
                        plot.background = element_rect(color = "black")
                  ) +
                  NoGrid() +
                  stat_smooth(
                        data = smooth_df,
                        method = lm,
                        linewidth = 0.5,
                        colour = "#61e4e8"
                  ) +
                  scale_color_gradientn(
                        colors = viridis_plasma_light_high,
                        limits = c(percent.mt.cutoff, NA),
                        na.value = "#4b4d53"
                  )
      }

      scatter_before <- make_scatter(metadata, "# of UMIs vs # of genes")
      scatter_after <- make_scatter(
            filtered_metadata,
            "# of UMIs vs # of genes after filtering"
      )
      scatter_combined <- plot_grid(scatter_before, scatter_after)

      if (!is.null(fig_dir)) {
            ggsave(
                  file.path(fig_dir, "ncells_before_and_after_filtering.png"),
                  plot = bar_combined,
                  dpi = 320,
                  width = 20,
                  height = 10
            )
            ggsave(
                  file.path(
                        fig_dir,
                        "ncounts_vs_ngenes_coloredby_percentmt_before_and_after_filtering.png"
                  ),
                  plot = scatter_combined,
                  dpi = 320,
                  width = 20,
                  height = 10
            )
            message(
                  "[plot_before_and_after] Saved 2 plots to: ",
                  fig_dir
            )
      }

      invisible(list(bar = bar_combined, scatter = scatter_combined))
}


# ============================================================
# Doublet QC violins ----
# ============================================================

#' QC-metric violin plots split by final doublet call.
#'
#' Expects metadata column \code{doublet_final} from
#' \code{run_doublet_detection()}. UMAP visualizations of the
#' round-1/round-2/combined calls are built directly in the qc
#' script via wrap_plots, not here.
plot_doublets <- function(seurat_obj, fig_dir = NULL) {
      metrics <- c(
            "nCount_RNA",
            "nFeature_RNA",
            "complexity",
            "percent.mt",
            "percent.hb"
      )
      metrics <- intersect(metrics, colnames(seurat_obj@meta.data))

      vln_sample <- VlnPlot(
            seurat_obj,
            group.by = "sample",
            split.by = "doublet_final",
            features = metrics,
            ncol = 3,
            pt.size = 0
      ) &
            theme(legend.position = "right")

      vln_condition <- VlnPlot(
            seurat_obj,
            group.by = "condition",
            split.by = "doublet_final",
            features = metrics,
            ncol = 3,
            pt.size = 0
      ) &
            theme(legend.position = "right")

      if (!is.null(fig_dir)) {
            ggsave(
                  file.path(fig_dir, "doublets_by_sample.png"),
                  plot = vln_sample,
                  dpi = 320,
                  width = 20,
                  height = 20
            )
            ggsave(
                  file.path(fig_dir, "doublets_by_condition.png"),
                  plot = vln_condition,
                  dpi = 320,
                  width = 20,
                  height = 20
            )
            message("[plot_doublets] Saved 2 plots to: ", fig_dir)
      }

      invisible(list(
            vln_sample = vln_sample,
            vln_condition = vln_condition
      ))
}
