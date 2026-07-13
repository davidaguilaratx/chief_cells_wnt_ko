library(Seurat)
library(sctransform)
library(tidyverse)
library(scuttle)
library(zeallot) # for unpacking multiple values from functions
library(qs2)
library(RhpcBLASctl)
library(ggraph) # required for clustree

seed = 281330800

set.seed(seed) # for reproducibility

# set up directories and parameters ----
base_dir <- "/home/aguilada/aguilada-colacino/scRNAseq/chief_cells_analysis/SingleCell2016"
data_dir <- file.path(base_dir, "data")
fig_dir <- file.path(base_dir, "clustering")
dir_umaps <- file.path(fig_dir, "umaps")
dir_markers <- file.path(dir_umaps, "cellmarkers")

dirs_to_make <- c(data_dir, fig_dir, dir_umaps, dir_markers)
walk(dirs_to_make, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# source utility functions
source(here::here(file.path(base_dir, "scripts/R/utils.R")))
source(here::here(file.path(base_dir, "scripts/R/plot_utils.R")))

# number of workers for parallel processes
# nthreads <- 2
c(nthreads, bp_param) %<-%
  configure_parallelism(rng_seed = seed, max_workers = 4)

# load data ----
filtered_seurat <- qs_read(
  file.path(data_dir, 'filtered_seurat.qs'),
  nthreads = nthreads
)

## Check doublet calls ----
# rerun clustering
filtered_seurat <- FindVariableFeatures(
  filtered_seurat,
  selection.method = "vst",
  nfeatures = 3000
)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat)

elbow = calc_pcs_elbow(filtered_seurat)
plot_pcs_elbow(elbow)

filtered_seurat <- FindNeighbors(
  filtered_seurat,
  reduction = 'pca',
  dims = 1:20
)

# leiden clustering more  sensitive to resolution parameter than Louvain,
# will extend lower range to capture more broad clusters
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5)
filtered_seurat <- FindClusters(
  filtered_seurat,
  algorithm = 4, # Leiden clustering
  resolution = resolutions
)
filtered_seurat = RunUMAP(filtered_seurat, reduction = 'pca', dims = 1:20)

# calcualte percent.mt quartiles again since we filtered cells
filtered_seurat@meta.data$mito.level = cut(
  filtered_seurat@meta.data$percent.mt,
  breaks = c(-Inf, 0.7, 1.4, 2.4, Inf),
  labels = c("1st", "2nd", "3rd", "4th"),
  ordered_result = TRUE
)

# add doublet round metadata column for visualizing doublet rounds in umaps and violin plots
filtered_seurat$doublet_round <- dplyr::case_when(
  filtered_seurat$scDblFinder.class_r1 == "doublet" ~ "round1_doublet",
  !is.na(filtered_seurat$scDblFinder.class_r2) &
    filtered_seurat$scDblFinder.class_r2 == "doublet" ~ "round2_doublet",
  TRUE ~ "singlet"
)

# qs_save(
#     filtered_seurat,
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )

# filtered_seurat <- qs_read(
#     file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )

## Panel: doublet UMAPs ----

# visualize clustering resolution stability with clustree
clustree::clustree(filtered_seurat, prefix = 'RNA_snn_res.')
ggsave(
  filename = file.path(fig_dir, 'clustree_before_doublet_filtering.png'),
  width = 10,
  height = 15
)

# clustree::clustree_overlay(
#     filtered_seurat,
#     prefix = 'RNA_snn_res.',
#     overlay = 'doublet_final',
#     node_colour = "doublet_final",
#     node_colour_palette = doublet_cols
# )

panel_resolutions <- wrap_plots(
  lapply(resolutions, function(res) {
    DimPlot(
      filtered_seurat,
      reduction = 'umap',
      group.by = paste0('RNA_snn_res.', res),
      shuffle = TRUE
    ) +
      labs(title = paste0('Resolution = ', res)) +
      theme_classic2()
  }),
  ncol = 3
)
save_panel(
  panel_resolutions,
  dir_umaps,
  'umap_resolutions_before_doublet_filtering.png',
  width = 18,
  height = 12
)

# counts for the combined panel subtitle
n_total_d <- ncol(filtered_seurat) # total number of cells
n_r1_d <- sum(
  filtered_seurat$scDblFinder.class_r1 == "doublet",
  na.rm = TRUE
)
n_r2_d <- sum(
  filtered_seurat$scDblFinder.class_r2 == "doublet",
  na.rm = TRUE
)
n_any_d <- n_r1_d + n_r2_d
pct_d <- function(n) sprintf("%.2f%%", 100 * n / n_total_d)
doublet_cols <- c(singlet = "grey80", doublet = "#F8766D")

panel_doublets <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = "umap",
    group.by = "RNA_snn_res.0.15",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) +
    NoLegend() +
    labs(title = "Clusters"),
  DimPlot(
    filtered_seurat,
    reduction = "umap",
    group.by = "scDblFinder.class_r1",
    cols = doublet_cols,
    shuffle = TRUE
  ) +
    labs(
      title = "Round 1 doublets",
      subtitle = sprintf("doublets = %d (%s)", n_r1_d, pct_d(n_r1_d))
    ),
  DimPlot(
    filtered_seurat,
    reduction = "umap",
    group.by = "scDblFinder.class_r2",
    cols = doublet_cols,
    na.value = "grey95",
    shuffle = TRUE
  ) +
    labs(
      title = "Round 2 doublets (R1 singlets only)",
      subtitle = sprintf("doublets = %d (%s)", n_r2_d, pct_d(n_r2_d))
    ),
  DimPlot(
    filtered_seurat,
    reduction = "umap",
    group.by = "doublet_final",
    cols = doublet_cols,
    shuffle = TRUE
  ) +
    labs(
      title = "Combined doublet calls",
      subtitle = sprintf(
        "total = %d | R1 = %d (%s) | R2 = %d (%s) | any = %d (%s)",
        n_total_d,
        n_r1_d,
        pct_d(n_r1_d),
        n_r2_d,
        pct_d(n_r2_d),
        n_any_d,
        pct_d(n_any_d)
      )
    ),
  ncol = 2
)
save_panel(
  panel_doublets,
  dir_umaps,
  'doublet_umaps.png',
  width = 16,
  height = 14
)

## QC violins split by doublet call and doublet round ----
plot_doublets(filtered_seurat, fig_dir = fig_dir)
plot_doublet_rounds(
  filtered_seurat,
  reduction = "umap",
  group.by = "RNA_snn_res.0.15",
  fig_dir = fig_dir,
  umap_dir = dir_umaps
)

panel_doublet_qc <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.15',
    shuffle = TRUE,
    label = TRUE
  ),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'doublet_final',
    shuffle = TRUE,
    cols = c("#00BFC4", "#F8766D")
  ),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'condition',
    shuffle = TRUE
  ),
  qc_feature_plot(
    filtered_seurat,
    'scDblFinder.score_r1',
    min_cutoff = 0,
    max_cutoff = 100
  ),
  qc_feature_plot(
    filtered_seurat,
    'scDblFinder.score_r2',
    min_cutoff = 0,
    max_cutoff = 100
  ),
  qc_feature_plot(
    filtered_seurat,
    'nCount_RNA',
    min_cutoff = 'q10',
    max_cutoff = 'q95'
  ),
  qc_feature_plot(
    filtered_seurat,
    'nFeature_RNA',
    min_cutoff = 'q10',
    max_cutoff = 'q90'
  ),
  qc_feature_plot(
    filtered_seurat,
    'EGFP',
    min_cutoff = 0,
    max_cutoff = 'q75'
  ),
  qc_feature_plot(
    filtered_seurat,
    'percent.mt',
    min_cutoff = 0,
    max_cutoff = 'q75'
  ),
  ncol = 3
)
save_panel(
  panel_doublet_qc,
  dir_umaps,
  'filtered_before_removing_doublets_umaps.png',
  width = 18,
  height = 12
)

panel_doublets_vln <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.15',
    shuffle = TRUE,
    label = TRUE
  ) +
    labs(title = "Clusters") +
    theme_classic2(),
  VlnPlot(
    filtered_seurat,
    features = 'scDblFinder.score_r1',
    group.by = 'RNA_snn_res.0.15',
    pt.size = 0.1,
    alpha = 0.05
  ) +
    labs(title = "Round 1 doublet score") +
    theme_classic2(),
  VlnPlot(
    filtered_seurat,
    features = 'nCount_RNA',
    group.by = 'RNA_snn_res.0.15',
    pt.size = 0.1,
    alpha = 0.05
  ) +
    labs(title = "UMI counts per cluster") +
    theme_classic2(),
  VlnPlot(
    filtered_seurat,
    features = 'nFeature_RNA',
    group.by = 'RNA_snn_res.0.15',
    pt.size = 0.1,
    alpha = 0.05
  ) +
    labs(title = "Gene counts per cluster") +
    theme_classic2(),
  ncol = 2
)
save_panel(
  panel_doublets_vln,
  dir_umaps,
  'vlns_filtered_before_removing_doublets.png',
  width = 18,
  height = 12
)

round_cols <- c(
  singlet = "grey80",
  round1_doublet = "#F8766D",
  round2_doublet = "#00BFC4"
)

panel_cc_f <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.15',
    shuffle = TRUE,
    repel = TRUE,
    label = TRUE
  ) +
    labs(title = "Cluster") +
    theme_classic2(),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'Phase',
    shuffle = TRUE
  ) +
    labs(title = "Cell Cycle Phase") +
    theme_classic2(),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'doublet_round',
    cols = round_cols,
    shuffle = TRUE
  ) +
    labs(title = "Doublet Calls") +
    theme_classic2(),
  qc_feature_plot(
    filtered_seurat,
    'scDblFinder.score_r1',
    min_cutoff = 0,
    max_cutoff = '100'
  ) +
    labs(title = "Round 1 Doublet Score"),
  qc_feature_plot(
    filtered_seurat,
    'scDblFinder.score_r2',
    min_cutoff = 0,
    max_cutoff = '100'
  ) +
    labs(title = "Round 2 Doublet Score"),
  ncol = 3
)
save_panel(
  panel_cc_f,
  dir_umaps,
  'cell_cycle_umaps_before_removing_doublets.png',
  width = 18,
  height = 12
)

## flag high doublet score cells on UMAP ----
filtered_seurat$doublet_score_r1_high <-
  !is.na(filtered_seurat$scDblFinder.score_r1) &
  filtered_seurat$scDblFinder.score_r1 > 0.9
filtered_seurat$doublet_score_r2_high <-
  !is.na(filtered_seurat$scDblFinder.score_r2) &
  filtered_seurat$scDblFinder.score_r2 > 0.9

high_score_cols <- c("FALSE" = "grey80", "TRUE" = "#F8766D")

panel_doublet_score_high <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.15',
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) +
    NoLegend() +
    labs(title = "Clusters (res 0.15)"),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'doublet_round',
    cols = round_cols,
    shuffle = TRUE
  ) +
    NoLegend() +
    labs(title = "Doublet Calls by Round"),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'doublet_score_r1_high',
    cols = high_score_cols,
    shuffle = TRUE
  ) +
    labs(
      title = "Round 1: doublet score > 0.9",
      subtitle = sprintf(
        "n = %d (%.2f%%)",
        sum(filtered_seurat$doublet_score_r1_high),
        100 * mean(filtered_seurat$doublet_score_r1_high)
      )
    ),
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'doublet_score_r2_high',
    cols = high_score_cols,
    shuffle = TRUE
  ) +
    labs(
      title = "Round 2: doublet score > 0.9",
      subtitle = sprintf(
        "n = %d (%.2f%%)",
        sum(filtered_seurat$doublet_score_r2_high),
        100 * mean(filtered_seurat$doublet_score_r2_high)
      )
    ),
  ncol = 2
)

save_panel(
  panel_doublet_score_high,
  dir_umaps,
  'doublet_score_high_umaps.png',
  width = 18,
  height = 6
)

## Panel: nUMI / nGene vs doublet score (per round) ----
# diagnostic for choosing a score-based threshold over the algorithm's
# class call. Colored by the algorithm's class; dashed line at 0.9.
score_scatter <- function(score_col, y_col, class_col, title) {
  md <- filtered_seurat@meta.data
  md <- md[!is.na(md[[score_col]]), ]
  ggplot(
    md,
    aes(
      x = .data[[score_col]],
      y = .data[[y_col]],
      color = .data[[class_col]]
    )
  ) +
    scattermore::geom_scattermore(pointsize = 1.4, alpha = 0.4) +
    scale_y_log10() +
    scale_color_manual(values = doublet_cols, na.value = "grey95") +
    geom_vline(xintercept = 0.9, linetype = "dashed", color = "grey30") +
    labs(title = title, x = score_col, y = y_col, color = class_col) +
    theme_classic2()
}

panel_score_scatter <- wrap_plots(
  score_scatter(
    "scDblFinder.score_r1",
    "nCount_RNA",
    "scDblFinder.class_r1",
    "Round 1: nUMI vs score"
  ),
  score_scatter(
    "scDblFinder.score_r1",
    "nFeature_RNA",
    "scDblFinder.class_r1",
    "Round 1: nGene vs score"
  ),
  score_scatter(
    "scDblFinder.score_r2",
    "nCount_RNA",
    "scDblFinder.class_r2",
    "Round 2: nUMI vs score (r1 singlets)"
  ),
  score_scatter(
    "scDblFinder.score_r2",
    "nFeature_RNA",
    "scDblFinder.class_r2",
    "Round 2: nGene vs score (r1 singlets)"
  ),
  ncol = 2
)
save_panel(
  panel_score_scatter,
  fig_dir,
  'doublet_score_vs_libsize_scatter.png',
  width = 14,
  height = 12
)

## Panel: per-cluster doublet score summary ----
# % cells with score > 0.9 per cluster, side-by-side r1 / r2,
# paired with the cluster UMAP for cross-reference. Clusters with
# high bars are likely doublet-driven and should be reviewed regardless
# of the chosen per-cell threshold.
cluster_doublet_summary <- filtered_seurat@meta.data %>%
  as_tibble() %>%
  group_by(cluster = RNA_snn_res.0.15) %>%
  summarise(
    n_cells = dplyr::n(),
    pct_r1_high = 100 *
      mean(scDblFinder.score_r1 > 0.9, na.rm = TRUE),
    pct_r2_high = 100 *
      mean(scDblFinder.score_r2 > 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cluster = forcats::fct_reorder(
      cluster,
      pmax(pct_r1_high, pct_r2_high, na.rm = TRUE),
      .desc = TRUE
    )
  )

bar_pct_high <- cluster_doublet_summary %>%
  pivot_longer(
    cols = c(pct_r1_high, pct_r2_high),
    names_to = "round",
    values_to = "pct"
  ) %>%
  mutate(
    round = recode(
      round,
      pct_r1_high = "Round 1",
      pct_r2_high = "Round 2"
    )
  ) %>%
  ggplot(aes(x = cluster, y = pct, fill = round)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  labs(
    title = "% cells with doublet score > 0.9 per cluster",
    x = "Cluster (res 0.15)",
    y = "% cells > 0.9",
    fill = NULL
  ) +
  theme_classic2()

panel_cluster_doublet <- wrap_plots(
  DimPlot(
    filtered_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.15',
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) +
    NoLegend() +
    labs(title = "Clusters (res 0.15)"),
  bar_pct_high,
  ncol = 2,
  widths = c(1, 1.3)
)
save_panel(
  panel_cluster_doublet,
  fig_dir,
  'doublet_cluster_pct_high.png',
  width = 18,
  height = 8
)


## check cluster markers for for doublet clusters ----
Idents(filtered_seurat) <- filtered_seurat$RNA_snn_res.0.15
all.pos.markers = FindAllMarkers(
  filtered_seurat,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf
)

# Combine markers with gene descriptions
annotations = read.csv(file.path(data_dir, 'annotationhub_mice_genes.csv'))
ann.pos.markers = left_join(
  x = all.pos.markers,
  y = annotations[, c("gene_name", "description")],
  by = c("gene" = "gene_name")
) %>%
  unique()

# add pct.diff col
ann.pos.markers$pct.diff = ann.pos.markers$pct.1 - ann.pos.markers$pct.2

# Rearrange the columns to be more intuitive
ann.pos.markers = ann.pos.markers[, c(6, 7, 2:4, 9, 1, 5, 8)]

# Order the rows by log2FC values
ann.pos.markers = ann.pos.markers %>%
  dplyr::arrange(cluster, desc(avg_log2FC))

# save markers
write.csv(
  ann.pos.markers,
  file.path(data_dir, 'annotated_pos_markers_res.0.15_filtered.csv'),
  row.names = FALSE
)


# Extract top 20 markers per cluster
top20 = ann.pos.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05, pct.diff > 0.25) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20)
write.csv(
  top20,
  file.path(data_dir, 'top20_markers_res.0.15_filtered_seurat.csv'),
  row.names = FALSE
)


## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster, assay = 'RNA', group.by = 'condition') {
  FindConservedMarkers(
    filtered_seurat, # need to specify exact seurat object here
    ident.1 = cluster,
    assay = assay,
    grouping.var = group.by, # condition variable
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    min.diff.pct = -Inf
  ) %>%
    rownames_to_column(var = "gene") %>%
    left_join(
      y = unique(annotations[, c("gene_name", "description")]),
      by = c("gene" = "gene_name")
    ) %>%
    cbind(cluster_id = cluster, .)
}


### get conserved markers across conditions
n = length(unique(filtered_seurat$RNA_snn_res.0.15)) # number of clusters at resolution 0.15
# get conserved markers # adjust 0:## to (number of clusters - 1)
conserved.pos.markers = map_dfr(c(1:n), get_conserved) # 17 clusters at res 0.15. Leiden clustering starts from index 1

# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
conserved.pos.markers = conserved.pos.markers %>%
  mutate(
    avg_log2FC = (ctrl_avg_log2FC + wnt_ko_avg_log2FC) / 2,
    max_adj_pval = ifelse(
      ctrl_p_val_adj > wnt_ko_p_val_adj,
      ctrl_p_val_adj,
      wnt_ko_p_val_adj
    ),
    avg.pct.diff = ((wnt_ko_pct.1 - wnt_ko_pct.2) +
      (ctrl_pct.1 - ctrl_pct.2)) /
      2,
    min.pct.1 = ifelse(ctrl_pct.1 < wnt_ko_pct.1, ctrl_pct.1, wnt_ko_pct.1),
    max.pct.2 = ifelse(ctrl_pct.2 > wnt_ko_pct.2, ctrl_pct.2, wnt_ko_pct.2)
  )

# Rearrange the columns to be more intuitive
conserved.pos.markers = conserved.pos.markers[, c(1:2, 16:20, 15, 3:14)]

# Order the rows by log2FC values
conserved.pos.markers = conserved.pos.markers %>%
  dplyr::arrange(cluster_id, desc(avg_log2FC))

# save markers
write.csv(
  conserved.pos.markers,
  file.path(
    data_dir,
    'annotated_conserved_pos_markers_res.0.15_filtered.csv'
  ),
  row.names = FALSE
)


# Extract top 20 conserved markers per cluster
top20.conserved = conserved.pos.markers %>%
  dplyr::group_by(cluster_id) %>%
  dplyr::filter(max_adj_pval < 0.05, avg.pct.diff > 0.25) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  # slice_max(order_by=tibble(avg.pct.diff),n=20)
  slice_max(order_by = tibble(avg_log2FC, avg.pct.diff), n = 20)
write.csv(
  top20.conserved,
  file.path(data_dir, 'top20_conserved_markers_res.0.15_filtered_seurat.csv'),
  row.names = FALSE
)

# top10.conserved = top10.conserved[!top10.conserved$cluster_id == 26,] # remove because only one condition

gc() # free up RAM

## gptcelltype annotation with top markers ----
Sys.setenv(OPENAI_API_KEY = '')
# get unique markers per cluster for annotation
res <- GPTCelltype::gptcelltype(
  ann.pos.markers,
  tissuename = 'mouse stomach with KO beta-catenin in chief cells',
  model = 'gpt-5.4',
  topgenenumber = 10
)
res

res <- GPTCelltype::gptcelltype(
  conserved.pos.markers,
  tissuename = 'mouse stomach with KO beta-catenin in chief cells',
  model = 'gpt-5.4',
  topgenenumber = 20
)
res

# harmonize dataset with doublets included ----
# keep doublets for now since whole clusters have been labeled as dobulets.
# # Remove later post-hoc once verified.

# slim down seurat object first to only what we need downstream
filtered_seurat <- DietSeurat(
  filtered_seurat,
  assays = 'RNA',
  layers = c('counts')
)

## split data by sample ----
tmp <- filtered_seurat
tmp[['RNA']] <- split(
  filtered_seurat[['RNA']],
  f = filtered_seurat$sample
)
gc()

options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
doublets_sct_seurat <- SCTransform(
  tmp,
  vars.to.regress = c('percent.mt', 'cc.difference'), # https://satijalab.org/seurat/articles/cell_cycle_vignette.html
  variable.features.n = 5000
)
rm(tmp) # free up memory
gc()


resolutions <- c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5)
doublets_sct_seurat <- RunPCA(doublets_sct_seurat, assay = 'SCT')
doublets_sct_seurat <- FindNeighbors(
  doublets_sct_seurat,
  reduction = 'pca',
  dims = 1:50
)
doublets_sct_seurat <- FindClusters(
  doublets_sct_seurat,
  algorithm = 4, # Leiden clustering
  resolution = resolutions
)
doublets_sct_seurat = RunUMAP(
  doublets_sct_seurat,
  reduction = 'pca',
  dims = 1:50
)

qs_save(
  doublets_sct_seurat,
  file = file.path(data_dir, 'doublets_sct_seurat.qs'),
  nthreads = nthreads
)

# doublets_sct_seurat <- qs_read(
#   file.path(data_dir, 'doublets_sct_seurat.qs'),
#   nthreads = nthreads
# )
gc()

## sample correction with Harmony ----
doublets_sct_seurat <- IntegrateLayers(
  doublets_sct_seurat,
  method = HarmonyIntegration,
  orig.reduction = 'pca',
  new.reduction = 'harmony',
  assay = 'SCT'
)
gc()
doublets_sct_seurat <- FindNeighbors(
  doublets_sct_seurat,
  reduction = 'harmony',
  dims = 1:50,
  graph.name = c('harmony_nn', 'harmony_snn')
)
doublets_sct_seurat <- FindClusters(
  doublets_sct_seurat,
  algorithm = 4, # Leiden clustering
  resolution = resolutions,
  graph.name = 'harmony_snn'
)
doublets_sct_seurat <- RunUMAP(
  doublets_sct_seurat,
  reduction = 'harmony',
  reduction.name = 'umap_harmony',
  dims = 1:50
)

qs_save(
  doublets_sct_seurat,
  file = file.path(data_dir, 'doublets_harmony_seurat.qs'),
  nthreads = nthreads
)


# subset to singlets data to produce doublet filtered dataset----
filtered_seurat <- subset(filtered_seurat, doublet_final == "singlet")

# SCTranform normalization ----
# split data by sample for sample-wise normalization
filtered_seurat[['RNA']] <- split(
  filtered_seurat[['RNA']],
  f = filtered_seurat$sample
)
gc()

# selecting 5000 features to account for the possible different rankings of highly variable genes
# between samples. Since we're normalizing with SCTtransform this # shouldn't change things much.
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
sct_seurat <- SCTransform(
  filtered_seurat,
  vars.to.regress = c('percent.mt', 'cc.difference'), # https://satijalab.org/seurat/articles/cell_cycle_vignette.html
  variable.features.n = 5000
)
rm(filtered_seurat) # free up memory
gc()

# qs_save(
#     sct_seurat,
#     file = file.path(data_dir, 'sct_seurat.qs'),
#     nthreads = nthreads
# )

# sct_seurat <- qs_read(
#     file.path(data_dir, 'sct_seurat.qs'),
#     nthreads = nthreads
# )

## SCT dimenion reduction and clustering ----

# leiden clustering more  sensitive to resolution parameter than Louvain,
# will extend lower range to capture more broad clusters
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0)

sct_seurat <- RunPCA(sct_seurat, assay = 'SCT')
sct_seurat <- FindNeighbors(sct_seurat, reduction = 'pca', dims = 1:50) # ok to use more PCs with SCT
sct_seurat <- FindClusters(
  sct_seurat,
  algorithm = 4, # Leiden clustering
  resolution = resolutions
)
sct_seurat = RunUMAP(sct_seurat, reduction = 'pca', dims = 1:50)

### visualize SCT UMAPs before batch correction ----
panel_sct_umap <- wrap_plots(
  DimPlot(
    sct_seurat,
    reduction = 'umap',
    group.by = 'SCT_snn_res.0.05',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap',
    group.by = 'sample',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap',
    group.by = 'condition',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap',
    group.by = 'Phase',
    shuffle = TRUE
  ),
  ncol = 2
)
save_panel(
  panel_sct_umap,
  dir_umaps,
  'sct_umap_overview.png',
  width = 14,
  height = 12
)
panel_sct_umap
# clusters don't mix well by sample or condition,
# looks like we do need harmony

## sample correction with Harmony ----
sct_seurat <- IntegrateLayers(
  sct_seurat,
  method = HarmonyIntegration,
  orig.reduction = 'pca',
  new.reduction = 'harmony',
  assay = 'SCT'
)
sct_seurat <- FindNeighbors(
  sct_seurat,
  reduction = 'harmony',
  dims = 1:50,
  graph.name = c('harmony_nn', 'harmony_snn')
)
sct_seurat <- FindClusters(
  sct_seurat,
  algorithm = 4, # Leiden clustering
  resolution = resolutions,
  graph.name = 'harmony_snn'
)
sct_seurat <- RunUMAP(
  sct_seurat,
  reduction = 'harmony',
  reduction.name = 'umap_harmony',
  dims = 1:50
)

### visualize Harmony UMAPs ----
panel_harmony_umap <- wrap_plots(
  DimPlot(
    sct_seurat,
    reduction = 'umap_harmony',
    group.by = 'harmony_snn_res.0.05',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap_harmony',
    group.by = 'sample',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap_harmony',
    group.by = 'condition',
    shuffle = TRUE
  ),
  DimPlot(
    sct_seurat,
    reduction = 'umap_harmony',
    group.by = 'Phase',
    shuffle = TRUE
  ),
  ncol = 2
)
save_panel(
  panel_harmony_umap,
  dir_umaps,
  'harmony_umap_overview.png',
  width = 14,
  height = 12
)

qs_save(
  sct_seurat,
  file = file.path(data_dir, 'harmony_seurat.qs'),
  nthreads = nthreads
)

#package versions----
sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

# Matrix products: default
# BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRblas.so
# LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
#  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# time zone: America/Detroit
# tzcode source: system (glibc)

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] future_1.70.0               scattermore_1.2             scCustomize_3.2.4
#  [4] scales_1.4.0                viridis_0.6.5               viridisLite_0.4.3
#  [7] ggpubr_0.6.3                ggthemes_5.2.0              patchwork_1.3.2
# [10] cowplot_1.2.0               scDblFinder_1.24.10         yaml_2.3.12
# [13] BiocParallel_1.44.0         RhpcBLASctl_0.23-42         qs2_0.1.7
# [16] zeallot_0.2.0.9000          scuttle_1.20.0              SingleCellExperiment_1.32.0
# [19] SummarizedExperiment_1.40.0 Biobase_2.70.0              GenomicRanges_1.62.1
# [22] Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.1
# [25] BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0
# [28] matrixStats_1.5.0           lubridate_1.9.4             forcats_1.0.1
# [31] stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.2
# [34] readr_2.1.5                 tidyr_1.3.2                 tibble_3.3.1
# [37] ggplot2_4.0.2               tidyverse_2.0.0             sctransform_0.4.3
# [40] Seurat_5.4.0                SeuratObject_5.4.0          sp_2.2-0

# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0     bitops_1.0-9              httr_1.4.8
#   [4] RColorBrewer_1.1-3        tools_4.5.1               backports_1.5.1
#   [7] R6_2.6.1                  lazyeval_0.2.3            uwot_0.2.4
#  [10] withr_3.0.2               gridExtra_2.3             progressr_0.19.0
#  [13] textshaping_1.0.3         cli_3.6.6                 spatstat.explore_3.5-3
#  [16] fastDummies_1.7.5         labeling_0.4.3            S7_0.2.1
#  [19] spatstat.data_3.1-8       ggridges_0.5.7            pbapply_1.7-4
#  [22] systemfonts_1.3.2         Rsamtools_2.26.0          harmony_1.2.4
#  [25] scater_1.38.1             sessioninfo_1.2.3         parallelly_1.46.1
#  [28] mcprogress_0.1.1          limma_3.66.0              shape_1.4.6.1
#  [31] BiocIO_1.20.0             ica_1.0-3                 spatstat.random_3.4-2
#  [34] car_3.1-5                 Matrix_1.7-3              ggbeeswarm_0.7.3
#  [37] abind_1.4-8               lifecycle_1.0.5           edgeR_4.8.2
#  [40] snakecase_0.11.1          carData_3.0-6             glmGamPoi_1.22.0
#  [43] SparseArray_1.10.10       Rtsne_0.17                paletteer_1.7.0
#  [46] grid_4.5.1                promises_1.5.0            dqrng_0.4.1
#  [49] crayon_1.5.3              miniUI_0.1.2              lattice_0.22-7
#  [52] beachmat_2.26.0           cigarillo_1.0.0           pillar_1.11.1
#  [55] metapod_1.18.0            rjson_0.2.23              xgboost_3.2.1.1
#  [58] future.apply_1.20.2       codetools_0.2-20          glue_1.8.0
#  [61] leidenbase_0.1.36         spatstat.univar_3.1-4     data.table_1.18.2.1
#  [64] vctrs_0.7.3               png_0.1-9                 spam_2.11-3
#  [67] gtable_0.3.6              rematch2_2.1.2            S4Arrays_1.10.1
#  [70] mime_0.13                 survival_3.8-3            statmod_1.5.1
#  [73] bluster_1.20.0            fitdistrplus_1.2-6        ROCR_1.0-12
#  [76] nlme_3.1-168              RcppAnnoy_0.0.23          GenomeInfoDb_1.46.2
#  [79] rprojroot_2.1.1           irlba_2.3.7               vipor_0.4.7
#  [82] KernSmooth_2.23-26        otel_0.2.0                colorspace_2.1-2
#  [85] ggrastr_1.0.2             tidyselect_1.2.1          compiler_4.5.1
#  [88] curl_7.0.0                BiocNeighbors_2.4.0       DelayedArray_0.36.1
#  [91] plotly_4.12.0             stringfish_0.18.0         rtracklayer_1.70.1
#  [94] lmtest_0.9-40             digest_0.6.39             goftest_1.2-3
#  [97] spatstat.utils_3.2-0      XVector_0.50.0            htmltools_0.5.9
# [100] pkgconfig_2.0.3           sparseMatrixStats_1.22.0  fastmap_1.2.0
# [103] rlang_1.2.0               GlobalOptions_0.1.4       htmlwidgets_1.6.4
# [106] UCSC.utils_1.6.1          DelayedMatrixStats_1.32.0 shiny_1.11.1
# [109] farver_2.1.2              zoo_1.8-15                jsonlite_2.0.0
# [112] BiocSingular_1.26.1       RCurl_1.98-1.18           magrittr_2.0.5
# [115] Formula_1.2-5             dotCall64_1.2             Rcpp_1.1.1
# [118] reticulate_1.46.0         stringi_1.8.7             MASS_7.3-65
# [121] plyr_1.8.9                listenv_0.10.1            ggrepel_0.9.8
# [124] deldir_2.0-4              Biostrings_2.78.0         splines_4.5.1
# [127] tensor_1.5.1              hms_1.1.3                 circlize_0.4.18
# [130] locfit_1.5-9.12           igraph_2.2.3              spatstat.geom_3.6-0
# [133] ggsignif_0.6.4            RcppHNSW_0.6.0            reshape2_1.4.5
# [136] ScaledMatrix_1.18.0       XML_3.99-0.23             RcppParallel_5.1.11-2
# [139] scran_1.38.1              ggprism_1.0.7             tzdb_0.5.0
# [142] httpuv_1.6.16             RANN_2.6.2                polyclip_1.10-7
# [145] janitor_2.2.1             rsvd_1.0.5                broom_1.0.10
# [148] xtable_1.8-4              restfulr_0.0.16           RSpectra_0.16-2
# [151] rstatix_0.7.3             later_1.4.8               ragg_1.5.0
# [154] beeswarm_0.4.0            GenomicAlignments_1.46.0  cluster_2.1.8.1
# [157] timechange_0.3.0          globals_0.19.1            here_1.0.2
