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
fig_dir <- file.path(base_dir, "qc")
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

#
# # Empty list to populate seurat object for each sample
# list_seurat <- list()
#
#
# for (sample in 1:8) {
#   # Path to data directory
#   data_dir <- paste0("~/SingleCell2016/per_sample_outs/15160-CH-",sample,"/sample_filtered_feature_bc_matrix.h5")
#
#   # Create a Seurat object for each sample
#   seurat_data <- Read10X_h5(filename = data_dir,)
#   seurat_obj <- CreateSeuratObject(counts = seurat_data,
#                                    project = 'wnt')
#
#   # Save seurat object to list
#   list_seurat[[sample]] <- seurat_obj
# }
#
# qs_save(list_seurat,file='~/SingleCell2016/data/seurat__list.qs',nthreads = nthreads)
#
# cell_ids = c('wnt_ko_1','ctrl_1','wnt_ko_2','ctrl_2','wnt_ko_3', 'ctrl_3','wnt_ko_4', 'ctrl_4')
#
#
# merged_seurat = merge(list_seurat[[1]], list_seurat[-1], add.cell.ids=cell_ids)
#
# add_metadata <- function(seurat_object) {
#   metadata <- seurat_object@meta.data
#
#   metadata$cells = colnames(seurat_object)
#
#   metadata$sample <- sub("^((?:wnt_ko|ctrl)_\\d+).*", "\\1", colnames(seurat_object))
#   metadata$condition <- sub("^(wnt_ko|ctrl).*", "\\1", colnames(seurat_object))
#   metadata$ko <- ifelse(metadata$condition =='wnt_ko',1,0)
#
#   seurat_object@meta.data = metadata
#
#   return(seurat_object)
# }
#
# merged_seurat = add_metadata(merged_seurat)

# # save merged seurat object to disk as qs file for faster loading in future steps
# qs_save(merged_seurat, file='~/SingleCell2016/data/merged_seurat.qs',nthreads = nthreads)

# # read merged seurat object from disk
# merged_seurat <- qs_read(
#     file.path(data_dir, "merged_seurat.qs"),
#     nthreads = nthreads
# )

# calculate QC metrics and add to metadata
merged_seurat <- calc_metrics(merged_seurat)
gc()

# visualize QC metrics and potential cutoffs for filtering ----
plot_qcs(merged_seurat, fig_dir = fig_dir)
gc()

# normalize, find variable features, and run PCA / UMAP for dimensionality reduction and visualization ----
## normlize ----
if (
  !("data" %in%
    slotNames(merged_seurat[["RNA"]]) &&
    ncol(merged_seurat[["RNA"]]@data) > 0)
) {
  merged_seurat <- NormalizeData(
    merged_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
}
## find variable features ----
merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 3000
)

### visualize variable features ----
top25 = head(VariableFeatures(merged_seurat), 25)
# plot variable features with and without labels
p = VariableFeaturePlot(merged_seurat)
LabelPoints(plot = p, points = top25, repel = TRUE) +
  labs(title = 'Top 25 variably expressed genes') +
  theme_classic2()
ggsave(
  filename = file.path(fig_dir, 'variable_features.png'),
  width = 10,
  height = 8
)

## scale data and run PCA ----
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)

# qs_save(
#     merged_seurat,
#     file = file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )

# merged_seurat <- qs_read(
#     file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )

### examine pcs with elbow plot ----
elbow <- calc_pcs_elbow(merged_seurat)
plot_pcs_elbow(elbow, fig_dir = fig_dir)
pcs <- elbow$pcs


## clustering and UMAP ----
resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0)

# going with 20 pcs based on percent variation captured
merged_seurat <- FindNeighbors(merged_seurat, reduction = 'pca', dims = 1:20)
merged_seurat = FindClusters(
  object = merged_seurat,
  algorithm = 4, # Leiden algorithm clustering. Need to install leidenalg package in Python for this to work, as well as install.packages('leidenbase')
  resolution = resolutions
)
merged_seurat <- RunUMAP(merged_seurat, reduction = 'pca', dims = 1:20)

# cell cycle and mitochondrial quartiles ----
# create mitochondrial quartiles for visualization
merged_seurat@meta.data$mito.level = cut(
  merged_seurat@meta.data$percent.mt,
  breaks = c(-Inf, 2, 4, 6, Inf),
  labels = c("1st", "2nd", "3rd", "4th"),
  ordered_result = TRUE
)

# cell cycle genes, formatted to match genes in matrix with first letter upper-cased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)

# Score cells for cell cycle
merged_seurat = CellCycleScoring(
  merged_seurat,
  g2m.features = g2m.genes,
  s.features = s.genes,
  set.ident = F
)

# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Seurat suggests regressing out the difference between S and G2M phase scores
# in tissues with differentiating cells,
merged_seurat$cc.difference = merged_seurat$S.Score - merged_seurat$G2M.Score


# qs_save(
#     merged_seurat,
#     file = file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )

# merged_seurat <- qs_read(
#     file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )

# QC diagnostic panels ----

## check PCA for cell cycle genes ----
cc_test = RunPCA(
  merged_seurat,
  features = c(s.genes, g2m.genes),
  reduction.name = 'pca'
)
Idents(cc_test) <- cc_test$Phase
DimPlot(cc_test, reduction = 'pca', group.by = 'Phase', shuffle = T) +
  labs(title = 'PCA of cell cycle genes') +
  theme_classic2()
ggsave(
  filename = file.path(fig_dir, 'pca_cell_cycle_genes.png'),
  width = 10,
  height = 8
)

DimHeatmap(cc_test, dims = 1, cells = 5000, balanced = TRUE)
# can see top features for PC1 include cell cycle genes, Top2a and Mki67
TopFeatures(cc_test, reduction = 'pca', dim = 1)[1:20]
# [1] "Top2a"  "Mki67"  "Nusap1" "Tpx2"   "Ube2c"  "Cdca3"  "Cdk1"
#  [8] "Birc5"  "Ect2"   "Hmmr"   "Kif11"  "Ccnb2"  "Anln"   "Cdc20"
# [15] "Tacc3"  "Cenpf"  "Ckap2"  "Gtse1"  "Ckap2l" "Nek2"
rm(cc_test) # free up memory
gc()

## Panel 1: sample / condition / clusters / mito.level ----
panel_embedding <- (DimPlot(
  merged_seurat,
  reduction = 'umap',
  group.by = 'RNA_snn_res.0.6',
  shuffle = TRUE
) +
  DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'sample',
    shuffle = TRUE
  )) /
  (DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'condition',
    shuffle = TRUE
  ) +
    DimPlot(
      merged_seurat,
      reduction = 'umap',
      group.by = 'mito.level',
      shuffle = TRUE
    ))
# save_panel is a ggsave wrapper in plot_utils.R with
# default height,width and dpi set for saving panels with multiple subplots
save_panel(panel_embedding, dir_umaps, 'umap_embedding_overview.png')

## Panel 2: cell-cycle diagnostics ----
panel_cc <- wrap_plots(
  DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'RNA_snn_res.0.6',
    shuffle = TRUE,
    label = TRUE
  ),
  qc_feature_plot(merged_seurat, 'cc.difference'),
  qc_feature_plot(merged_seurat, 'S.Score'),
  qc_feature_plot(merged_seurat, 'G2M.Score'),
  ncol = 2
)
save_panel(panel_cc, dir_umaps, 'cell_cycle_umaps.png')

## Panel 3: QC metrics on UMAP ----
panel_qc_metrics <- wrap_plots(
  qc_feature_plot(merged_seurat, 'nCount_RNA'),
  qc_feature_plot(merged_seurat, 'nFeature_RNA'),
  qc_feature_plot(merged_seurat, 'percent.mt'),
  qc_feature_plot(merged_seurat, 'percent.hb'),
  qc_feature_plot(merged_seurat, 'complexity'),
  ncol = 3
)
save_panel(
  panel_qc_metrics,
  dir_umaps,
  'qc_metrics_umap.png',
  width = 21,
  height = 14
)

# flag cells based on QC metrics ----
## label outliers ----
merged_seurat$outlier <- isOutlier(
  merged_seurat$nCount_RNA,
  nmads = 5,
  log = TRUE,
  batch = merged_seurat$sample
) |
  isOutlier(
    merged_seurat$nFeature_RNA,
    nmads = 5,
    log = TRUE,
    batch = merged_seurat$sample
  ) |
  isOutlier(
    merged_seurat$pct_counts_in_top_20_genes,
    nmads = 5,
    log = FALSE,
    batch = merged_seurat$sample
  ) |
  (merged_seurat$complexity < 0.8)

merged_seurat$mt_outlier <- isOutlier(
  merged_seurat$percent.mt,
  nmads = 3,
  batch = merged_seurat$sample
) |
  (merged_seurat$percent.mt > 8)

merged_seurat$flagged <- merged_seurat$outlier | merged_seurat$mt_outlier

# qs_save(
#     merged_seurat,
#     file = file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )
# merged_seurat <- qs_read(
#     file.path(data_dir, 'metrics_seurat.qs'),
#     nthreads = nthreads
# )

## visulize outliers ----
plot_outliers(merged_seurat, mt_hard_cutoff = 8, fig_dir = fig_dir)

## Panel 4: outlier flags on UMAP ----
# FALSE -> default teal, TRUE -> default salmon/red (swapped from ggplot default)
flag_cols <- c("FALSE" = "#00BFC4", "TRUE" = "#F8766D")

panel_outliers <- wrap_plots(
  DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'outlier',
    cols = flag_cols,
    order = "TRUE",
    shuffle = TRUE
  ),
  DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'mt_outlier',
    cols = flag_cols,
    order = "TRUE",
    shuffle = TRUE
  ),
  DimPlot(
    merged_seurat,
    reduction = 'umap',
    group.by = 'flagged',
    cols = flag_cols,
    order = "TRUE",
    shuffle = TRUE
  ),
  ncol = 2
)
save_panel(
  panel_outliers,
  dir_umaps,
  'outlier_umaps.png',
  width = 14,
  height = 12
)


# marker umaps, pre-filtering ----
# stomach markers, rough estimate from Perplexity
stomach_all_markers <- c(
  # --- CHIEF & SECRETORY ---
  "Chief (CRISPR KO)" = "EGFP", # Cldn18-EGFP; upregulated in Ctnnb1 KO
  "Chief (Endog)" = "Gif",
  "Chief (Digestive)" = "Pgc",
  "Mucous Neck" = "Muc6",

  # --- EPITHELIAL HIERARCHY ---
  "Stem Cells" = "Lgr5",
  "TA (Cycling)" = "Mki67",
  "Pit Cells (Antrum)" = "Muc5ac", # surface mucous pit cells
  "Pit Cells (Corpus)" = "Gkn1", # gastrokine-1, corpus/antrum pit marker
  "Parietal Cells" = "Atp4b", # H+/K+ ATPase beta subunit
  "Parietal (Alt)" = "Atp4a", # H+/K+ ATPase alpha subunit
  "Enteroendocrine" = "Chga",
  "Tuft Cells" = "Dclk1",

  # --- IMMUNE COMPARTMENT (CD45+) ---
  "Immune (Pan)" = "Ptprc",
  "T Cells" = "Cd3e",
  "B Cells" = "Cd79a",
  "Plasma Cells" = "Jchain",
  "Macrophages" = "Adgre1",
  "Dendritic Cells" = "Itgax",
  "Neutrophils" = "S100a8",
  "ILCs" = "Klrb1c",

  # --- STROMAL & WALL ---
  "Stromal (Pan)" = "Pdgfra",
  "Telocytes (Niche)" = "Foxl1",
  "Fibroblasts" = "Col1a1",
  "Myofibroblasts" = "Acta2",
  "Endothelial" = "Pecam1",
  "Lymphatic" = "Prox1",
  "Pericytes" = "Rgs5",
  "Enteric Glia" = "S100b",
  "Smooth Muscle" = "Myh11"
)

# only keep marker if in data
stomach_all_markers <- stomach_all_markers[
  stomach_all_markers %in% rownames(merged_seurat)
]

# mitochondrial genes
mt_genes <- grep(
  "^mt-",
  rownames(merged_seurat),
  ignore.case = TRUE,
  value = TRUE
)

## Panel 5: mitochondrial genes on UMAP ----
panel_mt <- qc_feature_plot(merged_seurat, mt_genes)
save_panel(panel_mt, dir_umaps, 'mt_genes_umap.png', width = 18, height = 18)

## umap stomach markers by condition ----
Iterate_FeaturePlot_scCustom(
  merged_seurat,
  features = stomach_all_markers,
  reduction = 'umap',
  split.by = 'condition',
  colors_use = viridis_inferno_dark_high,
  file_path = paste0(dir_markers, '/'),
  file_name = 'umap',
  file_type = '.png'
)

## umap mt_gene by condition ----
Iterate_FeaturePlot_scCustom(
  merged_seurat,
  features = mt_genes,
  reduction = 'umap',
  split.by = 'condition',
  colors_use = viridis_inferno_dark_high,
  file_path = paste0(dir_markers, '/'),
  file_name = 'umap',
  file_type = '.png'
)

# filter flagged cells ----
filtered_seurat <- subset(merged_seurat, flagged == FALSE)
# qs_save(
#     filtered_seurat,
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )
# filtered_seurat <- qs_read(
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )

## visualize before and after cell-level filtering ----
plot_before_and_after(
  merged_seurat,
  filtered_seurat,
  fig_dir = fig_dir,
)

# gene-level filtering ----
filtered_seurat <- filter_genes(filtered_seurat, min_cells = 10)
# [filter_genes] 19071 genes before | 16054 after | 3017 removed (min_cells = 10)

# doublet detection ----
# Two rounds of scDblFinder:
#   round 1: full object
#   round 2: only cells that were singlets in round 1
# Writes round-specific + combined doublet labels to metadata.
rm(merged_seurat) # to free up memory
gc()

# required joinedlayers object for scDblFinder
# will first slim down seurat object to only what we need
filtered_seurat <- JoinLayers(filtered_seurat)
gc()

# qs_save(
#     filtered_seurat,
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )
# filtered_seurat <- qs_read(
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )

# set threads to 1 to avoid oversubscribing cores when using BiocParallel
# also need to add export MALLOC_ARENA_MAX=2 in .bashrc to prevent memory fragmentation issues with multithreading in some packages (e.g. igraph)
RhpcBLASctl::omp_set_num_threads(nthreads) # covers xgboost, igraph, uwot, Matrix
RhpcBLASctl::blas_set_num_threads(nthreads) # harmless no-op on reference BLAS
data.table::setDTthreads(nthreads) # data.table has its own pool
RcppParallel::setThreadOptions(numThreads = nthreads) # TBB-based packages

threads_serial(n = nthreads)

filtered_seurat <- run_doublet_detection(
  filtered_seurat,
  # bpparam = bp_param,
  # bpparam = MulticoreParam(
  #     workers = 4,
  #     progressbar = TRUE,
  #     RNGseed = seed
  # ),
  bpparam = SerialParam(RNGseed = seed)
)
gc()

# qs_save(
#     filtered_seurat,
#     file = file.path(data_dir, 'filtered_seurat.qs'),
#     nthreads = nthreads
# )

filtered_seurat <- qs_read(
  file.path(data_dir, 'filtered_seurat.qs'),
  nthreads = nthreads
)


sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

# Matrix products: default
# BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRblas.so
# LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# time zone: America/Detroit
# tzcode source: system (glibc)

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] scattermore_1.2             scCustomize_3.2.4           scales_1.4.0                viridis_0.6.5               viridisLite_0.4.3           ggpubr_0.6.3
#  [7] ggthemes_5.2.0              patchwork_1.3.2             cowplot_1.2.0               scDblFinder_1.24.10         yaml_2.3.12                 BiocParallel_1.44.0
# [13] ggraph_2.2.2                RhpcBLASctl_0.23-42         qs2_0.2.0                   zeallot_0.2.0.9000          scuttle_1.20.0              SingleCellExperiment_1.32.0
# [19] SummarizedExperiment_1.40.0 Biobase_2.70.0              GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.1
# [25] BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0       matrixStats_1.5.0           lubridate_1.9.5             forcats_1.0.1
# [31] stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.2                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1
# [37] ggplot2_4.0.3               tidyverse_2.0.0             sctransform_0.4.3           Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1

# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0    bitops_1.0-9             httr_1.4.8               RColorBrewer_1.1-3       numDeriv_2016.8-1.1      tools_4.5.1              backports_1.5.1
#   [8] R6_2.6.1                 sn_2.1.3                 lazyeval_0.2.3           uwot_0.2.4               withr_3.0.2              gridExtra_2.3            progressr_0.19.0
#  [15] cli_3.6.6                spatstat.explore_3.8-0   fastDummies_1.7.6        sandwich_3.1-1           mvtnorm_1.3-7            S7_0.2.2                 spatstat.data_3.1-9
#  [22] ggridges_0.5.7           pbapply_1.7-4            Rsamtools_2.26.0         scater_1.38.1            parallelly_1.47.0        plotrix_3.8-14           mcprogress_0.1.1
#  [29] limma_3.66.0             shape_1.4.6.1            BiocIO_1.20.0            ica_1.0-3                spatstat.random_3.4-5    car_3.1-5                Matrix_1.7-5
#  [36] ggbeeswarm_0.7.3         abind_1.4-8              lifecycle_1.0.5          multcomp_1.4-30          edgeR_4.8.2              snakecase_0.11.1         carData_3.0-6
#  [43] mathjaxr_2.0-0           SparseArray_1.10.10      Rtsne_0.17               paletteer_1.7.0          grid_4.5.1               promises_1.5.0           dqrng_0.4.1
#  [50] crayon_1.5.3             miniUI_0.1.2             lattice_0.22-9           beachmat_2.26.0          cigarillo_1.0.0          pillar_1.11.1            metapod_1.18.0
#  [57] rjson_0.2.23             xgboost_3.2.1.1          future.apply_1.20.2      codetools_0.2-20         mutoss_0.1-14            glue_1.8.1               spatstat.univar_3.1-7
#  [64] data.table_1.18.2.1      Rdpack_2.6.6             vctrs_0.7.3              png_0.1-9                spam_2.11-3              gtable_0.3.6             rematch2_2.1.2
#  [71] cachem_1.1.0             rbibutils_2.4.1          S4Arrays_1.10.1          mime_0.13                tidygraph_1.3.1          survival_3.8-6           statmod_1.5.1
#  [78] bluster_1.20.0           TH.data_1.1-5            fitdistrplus_1.2-6       ROCR_1.0-12              nlme_3.1-169             RcppAnnoy_0.0.23         GenomeInfoDb_1.46.2
#  [85] rprojroot_2.1.1          irlba_2.3.7              vipor_0.4.7              KernSmooth_2.23-26       otel_0.2.0               colorspace_2.1-2         mnormt_2.1.2
#  [92] ggrastr_1.0.2            tidyselect_1.2.1         compiler_4.5.1           curl_7.1.0               BiocNeighbors_2.4.0      TFisher_0.2.0            DelayedArray_0.36.1
#  [99] plotly_4.12.0            stringfish_0.19.0        rtracklayer_1.70.1       lmtest_0.9-40            digest_0.6.39            goftest_1.2-3            presto_1.0.0
# [106] spatstat.utils_3.2-2     XVector_0.50.0           htmltools_0.5.9          pkgconfig_2.0.3          fastmap_1.2.0            GlobalOptions_0.1.4      rlang_1.2.0
# [113] htmlwidgets_1.6.4        UCSC.utils_1.6.1         shiny_1.13.0             farver_2.1.2             zoo_1.8-15               jsonlite_2.0.0           BiocSingular_1.26.1
# [120] RCurl_1.98-1.18          magrittr_2.0.5           Formula_1.2-5            dotCall64_1.2            Rcpp_1.1.1-1.1           reticulate_1.46.0        stringi_1.8.7
# [127] MASS_7.3-65              plyr_1.8.9               listenv_0.10.1           ggrepel_0.9.8            openai_0.4.1             deldir_2.0-4             Biostrings_2.78.0
# [134] graphlayouts_1.2.3       splines_4.5.1            multtest_2.66.0          tensor_1.5.1             circlize_0.4.18          hms_1.1.4                qqconf_1.3.2
# [141] locfit_1.5-9.12          igraph_2.3.0             spatstat.geom_3.7-3      ggsignif_0.6.4           RcppHNSW_0.6.0           reshape2_1.4.5           ScaledMatrix_1.18.0
# [148] XML_3.99-0.23            metap_1.13               RcppParallel_5.1.11-2    scran_1.38.1             renv_1.2.2               BiocManager_1.30.27      GPTCelltype_1.0.1
# [155] ggprism_1.0.7            tzdb_0.5.0               tweenr_2.0.3             httpuv_1.6.17            RANN_2.6.2               polyclip_1.10-7          future_1.70.0
# [162] ggforce_0.5.0            janitor_2.2.1            rsvd_1.0.5               broom_1.0.12             xtable_1.8-8             restfulr_0.0.16          RSpectra_0.16-2
# [169] rstatix_0.7.3            later_1.4.8              memoise_2.0.1            beeswarm_0.4.0           GenomicAlignments_1.46.0 cluster_2.1.8.2          timechange_0.4.0
# [176] globals_0.19.1           here_1.0.2
