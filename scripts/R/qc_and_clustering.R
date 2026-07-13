library(Seurat)
library(sctransform)
library(tidyverse)
library(scuttle)
library(zeallot) # for unpacking multiple values from functions
library(qs2)
library(RhpcBLASctl)

seed = 281330800

set.seed(seed) # for reproducibility

# set up directories and parameters ----
base_dir <- "/home/aguilada/chief_cells_analysis/SingleCell2016"
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
pct = merged_seurat[["pca"]]@stdev / sum(merged_seurat[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(
    which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
    decreasing = T
)[1] +
    1
co2
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs

# Create a dataframe to visualize percent variation captured by PCs
plot_df = data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

#### percent variation captured by PCs ----
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_classic2() +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    labs(
        title = 'Eblow Plot of percent variation explained by PCs',
        x = 'Cumulative percentage',
        y = 'Percent of variation'
    ) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 18)
    )
ggsave(
    filename = file.path(fig_dir, 'pca_variation_explained.png'),
    width = 10,
    height = 8
)

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
# helpers to keep each FeaturePlot styled consistently
qc_feature_plot <- function(feature, max_cutoff = 'q90') {
    FeaturePlot(
        merged_seurat,
        reduction = 'umap',
        features = feature,
        min.cutoff = 'q10',
        max.cutoff = max_cutoff,
        cols = viridis(256)
    ) +
        DarkTheme()
}

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
    qc_feature_plot('cc.difference'),
    qc_feature_plot('S.Score'),
    qc_feature_plot('G2M.Score'),
    ncol = 2
)
save_panel(panel_cc, dir_umaps, 'cell_cycle_umaps.png')

## Panel 3: QC metrics on UMAP ----
panel_qc_metrics <- wrap_plots(
    qc_feature_plot('nCount_RNA'),
    qc_feature_plot('nFeature_RNA'),
    qc_feature_plot('percent.mt'),
    qc_feature_plot('percent.hb'),
    qc_feature_plot('complexity'),
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
panel_mt <- FeaturePlot(
    merged_seurat,
    reduction = 'umap',
    features = mt_genes,
    min.cutoff = 'q10',
    max.cutoff = 'q90',
    cols = viridis(256)
) &
    DarkTheme()
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
RhpcBLASctl::omp_set_num_threads(1) # covers xgboost, igraph, uwot, Matrix
RhpcBLASctl::blas_set_num_threads(1) # harmless no-op on reference BLAS
data.table::setDTthreads(1) # data.table has its own pool
RcppParallel::setThreadOptions(numThreads = 1) # TBB-based packages

filtered_seurat <- run_doublet_detection(
    filtered_seurat,
    bpparam = bp_param,
    # bpparam = MulticoreParam(
    #     workers = 4,
    #     progressbar = TRUE,
    #     RNGseed = seed
    # ),
    # bpparam = SerialParam(progressbar = TRUE, RNGseed = seed),
    samples = 'sample'
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
        group.by = "seurat_clusters",
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
        order = "doublet",
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
        order = "doublet",
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
        order = "doublet",
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

## QC violins split by doublet call ----
plot_doublets(filtered_seurat, fig_dir = fig_dir)

# subset to singlets after visual inspection ----
# slim down seurat object first to ony what we need downstream
filtered_seurat <- DietSeurat(
    filtered_seurat,
    assays = 'RNA',
    layers = c('counts')
)

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
resolutiontest_s = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8, 1.0)

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
    assay = 'SCT',
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
