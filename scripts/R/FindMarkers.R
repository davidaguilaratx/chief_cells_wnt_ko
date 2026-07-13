library(Seurat)

Idents(filtered_seurat) <- filtered_seurat$RNA_snn_res.0.15
all.pos.markers = FindAllMarkers(
  filtered_seurat,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf
)

# Combine markers with gene descriptions
annotations = read.csv(
  '/home/aguilada/aguilada-colacino/scRNAseq/analysis/data/annotationhub_mice_genes.csv'
)
  
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
  group_by(cluster) %>%
  filter(p_val_adj < 0.05, pct.diff > 0.25) %>%
  slice_max(order_by = tibble(avg_log2FC), n = 20)
write.csv(
  top20,
  file.path(data_dir, 'top20_markers_res.0.15_filtered_seurat.csv'),
  row.names = FALSE
)


## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster, assay = 'RNA', group.by = 'condition') {
  FindConservedMarkers(
    integrated_seurat, # need to specify exact seurat object here
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
    avg.pct.diff = ((wnt_ko_pct.1 - wnt_ko_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2,
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
  paste0(
    annot_dir_integrated,
    'conserved_pos_markers_across_condition_res.0.15_filtered_seurat.csv'
  ),
  row.names = FALSE
)


# Extract top 20 conserved markers per cluster
top20.conserved = conserved.pos.markers %>%
  group_by(cluster_id) %>%
  filter(max_adj_pval < 0.05, avg.pct.diff > 0.25) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  # slice_max(order_by=tibble(avg.pct.diff),n=20)
  slice_max(order_by = tibble(avg_log2FC, avg.pct.diff), n = 20)
write.csv(
  top20.conserved,
  paste0(
    annot_dir_integrated,
    'top20.conserved_pos_markers_across_condition_harmonyres.0.25_df.csv'
  ),
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
  topgenenumber = 10 # pick anywhere from 5-25
)
res

# look at conserved markers across conditions
res <- GPTCelltype::gptcelltype(
  conserved.pos.markers,
  tissuename = 'mouse stomach with KO beta-catenin in chief cells',
  model = 'gpt-5.4',
  topgenenumber = 20 # pick anywhere from 5-25
)
res
