library(AnnotationHub)

dir <- "/home/aguilada/aguilada-colacino/scRNAseq/chief_cells_analysis/SingleCell2016/data"

# Connect to AnnotationHub
ah = AnnotationHub::AnnotationHub()
# snapshotDate(): 2025-10-29

# Access the Ensembl database for organism
ahDb = query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)

# Acquire the latest annotation files
id = ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb = ah[[id]]
# rdatadateadded: 2024-10-28

# Extract gene-level information from database
annotations = genes(edb, return.type = "data.frame")

# Select annotations of interest
annotations = annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

write.csv(
  annotations,
  file.path(dir, 'annotationhub_mice_genes.csv'),
  row.names = FALSE
)

rm(ah, ahDb, edb)
gc()


### sessionInfo() ----
sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 8.10 (Ootpa)

# Matrix products: default
# BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRblas.so
# LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.5.1/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8
#  [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

# time zone: America/Detroit
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] ensembldb_2.34.0        AnnotationFilter_1.34.0 GenomicFeatures_1.62.0  AnnotationDbi_1.72.0    Biobase_2.70.0          GenomicRanges_1.62.1    Seqinfo_1.0.0           IRanges_2.44.0
#  [9] S4Vectors_0.48.1        AnnotationHub_4.0.0     BiocFileCache_3.0.0     dbplyr_2.5.2            BiocGenerics_0.56.0     generics_0.1.4

# loaded via a namespace (and not attached):
#  [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 rjson_0.2.23                httr2_1.2.2                 lattice_0.22-9              vctrs_0.7.3
#  [7] tools_4.5.1                 bitops_1.0-9                curl_7.1.0                  parallel_4.5.1              tibble_3.3.1                RSQLite_2.4.6
# [13] blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-5                cigarillo_1.0.0             lifecycle_1.0.5             compiler_4.5.1
# [19] Rsamtools_2.26.0            Biostrings_2.78.0           codetools_0.2-20            GenomeInfoDb_1.46.2         lazyeval_0.2.3              RCurl_1.98-1.18
# [25] yaml_2.3.12                 pillar_1.11.1               crayon_1.5.3                BiocParallel_1.44.0         DelayedArray_0.36.1         cachem_1.1.0
# [31] abind_1.4-8                 tidyselect_1.2.1            dplyr_1.2.1                 purrr_1.2.2                 restfulr_0.0.16             BiocVersion_3.22.0
# [37] grid_4.5.1                  fastmap_1.2.0               SparseArray_1.10.10         cli_3.6.6                   magrittr_2.0.5              S4Arrays_1.10.1
# [43] XML_3.99-0.23               withr_3.0.2                 UCSC.utils_1.6.1            filelock_1.0.3              rappdirs_0.3.4              bit64_4.8.0
# [49] XVector_0.50.0              httr_1.4.8                  matrixStats_1.5.0           bit_4.6.0                   otel_0.2.0                  png_0.1-9
# [55] memoise_2.0.1               BiocIO_1.20.0               rtracklayer_1.70.1          rlang_1.2.0                 glue_1.8.1                  DBI_1.3.0
# [61] BiocManager_1.30.27         renv_1.2.2                  jsonlite_2.0.0              R6_2.6.1                    ProtGenerics_1.42.0         MatrixGenerics_1.22.0
# [67] GenomicAlignments_1.46.0
