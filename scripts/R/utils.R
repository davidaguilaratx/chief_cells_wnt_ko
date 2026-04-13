library(parallel)
library(BiocParallel)
library(Seurat)
library(yaml)
library(scDblFinder)
library(SingleCellExperiment)

# ============================================================
# Parallelism Configuration ----
# ============================================================

#' Configure parallelism based on OS and compute environment
#'
#' Detects whether the session is running on Windows, Linux/macOS,
#' or an HPC (via SLURM environment variables) and returns an
#' appropriate BiocParallel BPPARAM object and thread count.
#'
#' @return A named list with:
#'   - `nthreads`: integer number of workers
#'   - `bpparam`:  a BiocParallelParam object
configure_parallelism <- function(rng_seed = NULL, max_workers = 8L) {
  os <- Sys.info()[["sysname"]]
  is_hpc <- nchar(Sys.getenv("SLURM_JOB_ID")) > 0
  max_hpc_threads <- max_workers # diminishing returns beyond ~8 for BiocParallel / qs2 I/O?

  if (is_hpc) {
    # Prefer SLURM_CPUS_PER_TASK (--cpus-per-task), fall back to SLURM_NTASKS
    slurm_cpus <- suppressWarnings(as.integer(Sys.getenv(
      "SLURM_CPUS_ON_NODE",
      unset = NA
    )))
    slurm_tasks <- suppressWarnings(as.integer(Sys.getenv(
      "SLURM_NTASKS",
      unset = NA
    )))
    detected <- Filter(Negate(is.na), c(slurm_cpus, slurm_tasks, 1L))[[1]]
    nthreads <- min(max(1L, detected - 2L), max_hpc_threads)
    bpparam <- BiocParallel::MulticoreParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  } else if (os == "Windows") {
    # Forking is unavailable on Windows; use SNOW-based parallelism
    physical_cores <- parallel::detectCores(logical = FALSE)
    nthreads <- max(1L, floor(physical_cores / 4L))
    bpparam <- BiocParallel::SnowParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  } else {
    # Linux / macOS — fork-based, leave 2 cores free for the OS
    physical_cores <- parallel::detectCores(logical = FALSE)
    nthreads <- max(1L, physical_cores - 2L)
    bpparam <- BiocParallel::MulticoreParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  }

  message(sprintf(
    "[configure_parallelism] OS: %s | HPC: %s | workers: %d | backend: %s",
    os,
    is_hpc,
    nthreads,
    class(bpparam)
  ))

  list(nthreads = nthreads, bpparam = bpparam)
}

# ============================================================
# Get yaml parameter ----
# ============================================================

#' Retrieve a parameter from a nested YAML config with timepoint-level override
#'
#' Looks up \code{param_name} first under the timepoint-specific block, then
#' falls back to the \code{default} block. Errors loudly if neither exists.
#'
#' @param configs    Named list parsed from a YAML config file via \code{yaml::read_yaml()}
#' @param param_name Character string naming the parameter to retrieve
#' @param timepoint  Character string matching a top-level key in \code{configs}
#'
#' @return The value of \code{param_name} for the given timepoint, or the default
#'
#' @examples
#' lower <- get_param(configs, "emptydrops", "lower_umi_bound", "month7")
get_param <- function(
  configs,
  section,
  param_name,
  timepoint,
  fallback = NULL
) {
  value <- configs[[section]][[timepoint]][[param_name]] %||%
    configs[[section]]$default[[param_name]] %||%
    configs$default[[param_name]] %||%
    fallback
  if (is.null(value)) {
    stop(
      sprintf(
        "[get_param] Parameter '%s' not found in section '%s' ",
        param_name,
        section
      ),
      sprintf("for timepoint '%s' or in defaults.", timepoint)
    )
  }
  value
}


# ============================================================
# Quality metrics calculations ----
# ============================================================

calc_metrics <- function(seurat_object) {
  # ==========================================================
  # 1. Vectorized Metadata Calculations (MASSIVELY FASTER)
  # ==========================================================
  # Pull out metadata to avoid repeated S4 object copies
  meta <- seurat_object@meta.data
  original_cells <- rownames(meta)

  # Calculate Z-scores and complexity natively and vectorized
  meta <- meta |>
    dplyr::group_by(sample) |>
    dplyr::mutate(
      log10_genes = log10(nFeature_RNA),
      log10_counts = log10(nCount_RNA),
      log10_genes_zscored = (log10_genes - mean(log10_genes, na.rm = TRUE)) /
        sd(log10_genes, na.rm = TRUE),
      log10_counts_zscored = (log10_counts - mean(log10_counts, na.rm = TRUE)) /
        sd(log10_counts, na.rm = TRUE),
      complexity = log10_genes / log10_counts
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-log10_genes, -log10_counts) |>
    as.data.frame()

  # Restore rownames and inject back into the Seurat object
  rownames(meta) <- original_cells
  seurat_object@meta.data <- meta

  # Clean up memory
  rm(meta, original_cells)
  gc()

  # ==========================================================
  # 2. Percentage Feature Sets
  # ==========================================================
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^mt-",
    col.name = "percent.mt"
  )
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^Rp[sl]",
    col.name = "percent.rb"
  )
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^Hb[^(P|E|S)]",
    col.name = "percent.hb"
  )

  # ==========================================================
  # 3. Cumulative proportions of top genes (Memory-Optimized)
  # ==========================================================
  layer_names <- Layers(seurat_object, assay = "RNA", search = "counts")
  qc_list <- vector("list", length(layer_names))
  names(qc_list) <- layer_names

  for (lyr in layer_names) {
    # Extract just this layer
    counts_matrix <- LayerData(seurat_object, assay = "RNA", layer = lyr)

    # Calculate stats
    qc_stats <- scuttle::perCellQCMetrics(
      counts_matrix,
      percent_top = c(20, 50, 100, 200, 500)
    )
    qc_list[[lyr]] <- as.data.frame(qc_stats)

    # Free memory immediately
    rm(counts_matrix, qc_stats)
    gc()
  }

  # Bind and assign efficiently to Seurat metadata
  all_qc_stats <- do.call(rbind, unname(qc_list))

  seurat_object$pct_counts_in_top_20_genes <- all_qc_stats$percent.top_20
  seurat_object$pct_counts_in_top_50_genes <- all_qc_stats$percent.top_50
  seurat_object$pct_counts_in_top_100_genes <- all_qc_stats$percent.top_100
  seurat_object$pct_counts_in_top_200_genes <- all_qc_stats$percent.top_200
  seurat_object$pct_counts_in_top_500_genes <- all_qc_stats$percent.top_500
  seurat_object$sum_counts <- all_qc_stats$sum
  seurat_object$n_genes_detected <- all_qc_stats$detected

  rm(all_qc_stats, qc_list)
  gc()

  # ==========================================================
  # 4. MALAT1 Thresholding
  # ==========================================================
  if (any(grepl("^Malat1", rownames(seurat_object), ignore.case = TRUE))) {
    has_norm <- "data" %in%
      slotNames(seurat_object[["RNA"]]) &&
      ncol(seurat_object[["RNA"]]@data) > 0

    if (has_norm) {
      message(
        "[calc_metrics] RNA normalized data already present; using it for MALAT1 threshold."
      )
    } else {
      message(
        "[calc_metrics] RNA normalized data missing; calling NormalizeData() before MALAT1 threshold."
      )
      seurat_object <- NormalizeData(seurat_object)
    }

    seurat_object <- Add_MALAT1_Threshold(
      object = seurat_object,
      species = "mouse",
      sample_col = "sample",
      save_plots = TRUE,
      save_plot_path = fig_dir,
      save_plot_name = "MALAT1_Threshold_Plots"
    )
  }

  return(seurat_object)
}


# ============================================================
# Gene-level filtering ----
# ============================================================

# according to best practices, want to filter out genes that
# are not detected in at least 20 cells.
# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

filter_genes <- function(seurat_obj, min_cells = 3L) {
  lyrs <- Layers(seurat_obj, assay = "RNA", search = "counts")
  gene_ncells <- Reduce(
    "+",
    lapply(lyrs, function(lyr) {
      m <- LayerData(seurat_obj, assay = "RNA", layer = lyr)
      tabulate(m@i + 1L, nbins = nrow(m))
    })
  )
  keep <- gene_ncells >= min_cells
  message(sprintf(
    "[filter_genes] %d genes before | %d after | %d removed (min_cells = %d)",
    nrow(seurat_obj),
    sum(keep),
    sum(!keep),
    min_cells
  ))
  seurat_obj[keep, ]
}


# ============================================================
# Doublet detection ----
# ============================================================

#' Two-round scDblFinder doublet detection
#'
#' Runs scDblFinder once on the full object, then re-runs it on the
#' round-1 singlets. Writes round-specific score/class columns and a
#' combined \code{doublet_final} label to the Seurat metadata.
#'
#' For 10X Flex / multiplexed data where biological samples share one
#' physical capture, we do NOT pass \code{samples} to scDblFinder —
#' cross-sample doublets are real and should be detected. This also
#' avoids scDblFinder's per-sample worker fan-out, which serializes
#' large SCE chunks across workers and has been the source of OOM /
#' "error writing to connection" failures on this dataset.
#'
#' @param seurat_obj  Seurat object (post QC, gene-filtered)
#' @param bpparam     BiocParallelParam used for scDblFinder's internal
#'                    PCA / KNN steps only.
#' @return the input Seurat object with new metadata columns:
#'   scDblFinder.score_r1, scDblFinder.class_r1,
#'   scDblFinder.score_r2, scDblFinder.class_r2,
#'   doublet_final ("singlet" | "doublet")
run_doublet_detection <- function(seurat_obj, bpparam, samples = NULL) {
  run_once <- function(so, label) {
    message(sprintf("[doublet] round %s: %d cells", label, ncol(so)))
    sce <- as.SingleCellExperiment(so)
    # Ensure a clean backend for this call — avoids stale worker sockets
    # from the previous round exhausting file descriptors.
    if (BiocParallel::bpisup(bpparam)) {
      BiocParallel::bpstop(bpparam)
    }
    BiocParallel::bpstart(bpparam)
    on.exit(
      if (BiocParallel::bpisup(bpparam)) BiocParallel::bpstop(bpparam),
      add = TRUE
    )
    sce <- scDblFinder(
      sce,
      clusters = TRUE,
      samples = samples,
      dbr.sd = 1,
      BPPARAM = bpparam
    )
    out <- data.frame(
      score = sce$scDblFinder.score,
      class = as.character(sce$scDblFinder.class),
      row.names = colnames(sce),
      stringsAsFactors = FALSE
    )
    rm(sce)
    gc()
    out
  }

  # round 1 — all cells
  r1 <- run_once(seurat_obj, "1")
  seurat_obj$scDblFinder.score_r1 <- r1[colnames(seurat_obj), "score"]
  seurat_obj$scDblFinder.class_r1 <- r1[colnames(seurat_obj), "class"]

  # round 2 — only round-1 singlets
  singlet_cells <- rownames(r1)[r1$class == "singlet"]
  r2 <- run_once(seurat_obj[, singlet_cells], "2")

  seurat_obj$scDblFinder.score_r2 <- NA_real_
  seurat_obj$scDblFinder.class_r2 <- NA_character_
  seurat_obj$scDblFinder.score_r2[match(
    rownames(r2),
    colnames(seurat_obj)
  )] <- r2$score
  seurat_obj$scDblFinder.class_r2[match(
    rownames(r2),
    colnames(seurat_obj)
  )] <- r2$class

  # final label: doublet in either round
  is_doublet <- seurat_obj$scDblFinder.class_r1 == "doublet" |
    (!is.na(seurat_obj$scDblFinder.class_r2) &
      seurat_obj$scDblFinder.class_r2 == "doublet")
  seurat_obj$doublet_final <- factor(
    ifelse(is_doublet, "doublet", "singlet"),
    levels = c("singlet", "doublet")
  )

  message(sprintf(
    "[doublet] summary: singlet=%d | doublet=%d (R1=%d, R2=%d)",
    sum(!is_doublet),
    sum(is_doublet),
    sum(seurat_obj$scDblFinder.class_r1 == "doublet"),
    sum(
      seurat_obj$scDblFinder.class_r1 == "singlet" &
        !is.na(seurat_obj$scDblFinder.class_r2) &
        seurat_obj$scDblFinder.class_r2 == "doublet"
    )
  ))

  seurat_obj
}
