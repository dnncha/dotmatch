#' DotMatch R Interface
#'
#' Provides convenient wrappers around the Python `dotmatch` package using
#' reticulate. This enables Bioconductor / tidyverse users to benefit from
#' DotMatch's deterministic short-DNA assignment, explicit ambiguity handling,
#' and high performance from R.
#'
#' @docType package
#' @name dotmatch
NULL

#' @import reticulate
NULL

.onLoad <- function(libname, pkgname) {
  # Ensure python dotmatch is available
  # Users should do: reticulate::py_install("dotmatch", pip = TRUE)
  # or use a conda env with dotmatch installed.
  invisible()
}

#' Get the dotmatch Python module
#'
#' @return The imported Python module.
#' @export
dotmatch <- function() {
  reticulate::import("dotmatch", delay_load = TRUE)
}

#' Compute edit distance
#'
#' @param a,b Character strings.
#' @return Integer distance.
#' @export
distance <- function(a, b) {
  dm <- dotmatch()
  dm$distance(a, b)
}

#' Assign sequences to a known library (Python level)
#'
#' Thin wrapper. For large data, prefer calling the native `dotmatch` CLI
#' from R via `system2` or `processx` and then reading the outputs with
#' the anndata / data.table helpers.
#'
#' @param reads,targets Character vectors.
#' @param k Maximum distance.
#' @param policy "radius" or "best".
#' @return data.frame of results (via pandas -> R).
#' @export
assign <- function(reads, targets, k = 1L, policy = "radius") {
  dm <- dotmatch()
  res <- dm$assign(reads, targets, k = k, policy = policy)
  # results_to_dataframe returns pandas; convert
  reticulate::py_to_r(res)
}

#' Read a DotMatch counts TSV into a SingleCellExperiment or SummarizedExperiment
#'
#' Requires the `anndata` Python extra (or use `read.delim` + manual construction).
#' For full scverse interop from R, consider `zellkonverter` + AnnData.
#'
#' @param path Path to counts.mageck.tsv or similar.
#' @return A list with counts matrix, row/col data (minimal).
#' @export
read_counts <- function(path) {
  dm <- dotmatch()
  # This will only work if anndata extra is installed in the reticulate env
  ad <- tryCatch(dm$counts_tsv_to_anndata(path), error = function(e) stop(e))
  # Convert AnnData to R-friendly structure
  list(
    counts = reticulate::py_to_r(ad$X),
    features = reticulate::py_to_r(ad$var),
    samples = reticulate::py_to_r(ad$obs)
  )
}

#' Assign features from an in-memory data.frame (cell x seq)
#'
#' Convenience for users who have already extracted windows into a data.frame.
#' Returns a data.frame with assignment columns. For full SingleCellExperiment
#' workflows, use this then construct your object.
#'
#' @param df data.frame with at least a column of sequences.
#' @param seq_col Name of the column containing sequences.
#' @param library Library (path or data.frame of id/seq).
#' @param k,policy Assignment parameters.
#' @return The input df with added columns (assigned_feature, feature_distance, feature_status).
#' @export
assign_features_df <- function(df, seq_col = "feature_seq", library, k = 1L, policy = "radius") {
  dm <- dotmatch()
  seqs <- as.character(df[[seq_col]])
  # Use the tl helper if available
  tl <- tryCatch(dm$tl, error = function(e) NULL)
  if (!is.null(tl)) {
    # tl.assign_features expects an AnnData; we do a roundtrip for convenience
    # For pure df use, fall back to core assign
  }
  res <- dm$assign_dataframe(seqs, library, k = k, policy = policy)
  res_df <- reticulate::py_to_r(res)
  df[[paste0(seq_col, "_assigned")]] <- res_df$target_name
  df[[paste0(seq_col, "_dist")]] <- res_df$best_distance
  df[[paste0(seq_col, "_status")]] <- res_df$status_name
  df
}

#' Levenshtein distance <= k
#'
#' @param a,b Character strings.
#' @param k Max edits.
#' @return Logical.
#' @export
distance_leq <- function(a, b, k = 1L) {
  dm <- dotmatch()
  dm$distance_leq(a, b, k)
}

#' Create a reusable Matcher for repeated assignment
#'
#' @param targets Character vector of target sequences.
#' @return An environment with $assign method (wrapper around Python Matcher).
#' @export
Matcher <- function(targets) {
  dm <- dotmatch()
  py_matcher <- dm$Matcher(targets)
  structure(list(
    py = py_matcher,
    assign = function(reads, k = 1L, policy = "radius") {
      res <- py_matcher$assign(reads, k = k, policy = policy)
      reticulate::py_to_r(res)
    },
    close = function() py_matcher$close()
  ), class = "dotmatch_matcher")
}

#' @export
print.dotmatch_matcher <- function(x, ...) {
  cat("DotMatch Matcher (reticulate wrapper)\n")
  invisible(x)
}
  # Merge back (simple cbind by order)
  df$assigned_feature <- res_df$target_name
  df$feature_distance <- res_df$best_distance
  df$feature_status <- res_df$status_name
  df
}
