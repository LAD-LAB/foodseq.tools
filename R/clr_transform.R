#' @title Centered Log-Ratio (CLR) Transform a Phyloseq Object
#'
#' @description Applies the centered log-ratio (CLR) transform to a
#'   phyloseq object's count table, the standard compositional-data
#'   preprocessing step expected by \code{\link{pca_plot}} and
#'   \code{\link{paran_pc}} (both of which document "CLR-transformed and
#'   filtered data" as their input but, until now, left every caller to
#'   implement that step themselves). For each sample, CLR replaces the raw
#'   counts \eqn{x} with \eqn{\log(x + pseudocount) - mean(\log(x +
#'   pseudocount))}, where the mean is taken across taxa within that sample.
#'
#' @param ps A phyloseq object with a count-based otu_table (not yet
#'   transformed).
#' @param pseudocount Numeric pseudocount added before taking logs, to
#'   avoid \code{log(0)}. If \code{NULL} (default), uses half the smallest
#'   nonzero count in the table -- a common convention (e.g.
#'   \code{microbiome::transform(ps, "clr")}) that scales with the data
#'   rather than imposing an arbitrary constant like the traditional "add
#'   1". Pass a fixed value (e.g. \code{1}) to use that convention instead.
#' @param drop_zero_variance If \code{TRUE} (default), drop taxa that are
#'   constant across all samples (most commonly all-zero) before
#'   transforming -- these carry no information for downstream ordination
#'   and, left in, can produce a rank-deficient CLR matrix.
#' @param verbose If \code{TRUE} (default), report the pseudocount used and
#'   how many taxa were dropped for zero variance.
#'
#' @return The phyloseq object with its otu_table replaced by the
#'   CLR-transformed values (a matrix of the same dimensions, minus any
#'   zero-variance taxa dropped).
#'
#' @export
clr_transform <- function(ps, pseudocount = NULL, drop_zero_variance = TRUE,
                          verbose = TRUE) {

  stopifnot(inherits(ps, "phyloseq"))

  rows_mode <- phyloseq::taxa_are_rows(ps)
  otu       <- methods::as(phyloseq::otu_table(ps), "matrix")
  mat       <- if (rows_mode) t(otu) else otu   # samples x taxa, from here on

  if (drop_zero_variance) {
    keep <- apply(mat, 2, function(col) length(unique(col)) > 1)
    n_dropped <- sum(!keep)
    if (n_dropped > 0) {
      if (verbose) message(n_dropped, " zero-variance taxa dropped before CLR transform.")
      mat <- mat[, keep, drop = FALSE]
    }
  }

  if (!any(mat > 0)) {
    stop("clr_transform: otu_table has no nonzero values -- nothing to transform.")
  }

  if (is.null(pseudocount)) {
    pseudocount <- min(mat[mat > 0]) / 2
    if (verbose) message(sprintf("Using pseudocount = %.6g (half the smallest nonzero count).",
                                 pseudocount))
  } else if (verbose) {
    message(sprintf("Using user-supplied pseudocount = %.6g.", pseudocount))
  }

  n_empty_samples <- sum(rowSums(mat) == 0)
  if (n_empty_samples > 0) {
    warning(n_empty_samples, " sample(s) have zero total reads; their CLR values ",
            "will be a degenerate constant (all zeros) and should likely be filtered ",
            "out before ordination.")
  }

  log_mat <- log(mat + pseudocount)
  clr_mat <- log_mat - rowMeans(log_mat)

  out_mat <- if (rows_mode) t(clr_mat) else clr_mat
  phyloseq::otu_table(ps) <- phyloseq::otu_table(out_mat, taxa_are_rows = rows_mode)

  ps
}
