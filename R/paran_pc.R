#' @title Parallel Analysis for PC Retention
#'
#' @description Determines the number of significant principal components via
#'   permutation-based parallel analysis. For each of \code{B} permutations,
#'   each column of the OTU table is independently shuffled and eigenvalues
#'   are computed. A PC is considered significant when its observed eigenvalue
#'   exceeds the specified percentile of the permuted distribution. A QC plot
#'   compares internal eigenvalues to those stored in the \code{pca_plot()}
#'   output to verify consistency.
#'
#' @param ps Filtered and CLR-transformed phyloseq object on which
#'   \code{\link{pca_plot}} was run.
#' @param pca_output Output list returned by \code{pca_plot()}.
#' @param B Integer; number of permutations. Default \code{1000}.
#' @param centile Numeric; percentile threshold (0-100) for significance.
#'   Default \code{95}.
#' @param seed Integer; random seed for reproducibility. Default \code{123}.
#'
#' @return A named list with six elements:
#'   \describe{
#'     \item{cutoff_pc}{Integer; last PC whose eigenvalue exceeds the random
#'       threshold (i.e., retain PCs 1 through \code{cutoff_pc}).}
#'     \item{var_explained}{Numeric; proportion of total variance explained by
#'       PCs 1 through \code{cutoff_pc}.}
#'     \item{paran.df}{Long-format data frame used to build \code{paran.plot}.}
#'     \item{paran.plot}{ggplot2 parallel analysis scree plot with the
#'       significance cutoff marked.}
#'     \item{qc.df}{Data frame comparing internal vs. \code{pca_plot()}
#'       eigenvalues for quality control.}
#'     \item{qc.plot}{ggplot2 scatter plot of internal vs. \code{pca_plot()}
#'       eigenvalues (should fall on the diagonal if consistent).}
#'   }
#'
#' @export
paran_pc <- function(ps,
                     pca_output,
                     B         = 1000,
                     centile   = 95,
                     seed      = 123) {

  future::plan("multisession", workers = parallel::detectCores() - 1)
  set.seed(seed)

  seqtab <- as.matrix(phyloseq::otu_table(ps))
  if (phyloseq::taxa_are_rows(ps)) seqtab <- t(seqtab)

  n <- nrow(seqtab)
  p <- ncol(seqtab)

  # Observed eigenvalues
  obs_pca <- stats::prcomp(seqtab, center = TRUE, scale. = FALSE)
  obs_eig <- obs_pca$sdev^2

  # Permuted eigenvalues (column-wise shuffling)
  rand_eig <- future.apply::future_sapply(seq_len(B), function(b) {
    perm_data <- apply(seqtab, 2, sample)
    eigen(stats::cov(perm_data))$values
  })
  rand_eig <- t(rand_eig)  # rows = permutations, cols = PCs

  rand_means      <- apply(rand_eig, 2, mean)
  rand_percentile <- apply(rand_eig, 2, function(x) stats::quantile(x, probs = centile / 100))

  cutoff_pc        <- sum(obs_eig > rand_percentile)
  variance_explained <- sum(obs_eig[seq_len(cutoff_pc)]) / sum(obs_eig)

  paran.df <- data.frame(
    PC             = paste0("PC", seq_along(obs_eig)),
    PC_index       = seq_along(obs_eig),
    observed_eig   = obs_eig,
    rand_mean      = rand_means,
    rand_percentile = rand_percentile
  ) %>%
    tidyr::pivot_longer(
      cols      = c(observed_eig, rand_mean, rand_percentile),
      names_to  = "type",
      values_to = "eigenvalue"
    )

  paran.plot <- paran.df %>%
    ggplot2::ggplot(ggplot2::aes(x = PC_index, y = eigenvalue, color = type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::geom_vline(xintercept = cutoff_pc, color = "red", linetype = "dashed") +
    ggplot2::annotate("text",
             x     = cutoff_pc + 0.5,
             y     = max(obs_eig),
             label = paste0("Retain up to PC ", cutoff_pc,
                            " (", round(variance_explained * 100, 1), "% variance)"),
             hjust = 0, vjust = 1.2, color = "red") +
    ggplot2::labs(title = "Parallel Analysis (Permutation)",
         x     = "Principal Component",
         y     = "Eigenvalue")

  # QC: compare internal vs. pca_plot() eigenvalues
  qc.df <- paran.df %>%
    tidyr::pivot_wider(id_cols = PC, names_from = "type", values_from = "eigenvalue") %>%
    dplyr::left_join(pca_output$scree.table, by = "PC")

  qc.plot <- qc.df %>%
    ggplot2::ggplot(ggplot2::aes(x = observed_eig, y = Eigenvalue)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "QC: pca_plot() vs paran_pc() Eigenvalues",
         x = "paran_pc() observed eigenvalue",
         y = "pca_plot() Eigenvalue")

  print(paran.plot)

  return(list(
    cutoff_pc    = cutoff_pc,
    var_explained = variance_explained,
    paran.df     = paran.df,
    paran.plot   = paran.plot,
    qc.df        = qc.df,
    qc.plot      = qc.plot
  ))
}
