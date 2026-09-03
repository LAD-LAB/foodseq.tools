#' @title Scatter Plot of Read Counts for Two Taxa
#'
#' @description Plots the total read counts of two user-specified taxa against
#'   each other across all samples, with log10-scaled axes. When multiple
#'   ASVs match a taxon name (e.g., multiple species in a genus), their reads
#'   are summed per sample. Optionally colors points and adds 95\% ellipses
#'   by a metadata variable.
#'
#' @param ps Phyloseq object.
#' @param taxon_x Character; name of the taxon to plot on the x-axis.
#' @param taxon_y Character; name of the taxon to plot on the y-axis.
#' @param tax_rank Character; name of the tax_table column to search for
#'   \code{taxon_x} and \code{taxon_y} (e.g., \code{"species"},
#'   \code{"genus"}).
#' @param color_var Character; name of a sample_data column to use for
#'   coloring points. If \code{NULL} (default), all points are black.
#'
#' @return A ggplot2 scatter plot with log10-scaled axes.
#'
#' @export
plot_taxa_correlation <- function(ps,
                                  taxon_x   = "",
                                  taxon_y   = "",
                                  tax_rank  = "",
                                  color_var = NULL) {

  if (!nzchar(tax_rank)) stop("'tax_rank' must be specified.")
  if (!nzchar(taxon_x))  stop("'taxon_x' must be specified.")
  if (!nzchar(taxon_y))  stop("'taxon_y' must be specified.")

  get_asvs <- function(taxon_name) {
    tax <- as.data.frame(phyloseq::tax_table(ps))
    if (!tax_rank %in% colnames(tax)) {
      stop("'tax_rank' column '", tax_rank, "' not found in tax_table.")
    }
    rownames(tax)[!is.na(tax[[tax_rank]]) & tax[[tax_rank]] == taxon_name]
  }

  asv_x <- get_asvs(taxon_x)
  asv_y <- get_asvs(taxon_y)

  if (length(asv_x) == 0) stop("No ASVs found for taxon_x: '", taxon_x, "'")
  if (length(asv_y) == 0) stop("No ASVs found for taxon_y: '", taxon_y, "'")

  otu <- as.data.frame(phyloseq::otu_table(ps))
  if (phyloseq::taxa_are_rows(ps)) otu <- as.data.frame(t(otu))

  # NB: data.frame(), not as.data.frame() -- phyloseq's sample_data class
  # masquerades as a data.frame via an S4 .S3Class slot trick, so
  # as.data.frame() is a no-op that silently breaks rownames_to_column()
  # below (it ends up dropping all but a couple of rows).
  meta <- data.frame(phyloseq::sample_data(ps)) %>%
    tibble::rownames_to_column("Sample_ID")

  df <- otu %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::mutate(
      X_reads = rowSums(dplyr::select(., dplyr::all_of(asv_x))),
      Y_reads = rowSums(dplyr::select(., dplyr::all_of(asv_y)))
    ) %>%
    dplyr::select(Sample_ID, X_reads, Y_reads) %>%
    dplyr::left_join(meta, by = "Sample_ID")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = X_reads, y = Y_reads))

  if (!is.null(color_var) && color_var %in% colnames(df)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), size = 2, alpha = 0.75) +
      ggplot2::stat_ellipse(ggplot2::aes(color = .data[[color_var]]),
                   level = 0.95, type = "t", linetype = "dashed") +
      ggplot2::labs(color = color_var)
  } else {
    p <- p + ggplot2::geom_point(color = "black", size = 2, alpha = 0.75)
  }

  p +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x     = paste("Reads of", taxon_x),
      y     = paste("Reads of", taxon_y),
      title = paste("Correlation:", taxon_x, "vs.", taxon_y)
    ) +
    ggplot2::theme_minimal()
}
