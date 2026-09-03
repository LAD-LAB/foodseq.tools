#' @title Alpha Diversity
#'
#' @description Computes per-sample alpha diversity measures (richness,
#'   Shannon, Simpson, etc.) via \code{\link[phyloseq]{estimate_richness}},
#'   merges them with the phyloseq object's sample_data for immediate use,
#'   and optionally builds a grouped boxplot. Requires raw (or otherwise
#'   non-negative) counts -- alpha diversity indices are undefined on
#'   CLR-transformed data (see \code{\link{clr_transform}}), which contains
#'   negative values by construction; run this \emph{before} CLR-transforming
#'   a phyloseq object for ordination, not after.
#'
#' @param ps A phyloseq object with a non-negative (typically raw count)
#'   otu_table.
#' @param measures Character vector of diversity measures to compute; any
#'   subset of \code{"Observed"}, \code{"Chao1"}, \code{"ACE"},
#'   \code{"Shannon"}, \code{"Simpson"}, \code{"InvSimpson"}, \code{"Fisher"}
#'   (passed through to \code{\link[phyloseq]{estimate_richness}}). Default
#'   \code{c("Observed", "Shannon", "Simpson")}.
#' @param group_var Optional character; name of a sample_data column to
#'   group samples by in the returned plot (e.g. a treatment or timepoint
#'   variable). If \code{NULL} (default), no plot is built.
#' @param plot If \code{FALSE}, skip building a plot even when
#'   \code{group_var} is supplied. Default \code{TRUE}.
#'
#' @return A named list with two elements: \code{data} (a data frame with
#'   one row per sample, the requested diversity measures, and all
#'   sample_data columns), and \code{plot} (a ggplot2 object -- one panel
#'   per measure, boxplot + jittered points by \code{group_var} -- or
#'   \code{NULL} if \code{group_var} is \code{NULL} or \code{plot = FALSE}).
#'
#' @export
alpha_diversity <- function(ps, measures = c("Observed", "Shannon", "Simpson"),
                            group_var = NULL, plot = TRUE) {

  stopifnot(inherits(ps, "phyloseq"))

  otu <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (any(otu < 0, na.rm = TRUE)) {
    stop("alpha_diversity: otu_table contains negative values, which is not ",
         "valid for diversity indices. Did you pass CLR-transformed data ",
         "(see clr_transform())? Compute alpha diversity on raw counts instead.")
  }

  rich <- phyloseq::estimate_richness(ps, measures = measures)
  rich <- tibble::rownames_to_column(rich, var = "Sample")

  # NB: data.frame(), not as.data.frame() -- phyloseq's sample_data class
  # masquerades as a data.frame via an S4 .S3Class slot trick, so
  # as.data.frame() is effectively a no-op that leaves the S4 wrapper in
  # place and silently breaks downstream tidyverse calls (rownames_to_column
  # in particular ends up dropping all but a couple of rows).
  samdf <- data.frame(phyloseq::sample_data(ps, errorIfNULL = FALSE))
  if (!is.null(samdf) && ncol(samdf) > 0) {
    samdf <- tibble::rownames_to_column(samdf, var = "Sample")
    data <- dplyr::left_join(rich, samdf, by = "Sample")
  } else {
    data <- rich
  }

  p <- NULL
  if (plot && !is.null(group_var)) {
    if (!group_var %in% colnames(data)) {
      warning("alpha_diversity: group_var '", group_var,
              "' not found in sample_data; skipping plot.")
    } else {
      long <- tidyr::pivot_longer(data, dplyr::all_of(measures),
                                  names_to = "measure", values_to = "value")
      p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[group_var]], y = value)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.6) +
        ggplot2::facet_wrap(~measure, scales = "free_y") +
        ggplot2::labs(x = group_var, y = "Value", title = "Alpha diversity") +
        ggplot2::theme_minimal()
    }
  }

  list(data = data, plot = p)
}
