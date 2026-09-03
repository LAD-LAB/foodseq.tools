#' @title Broken Stick Method for PC Retention
#'
#' @description Applies the broken stick model (MacArthur 1957) to determine
#'   how many principal components to retain. A PC is retained when its
#'   eigenvalue exceeds the corresponding broken stick expectation under a
#'   random null model. The cutoff is overlaid on the scree plot returned by
#'   \code{\link{pca_plot}}.
#'
#' @param pca_output Output list returned by \code{pca_plot()}.
#'
#' @return A ggplot2 object: the scree plot from \code{pca_output} with a
#'   dashed red vertical line marking the suggested retention cutoff.
#'
#' @export
bstick_pc <- function(pca_output) {

  eigenvalues <- pca_output$scree.table %>%
    dplyr::pull(Eigenvalue)

  bs_values <- vegan::bstick(n = length(eigenvalues))
  bs_k      <- sum(eigenvalues > bs_values)

  pca_output$scree.plot +
    ggplot2::geom_vline(xintercept = bs_k, color = "red", linetype = "dashed") +
    ggplot2::annotate("text",
             x     = bs_k + 0.3,
             y     = max(eigenvalues) * 0.9,
             label = paste0("Retain PCs 1-", bs_k),
             hjust = 0, color = "red")
}
