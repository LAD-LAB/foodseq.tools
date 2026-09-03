#' @title Elbow Method for PC Retention
#'
#' @description Identifies the "elbow" in a PCA scree plot using the
#'   perpendicular-distance method: for each PC, the orthogonal distance from
#'   the point to the line connecting PC1 and the last PC is computed. The
#'   elbow is the PC with the maximum distance, i.e., the point of maximum
#'   curvature. The result is overlaid on the scree plot from
#'   \code{\link{pca_plot}}.
#'
#' @param pca_output Output list returned by \code{pca_plot()}.
#'
#' @return A ggplot2 object: the scree plot from \code{pca_output} with a
#'   dashed red vertical line at the elbow PC.
#'
#' @export
elbow_pc <- function(pca_output) {

  eigenvalues <- pca_output$scree.table %>%
    dplyr::pull(Eigenvalue)

  n      <- length(eigenvalues)
  points <- cbind(seq_len(n), eigenvalues)

  # Unit vector along the line from first to last point
  line_vec <- points[n, ] - points[1, ]
  line_vec <- line_vec / sqrt(sum(line_vec^2))

  # Perpendicular distance from each point to that line
  distances <- apply(points, 1, function(pt) {
    vec  <- pt - points[1, ]
    proj <- sum(vec * line_vec) * line_vec
    sqrt(sum((vec - proj)^2))
  })

  elbow <- which.max(distances)

  pca_output$scree.plot +
    ggplot2::geom_vline(xintercept = elbow, color = "red", linetype = "dashed") +
    ggplot2::annotate("text",
             x     = elbow + 0.3,
             y     = max(eigenvalues) * 0.9,
             label = paste0("Elbow = PC", elbow),
             hjust = 0, color = "red")
}
