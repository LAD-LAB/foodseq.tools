# Shared by pca_plot() (fits+plots a fresh PCA) and project_pca() (projects
# new data into an EXISTING PCA's fixed space) -- both end up needing the
# exact same thing: given a fitted PCA's rotation/sdev, a base ggplot
# (points already drawn), and a phyloseq object to pull taxon labels from,
# draw the top-nTaxa loading arrows and their quadrant-aware text labels.
# Pulled out of pca_plot() as a pure refactor (verified byte-identical
# output before/after) -- project_pca() never re-fits anything, so it
# always passes through the REFERENCE pca's rotation/sdev/labels here,
# never its own; the arrows/labels are the same fixed reference geometry
# either way.
.pca_biplot_layer <- function(base_plot, rotation, sdev, ps_labels,
                              xPC, yPC, nTaxa, bplab, colorName) {
  if (!is.null(bplab)) {
    # Relocate bplab to the end of tax table so it will be used for labeling
    phyloseq::tax_table(ps_labels) <- ps_labels@tax_table %>%
      data.frame() %>%
      dplyr::relocate(.data[[bplab]], .after = dplyr::everything()) %>%
      as.matrix() %>%
      phyloseq::tax_table()
  }

  # Calculate loadings
  V <- rotation             # Eigenvectors
  L <- diag(sdev)           # Diagonal matrix with square roots of eigenvalues
  loadings <- V %*% L
  colnames(loadings) <- colnames(V)  # Assign column names to loadings

  # Get loadings for specified PCs and format for plotting
  loadings.xy <- data.frame(loadings[, c(paste0('PC', xPC), paste0('PC', yPC))]) %>%
    dplyr::rename(PCx = paste0('PC', xPC), PCy = paste0('PC', yPC)) %>%
    dplyr::mutate(variable = row.names(loadings),
           length = sqrt(PCx^2 + PCy^2),
           ang = atan2(PCy, PCx) * (180 / pi))

  loadings.plot <- dplyr::top_n(loadings.xy, nTaxa, wt = length)

  # Adjust angles to keep labels upright
  loadings.plot <- loadings.plot %>%
    dplyr::mutate(adj_ang = ifelse(ang < -90, ang + 180,
                            ifelse(ang > 90, ang - 180, ang)))

  # Rename loadings with lowest taxonomic level
  loadings.taxtab <- phyloseq::tax_table(ps_labels)[row.names(loadings.plot)] %>%
    data.frame()
  loadings.taxtab <- loadings.taxtab[cbind(1:nrow(loadings.taxtab), max.col(!is.na(loadings.taxtab), ties.method = 'last'))] %>%
    data.frame()
  colnames(loadings.taxtab) <- c("name")
  loadings.taxtab$asv <- phyloseq::tax_table(ps_labels)[row.names(loadings.plot)] %>%
    data.frame() %>%
    rownames()

  loadings.plot <- loadings.taxtab %>%
    dplyr::select(asv, name) %>%
    dplyr::right_join(loadings.plot, by = c('asv' = 'variable'))

  # Determine the quadrant of each label
  q1 <- dplyr::filter(loadings.plot, PCx > 0 & PCy > 0)
  q2 <- dplyr::filter(loadings.plot, PCx < 0 & PCy > 0)
  q3 <- dplyr::filter(loadings.plot, PCx < 0 & PCy < 0)
  q4 <- dplyr::filter(loadings.plot, PCx > 0 & PCy < 0)

  pca.biplot <-
    base_plot +
    ggplot2::geom_segment(data = loadings.plot,
                 ggplot2::aes(x = 0, y = 0,
                     xend = PCx, yend = PCy),
                 color = 'black',
                 arrow = ggplot2::arrow(angle = 15,
                               length = ggplot2::unit(0.1, 'inches'))) +
    ggplot2::labs(color = colorName)

  # Add geom_text for each quadrant with adjusted angle and justification
  if (nrow(q1) != 0) {
    pca.biplot <- pca.biplot +
      ggplot2::geom_text(data = q1, ggplot2::aes(x = PCx, y = PCy, hjust = 0, vjust = 0, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q2) != 0) {
    pca.biplot <- pca.biplot +
      ggplot2::geom_text(data = q2, ggplot2::aes(x = PCx, y = PCy, hjust = 1, vjust = 0, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q3) != 0) {
    pca.biplot <- pca.biplot +
      ggplot2::geom_text(data = q3, ggplot2::aes(x = PCx, y = PCy, hjust = 1, vjust = 1, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }
  if (nrow(q4) != 0) {
    pca.biplot <- pca.biplot +
      ggplot2::geom_text(data = q4, ggplot2::aes(x = PCx, y = PCy, hjust = 0, vjust = 1, angle = adj_ang,
                               label = name,
                               fontface = 'bold'),
                color = 'black', show.legend = FALSE)
  }

  list(pca.biplot = pca.biplot, loadings = loadings)
}

#' @title PCA plot and biplot
#'
#' @description This function runs a Principal Component Analysis.
#'
#' @param ps clr transformed and filtered data
#' @param colorVar variable from samdf to color samples by
#' @param colorName what to display variable name as in legend
#' @param nTaxa number of taxa to display
#' @param customColors optional named vector of colors
#' @param customGradient optional df of low/high colors for gradient
#' @param mid for customGradient -- midpoint to median, mean, middle of range
#' @param xPC Principal Component for x-axis
#' @param yPC Principal Component for y-axis
#' @param ellipse optional add centroid ellipses
#' @param bplab name of column that you want to use for labeling the biplot
#'
#' @return A named list with six elements: \code{pca.df} (samdf with all PCs
#'   added as columns), \code{pca.biplot} (PCA biplot with \code{nTaxa}
#'   factor loadings displayed), \code{loadings} (data frame of loadings in
#'   PCA space), \code{pca.output} (raw output of \code{prcomp()}),
#'   \code{scree.table} (table of scree-plot numbers), and \code{scree.plot}
#'   (the scree plot).
#' @export
pca_plot <- function(ps, # clr transformed and filtered data
                    colorVar = NULL, # variable from samdf to color samples by
                    colorName = NULL, # what to display variable name as in legend
                    nTaxa = 10, # number of taxa to display
                    customColors = NULL, # optional named vector of colors
                    customGradient = NULL, # optional df of low/high colors for gradient
                    mid = "mean", # for customGradient -- midpoint to median, mean, middle of range
                    xPC = 1, # Principal Component for x-axis
                    yPC = 2,  # Principal Component for y-axis
                    ellipse = FALSE, # optional add centroid ellipses
                    bplab = NULL # optional variable name for biplot arrow labels
) {

  # Normalize orientation to samples-as-rows / taxa-as-columns, which
  # everything below assumes (prcomp() is run directly on ps@otu_table, and
  # pca$x's rownames are later joined against sample_data by sample name).
  # Without this, a taxa-as-rows phyloseq object -- an equally common
  # convention -- silently produces a transposed, nonsensical PCA.
  if (phyloseq::taxa_are_rows(ps)) {
    otu_samples_rows <- t(methods::as(phyloseq::otu_table(ps), "matrix"))
    phyloseq::otu_table(ps) <- phyloseq::otu_table(otu_samples_rows, taxa_are_rows = FALSE)
  }

  if ("name" %in% colnames(data.frame(ps@sam_data))) {
    # Prevent conflict with 'name' column
    phyloseq::sample_data(ps) <- ps@sam_data %>%
      data.frame() %>%
      dplyr::rename(name.x = name)
  }

  samdf <- data.frame(ps@sam_data) %>%
    tibble::rownames_to_column(var = 'name')

  # PCA
  pca <- stats::prcomp(ps@otu_table, center = TRUE, scale = FALSE)

  # % variance explained
  eigs <- pca$sdev^2
  varExplained <- 100 * eigs / sum(eigs)
  names(varExplained) <- paste0('PC', seq_along(varExplained))

  # Create a scree table with eigenvalues, variance explained, and cumulative variance
  scree.table <- data.frame(
    PC = paste0("PC", seq_along(varExplained)),
    Eigenvalue = eigs,
    VarianceExplained = varExplained,
    CumulativeVariance = cumsum(varExplained)
  )

  # Generate a scree plot using ggplot2
  scree.plot <- ggplot2::ggplot(scree.table, ggplot2::aes(x = as.numeric(gsub("PC", "", PC)), y = VarianceExplained)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
    ggplot2::theme_classic()

  # Extract variance explained for specified PCs
  ve.xPC <- as.character(round(varExplained[paste0('PC', xPC)], 3))
  ve.yPC <- as.character(round(varExplained[paste0('PC', yPC)], 3))

  # PCA scores
  pca.df <- data.frame(pca$x) %>%
    tibble::rownames_to_column(var = 'name')

  # Add back sample data
  pca.df <- dplyr::left_join(pca.df, samdf, by = "name")

  # Calculate plotting limits based on specified PCs
  limit <- max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))])) +
    0.05 * max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))]))

  # Initialize PCA plot
  if (!is.null(colorVar)) {
    pca.plot <- ggplot2::ggplot(pca.df, ggplot2::aes(x = .data[[paste0('PC', xPC)]],
                                                      y = .data[[paste0('PC', yPC)]],
                                                      color = .data[[colorVar]]))
  } else {
    pca.plot <- ggplot2::ggplot(pca.df, ggplot2::aes(x = .data[[paste0('PC', xPC)]],
                                                      y = .data[[paste0('PC', yPC)]]))
  }

  pca.plot <- pca.plot +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste0('PC', xPC, ' (', ve.xPC, '%)'),
         y = paste0('PC', yPC, ' (', ve.yPC, '%)')) +
    ggplot2::xlim(-limit, limit) + ggplot2::ylim(-limit, limit) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 1, color = 'black'),
          axis.ticks = ggplot2::element_line(color = 'black'),
          axis.title = ggplot2::element_text(size = 14, face = 'bold', color = 'black'))

  # Add custom color scale if provided
  if (!is.null(customColors)) {
    pca.plot <- pca.plot + ggplot2::scale_color_manual(values = customColors)
  }

  # Add custom color gradient if provided
  if (!is.null(customGradient)) {
    # Compute the midpoint only when it's actually needed -- avoids evaluating
    # pca.df[[colorVar]] (and erroring/warning on a NULL colorVar) otherwise.
    midpoint <- switch(mid,
                       "middle" = mean(range(pca.df[[colorVar]], na.rm = TRUE)),  # Mean of min and max
                       "median" = stats::median(pca.df[[colorVar]], na.rm = TRUE),       # Median
                       "mean" = mean(pca.df[[colorVar]], na.rm = TRUE),           # Mean
                       stop("'mid' must be one of 'middle', 'median', or 'mean'") # Error for invalid 'mid'
    )

    pca.plot <- pca.plot + ggplot2::scale_color_gradient2(low = paste0(customGradient$low),
                                                 mid = paste0(customGradient$mid),
                                                 high = paste0(customGradient$high),
                                                 midpoint = midpoint)
  }

  # Add optional ellipses
  if (ellipse) {
    if (!is.null(colorVar)) {
      pca.plot <- pca.plot + ggplot2::stat_ellipse(level = 0.95, ggplot2::aes(group = .data[[colorVar]]), linetype = "dashed")
    } else {
      pca.plot <- pca.plot + ggplot2::stat_ellipse(level = 0.95, linetype = "dashed")
    }
  }

  biplot_out <- .pca_biplot_layer(pca.plot, pca$rotation, pca$sdev, ps,
                                  xPC, yPC, nTaxa, bplab, colorName)

  return(list(pca.df = pca.df,
              pca.biplot = biplot_out$pca.biplot,
              loadings = biplot_out$loadings,
              pca.output = pca,
              scree.table = scree.table,
              scree.plot = scree.plot))
}
