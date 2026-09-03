#' @title Plot ASV Length Histograms
#'
#' @description This function plots a histogram of ASV lengths.
#'
#' @param marker The marker gene, either "trnL" or "12S"
#' @param asvtab The ASV table (samples x ASVs). Defaults to looking up
#'   \code{qiime.asvtab.<marker>} in the calling environment.
#' @param out_dir Path to the output directory for saving the plot. Defaults to
#'   looking up \code{qiime.dir.<marker>} in the calling environment.
#' @param proj The name of the project
#' @param bin_size The bin size of the histogram
#' @return A histogram of ASV lengths and a dataframe of lengths
#' @export
plot_asv_length_hist <- function(marker,
                                 asvtab = get(paste0("qiime.asvtab.", marker), envir = parent.frame()),
                                 out_dir = get(paste0("qiime.dir.", marker), envir = parent.frame()),
                                 proj = project,
                                 bin_size = 2) {

  # Calculate ASV lengths
  lengths <- data.frame(
    asv = colnames(asvtab),
    reads = colSums(asvtab)
  ) |>
    dplyr::mutate(length = nchar(asv))

  # # Guess marker based on median length
  # median_len <- median(lengths$length)

  if (marker == "trnL") {
    len_range <- c(27, 76)
  } else if (marker == "12S") {
    len_range <- c(91, 102)
  } else {}

  # ggplot_build() opens an internal device to compute layout; suppress font
  # warnings by using a null device that doesn't need the font to resolve
  gb <- suppressWarnings(ggplot2::ggplot_build(
    ggplot2::ggplot(lengths, ggplot2::aes(x = length)) +
      ggplot2::geom_histogram(binwidth = bin_size, boundary = 0)
  ))

  # geom_histogram layer is [[1]]
  max_y <- max(gb$data[[1]]$count, na.rm = TRUE)
  x_range <- gb$layout$panel_scales_x[[1]]$range$range
  y_range <- gb$layout$panel_scales_y[[1]]$range$range

  x_distance <- x_range[2] - x_range[1]

  # Make plot
  p <- ggplot2::ggplot(lengths, ggplot2::aes(x = length)) +
    ggplot2::annotate("rect", xmin = len_range[1], xmax = len_range[2], ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "red") +
    ggplot2::geom_histogram(binwidth = bin_size, boundary = 0) +
    ggplot2::geom_vline(xintercept = len_range,
               color = 'red', linetype = 'dashed') +
    ggplot2::annotate("text", x = len_range[1] - x_distance * 0.025, y = max_y * 0.5,
             label = paste0("95% of ", marker, " ASVs fall within this range"),
             angle = 90, color = "red", size = 3, hjust = 0.5) +
    ggplot2::labs(x = 'ASV length (bp)', y = 'Count',
         title = paste0(proj, ": Histogram of ", marker, " ASV lengths")) +
    ggplot2::scale_x_continuous(minor_breaks = seq(0, 250, 10),
                       breaks = seq(0, 250, 50))

  # Save output
  ggplot2::ggsave(file.path(out_dir, paste0("QC_seq-lengths-histogram_", marker, ".png")),
         plot = p, width = 8, height = 6, dpi = 300)

  # Return results
  list(plot = p,
       lengths = lengths)
}
