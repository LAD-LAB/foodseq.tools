#' @title Filter Phyloseq by Sample Metadata Field
#'
#' @description Subsets a phyloseq object to samples where a specified
#'   sample_data field equals a given value. A convenience wrapper around
#'   \code{phyloseq::prune_samples}.
#'
#' @param ps Phyloseq object.
#' @param field Character; name of the column in sample_data to filter on.
#' @param value The value to match; samples where \code{field == value} are
#'   retained.
#'
#' @return Phyloseq object containing only the matching samples.
#'
#' @export
filter_phyloseq <- function(ps, field, value) {

  if (!(field %in% colnames(phyloseq::sample_data(ps)))) {
    stop("Field '", field, "' not found in sample_data.")
  }

  keep_samples <- phyloseq::sample_data(ps)[[field]] == value
  phyloseq::prune_samples(keep_samples, ps)
}
