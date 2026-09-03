#' @title Percent of ASVs Assigned per Taxonomic Rank
#'
#' @description Computes the percentage of ASVs in a phyloseq object that
#'   have a non-missing, non-empty assignment at each taxonomic rank. Useful
#'   for assessing taxonomy completeness after assignment.
#'
#' @param ps Phyloseq object with a tax_table.
#'
#' @return A data frame with columns \code{Taxonomic_Level} and
#'   \code{Percent_Assigned} (numeric, 0–100).
#'
#' @export
percent_assigned_tax <- function(ps) {

  tax <- as.data.frame(phyloseq::tax_table(ps))

  pct <- sapply(tax, function(col) mean(!is.na(col) & col != "") * 100)

  data.frame(
    Taxonomic_Level  = names(pct),
    Percent_Assigned = unname(pct),
    row.names        = NULL
  )
}
