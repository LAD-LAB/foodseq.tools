#' @title Update Taxonomy
#'
#' @description This function updates the taxonomic assignments of a
#'   phyloseq object.
#'
#' @param ps A phyloseq object
#' @param reference A file path to a FASTA-formatted reference
#' @param marker The marker name as a string
#' @return Updated phyloseq object
#' @export
update_taxonomy <- function(ps, reference, marker) {
  marker <- tolower(marker)
  qiime.asvtab <- phyloseq::otu_table(ps)@.Data

  if (marker %in% c("plants", "plant", "trnl", "trnlgh")) {
    assignments <- assignment_trnL(qiime.asvtab, reference, NULL)$assignments
  } else if (marker %in% c("animals", "animal", "vertebrates", "vertebrate", "12s", "12sv5")) {
    assignments <- assignment_12S(qiime.asvtab, reference, NULL)$assignments
  } else {
    stop("The input received for `marker` is invalid.")
  }

  phyloseq::tax_table(ps) <- phyloseq::tax_table(assignments)

  return(ps)
}
