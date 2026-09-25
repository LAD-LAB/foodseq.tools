#' @title Extract the finest taxonomic level of ASVs from a phyloseq taxonomy
#'
#' @description Identifies the lowest taxonomic name available for each ASV in a
#'   phyloseq taxonomy table and aggregates them in a new column.
#'
#' @param taxtab The taxonomy table of the phyloseq object to be renamed --
#'   accepts \code{phyloseq::tax_table(ps)} directly (a \code{taxonomyTable}),
#'   a plain matrix, or a data frame; always coerced to a data frame
#'   internally before use.
#'
#' @details Takes the tax table's right-most non-\code{NA} column per row --
#'   this assumes columns are already ordered broad-to-specific (e.g.
#'   kingdom ... species, left to right), the standard \code{phyloseq}
#'   convention.
#'
#' @return A data frame version of \code{taxtab}, with an added column,
#'   'lowest_level', that contains the name of the lowest phylogenetic level
#'   to which that ASV is identified. Wrap in
#'   \code{phyloseq::tax_table(as.matrix(...))} to use as a phyloseq
#'   tax_table again.
#' @export
#'
#'

lowest_level <- function(taxtab){
     # Update taxa names from ASV sequence to identified taxon at the most
     # precise phylogenetic level possible

     # taxtab$<- assignment below requires a data.frame -- a phyloseq
     # tax_table()/taxonomyTable (or a plain matrix) is matrix-classed, and
     # `$<-` on a matrix silently coerces the WHOLE thing into a flat list
     # rather than adding a column, destroying its 2D structure. Coercing
     # up front means this works the same whether taxtab is tax_table(ps)
     # itself, a bare matrix, or an already-built data frame.
     if (!is.data.frame(taxtab)) {
          taxtab <- as.data.frame(taxtab, stringsAsFactors = FALSE)
     }

     # This gets the right-most, non-NA value
     lowest.index <- max.col(!is.na(taxtab), 'last')
     taxtab$lowest_level <- taxtab[cbind(seq_along(lowest.index),
                                         lowest.index)]

     taxtab
}
