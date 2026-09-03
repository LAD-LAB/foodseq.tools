#' @title Filter a trnL Phyloseq Object
#'
#' @description Standardized preprocessing for trnL (plant) phyloseq
#'   objects: computes a \code{lowest_level} taxonomy column and removes
#'   unassigned and (optionally) named control ASVs. Unlike
#'   \code{\link{filter_12S_taxa}}, does \emph{not} agglomerate by
#'   \code{lowest_level} -- ASVs are left at their original resolution.
#'   Safe to re-run on an already-filtered object -- the tax table is
#'   restricted to \code{tax_cols} first, so a stale \code{lowest_level}
#'   from a prior run is dropped and recomputed rather than silently
#'   reused or double-counted.
#'
#'   Any tax_table column \emph{not} listed in \code{tax_cols} (e.g.
#'   \code{taxa}/\code{common_name} from \code{\link{assign_common_names}})
#'   is set aside, kept out of the way while \code{lowest_level} is
#'   computed and filtering is applied (so it can't interfere), and
#'   reattached at the end for whichever ASVs survive filtering.
#'
#'   Counterpart to \code{\link{filter_12S_taxa}} for the trnL (plant) marker --
#'   no Homo sapiens concept applies here, and "unassigned" is determined
#'   solely by \code{NA} \code{superkingdom} (trnL's rank set has no
#'   \code{kingdom}/\code{family}/\code{order}-based signal analogous to
#'   12S's stricter definition).
#'
#' @param ps_trnL A trnL phyloseq object to filter.
#' @param controls Character vector of exact \code{species} values to treat
#'   as control sequences and remove (e.g. a synthetic spike-in). Matched
#'   against the \code{species} column, not \code{lowest_level}. If
#'   \code{NULL}, no control-species removal is performed. Default
#'   \code{"synthetic trnL ASV"}.
#' @param tax_cols Character vector of tax_table columns to use for
#'   filtering, \emph{in the order those columns should end up in the
#'   output}. Anything else (\code{taxa}/\code{common_name} from
#'   \code{\link{assign_common_names}}, etc.) is set aside and reattached
#'   at the end rather than dropped -- see Description. A stale
#'   \code{lowest_level} from a prior run is the one exception: always
#'   excluded and recomputed fresh, never set aside. Default (broadest to
#'   most specific, matching \code{\link{assignment_trnL}}'s own rank
#'   order) \code{c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "varietas", "forma")}.
#' @param tax_coal Character vector of tax_table columns to
#'   \code{dplyr::coalesce()}, in order, into \code{lowest_level} -- i.e.
#'   the most specific non-\code{NA} rank per ASV. Deliberately excludes
#'   \code{subspecies}/\code{varietas}/\code{forma} by default
#'   (\code{lowest_level} caps out at species); order matters here (most
#'   specific first) since it drives \code{coalesce()}. Default
#'   \code{c("species", "genus", "family", "order", "class", "phylum", "superkingdom")}.
#' @param export_NA_ASVs Optional file path to write a table of ASVs with
#'   \code{NA} \code{superkingdom} to, as a CSV, with columns \code{asv},
#'   \code{prevalence} (number of samples the ASV has nonzero reads in),
#'   \code{reads_total} (summed across all samples), and \code{reads_max}
#'   (its largest count in any single sample). If \code{NULL} (default),
#'   nothing is written.
#'
#' @return The filtered phyloseq object: unassigned (\code{NA}
#'   \code{superkingdom}) and (if \code{controls} is supplied) named
#'   control ASVs removed (not agglomerated -- ASVs are left at their
#'   original resolution). The tax table's column order runs broadest to
#'   most specific (superkingdom ... forma, \code{lowest_level} last),
#'   followed by any set-aside columns (see Description) reattached for
#'   the surviving ASVs.
#'
#' @export
filter_trnL_taxa <- function(ps_trnL = NULL,
                        controls = "synthetic trnL ASV",
                        tax_cols = c("superkingdom", "phylum", "class", "order", "family",
                                    "genus", "species", "subspecies", "varietas", "forma"),
                        tax_coal = c("species", "genus", "family", "order",
                                    "class", "phylum", "superkingdom"),
                        export_NA_ASVs = NULL) {

  stopifnot(inherits(ps_trnL, "phyloseq"))

  ps_out <- ps_trnL
  taxtab <- as.data.frame(phyloseq::tax_table(ps_out)@.Data, stringsAsFactors = FALSE)

  # -- Set aside (don't drop) any columns not in tax_cols, keyed by ASV -----
  # e.g. taxa/common_name from assign_common_names(). Kept out of the way
  # while lowest_level/filtering below run (so they can't interfere), then
  # reattached at the very end for whichever ASVs survive. `lowest_level`
  # is deliberately excluded from this set-aside even if a prior
  # filter_trnL() run left one behind -- it must always be recomputed
  # fresh, never reused (this is the only column requiring that special
  # handling: filter_trnL() has no tax_glom() step, so unlike filter_12S()
  # there's no glom-induced duplicate-lowest_level scenario to worry about
  # here -- every ASV keeps its own row regardless).
  aside_cols <- setdiff(colnames(taxtab), c(tax_cols, "lowest_level"))
  aside_df   <- taxtab[, aside_cols, drop = FALSE]

  # -- Restrict to tax_cols (idempotency; output column order follows
  # tax_cols' own order, broadest-to-most-specific by default) -------------
  missing_tax_cols <- setdiff(tax_cols, colnames(taxtab))
  if (length(missing_tax_cols) > 0) {
    message("filter_trnL: tax_cols not found in tax_table, skipping: ",
            paste(missing_tax_cols, collapse = ", "))
  }
  taxtab <- taxtab[, intersect(tax_cols, colnames(taxtab)), drop = FALSE]

  # -- lowest_level: most specific non-NA rank, per tax_coal's order -------
  coal_cols <- intersect(tax_coal, colnames(taxtab))
  if (length(coal_cols) == 0) {
    stop("filter_trnL: none of `tax_coal` are present among `tax_cols`/the tax table.")
  }
  taxtab$lowest_level <- do.call(dplyr::coalesce, unname(as.list(taxtab[coal_cols])))

  phyloseq::tax_table(ps_out) <- as.matrix(taxtab)

  # -- Optional: identify (and export) unassigned ASVs ----------------------
  # NA status for trnL is always just NA superkingdom -- no family/order
  # check like filter_12S(), which trnL's rank set doesn't support the same
  # way.
  has_sk <- "superkingdom" %in% colnames(taxtab)
  if (has_sk) {
    na_asvs <- rownames(taxtab)[is.na(taxtab$superkingdom)]
  } else {
    warning("filter_trnL: 'superkingdom' not present in tax_cols; ",
            "skipping unassigned-ASV detection and filtering.")
    na_asvs <- character(0)
  }
  if (!is.null(export_NA_ASVs)) {
    if (length(na_asvs) > 0) {
      otu_mat <- methods::as(phyloseq::otu_table(ps_out), "matrix")
      if (!phyloseq::taxa_are_rows(ps_out)) otu_mat <- t(otu_mat)
      na_mat <- otu_mat[na_asvs, , drop = FALSE]
      na_df <- data.frame(
        asv         = na_asvs,
        prevalence  = rowSums(na_mat > 0),
        reads_total = rowSums(na_mat),
        reads_max   = apply(na_mat, 1, max),
        stringsAsFactors = FALSE
      )
    } else {
      na_df <- data.frame(asv = character(0), prevalence = integer(0),
                          reads_total = numeric(0), reads_max = numeric(0),
                          stringsAsFactors = FALSE)
    }
    utils::write.csv(na_df, export_NA_ASVs, row.names = FALSE)
    message(nrow(na_df), " unassigned ASV(s) written to '", export_NA_ASVs, "'.")
  }

  # -- Filter: drop unassigned and named control ASVs -----------------------
  # NA-safe for the controls check -- an ASV with NA `species` (e.g.
  # resolved only to genus) is never dropped just because it doesn't match
  # a species-level comparison; only an ASV that's actually a named control
  # is removed.
  keep <- rep(TRUE, nrow(taxtab))
  if (has_sk) keep <- keep & !is.na(taxtab$superkingdom)
  if (!is.null(controls) && length(controls) > 0) {
    if ("species" %in% colnames(taxtab)) {
      keep <- keep & (is.na(taxtab$species) | !(taxtab$species %in% controls))
    } else {
      warning("filter_trnL: 'species' not present in tax_cols; skipping control-ASV filtering.")
    }
  }

  ps_out <- phyloseq::prune_taxa(rownames(taxtab)[keep], ps_out)

  # -- Reattach the set-aside columns for whichever ASVs survived ----------
  # No tax_glom() here, so unlike filter_12S() every surviving ASV keeps its
  # own row -- a straight lookup by rowname is all that's needed, no
  # representative-row ambiguity to resolve.
  if (ncol(aside_df) > 0) {
    final_taxtab <- as.data.frame(phyloseq::tax_table(ps_out)@.Data, stringsAsFactors = FALSE)
    aside_final  <- aside_df[rownames(final_taxtab), , drop = FALSE]
    final_taxtab <- cbind(final_taxtab, aside_final)
    phyloseq::tax_table(ps_out) <- as.matrix(final_taxtab)
  }

  ps_out
}
