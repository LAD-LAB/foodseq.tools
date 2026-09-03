#' @title Filter a 12SV5 Phyloseq Object
#'
#' @description Standardized preprocessing for 12SV5 (vertebrate) phyloseq
#'   objects: computes a \code{lowest_level} taxonomy column and removes
#'   Homo sapiens, unassigned, and (optionally) named control ASVs. Does
#'   \emph{not} agglomerate -- ASVs are left at their original resolution;
#'   glomming is reserved for a separate step elsewhere in the workflow.
#'   Safe to re-run on an already-filtered object -- the tax table is
#'   restricted to \code{tax_cols} first, so a stale \code{lowest_level}
#'   from a prior run is dropped and recomputed rather than silently reused
#'   or double-counted.
#'
#'   Any tax_table column \emph{not} listed in \code{tax_cols} (e.g.
#'   \code{taxa}/\code{common_name} from \code{\link{assign_common_names}})
#'   is set aside, kept out of the way while \code{lowest_level} is computed
#'   and filtering is applied (so it can't interfere -- e.g. accidentally
#'   get pulled into \code{lowest_level}), and reattached at the end for
#'   whichever ASVs survive filtering.
#'
#' @param ps_12SV5 A 12SV5 phyloseq object to filter.
#' @param controls Character vector of exact \code{species} values to treat
#'   as control sequences and remove (e.g. a synthetic spike-in). Matched
#'   against the \code{species} column, not \code{lowest_level}. If
#'   \code{NULL}, no control-species removal is performed. Default
#'   \code{"synthetic 12S ASV"}.
#' @param tax_cols Character vector of tax_table columns to use for
#'   filtering, \emph{in the order those columns should end up in
#'   the output}. Anything else (\code{taxa}/\code{common_name} from
#'   \code{\link{assign_common_names}}, etc.) is set aside and reattached
#'   at the end rather than dropped -- see Description. A stale
#'   \code{lowest_level} from a prior run is the one exception: always
#'   excluded and recomputed fresh, never set aside. Default (broadest to
#'   most specific, matching \code{\link{assignment_12S}}'s own rank order)
#'   \code{c("kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")}.
#' @param tax_coal Character vector of tax_table columns to
#'   \code{dplyr::coalesce()}, in order, into \code{lowest_level} --
#'   i.e. the most specific non-\code{NA} rank per ASV. Deliberately
#'   excludes \code{subspecies} by default (\code{lowest_level} caps out at
#'   species); order matters here (most specific first) since it drives
#'   \code{coalesce()}. Default
#'   \code{c("species", "genus", "family", "order", "class", "phylum", "kingdom")}.
#' @param export_human_ASVs Optional file path to write the list of ASVs
#'   identified as \code{Homo sapiens} (by \code{lowest_level}) to, as a
#'   CSV. If \code{NULL} (default), nothing is written.
#' @param export_NA_ASVs Optional file path to write a table of ASVs with
#'   \code{NA} \code{kingdom}, or \code{NA} at both \code{family} and
#'   \code{order}, to, as a CSV, with columns \code{asv}, \code{prevalence}
#'   (number of samples the ASV has nonzero reads in), \code{reads_total}
#'   (summed across all samples), and \code{reads_max} (its largest count
#'   in any single sample). If \code{NULL} (default), nothing is written.
#' @param calculate_human_reads_perc If \code{TRUE} (default), compute and
#'   report (via \code{message()}) the percentage of reads in
#'   \code{ps_12SV5} attributable to \code{Homo sapiens}, before filtering.
#'   Purely informational -- not returned, not required for the rest of
#'   the function, and skipped entirely (no computation) when \code{FALSE}.
#'
#' @return The filtered phyloseq object: \code{Homo sapiens}, unassigned
#'   (\code{NA} \code{kingdom}, or \code{NA} at both \code{family} and
#'   \code{order}), and (if \code{controls} is supplied) named control ASVs
#'   removed (not agglomerated -- ASVs are left at their original
#'   resolution). The tax table's column order runs broadest to most
#'   specific (kingdom ... subspecies, \code{lowest_level} last), followed
#'   by any set-aside columns (see Description) reattached for the
#'   surviving ASVs.
#'
#' @export
filter_12S_taxa <- function(ps_12SV5 = NULL,
                       controls = "synthetic 12S ASV",
                       tax_cols = c("kingdom", "phylum", "class", "order",
                                   "family", "genus", "species", "subspecies"),
                       tax_coal = c("species", "genus", "family", "order",
                                   "class", "phylum", "kingdom"),
                       export_human_ASVs = NULL,
                       export_NA_ASVs    = NULL,
                       calculate_human_reads_perc = TRUE) {

  stopifnot(inherits(ps_12SV5, "phyloseq"))

  ps_out <- ps_12SV5
  taxtab <- as.data.frame(phyloseq::tax_table(ps_out)@.Data, stringsAsFactors = FALSE)

  # -- Set aside (don't drop) any columns not in tax_cols, keyed by ASV -----
  # e.g. taxa/common_name from assign_common_names(). Kept out of the way
  # while lowest_level/filtering below run (so they can't interfere), then
  # reattached at the very end for whichever ASVs survive. `lowest_level`
  # is deliberately excluded from this set-aside even if a prior filter_12S()
  # run left one behind -- it must always be recomputed fresh, never reused.
  aside_cols <- setdiff(colnames(taxtab), c(tax_cols, "lowest_level"))
  aside_df   <- taxtab[, aside_cols, drop = FALSE]

  # -- Restrict to tax_cols (idempotency: drops a stale lowest_level from a
  # prior run). The output column order follows tax_cols' own order
  # (intersect(tax_cols, colnames(taxtab)), tax_cols first) -- so it stays
  # broadest-to-most-specific because tax_cols' default is listed that way,
  # not because of whatever order the input tax_table happened to arrive in.
  missing_tax_cols <- setdiff(tax_cols, colnames(taxtab))
  if (length(missing_tax_cols) > 0) {
    message("filter_12S: tax_cols not found in tax_table, skipping: ",
            paste(missing_tax_cols, collapse = ", "))
  }
  taxtab <- taxtab[, intersect(tax_cols, colnames(taxtab)), drop = FALSE]

  # -- Clean stray "NA" strings in subspecies (some assignment_12S() output
  # versions store the literal string rather than a real NA) --------------
  if ("subspecies" %in% colnames(taxtab)) {
    taxtab$subspecies <- ifelse(taxtab$subspecies == "NA", NA, taxtab$subspecies)
  }

  # -- lowest_level: most specific non-NA rank, per tax_coal's order -------
  coal_cols <- intersect(tax_coal, colnames(taxtab))
  if (length(coal_cols) == 0) {
    stop("filter_12S: none of `tax_coal` are present among `tax_cols`/the tax table.")
  }
  taxtab$lowest_level <- do.call(dplyr::coalesce, unname(as.list(taxtab[coal_cols])))

  phyloseq::tax_table(ps_out) <- as.matrix(taxtab)

  # -- Optional: identify (and export) Homo sapiens ASVs -------------------
  human_asvs <- rownames(taxtab)[!is.na(taxtab$lowest_level) & taxtab$lowest_level == "Homo sapiens"]
  if (!is.null(export_human_ASVs)) {
    utils::write.csv(data.frame(asv = human_asvs), export_human_ASVs, row.names = FALSE)
    message(length(human_asvs), " Homo sapiens ASV(s) written to '", export_human_ASVs, "'.")
  }

  # -- Optional: identify (and export) unassigned ASVs ----------------------
  has_kfo <- all(c("kingdom", "family", "order") %in% colnames(taxtab))
  if (has_kfo) {
    na_asvs <- rownames(taxtab)[is.na(taxtab$kingdom) | (is.na(taxtab$family) & is.na(taxtab$order))]
  } else {
    warning("filter_12S: 'kingdom'/'family'/'order' not all present in tax_cols; ",
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

  # -- Optional: report (not return) the % of reads from Homo sapiens ------
  if (calculate_human_reads_perc) {
    t_sums      <- phyloseq::taxa_sums(ps_out)
    total_reads <- sum(t_sums)
    human_reads <- sum(t_sums[intersect(human_asvs, names(t_sums))])
    pct         <- if (total_reads > 0) 100 * human_reads / total_reads else NA_real_
    message(sprintf("%.2f%% of reads in this object were from Homo sapiens.", pct))
  }

  # -- Filter: drop unassigned, Homo sapiens, and named control ASVs -------
  # NA-safe throughout -- an ASV with NA `species` (e.g. resolved only to
  # genus) is never dropped just because it doesn't match a species-level
  # comparison; only ASVs that are actually unassigned, actually human, or
  # actually a named control are removed.
  keep <- rep(TRUE, nrow(taxtab))
  if (has_kfo) {
    keep <- keep & !(is.na(taxtab$kingdom) | (is.na(taxtab$family) & is.na(taxtab$order)))
  }
  keep <- keep & (is.na(taxtab$lowest_level) | taxtab$lowest_level != "Homo sapiens")
  if (!is.null(controls) && length(controls) > 0) {
    if ("species" %in% colnames(taxtab)) {
      keep <- keep & (is.na(taxtab$species) | !(taxtab$species %in% controls))
    } else {
      warning("filter_12S: 'species' not present in tax_cols; skipping control-ASV filtering.")
    }
  }

  ps_out <- phyloseq::prune_taxa(rownames(taxtab)[keep], ps_out)

  # -- Reattach the set-aside columns for whichever ASVs survived ----------
  # No glomming here, so every surviving ASV keeps its own row -- a
  # straight lookup by rowname is all that's needed, no representative-row
  # ambiguity to resolve.
  if (ncol(aside_df) > 0) {
    final_taxtab  <- as.data.frame(phyloseq::tax_table(ps_out)@.Data, stringsAsFactors = FALSE)
    aside_final   <- aside_df[rownames(final_taxtab), , drop = FALSE]
    final_taxtab  <- cbind(final_taxtab, aside_final)
    phyloseq::tax_table(ps_out) <- as.matrix(final_taxtab)
  }

  ps_out
}
