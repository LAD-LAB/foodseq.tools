# R/apply_projection.R
#
# Applies a plan_projection() plan (with $decisions resolved) to build the
# projected phyloseq object -- the apply_harmonization() analog for
# projection. Never re-runs compare_asvs() or any other detection; trusts
# plan_projection()'s own stored state (ps_ref_asvs, taxtab, auto_map,
# needs_review, etc.) entirely, and re-derives keep_asvs/ambiguity from
# whatever's currently in plan$decisions (which the caller may have
# hand-edited since planning) via the SAME .resolve_projection_pairs()
# helper plan_projection() itself uses -- so a stale or inconsistent
# $decisions is caught here too, not just trusted blindly (the same
# defense-in-depth apply_harmonization() applies to its own `decisions`).
#
# Deliberately all-or-nothing, like the single-function design this
# replaces: if anything is still unresolved or ambiguous, $ps_projected is
# never built at all (returned as $status = "needs_review" instead) --
# unlike apply_harmonization(), where an unmerged ASV is a safe no-op
# default, an "unprojected" ASV has no safe default (it would either be
# silently missing from the projected space, or silently double-counted if
# guessed wrong), so a partial projection is never produced.
#
# Only takes `plan` and `ps_query` -- NOT ps_ref_outcomes again. This means
# apply_projection() can only ever detect ps_query-side drift since
# planning (checked below), not reference-side drift -- that protection now
# lives entirely in plan_projection()'s own prior_decisions checkpoint (see
# there), which does have a fresh ps_ref_outcomes to compare against. This
# matches apply_harmonization()'s own precedent (it has no reference-side
# staleness check either) at the cost of that specific protection if a
# `plan` is persisted and applied much later, after the reference's own
# harmonization has since changed -- re-run plan_projection() fresh in that
# case rather than trusting an old plan.
#
# One deliberate exception to "never touch samples": a ps_query sample left
# with zero reads across every ps_ref ASV (none of its ASVs correspond to
# anything in the reference space) is dropped from $ps_projected outright,
# not just flagged -- see the "Zero-read samples" section below for why.

# Checks ps_query's CURRENT taxa/samples against the fingerprint
# plan_projection() embedded on `plan` at planning time. ps_query-side only
# (see file header for why) -- never blocks, only ever warns, since
# decisions apply per-ASV and are otherwise unaffected by which ps_query
# samples/taxa are present now versus at planning time.
.apply_projection_checkpoint <- function(plan, ps_query) {
  hash_asvs <- function(asvs) vapply(asvs, function(x)
    digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))

  if (!is.null(plan$ps_query_taxa_hashes)) {
    current <- unname(hash_asvs(phyloseq::taxa_names(ps_query)))
    if (!setequal(plan$ps_query_taxa_hashes, current)) {
      warning("apply_projection: ps_query's taxa differ from when `plan` was built -- proceeding, ",
              "since decisions apply per-ASV: a newly-introduced ASV simply won't have a decision ",
              "yet (it'll surface if you re-run plan_projection()), and a removed one's decision is ",
              "simply unused.", call. = FALSE)
    }
  }
  if (!is.null(plan$ps_query_samples)) {
    current <- phyloseq::sample_names(ps_query)
    if (!setequal(plan$ps_query_samples, current)) {
      warning("apply_projection: ps_query's samples differ from when `plan` was built -- proceeding, ",
              "since decisions apply per-ASV and are otherwise unaffected by which samples are ",
              "present.", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' @title Apply a Query-to-Reference ASV Projection
#'
#' @description Applies a \code{\link{plan_projection}} plan (with
#'   \code{$decisions} resolved) to build the projected phyloseq object --
#'   the \code{\link{apply_harmonization}} analog for projection. Never
#'   re-runs detection; trusts \code{plan}'s own stored state entirely, and
#'   re-derives the merge map and ambiguity check from whatever's currently
#'   in \code{plan$decisions} (in case it was hand-edited since planning).
#'
#'   All-or-nothing: if anything is still unresolved or ambiguous,
#'   \code{$ps_projected} is never built at all -- \code{$status} comes
#'   back \code{"needs_review"} with an updated \code{$decisions} instead,
#'   the same shape \code{\link{plan_projection}}'s own \code{$decisions}
#'   has, so hand-editing it directly (or a \code{save_path}/spreadsheet
#'   round trip) and calling \code{apply_projection()} again keeps working.
#'   \code{apply_projection()} itself never re-opens the interactive review
#'   gadget (mirroring \code{\link{apply_harmonization}}, which doesn't
#'   either) -- for another interactive pass, re-run
#'   \code{\link{plan_projection}} with \code{prior_decisions} set to this
#'   result instead, which only re-shows whatever's still blank.
#'
#' @param plan The list returned by \code{\link{plan_projection}} (or by a
#'   previous \code{apply_projection()} call whose status was
#'   \code{"needs_review"}) -- must carry the internal state those
#'   functions attach, not just a bare decisions table (there's no
#'   detection state to rebuild one from). Safe to have been
#'   \code{saveRDS()}/reloaded in between.
#' @param ps_query The same phyloseq object \code{plan_projection()} was
#'   given -- re-supplied here since its real read counts are what's
#'   actually being projected.
#' @param verbose If \code{TRUE} (default), report progress and a summary.
#'
#' @section Query drift checkpoint:
#' \code{ps_query}'s current taxa and samples are checked against a
#' fingerprint \code{plan_projection()} embedded on \code{plan} -- always a
#' \code{warning()}, never blocks, since decisions apply per-ASV. There is
#' no reference-side check here at all (see this file's own header comment
#' for why) -- that protection lives in \code{\link{plan_projection}}'s own
#' \code{prior_decisions} checkpoint instead.
#'
#' @section Zero-read samples are dropped, not just flagged:
#' Unlike every other sample-level concern in this package (always left as
#' the caller's own explicit step -- see the file header), a \code{ps_query}
#' sample whose ASVs correspond to nothing in \code{ps_ref} ends up with
#' zero reads across every reference ASV, and is removed from
#' \code{$ps_projected} outright, with a \code{warning()} naming which. This
#' is the one exception to that stance: such a sample has no valid
#' composition at all, so there's no "keep it anyway" reading of it the way
#' there is for, say, a low-depth sample -- left in, it would silently
#' produce \code{NaN}/\code{Inf} (or a meaningless all-equal row) under a
#' CLR transform and distort a downstream PCA for every other sample too.
#'

#' @return A named list (class \code{"projection_result"}). If pairs still
#'   need review, \code{$status} is \code{"needs_review"} and
#'   \code{$decisions} holds them, in the same shape
#'   \code{\link{plan_projection}} itself returns (also carrying forward
#'   everything \code{apply_projection()} needs, so the result of a
#'   \code{"needs_review"} call can be fed straight back into
#'   \code{apply_projection()} again after further edits). Once every pair
#'   is resolved, \code{$status} is \code{"complete"} and the list also
#'   holds \code{$ps_projected} (\code{ps_query}'s data expressed in the
#'   reference's exact ASV/tax_table space, zero-filled for reference ASVs
#'   absent from \code{ps_query} -- \strong{except} any \code{ps_query}
#'   sample left with zero reads across every reference ASV, which is
#'   dropped entirely rather than kept as an all-zero row, with a
#'   \code{warning()} naming which; see the section below) and
#'   \code{$match_report} (one row per
#'   original \code{ps_query} ASV: its \code{fate} -- \code{exact_match},
#'   \code{ghost_redirected}, \code{auto_merged}, \code{reviewed_merged},
#'   \code{reviewed_distinct}, or \code{no_correspondence_found} -- and,
#'   where applicable, the reference ASV it was mapped to).
#'   \code{$compare_result} is always included.
#'
#' @examples
#' \dontrun{
#' # Interactive review (if anything needs it) happens inside plan_projection().
#' plan <- plan_projection(ref_harmonized, ps_query)
#' result <- apply_projection(plan, ps_query)
#' result$ps_projected
#' }
#'
#' @export
apply_projection <- function(plan, ps_query, verbose = TRUE) {

  stopifnot(inherits(plan, c("projection_plan", "projection_result")))
  stopifnot(inherits(ps_query, "phyloseq"))
  required <- c("ps_ref_asvs", "taxtab", "ghost_target", "exact", "auto_map",
                "needs_review", "is_s1", "cross", "ps_query_asvs_all")
  missing_fields <- setdiff(required, names(plan))
  if (length(missing_fields) > 0) {
    stop("apply_projection: `plan` is missing required field(s): ",
         paste(missing_fields, collapse = ", "), " -- it must be the list returned by ",
         "plan_projection() (or a prior apply_projection() 'needs_review' result), not a bare ",
         "decisions table (there's no detection state to rebuild one from).")
  }

  .apply_projection_checkpoint(plan, ps_query)

  res <- .resolve_projection_pairs(plan$auto_map, plan$needs_review, plan$decisions, plan$ps_ref_asvs)

  if (length(res$ambiguous_y) > 0 && verbose) {
    message(length(res$ambiguous_y), " query ASV(s) have CONFLICTING merge targets among ps_ref's ",
            "ASVs -- forcing these back to review rather than guessing: ",
            paste(res$ambiguous_y, collapse = ", "))
  }

  decisions_out <- .build_projection_decisions_out(
    plan$needs_review, plan$decisions, res$covered, res$ambiguous_y, plan$is_s1, plan$cross,
    plan$ps_ref_asvs, plan$ps_query_asvs_all, ps_query, plan$src_x, plan$src_y, verbose
  )

  # decisions_out is now the FULL table (not just what's still pending), so
  # "anything left to resolve" has to be checked directly via res$covered/
  # ambiguity rather than by decisions_out being non-NULL (it's non-NULL
  # whenever there was ever anything to review this round, resolved or not).
  s1_ambiguous_n <- sum(plan$is_s1 & plan$cross$asv_j %in% res$ambiguous_y)
  all_resolved   <- all(res$covered) && s1_ambiguous_n == 0

  if (!all_resolved) {
    return(structure(
      utils::modifyList(plan, list(status = "needs_review", decisions = decisions_out)),
      class = "projection_result"
    ))
  }

  # -- All decided (and no ambiguity remains): build the full query_asv -> ps_ref_asv map --
  map <- res$auto_map_applied
  if (length(res$keep_asvs_applied) > 0) map[names(res$keep_asvs_applied)] <- unname(res$keep_asvs_applied)
  for (e in plan$exact) map[e] <- e
  if (length(plan$ghost_target) > 0) map[names(plan$ghost_target)] <- unname(plan$ghost_target)

  # -- Build the projected OTU matrix (ps_query samples x ps_ref taxa) --------
  otu_q       <- methods::as(phyloseq::otu_table(ps_query), "matrix")
  rows_mode_q <- phyloseq::taxa_are_rows(ps_query)
  if (rows_mode_q) otu_q <- t(otu_q)   # normalize to samples x taxa
  samples_q <- rownames(otu_q)

  out_mat <- matrix(0, nrow = length(samples_q), ncol = length(plan$ps_ref_asvs),
                    dimnames = list(samples_q, plan$ps_ref_asvs))
  matched_q <- intersect(colnames(otu_q), names(map))
  for (q_asv in matched_q) {
    ref_asv <- map[[q_asv]]
    out_mat[, ref_asv] <- out_mat[, ref_asv] + otu_q[, q_asv]
  }

  empty_samples   <- samples_q[rowSums(out_mat) == 0]
  n_empty_samples <- length(empty_samples)
  if (n_empty_samples > 0) {
    # Automatically removed, not just flagged -- unlike apply_projection()'s
    # usual "never touch samples" stance (see file header: sample-level
    # filtering is always left as the caller's own explicit step), a sample
    # with zero reads across every ps_ref ASV isn't a filtering judgment
    # call the way e.g. dropping low-depth samples is -- none of its own
    # ASVs correspond to anything in the reference space at all, so it has
    # no valid composition, full stop. Left in, it would silently ride
    # along into $ps_projected and corrupt a downstream CLR transform
    # (NaN/Inf, or a meaningless all-equal row) and, in turn, the PCA --
    # there is no meaningful "keep it anyway" option, so it's dropped here
    # rather than deferred to the caller.
    warning(n_empty_samples, " ps_query sample(s) had zero total reads after projection ",
            "(none of their ASVs correspond to anything in ps_ref) and ",
            "have been removed from $ps_projected: ",
            paste(utils::head(empty_samples, 5), collapse = ", "),
            if (n_empty_samples > 5) ", ..." else "", ".")
    out_mat   <- out_mat[rowSums(out_mat) != 0, , drop = FALSE]
    samples_q <- rownames(out_mat)
  }

  ps_projected <- phyloseq::phyloseq(
    phyloseq::otu_table(out_mat, taxa_are_rows = FALSE),
    plan$taxtab
  )
  ps_query_samdf <- phyloseq::sample_data(ps_query, errorIfNULL = FALSE)
  if (!is.null(ps_query_samdf)) {
    ps_query_samdf <- ps_query_samdf[samples_q, , drop = FALSE]
    phyloseq::sample_data(ps_projected) <- ps_query_samdf
  }

  # -- match_report: one row per original ps_query ASV -----------------------
  ps_query_asvs_all <- plan$ps_query_asvs_all
  fate           <- stats::setNames(rep("no_correspondence_found", length(ps_query_asvs_all)), ps_query_asvs_all)
  matched_ps_ref <- stats::setNames(rep(NA_character_, length(ps_query_asvs_all)), ps_query_asvs_all)

  fate[intersect(ps_query_asvs_all, plan$exact)]        <- "exact_match"
  matched_ps_ref[intersect(ps_query_asvs_all, plan$exact)] <- intersect(ps_query_asvs_all, plan$exact)
  if (length(plan$ghost_target) > 0) {
    fate[names(plan$ghost_target)]        <- "ghost_redirected"
    matched_ps_ref[names(plan$ghost_target)] <- unname(plan$ghost_target)
  }
  if (length(res$auto_map_applied) > 0) {
    fate[names(res$auto_map_applied)]        <- "auto_merged"
    matched_ps_ref[names(res$auto_map_applied)] <- unname(res$auto_map_applied)
  }
  if (length(res$keep_asvs_applied) > 0) {
    fate[names(res$keep_asvs_applied)]        <- "reviewed_merged"
    matched_ps_ref[names(res$keep_asvs_applied)] <- unname(res$keep_asvs_applied)
  }
  declined <- setdiff(plan$needs_review$asv_j, names(res$keep_asvs_applied))
  fate[intersect(names(fate), declined)] <- "reviewed_distinct"

  reads_q_per_asv <- stats::setNames(rep(0, length(ps_query_asvs_all)), ps_query_asvs_all)
  present         <- intersect(ps_query_asvs_all, colnames(otu_q))
  reads_q_per_asv[present] <- colSums(otu_q[, present, drop = FALSE])

  match_report <- data.frame(
    asv                = ps_query_asvs_all,
    fate               = unname(fate[ps_query_asvs_all]),
    matched_ps_ref_asv = unname(matched_ps_ref[ps_query_asvs_all]),
    reads_total        = unname(reads_q_per_asv[ps_query_asvs_all]),
    stringsAsFactors   = FALSE
  )

  if (verbose) {
    n_mapped <- sum(match_report$fate %in%
                     c("exact_match", "ghost_redirected", "auto_merged", "reviewed_merged"))
    message(sprintf(
      "apply_projection: %d/%d ps_query ASV(s) mapped into ps_ref's space (%d exact, %d ghost-redirected, %d auto-merged, %d reviewed); %d dropped (no correspondence or reviewed as distinct).",
      n_mapped, length(ps_query_asvs_all),
      sum(match_report$fate == "exact_match"),
      sum(match_report$fate == "ghost_redirected"),
      sum(match_report$fate == "auto_merged"),
      sum(match_report$fate == "reviewed_merged"),
      length(ps_query_asvs_all) - n_mapped
    ))
  }

  structure(
    utils::modifyList(plan, list(status = "complete", ps_projected = ps_projected,
                          match_report = match_report)),
    class = "projection_result"
  )
}
