# R/apply_harmonization.R
#
# Applies the merge/rename decisions detected by compare_asvs() (see
# plan_harmonization.R for the full scenario glossary S1-S5/conflict/
# unassigned and the plan_harmonization()/apply_harmonization() usage
# pattern) and optionally reviewed via resolve_conflicts_interactive() or by
# hand-editing plan_harmonization()'s $decisions table.

# Validates and derives keep_asvs/manual_updates/s1_skip_keys/rationale from
# a decisions data.frame (plan_harmonization()'s $decisions table, or an
# equivalent hand-edited/re-imported one). Keyed directly by the table's own
# asv_i/asv_j columns -- no separate ID-lookup step needed.
.decisions_to_vectors <- function(decisions_df) {
  keywords    <- c("merge_keep_i", "merge_keep_j", "keep_distinct")
  known_asvs  <- unique(c(decisions_df$asv_i, decisions_df$asv_j))
  # A decision is valid if it's one of the 3 keywords, OR a literal ASV
  # sequence that itself appears somewhere in `decisions` -- i.e. a real
  # merge target, used when both asv_i and asv_j were merged into a third,
  # external representative (e.g. a multi-ASV group resolved via
  # resolve_conflicts_interactive()'s hub topology, where every
  # non-representative member points straight at one shared ASV rather
  # than at each other). Anything else is a typo/invalid value.
  bad <- unique(decisions_df$decision[!is.na(decisions_df$decision) &
                                        !decisions_df$decision %in% keywords &
                                        !decisions_df$decision %in% known_asvs])
  if (length(bad) > 0) {
    stop("apply_harmonization: invalid decision value(s) in `decisions`: ",
         paste(bad, collapse = ", "), ". Must be one of: ",
         paste(keywords, collapse = ", "), ", a literal ASV sequence ",
         "appearing elsewhere in `decisions` (the external merge target), or blank.")
  }

  # Guarantee these columns exist (as real NA, not just absent) before the
  # loop below reads them -- a decisions_df missing either one (e.g. a
  # minimal hand-built table) would otherwise make row$chosen_name/
  # row$rationale NULL, and is.na(NULL) is logical(0), which errors under
  # R's stricter `&&` (R >= 4.3) rather than just silently being falsy.
  if (!"chosen_name" %in% colnames(decisions_df)) decisions_df$chosen_name <- NA_character_
  if (!"rationale"   %in% colnames(decisions_df)) decisions_df$rationale   <- NA_character_

  s1_skip_keys        <- character(0)
  rationale_out        <- character(0)
  manual_updates_list <- list()  # asv_seq -> all chosen_name values seen, for a consistency check
  keep_asvs_list      <- list()  # discard_asv -> all keep_asv target(s) seen, for conflict detection

  for (i in seq_len(nrow(decisions_df))) {
    row      <- decisions_df[i, ]
    asv_i    <- row$asv_i
    asv_j    <- row$asv_j
    pair_key <- paste0(asv_i, "|", asv_j)
    is_auto  <- identical(row$action, "auto_merge_keep_shorter")
    dec      <- row$decision

    if (!is.na(row$rationale) && nzchar(row$rationale)) rationale_out[pair_key] <- row$rationale

    if (is.na(dec)) next                          # S1: default auto-merge proceeds; else: still unresolved
    if (identical(dec, "keep_distinct")) {
      if (is_auto) s1_skip_keys <- c(s1_skip_keys, pair_key)
      next
    }

    # dec is "merge_keep_i", "merge_keep_j", or a literal ASV sequence (the
    # external representative both asv_i and asv_j were merged into).
    if (is_auto) s1_skip_keys <- c(s1_skip_keys, pair_key)  # take explicit control of direction
    if (identical(dec, "merge_keep_i")) {
      keep_asvs_list[[asv_j]] <- c(keep_asvs_list[[asv_j]], asv_i)
    } else if (identical(dec, "merge_keep_j")) {
      keep_asvs_list[[asv_i]] <- c(keep_asvs_list[[asv_i]], asv_j)
    } else if (identical(dec, asv_i)) {
      # Literal form of "merge_keep_i" (e.g. hand-edited CSV).
      keep_asvs_list[[asv_j]] <- c(keep_asvs_list[[asv_j]], asv_i)
    } else if (identical(dec, asv_j)) {
      # Literal form of "merge_keep_j" (e.g. hand-edited CSV).
      keep_asvs_list[[asv_i]] <- c(keep_asvs_list[[asv_i]], asv_j)
    } else {
      # External target: both asv_i and asv_j discard into it. Downstream
      # apply_harmonization() already resolves multi-hop chains
      # (resolve_final_keep()), so it's fine if `dec` is itself an
      # intermediate hop rather than the ultimate final representative.
      keep_asvs_list[[asv_i]] <- c(keep_asvs_list[[asv_i]], dec)
      keep_asvs_list[[asv_j]] <- c(keep_asvs_list[[asv_j]], dec)
    }

    name <- if (!is.na(row$chosen_name) && nzchar(row$chosen_name)) row$chosen_name else NA_character_
    if (!is.na(name)) {
      manual_updates_list[[asv_i]] <- c(manual_updates_list[[asv_i]], name)
      manual_updates_list[[asv_j]] <- c(manual_updates_list[[asv_j]], name)
    }
  }

  inconsistent <- vapply(manual_updates_list, function(x) length(unique(x)) > 1, logical(1))
  if (any(inconsistent)) {
    warning("apply_harmonization: chosen_name is inconsistent across rows for ASV(s): ",
            paste(names(manual_updates_list)[inconsistent], collapse = ", "),
            ". Using the last value seen for each; consider making these consistent in `decisions`.")
  }
  manual_updates <- if (length(manual_updates_list) > 0) {
    stats::setNames(vapply(manual_updates_list, function(x) x[length(x)], character(1)),
                    names(manual_updates_list))
  } else character(0)

  # A discard ASV can only have ONE ultimate fate -- if different rows
  # decided it should merge into two different, distinct targets, that's a
  # genuine contradiction (not a stylistic inconsistency like chosen_name
  # above), and silently picking one (last-value-wins) would risk
  # misattributing reads with no visible trace of the conflict. Refuse to
  # produce an ambiguous keep_asvs at all; make the caller fix `decisions`.
  conflicting <- vapply(keep_asvs_list, function(x) length(unique(x)) > 1, logical(1))
  if (any(conflicting)) {
    conflict_msgs <- vapply(names(keep_asvs_list)[conflicting], function(k) {
      sprintf("%s -> {%s}", k, paste(unique(keep_asvs_list[[k]]), collapse = ", "))
    }, character(1))
    stop("apply_harmonization: conflicting merge decisions for the same ASV in `decisions` -- ",
         "it was decided to merge into more than one different, distinct target across ",
         "different rows. Make the decisions consistent (or 'keep_distinct' for the extra ",
         "row(s)) before re-applying: ", paste(conflict_msgs, collapse = "; "))
  }
  keep_asvs <- if (length(keep_asvs_list) > 0) {
    stats::setNames(vapply(keep_asvs_list, `[[`, character(1), 1), names(keep_asvs_list))
  } else character(0)

  list(keep_asvs = keep_asvs, manual_updates = manual_updates,
      s1_skip_keys = unique(s1_skip_keys), rationale = rationale_out)
}

# Structural validator for an ASV ledger (apply_harmonization()'s own
# $asv_outcomes, enriched with full taxonomy and cumulative across rounds --
# see the `prior_outcomes` argument -- and consumed by apply_projection() as its
# replacement for a live `psx` object). Shared by both: apply_asv_
# corrections() runs it defensively on its own output and on any incoming
# `prior_outcomes`; apply_projection() runs it on whatever ledger it's handed,
# since that may be a hand-edited CSV. Every check here corresponds to a way
# the ledger could otherwise silently corrupt a downstream projection:
#   - duplicate rows for the same ASV -- which one wins is undefined
#   - an invalid asv_status -- typo, not one of the 4 known values
#   - a merged_into chain that cycles, dangles (points at an ASV absent from
#     the ledger entirely), or resolves to something that isn't itself a
#     currently-final (representative/unchanged) ASV
#   - a representative/unchanged ASV with no taxonomy at all -- detection
#     against it would be silently meaningless
.validate_asv_ledger <- function(df, tax_rank_cols = NULL) {
  stopifnot(is.data.frame(df))
  errs <- character(0)

  if (!all(c("asv", "asv_status", "merged_into") %in% colnames(df))) {
    stop("Invalid ASV ledger: missing required column(s) -- must have at least ",
         "'asv', 'asv_status', and 'merged_into'.")
  }

  dup <- unique(df$asv[duplicated(df$asv)])
  if (length(dup) > 0) {
    errs <- c(errs, sprintf("duplicate row(s) for the same ASV: %s",
                            paste(utils::head(dup, 5), collapse = ", ")))
  }

  valid_status <- c("unchanged", "representative", "dropped after merge", "unresolved")
  bad_status   <- unique(df$asv_status[!df$asv_status %in% valid_status])
  if (length(bad_status) > 0) {
    errs <- c(errs, sprintf("invalid asv_status value(s): %s (must be one of: %s)",
                            paste(bad_status, collapse = ", "), paste(valid_status, collapse = ", ")))
  }

  drop_rows <- df[!is.na(df$asv_status) & df$asv_status == "dropped after merge", ]
  if (nrow(drop_rows) > 0) {
    drop_map <- stats::setNames(drop_rows$merged_into, drop_rows$asv)
    status_by_asv <- stats::setNames(df$asv_status, df$asv)
    for (a in names(drop_map)) {
      seen <- character(0); cur <- a; broken <- NA_character_
      while (cur %in% names(drop_map)) {
        if (cur %in% seen) { broken <- "cycle"; break }
        seen <- c(seen, cur); cur <- unname(drop_map[[cur]])
        if (is.na(cur) || !nzchar(cur)) { broken <- "blank merged_into"; break }
      }
      if (!is.na(broken)) {
        errs <- c(errs, sprintf("'%s' has a %s in its merge chain", a, broken)); next
      }
      if (!cur %in% names(status_by_asv)) {
        errs <- c(errs, sprintf("'%s' resolves to '%s', which doesn't appear in the ledger at all", a, cur)); next
      }
      final_status <- unname(status_by_asv[[cur]])
      if (!identical(final_status, "representative") && !identical(final_status, "unchanged")) {
        errs <- c(errs, sprintf("'%s' resolves to '%s', whose own status is '%s' (not representative/unchanged)",
                                a, cur, final_status))
      }
    }
  }

  # A representative/unchanged ASV with no taxonomy at all is NOT treated as
  # corruption -- an entirely unassigned ASV that was never flagged (nothing
  # to compare it against) is a completely normal, expected state for real
  # data (taxonomy-based cleanup is a separate, later step in the workflow,
  # so unassigned ASVs are routinely still present when this runs). It just
  # means detection will never match
  # against it -- still a perfectly valid target for an exact-sequence match
  # or a manual decision. A hard stop() here would break the ordinary
  # apply_harmonization() workflow (which validates its own output
  # unconditionally) for any dataset with unassigned ASVs at all, which is
  # the common case, not the exception. Reported as a warning instead.
  if (!is.null(tax_rank_cols)) {
    rank_cols_present <- intersect(tax_rank_cols, colnames(df))
    final_rows <- df[!is.na(df$asv_status) & df$asv_status %in% c("representative", "unchanged"), ]
    if (nrow(final_rows) > 0 && length(rank_cols_present) > 0) {
      has_any_tax <- apply(final_rows[, rank_cols_present, drop = FALSE], 1, function(r)
        any(!is.na(r) & nzchar(as.character(r))))
      blank_asvs <- final_rows$asv[!has_any_tax]
      if (length(blank_asvs) > 0) {
        warning(sprintf("ASV ledger: %d representative/unchanged ASV(s) have no taxonomy at all ",
                        length(blank_asvs)), "(fine for an exact-sequence match; will never match ",
                "anything via substring/species-set detection): ",
                paste(utils::head(blank_asvs, 5), collapse = ", "),
                if (length(blank_asvs) > 5) ", ..." else "", call. = FALSE)
      }
    }
  }

  if (length(errs) > 0) {
    stop("Invalid ASV ledger:\n  - ", paste(errs, collapse = "\n  - "))
  }
  invisible(TRUE)
}

#' @title Apply ASV Harmonization Corrections
#'
#' @description Applies the merge/rename decisions in a \code{decisions}
#'   table (typically \code{\link{plan_harmonization}}'s \code{$decisions},
#'   hand-edited as needed). Adds a \code{glom_name} column to the tax
#'   table, applies S2 rank alignment for merged pairs, physically merges
#'   ASVs where a merge was decided (transferring read counts and pruning
#'   the discarded ASV -- including through transitive merge chains, e.g. A
#'   merged into B which is itself merged into C, without losing reads
#'   regardless of processing order), and applies all chosen-name updates.
#'
#' @param ps A phyloseq object -- the same object (or the merged result of
#'   the same \code{ps_list}) that was passed to
#'   \code{\link{plan_harmonization}}.
#' @param decisions One of: the list returned by
#'   \code{\link{plan_harmonization}}; a decisions data frame (e.g.
#'   \code{plan$decisions}, possibly hand-edited); or a file path to a
#'   decisions table saved as a CSV (e.g. via
#'   \code{write.csv(plan$decisions, "decisions.csv", row.names = FALSE)}).
#'   Each row's \code{decision} column must be one of \code{"merge_keep_i"},
#'   \code{"merge_keep_j"}, \code{"keep_distinct"}, a literal ASV sequence
#'   (used when both \code{asv_i} and \code{asv_j} were merged into a
#'   third, external representative ASV -- e.g. a multi-ASV group resolved
#'   interactively, where every non-representative member points straight
#'   at one shared ASV rather than at each other), or blank (still
#'   unresolved); \code{chosen_name} and \code{rationale} are optional free
#'   text.
#' @param tax_rank_cols Character vector of formal taxonomic rank column
#'   names; inferred from \code{ps}'s tax table when \code{NULL}.
#' @param prior_outcomes Optional prior \code{apply_harmonization()} result
#'   (or its \code{$asv_outcomes} data frame, or a file path to one saved as
#'   a CSV) from an earlier round on a related object -- makes
#'   \code{$asv_outcomes} cumulative across rounds rather than starting
#'   fresh each call. Needed because \code{ps} physically loses an ASV the
#'   moment it's pruned; without this, \code{$asv_outcomes} would only ever
#'   remember the most recent round's merges, silently forgetting anything
#'   discarded in an earlier one. When supplied, this round's fresh rows are
#'   upserted onto the prior ones (this round wins for any ASV appearing in
#'   both), and every \code{merged_into} chain across the *combined* table is
#'   re-resolved so it always points at the current, ultimate representative
#'   -- e.g. if round 1 recorded A -> B and round 2 later discards B into C,
#'   A's entry is updated to A -> C rather than left pointing at the
#'   now-defunct B. See \code{\link{apply_projection}}, which consumes this
#'   cumulative ledger in place of a live reference object.
#'
#' @section Sample/taxa checkpoints:
#' If \code{decisions} carries the \code{plan_samples}/\code{plan_taxa_hashes}
#' manifest columns written by \code{\link{plan_harmonization}}, \code{ps} is
#' checked against them before anything is applied. Neither check ever
#' blocks execution -- there is no override parameter -- they only ever
#' \code{warning()}/\code{message()}: an ASV in \code{ps} the plan never saw
#' is reported (and left unchanged, since no decision exists for it; consider
#' re-running \code{plan_harmonization()} with \code{prior_decisions}); a
#' sample in \code{ps} the plan never saw is reported too, though lower risk
#' since every decision applies per-ASV regardless of which samples are
#' present (a \code{ps} with a pure *subset* of the plan's original samples
#' passes silently, just an informational message); if both checks fail at
#' once the warning escalates, recommending against trusting the output
#' without double-checking \code{ps}. Older/hand-built \code{decisions}
#' without the manifest columns, or a plan with zero flagged pairs (nothing
#' to carry the manifest on), skip both checks silently.
#'
#' @return A list with \code{ps} (the corrected phyloseq object),
#'   \code{taxtab} (its tax table as a data frame), \code{unresolved}
#'   (decision rows still blank, needing further review),
#'   \code{asv_outcomes} -- the \strong{ASV ledger}: one row per original
#'   ASV (as of \code{ps} at this call, plus anything carried forward via
#'   \code{prior_outcomes}), with its fate (\code{asv_status}:
#'   \code{"unchanged"} / \code{"representative"} / \code{"dropped after
#'   merge"} / \code{"unresolved"} -- the last meaning this ASV is part of
#'   a still-blank flagged decision, genuinely unsettled, not the same as
#'   never having been flagged at all), \code{merged_with}/\code{merged_into}
#'   (the latter always resolved to the ultimate surviving representative,
#'   through arbitrarily long chains), the \strong{full rank-by-rank
#'   taxonomy} (one column per \code{tax_rank_cols} entry -- enough for this
#'   table alone to stand in for a live reference object, e.g. in
#'   \code{\link{apply_projection}}), \code{batch_source} (derived from
#'   \code{decisions}' own \code{source_i}/\code{source_j} columns --
#'   populated only for ASVs that appear in at least one flagged pair; an
#'   ASV never flagged has no row to derive it from and stays \code{NA}),
#'   read counts, and rationale; \code{n_s2_corrected}; and \code{n_manual}.
#'   To resolve unresolved pairs manually: edit \code{decisions}'
#'   \code{decision} column and re-run.
#'
#' @export
apply_harmonization <- function(ps, decisions, tax_rank_cols = NULL,
                                  prior_outcomes = NULL) {

  stopifnot(requireNamespace("phyloseq", quietly = TRUE))
  stopifnot(inherits(ps, "phyloseq"))

  # -- Normalize prior_outcomes (same 3-input-type pattern used everywhere
  # else in this package) and validate it immediately, before merging, so
  # corruption in an incoming ledger is caught with a clear message rather
  # than silently propagated. --------------------------------------------
  prior_ledger <- NULL
  if (!is.null(prior_outcomes)) {
    prior_ledger <- if (is.list(prior_outcomes) && !is.data.frame(prior_outcomes) &&
                        "asv_outcomes" %in% names(prior_outcomes)) {
      prior_outcomes$asv_outcomes
    } else if (is.data.frame(prior_outcomes)) {
      prior_outcomes
    } else if (is.character(prior_outcomes) && length(prior_outcomes) == 1) {
      utils::read.csv(prior_outcomes, stringsAsFactors = FALSE, na.strings = "")
    } else {
      stop("apply_harmonization: `prior_outcomes` must be a previous apply_harmonization() ",
           "result, its $asv_outcomes data.frame, or a file path to one saved as a CSV.")
    }
    stopifnot(is.data.frame(prior_ledger))
    .validate_asv_ledger(prior_ledger)
  }

  flagged <- if (inherits(decisions, "harmonization_plan")) {
    decisions$decisions
  } else if (is.data.frame(decisions)) {
    decisions
  } else if (is.character(decisions) && length(decisions) == 1) {
    utils::read.csv(decisions, stringsAsFactors = FALSE, na.strings = "")
  } else {
    stop("`decisions` must be the list returned by plan_harmonization(), a decisions ",
         "data.frame, or a file path to a decisions CSV.")
  }
  stopifnot(is.data.frame(flagged))
  flagged <- as.data.frame(
    lapply(flagged, function(x) {
      if (is.character(x)) return(as.character(x))
      if (is.integer(x))   return(as.integer(x))
      if (is.numeric(x))   return(as.numeric(x))
      if (is.logical(x))   return(as.logical(x))
      if (is.factor(x))    return(as.factor(x))
      unclass(x)
    }),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  rownames(flagged) <- NULL

  # -- Taxa/sample checkpoints -----------------------------------------------
  # Both are purely informational -- neither ever blocks execution, and there
  # is no override parameter. Every merge/prune below is keyed entirely off
  # the tax_table dimension and is otherwise unaffected by which samples are
  # present, so a sample mismatch alone is low-risk; a taxa mismatch means
  # some ASVs simply never got a decision made for them. plan_samples/
  # plan_taxa_hashes are written by plan_harmonization() (see there) --
  # repeated on every decisions row so the checkpoint survives hand-editing/
  # row-filtering. A `decisions` predating this feature (or hand-built with
  # no such columns), or one whose plan happened to have zero flagged pairs
  # (no row to carry the manifest on), has nothing to check against -- both
  # checkpoints are silently skipped in that case rather than erroring.
  foreign_taxa <- character(0)
  taxa_checked <- FALSE
  if ("plan_taxa_hashes" %in% colnames(flagged) && nrow(flagged) > 0) {
    taxa_checked  <- TRUE
    plan_hash_str <- flagged$plan_taxa_hashes[[1]]
    plan_hashes   <- if (!is.na(plan_hash_str) && nzchar(plan_hash_str))
      strsplit(plan_hash_str, ";", fixed = TRUE)[[1]] else character(0)
    ps_asvs   <- phyloseq::taxa_names(ps)
    ps_hashes <- vapply(ps_asvs, function(x)
      digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))
    foreign_taxa <- ps_asvs[!ps_hashes %in% plan_hashes]
    if (length(foreign_taxa) > 0) {
      warning("apply_harmonization: ", length(foreign_taxa), " ASV(s) in `ps` were not ",
              "present in the plan used to generate these decisions and were left ",
              "unchanged. Consider re-running plan_harmonization() with ",
              "prior_decisions = <this plan> to properly harmonize them before merging.",
              call. = FALSE)
    }
  }
  taxa_ok <- !taxa_checked || length(foreign_taxa) == 0

  if ("plan_samples" %in% colnames(flagged) && nrow(flagged) > 0) {
    plan_sample_str <- flagged$plan_samples[[1]]
    plan_samples    <- if (!is.na(plan_sample_str) && nzchar(plan_sample_str))
      strsplit(plan_sample_str, ";", fixed = TRUE)[[1]] else character(0)
    ps_samples      <- phyloseq::sample_names(ps)
    foreign_samples <- setdiff(ps_samples, plan_samples)
    missing_samples <- setdiff(plan_samples, ps_samples)

    if (length(foreign_samples) == 0) {
      if (length(missing_samples) > 0) {
        message("apply_harmonization: ", length(missing_samples), " sample(s) from the ",
                "plan are not present in `ps` -- proceeding, since decisions apply per-ASV ",
                "and this looks like a subset of the plan's original samples.")
      }
    } else if (taxa_ok) {
      warning("apply_harmonization: ", length(foreign_samples), " sample(s) in `ps` were ",
              "not present when these decisions were made. Decisions apply per-ASV and are ",
              "unaffected by which samples are present, so this is likely fine if you're ",
              "intentionally applying the same decisions to a related-but-different sample ",
              "set -- but double check `ps` is the object you meant to pass.", call. = FALSE)
    } else {
      warning("apply_harmonization: BOTH the taxa and sample checkpoints failed for this ",
              "`ps` -- this usually means `ps` is not derived from the same data ",
              "plan_harmonization() ran on. Recommend verifying `ps` before trusting this ",
              "output; ", length(foreign_taxa), " ASV(s) were ignored (see above) and results ",
              "may be incomplete or wrong.", call. = FALSE)
    }
  }

  # Per-ASV source/batch label, derived directly from `decisions`' own
  # source_i/source_j columns -- no separate source-map object needed. Only
  # covers ASVs that appear in at least one flagged pair; an ASV never
  # flagged has no row to derive it from and stays unlabeled.
  source_map <- character(0)
  if (all(c("source_i", "source_j") %in% colnames(flagged))) {
    source_map <- stats::setNames(c(flagged$source_i, flagged$source_j),
                                  c(flagged$asv_i, flagged$asv_j))
    source_map <- source_map[!duplicated(names(source_map))]
  }

  vecs                <- .decisions_to_vectors(flagged)
  keep_asvs           <- vecs$keep_asvs
  manual_updates      <- vecs$manual_updates
  s1_skip_keys        <- vecs$s1_skip_keys
  rationale           <- vecs$rationale

  all_tax_cols    <- colnames(phyloseq::tax_table(ps))
  taxon_set_col   <- if ("taxa"        %in% all_tax_cols) "taxa"
                     else utils::tail(setdiff(all_tax_cols, "glom_name"), 1)
  common_name_col <- if ("common_name" %in% all_tax_cols) "common_name" else NULL

  if (is.null(tax_rank_cols)) {
    exclude       <- c(taxon_set_col, common_name_col, "glom_name")
    tax_rank_cols <- setdiff(all_tax_cols, exclude)
  }

  taxtab <- as.data.frame(phyloseq::tax_table(ps)@.Data, stringsAsFactors = FALSE)
  all_original_asvs <- rownames(taxtab)

  # Pre-transfer read counts (per ASV, summed across all samples)
  otu_mat_pre <- as.matrix(phyloseq::otu_table(ps))
  reads_pre   <- if (phyloseq::taxa_are_rows(ps)) rowSums(otu_mat_pre)
                 else                              colSums(otu_mat_pre)

  # -- Add `glom_name` column --------------------------------------------------------------------------
  # Default: each ASV gets its own sequence as a unique glom_name so that
  # no two ASVs can be silently merged unless they were explicitly decided on
  # (S1 auto-merge or user merge decision). Only those ASVs will share a name.
  taxtab$glom_name <- rownames(taxtab)

  # -- Apply S2 rank alignment when a manual merge decision was made ------
  # When the user merges an S2 pair, the higher-res ASV's rank columns are
  # overwritten with the lower-res ASV's values so the surviving representative
  # row is internally consistent regardless of which ASV has higher abundance.
  s2_rows     <- flagged[!is.na(flagged$scenario) & flagged$scenario == 2L, ]
  n_corrected <- 0L

  if (nrow(s2_rows) > 0 && length(manual_updates) > 0) {
    for (k in seq_len(nrow(s2_rows))) {
      row      <- s2_rows[k, ]
      keep_asv <- row$keep_asv   # lower-res ASV (hint stored during flagging)
      if (is.na(keep_asv)) next
      to_rename <- if (keep_asv == row$asv_i) row$asv_j else row$asv_i
      # Only apply alignment if both ASVs were merged to the same name
      nm_i <- manual_updates[row$asv_i]; nm_j <- manual_updates[row$asv_j]
      if (is.na(nm_i) || is.na(nm_j) || nm_i != nm_j) next
      if (!to_rename %in% rownames(taxtab)) {
        warning("S2 ASV not found in tax table rownames: ", to_rename); next
      }
      if (!keep_asv %in% rownames(taxtab)) {
        warning("S2 keep_asv not found in tax table rownames: ", keep_asv); next
      }
      taxtab[to_rename, tax_rank_cols]  <- taxtab[keep_asv, tax_rank_cols]
      if (!is.null(common_name_col) && common_name_col %in% colnames(taxtab))
        taxtab[to_rename, common_name_col] <- taxtab[keep_asv, common_name_col]
      if (taxon_set_col %in% colnames(taxtab))
        taxtab[to_rename, taxon_set_col] <- taxtab[keep_asv, taxon_set_col]
      n_corrected <- n_corrected + 1L
    }
  }

  # -- Apply manual updates ---------------------------------------------------
  n_manual <- 0L
  if (length(manual_updates) > 0) {
    for (asv_seq in names(manual_updates)) {
      if (!asv_seq %in% rownames(taxtab)) {
        warning("manual_updates: ASV not found in tax table: ", asv_seq); next
      }
      taxtab[asv_seq, "glom_name"] <- manual_updates[[asv_seq]]
      n_manual <- n_manual + 1L
    }
  }

  methods::slot(ps, "tax_table", check = FALSE) <- phyloseq::tax_table(as.matrix(taxtab))

  # -- Transfer counts and prune discarded ASVs ----------------------------
  # Build a combined discard->keep map from two sources:
  #   (a) S1 auto-merges: the longer ASV is discarded; keep_asv is the shorter.
  #   (b) Decisions from the `decisions` table's `decision` column.
  # All transfers and pruning happen in one pass.
  combined_keep      <- character(0)  # named: discard_asv -> keep_asv
  combined_keep_type <- character(0)  # named: discard_asv -> decision_type

  s1_rows <- flagged[!is.na(flagged$action) &
                       flagged$action == "auto_merge_keep_shorter", ]
  if (nrow(s1_rows) > 0 && "keep_asv" %in% colnames(s1_rows)) {
    for (k in seq_len(nrow(s1_rows))) {
      row_k      <- s1_rows[k, ]
      pair_key_k <- paste0(row_k$asv_i, "|", row_k$asv_j)
      # Skip pairs where the decision overrode the default or chose to keep distinct
      if (pair_key_k %in% s1_skip_keys) next
      keep_k    <- row_k$keep_asv
      if (is.na(keep_k) || !nzchar(keep_k)) next
      discard_k <- if (keep_k == row_k$asv_i) row_k$asv_j else row_k$asv_i
      combined_keep[discard_k]      <- keep_k
      combined_keep_type[discard_k] <- "s1_auto"
    }
  }

  if (length(keep_asvs) > 0) {
    for (d in names(keep_asvs)) {
      k_asv <- keep_asvs[[d]]
      combined_keep[d] <- k_asv
      is_s2 <- nrow(s2_rows) > 0 &&
                any((s2_rows$asv_i == d & s2_rows$asv_j == k_asv) |
                    (s2_rows$asv_i == k_asv & s2_rows$asv_j == d))
      combined_keep_type[d] <- if (is_s2) "s2_rank_align" else "user_merge"
    }
  }

  taxtab_pre_prune <- taxtab  # capture before any pruning -- used for discarded ASV metadata
  n_pruned <- 0L

  # Resolve each discard to its ULTIMATE surviving target up front.
  # combined_keep can contain merge chains (A -> B, B -> C), and summing
  # pairwise in insertion order is order-dependent: if B is processed as a
  # discard before A's counts have been folded into it, A's reads are added
  # to B right before B itself is pruned, and are lost. Resolving the full
  # chain first (union-find style) makes every transfer land on a surviving
  # ASV exactly once, regardless of order -- and lets the decisions table
  # below report the true final representative rather than an intermediate
  # ASV that was itself pruned.
  resolve_final_keep <- function(discard, map) {
    seen <- character(0)
    cur  <- discard
    while (cur %in% names(map)) {
      if (cur %in% seen) {
        warning("Cycle detected in ASV merge chain at: ", cur); return(NA_character_)
      }
      seen <- c(seen, cur)
      cur  <- map[[cur]]
    }
    cur
  }
  final_keep <- if (length(combined_keep) > 0) {
    stats::setNames(vapply(names(combined_keep), resolve_final_keep,
                           character(1), map = combined_keep),
                    names(combined_keep))
  } else character(0)

  if (length(combined_keep) > 0) {
    all_taxa  <- phyloseq::taxa_names(ps)
    rows_mode <- phyloseq::taxa_are_rows(ps)
    otu_mat   <- methods::as(phyloseq::otu_table(ps), "matrix")

    for (discard in names(combined_keep)) {
      if (!discard %in% all_taxa) {
        warning("discard ASV not found in phyloseq: ", discard); next
      }
      keep <- final_keep[[discard]]
      if (is.na(keep) || !keep %in% all_taxa) {
        warning("keep ASV not found in phyloseq (or a chain cycle) for discard: ",
                discard); next
      }
      if (rows_mode) {
        otu_mat[keep, ] <- otu_mat[keep, ] + otu_mat[discard, ]
      } else {
        otu_mat[, keep] <- otu_mat[, keep] + otu_mat[, discard]
      }
    }

    # Only prune discards that were actually transferred above (valid discard
    # AND a resolvable, valid final keep target) -- otherwise their reads
    # would vanish with no representative ASV left to hold them.
    valid_discard <- Filter(function(discard) {
      discard %in% all_taxa &&
        !is.na(final_keep[[discard]]) &&
        final_keep[[discard]] %in% all_taxa
    }, names(combined_keep))
    if (length(valid_discard) > 0) {
      surviving <- setdiff(all_taxa, valid_discard)
      if (rows_mode) {
        otu_mat <- otu_mat[surviving, , drop = FALSE]
      } else {
        otu_mat <- otu_mat[, surviving, drop = FALSE]
      }
      n_pruned <- length(valid_discard)
    }

    new_otu <- phyloseq::otu_table(otu_mat, taxa_are_rows = rows_mode)
    methods::slot(ps, "otu_table", check = FALSE) <- new_otu

    # Sync taxtab to surviving taxa
    surviving_taxa <- if (rows_mode) rownames(otu_mat) else colnames(otu_mat)
    taxtab <- taxtab[rownames(taxtab) %in% surviving_taxa, , drop = FALSE]
  }

  # -- Unresolved: decision rows still blank ------------------------------
  review_actions <- c("flag_lower_resolution", "flag_same_taxon",
                      "flag_split_distinct_names", "flag_taxon_subset",
                      "flag_taxon_overlap", "flag_conflict", "flag_unassigned")
  needs_manual <- flagged$action %in% review_actions
  unresolved   <- flagged[needs_manual & is.na(flagged$decision), ]

  # Counted directly as flagged (review-needing) rows with a non-NA
  # decision -- the exact complement of `unresolved` above. n_manual (the
  # glom_name-update count) isn't a reliable stand-in: a multi-ASV merge
  # group's chosen_name is applied to every member, not just one per pair,
  # so n_manual / 2 over- or under-counts (and can be fractional) whenever
  # a group has more than 2 members.
  n_resolved <- sum(needs_manual & !is.na(flagged$decision), na.rm = TRUE)

  n_auto <- sum(flagged$action == "auto_merge_keep_shorter", na.rm = TRUE)
  if (n_corrected > 0)
    message(n_corrected, " S2 rank alignment(s) applied (merged pairs).")
  if (n_auto > 0)
    message(n_auto, " S1 pair(s) auto-merged: counts transferred, longer ASV pruned.")
  if (n_pruned > 0)
    message(n_pruned, " ASV(s) physically removed (counts transferred to representative).")
  if (n_resolved > 0)
    message(n_resolved, " pair(s) resolved.")
  if (nrow(unresolved) > 0)
    message(nrow(unresolved), " pair(s) still unresolved (see $unresolved).")
  message("Corrections complete. corrected$ps is ready for downstream analysis.")

  # -- Remove glom_name from tax table --------------------------------------
  if ("glom_name" %in% colnames(taxtab)) {
    taxtab$glom_name <- NULL
    methods::slot(ps, "tax_table", check = FALSE) <- phyloseq::tax_table(as.matrix(taxtab))
  }

  # -- Build per-ASV outcomes table -------------------------------------------
  # NOTE: the "rep" (representative) bookkeeping below is built from
  # final_keep -- the fully-resolved end of each merge chain -- rather than
  # combined_keep's immediate targets, so that e.g. a 3-ASV chain A -> B -> C
  # correctly reports C (the surviving ASV) as the representative for BOTH
  # A and B, instead of reporting B (which was itself pruned) as A's
  # representative.
  discarded_asvs   <- names(combined_keep)
  reps_with_merges <- unique(unname(final_keep))

  # Rep -> semicolon-delimited list of ASVs merged into it (direct + chained)
  if (length(final_keep) > 0) {
    rep_to_merged <- tapply(names(final_keep), unname(final_keep),
                            function(x) paste(sort(x), collapse = "; "))
  } else {
    rep_to_merged <- character(0)
  }

  # Rep -> highest-priority decision_type among all contributing merges
  # (each discard's OWN immediate decision type still applies, even though
  # it now rolls up to the chain's final representative).
  priority_map      <- c(s2_rank_align = 3L, user_merge = 2L, s1_auto = 1L)
  rep_decision_type <- character(0)
  for (d in names(final_keep)) {
    k  <- final_keep[d]
    dt <- combined_keep_type[d]
    if (is.na(k)) next
    if (!k %in% names(rep_decision_type) ||
        priority_map[dt] > priority_map[rep_decision_type[k]])
      rep_decision_type[k] <- dt
  }

  # Pre-allocate outcome columns as plain vectors
  n_out             <- length(all_original_asvs)
  out_asv_status    <- character(n_out)
  out_merged_with   <- rep(NA_character_, n_out)
  out_merged_into   <- rep(NA_character_, n_out)
  out_dtype         <- character(n_out)
  out_scenario      <- rep(NA_character_, n_out)
  out_common_name   <- rep(NA_character_, n_out)
  out_taxa          <- rep(NA_character_, n_out)
  out_lowest_rank   <- rep(NA_character_, n_out)
  out_lowest_val    <- rep(NA_character_, n_out)
  out_batch         <- rep(NA_character_, n_out)
  out_reads         <- rep(NA_real_,      n_out)
  out_rationale     <- rep(NA_character_, n_out)
  # Full rank-by-rank taxonomy, not just the deepest resolved rank -- needed
  # so this table alone (the "ASV ledger") can stand in for a live reference
  # object's tax_table, e.g. in apply_projection(). Captured for every ASV
  # regardless of status; only representative/unchanged rows are required
  # (by .validate_asv_ledger()) to actually have any of it non-blank.
  out_taxonomy <- matrix(NA_character_, nrow = n_out, ncol = length(tax_rank_cols),
                         dimnames = list(NULL, tax_rank_cols))

  has_action_col <- "action" %in% colnames(flagged)
  # An ASV that's part of at least one flagged, review-needing pair whose
  # decision is still blank is NOT actually "unchanged" -- its fate is
  # genuinely unsettled, and treating it the same as a never-flagged ASV
  # would hide that from anything (e.g. apply_projection()) trusting this table
  # as a complete record of what's actually resolved.
  unresolved_asvs <- unique(c(unresolved$asv_i, unresolved$asv_j))

  for (idx in seq_along(all_original_asvs)) {
    asv_seq <- all_original_asvs[[idx]]

    # -- Fate / merge columns -----------------------------------------------
    if (asv_seq %in% discarded_asvs) {
      out_asv_status[idx]  <- "dropped after merge"
      # merged_into reports the chain's true final representative (where the
      # reads actually ended up); the flagged-pair lookup below stays on the
      # immediate partner, since `flagged` only records direct pairwise rows.
      out_merged_into[idx] <- unname(final_keep[[asv_seq]])
      out_dtype[idx]       <- unname(combined_keep_type[[asv_seq]])
      rep_asv <- unname(combined_keep[[asv_seq]])
      fl <- flagged[((flagged$asv_i == asv_seq & flagged$asv_j == rep_asv) |
                     (flagged$asv_i == rep_asv  & flagged$asv_j == asv_seq)), ]
      if (nrow(fl) > 0) {
        sc <- fl$scenario[[1]]
        out_scenario[idx] <- if (!is.na(sc)) paste0("S", sc)
                             else if (has_action_col && grepl("conflict",   fl$action[[1]], fixed = TRUE)) "conflict"
                             else if (has_action_col && grepl("unassigned", fl$action[[1]], fixed = TRUE)) "unassigned"
                             else NA_character_
        if (length(rationale) > 0) {
          k1 <- paste0(fl$asv_i[[1]], "|", fl$asv_j[[1]])
          k2 <- paste0(fl$asv_j[[1]], "|", fl$asv_i[[1]])
          rv <- rationale[k1]; if (is.na(rv)) rv <- rationale[k2]
          if (!is.na(rv)) out_rationale[idx] <- unname(rv)
        }
      }

    } else if (asv_seq %in% reps_with_merges) {
      out_asv_status[idx] <- "representative"
      out_merged_with[idx] <- if (asv_seq %in% names(rep_to_merged))
                                unname(rep_to_merged[[asv_seq]])
                              else NA_character_
      out_dtype[idx]     <- if (asv_seq %in% names(rep_decision_type))
                               unname(rep_decision_type[[asv_seq]])
                            else "unchanged"
      discs <- names(combined_keep)[combined_keep == asv_seq]
      sc_vals  <- character(0)
      rat_vals <- character(0)
      for (d in discs) {
        fl <- flagged[((flagged$asv_i == d & flagged$asv_j == asv_seq) |
                       (flagged$asv_i == asv_seq & flagged$asv_j == d)), ]
        if (nrow(fl) > 0) {
          sc <- fl$scenario[[1]]
          sc_label <- if (!is.na(sc)) paste0("S", sc)
                      else if (has_action_col && grepl("conflict",   fl$action[[1]], fixed = TRUE)) "conflict"
                      else if (has_action_col && grepl("unassigned", fl$action[[1]], fixed = TRUE)) "unassigned"
                      else NA_character_
          if (!is.na(sc_label)) sc_vals <- unique(c(sc_vals, sc_label))
          if (length(rationale) > 0) {
            k1 <- paste0(fl$asv_i[[1]], "|", fl$asv_j[[1]])
            k2 <- paste0(fl$asv_j[[1]], "|", fl$asv_i[[1]])
            rv <- rationale[k1]; if (is.na(rv)) rv <- rationale[k2]
            if (!is.na(rv)) rat_vals <- c(rat_vals, unname(rv))
          }
        }
      }
      if (length(sc_vals)  > 0) out_scenario[idx]  <- paste(sc_vals,  collapse = ";")
      if (length(rat_vals) > 0) out_rationale[idx] <- paste(rat_vals, collapse = " | ")

    } else if (asv_seq %in% unresolved_asvs) {
      out_asv_status[idx] <- "unresolved"
      out_dtype[idx] <- "unresolved"
    } else {
      out_asv_status[idx] <- "unchanged"
      out_dtype[idx] <- "unchanged"
    }

    # -- Taxonomy columns ---------------------------------------------------
    tab <- if (asv_seq %in% discarded_asvs) taxtab_pre_prune else taxtab
    if (asv_seq %in% rownames(tab)) {
      if (!is.null(common_name_col) && common_name_col %in% colnames(tab)) {
        v <- tab[[common_name_col]][[which(rownames(tab) == asv_seq)[[1]]]]
        if (!is.na(v) && nzchar(v)) out_common_name[idx] <- as.character(v)
      }
      if (taxon_set_col %in% colnames(tab)) {
        v <- tab[[taxon_set_col]][[which(rownames(tab) == asv_seq)[[1]]]]
        if (!is.na(v) && nzchar(v)) out_taxa[idx] <- as.character(v)
      }
      row_vals <- unlist(tab[asv_seq, tax_rank_cols, drop = FALSE])
      row_vals <- as.character(row_vals)
      out_taxonomy[idx, ] <- row_vals
      non_na   <- which(!is.na(row_vals) & nzchar(row_vals))
      if (length(non_na) > 0) {
        out_lowest_rank[idx] <- tax_rank_cols[[max(non_na)]]
        out_lowest_val[idx]  <- row_vals[[max(non_na)]]
      }
    }

    if (asv_seq %in% names(reads_pre)) out_reads[idx] <- reads_pre[[asv_seq]]
    if (!is.null(source_map) && asv_seq %in% names(source_map)) {
      sv <- source_map[[asv_seq]]
      if (!is.na(sv)) out_batch[idx] <- as.character(sv)
    }
  }

  asv_outcomes <- data.frame(
    asv               = all_original_asvs,
    asv_status        = out_asv_status,
    merged_with       = out_merged_with,
    merged_into       = out_merged_into,
    decision_type     = out_dtype,
    scenario          = out_scenario,
    common_name       = out_common_name,
    taxa              = out_taxa,
    lowest_rank       = out_lowest_rank,
    lowest_rank_value = out_lowest_val,
    batch_source      = out_batch,
    total_reads       = out_reads,
    rationale         = out_rationale,
    stringsAsFactors  = FALSE
  )
  asv_outcomes <- cbind(asv_outcomes, as.data.frame(out_taxonomy, stringsAsFactors = FALSE))

  # -- Fold in prior_outcomes, making the ledger cumulative across rounds --
  # Upsert by ASV (this round's fresh row wins for anything appearing in
  # both), then re-resolve every merged_into chain across the COMBINED
  # table -- a discard recorded in an earlier round may point at an ASV
  # this round has itself since discarded again, and that pointer must
  # always end up at the current, ultimate representative, never a
  # now-defunct intermediate.
  if (!is.null(prior_ledger)) {
    carried_over <- prior_ledger[!prior_ledger$asv %in% asv_outcomes$asv, , drop = FALSE]
    all_cols     <- union(colnames(asv_outcomes), colnames(carried_over))
    for (col in setdiff(all_cols, colnames(asv_outcomes))) asv_outcomes[[col]] <- NA
    for (col in setdiff(all_cols, colnames(carried_over))) carried_over[[col]] <- NA
    asv_outcomes <- rbind(asv_outcomes[, all_cols, drop = FALSE], carried_over[, all_cols, drop = FALSE])

    drop_idx <- which(asv_outcomes$asv_status == "dropped after merge")
    if (length(drop_idx) > 0) {
      drop_map <- stats::setNames(asv_outcomes$merged_into[drop_idx], asv_outcomes$asv[drop_idx])
      resolve_ledger_chain <- function(a) {
        seen <- character(0); cur <- a
        while (cur %in% names(drop_map)) {
          if (cur %in% seen) {
            warning("apply_harmonization: cycle detected in cumulative ledger merge ",
                    "chain at: ", cur); return(NA_character_)
          }
          seen <- c(seen, cur); cur <- unname(drop_map[[cur]])
        }
        cur
      }
      resolved <- vapply(names(drop_map), resolve_ledger_chain, character(1))
      asv_outcomes$merged_into[drop_idx] <- unname(resolved[asv_outcomes$asv[drop_idx]])
    }
  }
  .validate_asv_ledger(asv_outcomes, tax_rank_cols)

  list(
    ps             = ps,
    taxtab         = taxtab,
    unresolved     = unresolved,
    asv_outcomes   = asv_outcomes,
    n_s2_corrected = n_corrected,
    n_manual       = n_resolved
  )
}
