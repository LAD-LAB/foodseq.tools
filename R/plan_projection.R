# R/plan_projection.R
#
# Detects correspondences between a query dataset's ASVs and a reference's
# fixed ASV space -- the apply_projection() analog of plan_harmonization(): pure
# detection, no application. Paired with apply_projection() (see there),
# mirroring plan_harmonization()/apply_harmonization()'s own split -- kept
# separate for the same reason: detection and execution are different
# enough operations (very different data flow, very different risk
# profile) that one function doing both invites bugs.
#
# The reference is never a live phyloseq object -- it's ps_ref_outcomes,
# the cumulative ASV ledger apply_harmonization() produces (its
# $asv_outcomes, enriched with full taxonomy and made cumulative across
# prior_outcomes rounds -- see there). This closes a real gap a live-object
# design has no way to detect: a reference object can silently contain an
# ASV that's actually a stale, not-yet-pruned discard from an unfinished
# harmonization round, and a query could end up matched against it with no
# warning. The ledger makes the reference's own harmonization history an
# explicit, checkable input instead of an invisible assumption --
# plan_projection() refuses to run at all if the ledger says anything is
# still unresolved, and silently redirects any exact match against a known
# historical discard to its true current representative rather than
# treating it as a fresh, ambiguous relationship.
#
# Reuses compare_asvs() (plan_harmonization.R) for detection against the
# reference's reconstructed, CURRENT final ASV set -- it never mutates its
# inputs.

# Swaps every paired _i/_j column for the given row mask, so asv_i always
# ends up on whichever side the caller designates -- used to normalize
# compare_asvs()'s output (which orders i/j by internal sort order, with no
# notion of "reference" vs "new dataset") into a consistent x-is-i,
# y-is-j convention. keep_asv/keep_name/keep_taxonomy hold an actual ASV
# value rather than an i/j-relative reference, so they need no swap.
.swap_ij <- function(df, rows) {
  if (!any(rows)) return(df)
  paired <- c("asv", "len", "source", "reads_total", "reads_max", "reads_pct",
             "reads_prev", "rank", "deepest_name", "taxon", "common_name")
  for (p in paired) {
    ci <- paste0(p, "_i"); cj <- paste0(p, "_j")
    if (ci %in% colnames(df) && cj %in% colnames(df)) {
      tmp        <- df[rows, ci]
      df[rows, ci] <- df[rows, cj]
      df[rows, cj] <- tmp
    }
  }
  df
}

# Relabels compare_asvs()'s generic per-ASV `source_i`/`source_j` values --
# comma-joined tokens drawn from whatever `source_names` was actually passed
# in (default "x"/"y") -- into "ps_ref"/"ps_query" for anything user-facing.
# asv_i/asv_j's roles are a fixed invariant here (i is always the
# reference-side candidate, j is always the query ASV in question) no
# matter what source_names was set to, so $decisions always reads that way
# too. Token-exact (split on "," and match a whole token), never a plain
# string replace -- "ps_query" would otherwise collide with any other
# token containing it as a literal substring.
.relabel_sources <- function(x, src_x, src_y) {
  if (is.null(x)) return(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    toks <- strsplit(s, ",", fixed = TRUE)[[1]]
    toks <- ifelse(toks == src_x, "ps_ref", ifelse(toks == src_y, "ps_query", toks))
    paste(toks, collapse = ",")
  }, character(1), USE.NAMES = FALSE)
}

# Columns apply_harmonization() always writes onto $asv_outcomes that are
# NOT taxonomy -- used to infer tax_rank_cols from whatever's left over when
# the caller doesn't specify them explicitly.
.ledger_fixed_cols <- c("asv", "asv_status", "merged_with", "merged_into",
                        "decision_type", "scenario", "common_name", "taxa",
                        "lowest_rank", "lowest_rank_value", "batch_source",
                        "total_reads", "rationale")

# Turns a validated ASV ledger into everything plan_projection() actually
# needs to run detection: a reconstructed tax_table of the reference's
# CURRENT final ASVs (representative/unchanged rows only), and an
# exact-match ghost map (discard sequence -> its final representative,
# already fully resolved by apply_harmonization() itself -- see its
# prior_outcomes handling) for the rare case a new batch reproduces a
# sequence that used to exist verbatim in the reference before its own
# harmonization merged it away.
#
# Deliberately exact-match only, not a second full compare_asvs() detection
# pass against discards -- extending detection to ghosts would need full
# taxonomy on every discard too, plus a post-hoc redirect step for anything
# that classifies against one. A query ASV merely *similar* to a pure
# historical discard (never identical) just falls through to "no
# correspondence found" instead -- a graceful miss, not a wrong answer, for
# a case rare enough not to justify the added complexity.
.ledger_to_reference <- function(ledger_df, tax_rank_cols) {
  final_rows <- ledger_df[!is.na(ledger_df$asv_status) &
                            ledger_df$asv_status %in% c("representative", "unchanged"), ]
  if (nrow(final_rows) == 0) {
    stop("plan_projection: ps_ref_outcomes has no representative/unchanged ASVs at all -- ",
         "nothing to project ps_query onto.")
  }
  ps_ref_asvs <- final_rows$asv

  # Extra, non-rank columns worth carrying into $ps_projected's tax_table if
  # present and not all-blank -- mirrors what a live reference object's own
  # tax_table would have had.
  extra_cols <- intersect(c("common_name", "taxa"), colnames(final_rows))
  extra_cols <- extra_cols[vapply(extra_cols, function(cn)
    any(!is.na(final_rows[[cn]]) & nzchar(final_rows[[cn]])), logical(1))]

  taxtab_cols <- c(tax_rank_cols, extra_cols)
  taxtab_mat  <- as.matrix(final_rows[, taxtab_cols, drop = FALSE])
  rownames(taxtab_mat) <- ps_ref_asvs

  drop_rows <- ledger_df[!is.na(ledger_df$asv_status) &
                          ledger_df$asv_status == "dropped after merge", ]
  ghost_map <- stats::setNames(drop_rows$merged_into, drop_rows$asv)

  asvs      <- ps_ref_asvs
  dummy_otu <- matrix(0, nrow = 1, ncol = length(asvs),
                      dimnames = list("_ps_ref_no_sample_data_", asvs))
  shell <- phyloseq::phyloseq(phyloseq::otu_table(dummy_otu, taxa_are_rows = FALSE),
                              phyloseq::tax_table(taxtab_mat))

  list(ps_ref_asvs = ps_ref_asvs, taxtab = phyloseq::tax_table(taxtab_mat),
      ghost_map = ghost_map, shell = shell)
}

# Translates a plan_projection()-shaped decisions data.frame into the
# c(ps_query_asv = ps_ref_asv) `keep_asvs` map plan_projection()/
# apply_projection() work with internally. Deliberately narrower than
# apply_harmonization.R's .decisions_to_vectors() -- and NOT a reuse of it
# -- for two reasons specific to projection:
#
#   1. Vocabulary: ps_ref must always be the side kept (it's the fixed
#      coordinate system), so only "merge_keep_i" (or, for a hand-edited
#      CSV, the literal ps_ref ASV sequence itself) and "keep_distinct" are
#      meaningful; "merge_keep_j" (or any other literal value) is rejected
#      outright, with a specific, row-identifying error.
#   2. Conflict handling: .decisions_to_vectors() treats the same discard
#      ASV mapping to two different keep targets across rows as a hard,
#      unrecoverable stop(). Projection deliberately does NOT do that here
#      -- a query ASV independently flagged against two different ps_ref
#      candidates is a genuine, recoverable ambiguity, and
#      .resolve_projection_pairs()'s own ambiguity logic already forces
#      such pairs back to review rather than erroring the whole call. This
#      function stays a dumb per-row translator so that logic keeps
#      working exactly as it did.
#
# Also enforces, per row, that the ps_ref ASV a decision points to still
# exists among ps_ref's CURRENT final ASVs (per the ledger) -- a decisions
# table can outlive changes to ps_ref's own harmonization since it was made,
# and a stale asv_i must never be silently kept.
.projection_decisions_to_keep_asvs <- function(decisions_df, ps_ref_asvs) {
  keep_asvs <- character(0)
  if (is.null(decisions_df) || nrow(decisions_df) == 0) return(keep_asvs)
  stopifnot(all(c("asv_i", "asv_j") %in% colnames(decisions_df)))
  dec_col <- if ("decision" %in% colnames(decisions_df)) decisions_df$decision
             else rep(NA_character_, nrow(decisions_df))

  for (i in seq_len(nrow(decisions_df))) {
    dec <- dec_col[i]
    if (is.na(dec) || !nzchar(dec)) next   # not yet decided
    asv_i <- decisions_df$asv_i[i]; asv_j <- decisions_df$asv_j[i]

    if (identical(dec, "merge_keep_i") || identical(dec, asv_i)) {
      if (!asv_i %in% ps_ref_asvs) {
        stop("plan_projection: decision for query ASV '", asv_j, "' points to ps_ref ASV '", asv_i,
             "', which is no longer among ps_ref's current final ASVs (ps_ref's own harmonization ",
             "may have changed since these decisions were made). Re-run plan_projection() fresh ",
             "against the current ps_ref_outcomes.")
      }
      keep_asvs[asv_j] <- asv_i
    } else if (identical(dec, "keep_distinct")) {
      next   # explicit "no" -- nothing to add
    } else {
      stop("plan_projection: invalid decision '", dec, "' for pair (ps_ref ASV '", asv_i,
           "', ps_query ASV '", asv_j, "') -- only 'merge_keep_i' (or, equivalently, the literal ",
           "ps_ref ASV sequence itself) and 'keep_distinct' are meaningful decisions here; ps_ref's ",
           "ASV must always be the one kept, so 'merge_keep_j' (or any other value) is rejected.")
    }
  }
  keep_asvs
}

# Combines auto_map (S1, from detection -- fixed once, computed by
# plan_projection()) with whatever's in `decisions_df` (prior_decisions at
# planning time, or plan$decisions at apply time) to determine which pairs
# are still unresolved. Shared by plan_projection() and apply_projection()
# so this logic -- especially ambiguity detection, easy to get subtly wrong
# -- lives in exactly one place rather than being duplicated across the
# plan/apply boundary.
.resolve_projection_pairs <- function(auto_map, needs_review, decisions_df, ps_ref_asvs) {
  keep_asvs <- .projection_decisions_to_keep_asvs(decisions_df, ps_ref_asvs)

  # A query ASV assigned to more than one DISTINCT ps_ref ASV -- across
  # auto-resolved (S1) and decisions-derived keep_asvs combined -- is a
  # genuine contradiction. Never silently resolved by picking one; ambiguous
  # ASVs are always forced back to review instead.
  all_decided <- c(auto_map, keep_asvs)
  ambiguous_y <- if (length(all_decided) > 0) {
    by_y <- split(unname(all_decided), names(all_decided))
    names(by_y)[vapply(by_y, function(v) length(unique(v)) > 1, logical(1))]
  } else character(0)

  auto_map_applied  <- auto_map[!names(auto_map) %in% ambiguous_y]
  keep_asvs_applied <- keep_asvs[!names(keep_asvs) %in% ambiguous_y]

  # A pair is resolved either by merging (covered by keep_asvs_applied) OR by
  # an explicit "keep_distinct" decision -- the latter is a genuine, final
  # answer ("no, these are not the same") and must count as covered too, not
  # just "didn't get merged" (which would otherwise leave every
  # keep_distinct decision reappearing in $decisions forever).
  kd_pair_keys <- character(0)
  if (!is.null(decisions_df) && all(c("asv_i", "asv_j", "decision") %in% colnames(decisions_df))) {
    kd_rows <- decisions_df[!is.na(decisions_df$decision) & decisions_df$decision == "keep_distinct", ]
    if (nrow(kd_rows) > 0) kd_pair_keys <- paste0(kd_rows$asv_i, "|", kd_rows$asv_j)
  }
  needs_review_keys <- paste0(needs_review$asv_i, "|", needs_review$asv_j)
  covered <- (needs_review$asv_j %in% names(keep_asvs_applied) | needs_review_keys %in% kd_pair_keys) &
             !needs_review$asv_j %in% ambiguous_y

  list(keep_asvs_applied = keep_asvs_applied, auto_map_applied = auto_map_applied,
      ambiguous_y = ambiguous_y, covered = covered)
}

# Builds the $decisions table stored on `plan`/`result` -- ALWAYS the FULL
# set of this round's cross-pairs needing review (not just the still-blank
# ones), with decision/rationale populated from `decisions_df` (prior_decisions
# at planning time, or plan$decisions at apply time) wherever a value is
# already known, blank otherwise. Returning only the still-unresolved subset
# would silently discard already-resolved decisions the moment nothing is
# left pending -- exactly the bug that motivated this: apply_projection()
# needs the FULL decision set to build keep_asvs correctly, not just
# whatever's still outstanding. resolve_projection_interactive() is
# responsible for only stepping the user through the genuinely-blank rows
# (mirroring plan_harmonization()'s own "an exact-repeat of an
# already-resolved pair never needs review again" pattern) -- it still
# receives, and returns, this same full table.
#
# Returns NULL only when there is truly nothing to review this round at all
# (needs_review itself has zero rows) -- not merely when everything
# happens to already be resolved.
.build_projection_decisions_out <- function(needs_review, decisions_df, covered, ambiguous_y,
                                            is_s1, cross, ps_ref_asvs, ps_query_asvs_all, ps_query,
                                            src_x, src_y, verbose) {
  if (nrow(needs_review) == 0) return(NULL)

  decisions_out <- needs_review
  decisions_out$decision  <- NA_character_
  decisions_out$rationale <- NA_character_

  if (!is.null(decisions_df) && all(c("asv_i", "asv_j") %in% colnames(decisions_df))) {
    prior_keys <- paste0(decisions_df$asv_i, "|", decisions_df$asv_j)
    out_keys   <- paste0(decisions_out$asv_i, "|", decisions_out$asv_j)
    idx        <- match(out_keys, prior_keys)
    matched    <- !is.na(idx)
    if ("decision" %in% colnames(decisions_df))
      decisions_out$decision[matched]  <- decisions_df$decision[idx[matched]]
    if ("rationale" %in% colnames(decisions_df))
      decisions_out$rationale[matched] <- decisions_df$rationale[idx[matched]]
  }

  is_ambig_row <- decisions_out$asv_j %in% ambiguous_y
  if (any(is_ambig_row)) {
    decisions_out$note[is_ambig_row] <- paste0(
      "AMBIGUOUS: this ps_query ASV is also flagged toward a different, unmerged ps_ref ASV -- ",
      "resolve consistently. ", decisions_out$note[is_ambig_row]
    )
  }

  s1_ambiguous_rows <- cross[is_s1 & cross$asv_j %in% ambiguous_y, , drop = FALSE]
  n_pending <- sum(!covered) + nrow(s1_ambiguous_rows)
  if (verbose && n_pending > 0) {
    message(n_pending, " pair(s) need review before projecting. asv_i is ",
            "always the reference (ps_ref) ASV and asv_j is always the ps_query ASV in question -- ",
            "only 'merge_keep_i' is a meaningful decision; 'merge_keep_j' will be rejected. Edit ",
            "$decisions$decision (in R, or via write.csv()/the save_path argument + a spreadsheet ",
            "editor), then apply_projection(plan, ps_query), or re-run plan_projection() with ",
            "prior_decisions = <the edited result> for a fresh ambiguity check.")
  } else if (verbose) {
    message("0 pair(s) need review -- every cross-pair this round was already resolved by ",
            "prior_decisions.")
  }

  for (col in intersect(c("source_i", "source_j"), colnames(decisions_out))) {
    decisions_out[[col]] <- .relabel_sources(decisions_out[[col]], src_x, src_y)
  }

  hash_asvs <- function(asvs) vapply(asvs, function(x)
    digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))
  n_out <- nrow(decisions_out)
  decisions_out$ps_ref_taxa_hashes   <- rep(paste(unname(hash_asvs(ps_ref_asvs)), collapse = ";"), n_out)
  decisions_out$ps_query_taxa_hashes <- rep(paste(unname(hash_asvs(ps_query_asvs_all)), collapse = ";"), n_out)
  decisions_out$ps_query_samples     <- rep(paste(phyloseq::sample_names(ps_query), collapse = ";"), n_out)

  decisions_out
}

# Refuses to run at all when ps_ref and ps_query share any sample names --
# projection assumes the two are always cleanly distinguishable by name
# (project_pca()'s together = TRUE view, and any downstream join, depend on
# this). Only checkable when a LIVE reference object is actually available
# -- true for the common case (the full apply_harmonization() result passed
# as ps_ref_outcomes), but not when ps_ref_outcomes is a bare $asv_outcomes
# data.frame or CSV path (no live reference samples to check against at
# all). Skipped with a message() in that case, not silently -- so the gap
# is visible rather than assumed away.
.check_sample_name_uniqueness <- function(ps_ref_outcomes, ps_query, verbose) {
  live_ref <- if (is.list(ps_ref_outcomes) && !is.data.frame(ps_ref_outcomes) &&
                  "ps" %in% names(ps_ref_outcomes) && inherits(ps_ref_outcomes$ps, "phyloseq")) {
    ps_ref_outcomes$ps
  } else NULL

  if (is.null(live_ref)) {
    if (verbose) message("plan_projection: ps_ref_outcomes has no live reference phyloseq object ",
                         "(a bare $asv_outcomes data.frame or CSV path was supplied) -- skipping the ",
                         "sample-name uniqueness check. Pass apply_harmonization()'s full result as ",
                         "ps_ref_outcomes to enable it.")
    return(invisible(NULL))
  }

  dupes <- intersect(phyloseq::sample_names(live_ref), phyloseq::sample_names(ps_query))
  if (length(dupes) > 0) {
    stop("plan_projection: ps_ref and ps_query share ", length(dupes), " sample name(s) -- ",
         "e.g. '", dupes[1], "'. Reference and query sample names must be disjoint (rename before ",
         "retrying) -- there is no override, since project_pca() and any downstream join by sample ",
         "name assume this.")
  }
  invisible(NULL)
}

# Checks `prior_decisions`' embedded checkpoint columns (written by
# plan_projection() itself, see .build_projection_decisions_out()) against
# the CURRENT ps_ref_asvs and ps_query. Asymmetric by design: ps_ref's
# final ASV set is the fixed coordinate system every decision refers to, so
# drift there stop()s outright; ps_query-side drift (its own taxa, or its
# samples) only ever warns and proceeds, since decisions apply per-ASV. A
# `prior_decisions` predating this feature (no hash/sample columns at all),
# or NULL entirely, skips the checkpoint.
.project_query_checkpoint <- function(prior_df, ps_ref_asvs, ps_query) {
  if (is.null(prior_df) || nrow(prior_df) == 0) return(invisible(NULL))
  has_ps_ref_hash   <- "ps_ref_taxa_hashes" %in% colnames(prior_df)
  has_ps_query_hash <- "ps_query_taxa_hashes" %in% colnames(prior_df)
  has_ps_query_samp <- "ps_query_samples"     %in% colnames(prior_df)
  if (!has_ps_ref_hash && !has_ps_query_hash && !has_ps_query_samp) return(invisible(NULL))

  hash_asvs <- function(asvs) vapply(asvs, function(x)
    digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))

  if (has_ps_ref_hash) {
    recorded <- strsplit(prior_df$ps_ref_taxa_hashes[[1]], ";", fixed = TRUE)[[1]]
    current  <- unname(hash_asvs(ps_ref_asvs))
    if (!setequal(recorded, current)) {
      stop("plan_projection: ps_ref's final ASV set has changed since `prior_decisions` was made ",
           "(ps_ref_taxa_hashes no longer matches ps_ref_outcomes' current representative/unchanged ",
           "ASVs) -- ps_ref is the fixed coordinate system every decision here refers to, so this ",
           "risks silently mapping ps_query's reads onto the wrong ASVs. Re-run plan_projection() ",
           "fresh, without a stale prior_decisions, against the current ps_ref_outcomes.")
    }
  }
  if (has_ps_query_hash) {
    recorded <- strsplit(prior_df$ps_query_taxa_hashes[[1]], ";", fixed = TRUE)[[1]]
    current  <- unname(hash_asvs(phyloseq::taxa_names(ps_query)))
    if (!setequal(recorded, current)) {
      warning("plan_projection: ps_query's taxa differ from when `prior_decisions` was made -- ",
              "proceeding, since decisions apply per-ASV: a newly-introduced ASV simply won't have ",
              "a decision yet, and a removed one's decision is simply unused.", call. = FALSE)
    }
  }
  if (has_ps_query_samp) {
    recorded <- strsplit(prior_df$ps_query_samples[[1]], ";", fixed = TRUE)[[1]]
    current  <- phyloseq::sample_names(ps_query)
    if (!setequal(recorded, current)) {
      warning("plan_projection: ps_query's samples differ from when `prior_decisions` was made -- ",
              "proceeding, since decisions apply per-ASV and are otherwise unaffected by which ",
              "samples are present.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# resolve_projection_interactive -----------------------------------------------
#
# Reviews plan_projection()'s in-progress $decisions pairs via a Shiny
# gadget (RStudio Viewer pane) -- the projection analog of
# resolve_conflicts_interactive() for plan_harmonization(). Same relationship
# between the two, and lives here for the same reason resolve_conflicts_
# interactive() lives inside plan_harmonization.R rather than its own file:
# this is a purely internal helper, called automatically from within
# plan_projection() itself (never exported, never called directly by users)
# -- see plan_projection()'s own Step 2 below. Not a separate step the
# caller invokes by hand.
#
# Deliberately much simpler than resolve_conflicts_interactive()'s gadget:
# projection only ever asks one question per pair -- does this ps_query ASV
# belong with ps_ref's ASV (asv_i), or is it genuinely distinct? There's no
# naming, grouping, or sequence-selection phase here, since ps_ref's ASV and
# identity are always kept fixed by projection's own design -- a merge never
# renames or reselects anything, it only decides whether ps_query's reads
# land on ps_ref's ASV or not. No scenario-filter phase either, for the same
# reason: with only two meaningful answers per pair, there's nothing worth
# bulk-skipping by scenario type the way harmonization's much larger
# decision space needs.
#
# Falls back to a warning if shiny/miniUI/DT are not installed, to a message
# if there's nothing to review, or to returning `decisions_df` unchanged if
# the session isn't interactive -- in every one of those cases,
# plan_projection() just returns with those rows still blank, exactly as if
# this step didn't exist; there is no error path here that could abort
# planning altogether.
#
# Called as its own statement inside plan_projection() (rather than nested
# as an argument) purely for readability -- unlike an earlier theory
# recorded here, nesting was never actually the cause of the "object '...'
# not found" errors seen during development. The real cause was
# shiny::runGadget()'s viewer: shiny::dialogViewer() has a confirmed,
# unresolved RStudio bug that corrupts the calling function's own
# execution context on close (https://github.com/rstudio/rstudio/issues/11714),
# surfacing as whatever code happens to run next throwing a bizarre,
# effectively-random error unrelated to its own logic. Fixed below by using
# shiny::paneViewer() instead, which that issue confirms is NOT affected.
#
# Takes/returns a plain decisions data.frame -- no multi-input-type
# dispatch, since plan_projection() is the only caller and always has one
# already built. `save_path` (writing the final, reviewed decisions to a
# CSV when Finish/Exit is clicked, independent of the R return value) is a
# safety net for a real observed issue: some execution contexts (a Shiny
# gadget launched via RStudio's per-chunk notebook execution, as opposed to
# running the same code at the console) can lose a blocking gadget's return
# value on the way back out even though the review itself completed
# correctly inside the gadget. Anything left blank when Exit is clicked
# before reaching the end (or skipped) keeps decision = NA, exactly as if
# it had never been reviewed, so it simply reappears in $decisions the next
# time plan_projection()/apply_projection() is called.
resolve_projection_interactive <- function(decisions_df, save_path = NULL) {

  stopifnot(is.data.frame(decisions_df))
  stopifnot(all(c("asv_i", "asv_j") %in% colnames(decisions_df)))
  rownames(decisions_df) <- NULL

  if (nrow(decisions_df) == 0) {
    message("No pairs require review.")
    return(decisions_df)
  }

  if (!interactive()) {
    warning("resolve_projection_interactive() called in a non-interactive ",
            "session. Returning `decisions_df` unchanged.")
    return(decisions_df)
  }

  if (is.null(save_path)) {
    save_path <- file.path(tempdir(), paste0(
      "projection_decisions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  }

  if (!requireNamespace("shiny",  quietly = TRUE) ||
      !requireNamespace("miniUI", quietly = TRUE) ||
      !requireNamespace("DT",     quietly = TRUE)) {
    message("Install required packages for the interactive gadget:\n",
            "  install.packages(c('shiny','miniUI','DT'))")
    return(decisions_df)
  }

  if (!"decision"  %in% colnames(decisions_df)) decisions_df$decision  <- NA_character_
  if (!"rationale" %in% colnames(decisions_df)) decisions_df$rationale <- NA_character_

  # -- Scenario metadata (same vocabulary as resolve_conflicts_interactive()) --
  action_to_scn <- c(
    flag_lower_resolution     = "S2",  flag_same_taxon      = "S3",
    flag_split_distinct_names = "S1b", flag_taxon_subset    = "S4",
    flag_taxon_overlap        = "S5",  flag_conflict        = "conflict",
    flag_unassigned           = "unassigned"
  )
  scn_desc <- c(
    S1b = "Substring -- same LCA, different species sets",
    S2 = "Substring -- different resolution, shared branch",
    S3 = "Non-substring -- identical species sets",
    S4 = "Non-substring -- species set subset",
    S5 = "Non-substring -- overlapping species sets",
    conflict = "Substring -- conflicting taxonomy at same rank",
    unassigned = "Substring -- missing taxonomy"
  )
  scn_colour <- c(
    S1b = "#fcf8e3", S2 = "#fcf8e3", S3 = "#d9edf7", S4 = "#f2dede",
    S5 = "#f5f5f5", conflict = "#f2dede", unassigned = "#f5f5f5"
  )

  pair_short_desc <- function(row) {
    sub_txt <- if (isTRUE(row$is_substring)) "asv_j is a substring of asv_i (or vice versa)"
               else "non-nested sequences"
    switch(row$action,
      flag_split_distinct_names = paste0(sub_txt, " -- same formal taxon, different species sets"),
      flag_lower_resolution     = {
        ri <- if (!is.null(row$rank_i) && !is.na(row$rank_i) && nzchar(row$rank_i)) row$rank_i else "?"
        rj <- if (!is.null(row$rank_j) && !is.na(row$rank_j) && nzchar(row$rank_j)) row$rank_j else "?"
        paste0(sub_txt, sprintf(" -- ps_ref at %s level, ps_query at %s level", ri, rj))
      },
      flag_same_taxon           = paste0(sub_txt, " -- identical species sets"),
      flag_taxon_subset         = {
        si <- trimws(strsplit(if (!is.null(row$taxon_i) && !is.na(row$taxon_i) && nzchar(row$taxon_i))
                                row$taxon_i else "", ";")[[1]])
        sj <- trimws(strsplit(if (!is.null(row$taxon_j) && !is.na(row$taxon_j) && nzchar(row$taxon_j))
                                row$taxon_j else "", ";")[[1]])
        if (length(si) > 0 && length(sj) > 0 && all(si %in% sj))
          "Non-nested -- ps_ref's species set is a subset of ps_query's"
        else
          "Non-nested -- ps_query's species set is a subset of ps_ref's"
      },
      flag_taxon_overlap        = {
        oc <- if (!is.null(row$overlap_coef) && !is.na(row$overlap_coef))
                sprintf(" (overlap = %.2f)", row$overlap_coef) else ""
        paste0(sub_txt, " -- partially overlapping species sets", oc)
      },
      flag_conflict             = {
        rk <- if (!is.null(row$rank_i) && !is.na(row$rank_i) && nzchar(row$rank_i)) row$rank_i else "same rank"
        paste0(sub_txt, sprintf(" -- conflicting %s assignments", rk))
      },
      flag_unassigned           = paste0(sub_txt, " -- one or both ASVs lack taxonomy"),
      sub_txt
    )
  }

  # -- ASV metadata lookup (mirrors resolve_conflicts_interactive()'s) --------
  has_reads <- "reads_total_i" %in% colnames(decisions_df)
  meta_cols_i <- c("asv_i", "len_i",
                   if (has_reads) c("reads_total_i", "reads_max_i", "reads_pct_i", "reads_prev_i"),
                   "taxon_i", "common_name_i", "rank_i", "deepest_name_i")
  meta_cols_j <- c("asv_j", "len_j",
                   if (has_reads) c("reads_total_j", "reads_max_j", "reads_pct_j", "reads_prev_j"),
                   "taxon_j", "common_name_j", "rank_j", "deepest_name_j")
  base_names <- c("asv", "len",
                  if (has_reads) c("reads_total", "reads_max", "reads_pct", "reads_prev"),
                  "taxon", "common_name", "rank", "deepest_name")

  mi <- decisions_df[, intersect(meta_cols_i, colnames(decisions_df)), drop = FALSE]
  mj <- decisions_df[, intersect(meta_cols_j, colnames(decisions_df)), drop = FALSE]
  colnames(mi) <- base_names[seq_len(ncol(mi))]
  colnames(mj) <- base_names[seq_len(ncol(mj))]
  asv_meta <- unique(rbind(mi, mj))
  asv_meta <- asv_meta[!duplicated(asv_meta$asv), ]
  rownames(asv_meta) <- asv_meta$asv

  # -- CSS (same visual language as resolve_conflicts_interactive()) ---------
  gadget_css <- shiny::tags$style(shiny::HTML(
    'body,.container-fluid{padding:0 8px;font-size:13px;}
     .hdr{border-bottom:2px solid #337ab7;padding-bottom:5px;margin-bottom:10px;padding-right:115px;}
     .hdr h4{margin:0;display:inline;} .prog{float:right;font-size:11px;color:#888;margin-top:4px;}
     .apanel{border:1px solid #ddd;border-radius:4px;padding:8px;background:#fafafa;
              height:260px;overflow-y:auto;}
     .apanel.ref{border-color:#337ab7;background:#f0f7ff;}
     .seq{font-family:monospace;font-size:10px;color:#666;word-break:break-all;margin:2px 0 4px;}
     .splist{font-size:11px;padding-left:16px;margin:3px 0;}
     .role{background:#d9edf7;border:1px solid #bce8f1;border-radius:3px;
            padding:2px 5px;font-size:10px;margin-bottom:4px;display:inline-block;}
     .scnbar{border-radius:3px;padding:5px 8px;margin-bottom:6px;font-size:12px;font-weight:bold;}
     .note{font-size:11px;color:#666;font-style:italic;margin-bottom:5px;}
     .pair-subdesc{font-size:12px;color:#333;margin-bottom:6px;padding:2px 0;}
     .btn-choice{display:block; width:fit-content; margin:2px 0;} hr.thin{margin:6px 0;}
     .pair-wrap{display:flex; flex-direction:column;}
     .pair-context{flex:0 0 auto;}
     .pair-actions{flex:0 0 auto; border-top:2px solid #e0e0e0; padding-top:8px; background:#fff;}
     .global-btn{position:fixed; top:8px; right:12px; z-index:9999;}
     #decisions_dt table.dataTable td{max-width:220px !important; overflow:hidden !important;
       text-overflow:ellipsis !important; white-space:nowrap !important; cursor:text;}'
  ))

  asv_panel <- function(asv_seq, label, role_cls = NULL) {
    if (!asv_seq %in% rownames(asv_meta))
      return(shiny::div(class = "apanel", label, ": not found"))
    a <- asv_meta[asv_seq, ]
    sp_ul <- NULL
    if (!is.na(a$taxon) && nzchar(a$taxon)) {
      sp    <- trimws(strsplit(a$taxon, ";")[[1]])
      sp_ul <- shiny::tagList(
        shiny::tags$strong(sprintf("Species (%d):", length(sp))),
        shiny::tags$ul(class = "splist", lapply(sp, shiny::tags$li)))
    }
    reads_div <- NULL
    if (has_reads && "reads_total" %in% colnames(asv_meta) && !is.na(a$reads_total))
      reads_div <- shiny::div(class = "text-muted",
        sprintf("Reads: %s total | %s max | %.1f%% | %.1f%% prev",
                format(a$reads_total, big.mark = ","), format(a$reads_max, big.mark = ","),
                a$reads_pct, a$reads_prev))
    shiny::div(class = paste("apanel", if (!is.null(role_cls)) role_cls else ""),
      shiny::tags$h6(style = "margin-top:0;", label),
      shiny::div(class = "seq", substr(asv_seq, 1, 55),
                 if (nchar(asv_seq) > 55) "..." else ""),
      shiny::div(shiny::tags$strong(sprintf("%d bp", a$len))),
      reads_div, shiny::tags$hr(class = "thin"),
      shiny::div(sprintf("%s: %s",
                         if (!is.na(a$rank)) a$rank else "unassigned",
                         if (!is.na(a$deepest_name)) a$deepest_name else "(unassigned)")),
      shiny::div(shiny::tags$em("Common name: "),
                 if (!is.na(a$common_name)) a$common_name else "(none)"),
      sp_ul)
  }

  # `decisions_df` may now include rows already resolved by an earlier
  # prior_decisions round (plan_projection() carries the FULL table forward
  # deliberately -- see its own source comments) -- only step the user
  # through the genuinely-blank ones, mirroring plan_harmonization()'s own
  # "an exact-repeat of an already-resolved pair never needs review again"
  # pattern. The full table (already-resolved rows untouched) is still what
  # gets returned at the end.
  pending_idx <- which(is.na(decisions_df$decision))
  n <- length(pending_idx)
  if (n == 0) {
    message("No pairs require review -- every row already has a decision.")
    return(decisions_df)
  }
  all_cols <- colnames(decisions_df)
  default_cols <- intersect(
    c("asv_i", "asv_j", "action", "scenario", "decision", "rationale",
      "taxon_i", "taxon_j", "common_name_i", "common_name_j", "overlap_coef"),
    all_cols)

  # -- UI -----------------------------------------------------------------
  ui <- shiny::fluidPage(
    gadget_css,
    shiny::div(class = "global-btn",
               shiny::actionButton("show_table", "Decisions table",
                                   class = "btn-default btn-sm",
                                   icon  = shiny::icon("table"))),
    shiny::uiOutput("main_ui")
  )

  # -- Server ---------------------------------------------------------------
  server <- function(input, output, session) {

    rv <- shiny::reactiveValues(phase = "summary", k = 1L,
                                rows = decisions_df, snaps = list())
    session$userData$rv <- rv   # exposed for shiny::testServer(); not otherwise used

    output$main_ui <- shiny::renderUI({
      switch(rv$phase,
        summary = shiny::tagList(
          shiny::div(class = "hdr", shiny::tags$h4("Projection -- Decisions")),
          shiny::p(class = "text-muted small",
                   sprintf("%d pair(s) need review. asv_i is always the reference (ps_ref) ASV; asv_j is always the ps_query ASV in question.", n)),
          shiny::actionButton("summary_view_table", "View full table", class = "btn-default"),
          shiny::actionButton("summary_start", "Start reviewing ->",
                              class = "btn-primary", style = "margin-left:8px;")
        ),
        review = {
          k <- rv$k
          if (k > n) {
            decided <- sum(!is.na(rv$rows$decision[pending_idx]))
            return(shiny::tagList(
              shiny::div(class = "alert alert-success",
                         shiny::tags$strong("Review complete."),
                         sprintf(" %d / %d pair(s) decided.", decided, n)),
              shiny::actionButton("finish", "Finish -- return results", class = "btn-primary btn-lg")))
          }
          row <- rv$rows[pending_idx[k], ]
          scn <- action_to_scn[row$action]; if (is.na(scn)) scn <- if (!is.na(row$action)) row$action else "?"
          bg  <- if (!is.na(scn_colour[scn])) scn_colour[scn] else "#f5f5f5"
          desc <- if (!is.na(scn_desc[scn])) scn_desc[scn] else scn
          can_back <- k > 1 && !is.null(tryCatch(rv$snaps[[k - 1]], error = function(e) NULL))
          shiny::tagList(
            shiny::div(class = "hdr", shiny::tags$h4(sprintf("Pair %d / %d", k, n)),
                       shiny::div(class = "prog", desc)),
            shiny::div(class = "pair-wrap",
              shiny::div(class = "pair-context",
                shiny::div(class = "scnbar", style = sprintf("background:%s;", bg), desc),
                shiny::div(class = "pair-subdesc", pair_short_desc(row)),
                if (!is.null(row$note) && !is.na(row$note) && nzchar(row$note))
                  shiny::div(class = "note", row$note),
                shiny::fluidRow(
                  shiny::column(6, asv_panel(row$asv_i, "ps_ref ASV (fixed -- always kept)", "ref")),
                  shiny::column(6, asv_panel(row$asv_j, "ps_query ASV (the one in question)")))),
              shiny::div(class = "pair-actions",
                shiny::actionButton("rv_merge", "Merge -- same as ps_ref's ASV",
                                    class = "btn-success btn-choice"),
                shiny::actionButton("rv_distinct", "Keep distinct -- not the same",
                                    class = "btn-default btn-choice"),
                shiny::textInput("rv_rat", "Rationale (optional):", value = "", width = "100%"),
                if (can_back) shiny::actionButton("rv_back", "<- Back", class = "btn-default btn-sm"),
                shiny::actionButton("rv_skip", "Skip", class = "btn-link btn-sm"),
                shiny::actionButton("rv_exit", "Exit (keep progress so far)",
                                    class = "btn-link btn-sm", style = "margin-left:6px;")))
          )
        }
      )
    })

    shiny::observeEvent(input$summary_start, { rv$phase <- "review" })

    open_table_modal <- function() {
      shiny::showModal(shiny::modalDialog(
        title = "Decisions table", size = "l",
        shiny::checkboxGroupInput("tbl_cols", "Columns to display:",
                                  choices = all_cols, selected = default_cols, inline = TRUE),
        shiny::tags$hr(class = "thin"),
        DT::DTOutput("decisions_dt"),
        shiny::tags$hr(class = "thin"),
        shiny::fluidRow(
          shiny::column(8, shiny::textInput("csv_path", NULL, value = "",
                                            placeholder = "path/to/decisions.csv", width = "100%")),
          shiny::column(4, shiny::actionButton("export_csv", "Export CSV",
                                               class = "btn-default btn-sm", style = "margin-top:1px;"))),
        shiny::p(class = "text-muted small", style = "margin-top:4px;",
                 "Reflects progress made so far in this session."),
        footer = shiny::modalButton("Close"), easyClose = TRUE
      ))
    }
    shiny::observeEvent(input$summary_view_table, { open_table_modal() })
    shiny::observeEvent(input$show_table,         { open_table_modal() })

    output$decisions_dt <- DT::renderDT({
      cols <- if (!is.null(input$tbl_cols) && length(input$tbl_cols) > 0)
                intersect(input$tbl_cols, colnames(rv$rows))
              else default_cols
      df <- rv$rows[, cols, drop = FALSE]
      for (col in intersect(c("asv_i", "asv_j"), names(df)))
        df[[col]] <- paste0(substr(df[[col]], 1, 25), "...")
      DT::datatable(df, filter = "top", rownames = FALSE, class = "nowrap",
                    options = list(scrollX = TRUE, pageLength = 15, autoWidth = TRUE))
    }, server = TRUE)

    shiny::observeEvent(input$export_csv, {
      path <- if (!is.null(input$csv_path)) trimws(input$csv_path) else ""
      if (nzchar(path)) {
        utils::write.csv(rv$rows, path, row.names = FALSE, na = "")
        shiny::showNotification(paste("Written to", path), type = "message")
      }
    })

    record_and_advance <- function(dec) {
      k <- rv$k
      rv$snaps[[k]] <- rv$rows
      rat      <- if (!is.null(input$rv_rat)) trimws(input$rv_rat) else ""
      rows     <- rv$rows
      row_idx  <- pending_idx[k]
      if (!is.na(dec)) rows$decision[row_idx] <- dec
      if (nzchar(rat)) rows$rationale[row_idx] <- rat
      rv$rows <- rows
      rv$k    <- k + 1L
      shiny::updateTextInput(session, "rv_rat", value = "")
    }
    shiny::observeEvent(input$rv_merge,    { record_and_advance("merge_keep_i") })
    shiny::observeEvent(input$rv_distinct, { record_and_advance("keep_distinct") })
    shiny::observeEvent(input$rv_skip,     { record_and_advance(NA_character_) })
    shiny::observeEvent(input$rv_back, {
      k    <- rv$k
      snap <- tryCatch(rv$snaps[[k - 1]], error = function(e) NULL)
      if (is.null(snap)) return()
      rv$rows        <- snap
      rv$snaps[[k - 1]] <- NULL
      rv$k           <- k - 1L
    })

    # Unconditional save-to-disk on the way out -- see save_path's docs for
    # why: a blocking gadget's R-level return value can be lost in some
    # execution contexts even when the review inside it completed
    # correctly, and this is the fallback that keeps that work recoverable
    # regardless of what happens to the return value itself.
    finish_and_save <- function() {
      tryCatch({
        utils::write.csv(rv$rows, save_path, row.names = FALSE, na = "")
        message("Reviewed decisions also written to '", save_path, "' -- ",
                "use this path (decisions = \"", save_path, "\") if the ",
                "return value of this call doesn't come through as expected.")
      }, error = function(e) {
        warning("resolve_projection_interactive: could not write save_path ",
                "'", save_path, "': ", conditionMessage(e))
      })
      shiny::stopApp(invisible(rv$rows))
    }
    shiny::observeEvent(input$finish,  { finish_and_save() })
    shiny::observeEvent(input$rv_exit, { finish_and_save() })
  }

  # Deliberately shiny::paneViewer(), NOT shiny::dialogViewer() -- the
  # latter has a confirmed, unresolved RStudio bug where it corrupts the
  # calling function's own execution context on close, surfacing as the
  # NEXT line(s) of code (whatever they happen to be -- unrelated to this
  # gadget or its actual logic) throwing a bizarre, effectively-random
  # "object '...' not found" error. See
  # https://github.com/rstudio/rstudio/issues/11714, which explicitly
  # confirms paneViewer/browserViewer are NOT affected, only dialogViewer.
  # This gadget's own UI already renders its own header/title internally,
  # so nothing is lost by paneViewer() not taking a window title.
  shiny::runGadget(ui, server, viewer = shiny::paneViewer(minHeight = 500))
}

#' @title Plan a Query-to-Reference ASV Projection
#'
#' @description Detects correspondences between \code{ps_query}'s ASVs and
#'   a reference ASV space, then runs interactive review (a Shiny gadget,
#'   RStudio Viewer pane) for anything not resolved automatically or by
#'   \code{prior_decisions} -- the \code{plan_harmonization()} analog for
#'   projection, including that same automatic-review step (there is no
#'   separate function to call for it, mirroring how
#'   \code{plan_harmonization()} calls its own interactive review
#'   internally rather than exposing it as its own step). Paired with
#'   \code{\link{apply_projection}}, which actually builds the projected
#'   phyloseq object once every pair is resolved.
#'
#'   The reference is never a live phyloseq object -- it's
#'   \code{ps_ref_outcomes}, the cumulative ASV ledger
#'   \code{\link{apply_harmonization}} produces (its \code{$asv_outcomes},
#'   enriched with full taxonomy and, via that function's
#'   \code{prior_outcomes} argument, cumulative across as many rounds of
#'   harmonization as have happened so far). The reference is never
#'   modified, renamed, or pruned by anything here.
#'
#'   Detection reuses \code{compare_asvs()} (exact-sequence matches, and
#'   substring/species-set relationships) against the reference's
#'   reconstructed, current final ASVs. Algorithmically confident
#'   correspondences (scenario S1: substring + identical species set) are
#'   resolved automatically, with the reference's ASV always kept
#'   regardless of which sequence happens to be shorter -- a deliberate
#'   departure from \code{\link{apply_harmonization}}'s ordinary "keep the
#'   shorter" rule, since the point here is a fixed coordinate system, not
#'   the more parsimonious representative. Every other relationship is left
#'   for human review via the returned \code{$decisions} table.
#'
#' @param ps_ref_outcomes The cumulative ASV ledger describing the
#'   reference -- \code{\link{apply_harmonization}}'s return value, its
#'   \code{$asv_outcomes} data frame directly, or a file path to one saved
#'   as a CSV. Must have at least one \code{"representative"}/
#'   \code{"unchanged"} row and zero \code{"unresolved"} ones (see the
#'   unresolved gate below). Passing the full \code{apply_harmonization()}
#'   result (rather than a bare data frame or CSV) also enables the
#'   sample-name uniqueness check below.
#' @param ps_query New phyloseq object to project onto the reference's ASV
#'   space.
#' @param prior_decisions Optional decisions from an earlier
#'   \code{plan_projection()} round on this (or a closely related)
#'   reference/\code{ps_query} pair -- the list it returned, its
#'   \code{$decisions} data frame, or a path to one saved as a CSV. Folded
#'   into a fresh ambiguity check against this round's own detection, so
#'   already-decided pairs don't reappear and newly-introduced conflicts are
#'   still caught. If \code{NULL} (default), every pair is planned fresh.
#' @param source_names Length-2 character vector naming the reference and
#'   \code{ps_query} respectively (passed through to \code{compare_asvs()}).
#'   Default \code{c("x", "y")}.
#' @param min_overlap Numeric; S5 species-set overlap threshold, passed
#'   through to \code{compare_asvs()}. Default \code{0.01}.
#' @param tax_rank_cols Character vector of formal taxonomic rank column
#'   names present in \code{ps_ref_outcomes}; inferred (everything not one
#'   of the ledger's own fixed columns) when \code{NULL}.
#' @param save_path Optional file path to write \code{$decisions} to as a
#'   CSV. If \code{NULL} (default), nothing is written.
#' @param verbose If \code{TRUE} (default), report progress and a summary.
#'
#' @section The unresolved gate:
#' Before anything else runs, every row of \code{ps_ref_outcomes} is checked
#' -- if even one ASV is still \code{"unresolved"} (part of a flagged pair
#' from the reference's own harmonization that was never decided), the
#' whole call \code{stop()}s, naming which ASVs are pending. There is no
#' override.
#'
#' @section Sample-name uniqueness:
#' \code{ps_ref} and \code{ps_query} must never share a sample name --
#' checked with a hard \code{stop()} (no override) whenever a live
#' reference object is available (i.e. \code{ps_ref_outcomes} is the full
#' \code{apply_harmonization()} result, not a bare ledger). Skipped with a
#' \code{message()}, not silently, when only a bare \code{$asv_outcomes}
#' data frame or CSV was supplied, since there's no live reference object
#' to check names against at all.
#'
#' @section Reference/query drift checkpoint:
#' When \code{prior_decisions} carries the \code{ps_ref_taxa_hashes}/
#' \code{ps_query_taxa_hashes}/\code{ps_query_samples} columns
#' \code{plan_projection()} itself writes, the reference's current final
#' ASV set and \code{ps_query} are checked against them before detection
#' runs. Asymmetric: a mismatch on the reference's side \code{stop()}s
#' outright (a drifted reference risks silently mapping \code{ps_query}'s
#' reads onto the wrong ASVs); a \code{ps_query}-taxa or \code{ps_query}
#' -sample mismatch only ever \code{warning()}s and proceeds.
#'
#' @section Filtering the reference -- before harmonizing, not after:
#' If the reference needs taxonomy-based filtering (removing unassigned
#' ASVs, known control/spike-in sequences, host reads, etc.), apply it
#' \strong{before} running \code{\link{plan_harmonization}}/
#' \code{\link{apply_harmonization}} on it -- never after, and never by
#' filtering \code{ps_query} or \code{$ps_projected} independently
#' afterward. Harmonization physically merges some ASVs away, so filtering
#' the reference only after harmonizing filters a different ASV universe
#' than the one \code{ps_ref_outcomes} was actually built from. Filter
#' first, and \code{$ps_projected} will already reflect that filtering.
#'
#' @return A named list (class \code{"projection_plan"}), meant to be passed
#'   as-is to \code{\link{apply_projection}} (and safe to \code{saveRDS()}
#'   for later reuse, e.g. if you don't want to redo the review process).
#'   \code{$decisions} holds the pairs needing review (\code{asv_i} is
#'   always the reference-side candidate, \code{asv_j} is always the
#'   \code{ps_query} ASV in question, \code{source_i}/\code{source_j}
#'   always read \code{"ps_ref"}/\code{"ps_query"} regardless of
#'   \code{source_names}, plus \code{decision}/\code{rationale} columns --
#'   already filled in for anything resolved automatically, by
#'   \code{prior_decisions}, or by the interactive review that just ran,
#'   and still blank for anything skipped or left over when the gadget's
#'   \strong{Exit} was clicked early -- and the checkpoint columns
#'   described above). Any row still blank can be edited by hand directly
#'   on \code{$decisions} (or via a \code{save_path}/CSV round trip)
#'   before calling \code{apply_projection()}; alternatively, re-run
#'   \code{plan_projection()} with \code{prior_decisions} set to this same
#'   result to pick the interactive review back up where it left off
#'   (only the still-blank rows are shown again). \code{NULL} when nothing
#'   ever needed review (every pair auto-resolved). \code{$compare_result}
#'   (the full \code{compare_asvs()} output) is always included.
#'
#' @examples
#' \dontrun{
#' # Interactive review (if anything needs it) happens automatically here.
#' plan <- plan_projection(ref_harmonized, ps_query)
#' result <- apply_projection(plan, ps_query)
#' result$ps_projected
#' }
#'
#' @export
plan_projection <- function(ps_ref_outcomes, ps_query, prior_decisions = NULL,
                            source_names = c("x", "y"), min_overlap = 0.01,
                            tax_rank_cols = NULL, save_path = NULL, verbose = TRUE) {

  stopifnot(inherits(ps_query, "phyloseq"))
  stopifnot(length(source_names) == 2)
  src_x <- source_names[[1]]; src_y <- source_names[[2]]

  # -- Normalize ps_ref_outcomes (same 3-input-type pattern used everywhere
  # else in this package) ----------------------------------------------------
  ledger_df <- if (is.list(ps_ref_outcomes) && !is.data.frame(ps_ref_outcomes) &&
                   "asv_outcomes" %in% names(ps_ref_outcomes)) {
    ps_ref_outcomes$asv_outcomes
  } else if (is.data.frame(ps_ref_outcomes)) {
    ps_ref_outcomes
  } else if (is.character(ps_ref_outcomes) && length(ps_ref_outcomes) == 1) {
    utils::read.csv(ps_ref_outcomes, stringsAsFactors = FALSE, na.strings = "")
  } else {
    stop("plan_projection: `ps_ref_outcomes` must be an apply_harmonization() result, its ",
         "$asv_outcomes data.frame, or a file path to one saved as a CSV.")
  }
  stopifnot(is.data.frame(ledger_df))

  if (is.null(tax_rank_cols)) {
    tax_rank_cols <- setdiff(colnames(ledger_df), .ledger_fixed_cols)
  }
  .validate_asv_ledger(ledger_df, tax_rank_cols)

  # -- The unresolved gate: refuse to run at all if ps_ref isn't settled -------
  pending <- ledger_df$asv[!is.na(ledger_df$asv_status) & ledger_df$asv_status == "unresolved"]
  if (length(pending) > 0) {
    stop("plan_projection: ps_ref_outcomes has ", length(pending), " ASV(s) still unresolved from ",
         "the reference's own harmonization -- finish deciding those first (there is no ",
         "override): ", paste(utils::head(pending, 5), collapse = ", "),
         if (length(pending) > 5) ", ..." else "")
  }

  .check_sample_name_uniqueness(ps_ref_outcomes, ps_query, verbose)

  ref     <- .ledger_to_reference(ledger_df, tax_rank_cols)
  ps_ref_asvs <- ref$ps_ref_asvs

  # -- Normalize `prior_decisions` up front (same 3-input-type pattern used
  # elsewhere) -- needed before the checkpoint below. ------------------------
  prior_df <- NULL
  if (!is.null(prior_decisions)) {
    prior_df <- if (inherits(prior_decisions, c("projection_plan", "projection_result"))) {
      prior_decisions$decisions
    } else if (is.data.frame(prior_decisions)) {
      prior_decisions
    } else if (is.character(prior_decisions) && length(prior_decisions) == 1) {
      utils::read.csv(prior_decisions, stringsAsFactors = FALSE, na.strings = "")
    } else {
      stop("plan_projection: `prior_decisions` must be the list returned by plan_projection() or ",
           "apply_projection(), a decisions data.frame, or a file path to a decisions CSV.")
    }
    if (is.null(prior_df)) {
      stop("plan_projection: the supplied `prior_decisions` has no $decisions element to use -- ",
           "was everything already resolved (nothing was pending)?")
    }
    stopifnot(is.data.frame(prior_df))
    prior_df <- as.data.frame(
      lapply(prior_df, function(x) {
        if (is.character(x)) return(as.character(x))
        if (is.integer(x))   return(as.integer(x))
        if (is.numeric(x))   return(as.numeric(x))
        if (is.logical(x))   return(as.logical(x))
        if (is.factor(x))    return(as.factor(x))
        unclass(x)
      }),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    rownames(prior_df) <- NULL
  }

  .project_query_checkpoint(prior_df, ps_ref_asvs, ps_query)

  # -- Exact-match ghost redirect, before any detection at all --------------
  # A query ASV whose sequence is IDENTICAL to a known historical discard is
  # resolved immediately and definitively -- excluded from compare_asvs()
  # entirely, since there's nothing to detect: its fate is already known.
  ps_query_asvs_all <- phyloseq::taxa_names(ps_query)
  ghost_hits   <- intersect(ps_query_asvs_all, names(ref$ghost_map))
  ghost_target <- ref$ghost_map[ghost_hits]   # already fully resolved by apply_harmonization()

  ps_query_for_detect <- if (length(ghost_hits) > 0) {
    phyloseq::prune_taxa(setdiff(ps_query_asvs_all, ghost_hits), ps_query)
  } else ps_query

  cr <- if (length(setdiff(ps_query_asvs_all, ghost_hits)) > 0) {
    compare_asvs(list(ref$shell, ps_query_for_detect), source_names = c(src_x, src_y),
                min_overlap = min_overlap, tax_rank_cols = tax_rank_cols, verbose = verbose)
  } else {
    # Every query ASV was a ghost hit -- nothing left to run detection on.
    compare_asvs(list(ref$shell), source_names = src_x, min_overlap = min_overlap,
                tax_rank_cols = tax_rank_cols, verbose = FALSE)
  }
  flagged <- cr$flagged
  exact   <- intersect(ps_query_asvs_all, ps_ref_asvs)

  # -- Narrow $flagged to rows that are a genuine ps_ref<->ps_query CROSS pair --
  # A row only bears on projection when exactly one side is a genuine,
  # standalone query ASV (touches y, and NOT x -- "pure_y") and the other
  # side is a valid ps_ref candidate (touches x at all -- "ref_side").
  # ref_side deliberately does NOT require "not y" too: an ASV can
  # legitimately be BOTH a ps_ref final representative AND an exact-sequence
  # match already present in ps_query (dedup pools identical sequences
  # across inputs into one row with source "x,y" combined) -- that dual
  # role doesn't disqualify its OWN separate species-set relationship to
  # some genuinely different query ASV; it's still unambiguously ps_ref's
  # identity in that relationship.
  if (nrow(flagged) > 0) {
    i_touches_x  <- grepl(src_x, flagged$source_i, fixed = TRUE)
    i_touches_y  <- grepl(src_y, flagged$source_i, fixed = TRUE)
    j_touches_x  <- grepl(src_x, flagged$source_j, fixed = TRUE)
    j_touches_y  <- grepl(src_y, flagged$source_j, fixed = TRUE)
    i_is_ref_side <- i_touches_x
    i_is_pure_y   <- i_touches_y & !i_touches_x
    j_is_ref_side <- j_touches_x
    j_is_pure_y   <- j_touches_y & !j_touches_x
    relevant      <- (i_is_ref_side & j_is_pure_y) | (i_is_pure_y & j_is_ref_side)

    if (any(!relevant) && verbose) {
      message(sum(!relevant), " flagged pair(s) do not involve a genuine ps_ref<->ps_query ",
              "relationship (within-ps_ref or within-ps_query redundancy) and are out of scope ",
              "for projection.")
    }

    cross <- flagged[relevant, , drop = FALSE]
    if (nrow(cross) > 0) {
      # Swap whenever asv_i is currently the pure-query side, so asv_i always ends up ps_ref.
      cross <- .swap_ij(cross, i_is_pure_y[relevant])
    }
  } else {
    cross <- flagged[FALSE, , drop = FALSE]
  }

  is_s1    <- !is.na(cross$action) & cross$action == "auto_merge_keep_shorter"
  auto_map <- stats::setNames(cross$asv_i[is_s1], cross$asv_j[is_s1])  # query_asv -> ps_ref_asv
  needs_review <- cross[!is_s1, , drop = FALSE]

  res <- .resolve_projection_pairs(auto_map, needs_review, prior_df, ps_ref_asvs)

  if (length(res$ambiguous_y) > 0 && verbose) {
    message(length(res$ambiguous_y), " query ASV(s) have CONFLICTING merge targets among ps_ref's ",
            "ASVs -- forcing these back to review rather than guessing: ",
            paste(res$ambiguous_y, collapse = ", "))
  }

  decisions_out <- .build_projection_decisions_out(
    needs_review, prior_df, res$covered, res$ambiguous_y, is_s1, cross,
    ps_ref_asvs, ps_query_asvs_all, ps_query, src_x, src_y, verbose
  )

  # Step 2: interactive review (Shiny gadget) -- same relationship
  # plan_harmonization() has to resolve_conflicts_interactive(): called
  # automatically here, not a separate step the caller invokes by hand.
  # Only touches rows still blank after the prior_decisions carry-forward
  # above (see resolve_projection_interactive()'s own pending_idx logic);
  # falls back to leaving decisions_out unchanged (with a warning/message)
  # when the session isn't interactive, shiny/miniUI/DT aren't installed,
  # or there's nothing left to review.
  if (!is.null(decisions_out) && any(is.na(decisions_out$decision))) {
    decisions_out <- resolve_projection_interactive(decisions_out)
  }

  if (!is.null(decisions_out) && !is.null(save_path)) {
    utils::write.csv(decisions_out, save_path, row.names = FALSE, na = "")
    if (verbose) message("Decisions written to '", save_path, "'.")
  }

  hash_asvs <- function(asvs) vapply(asvs, function(x)
    digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))

  structure(
    list(
      decisions          = decisions_out,
      ps_ref_asvs        = ps_ref_asvs,
      taxtab             = ref$taxtab,
      ghost_target       = ghost_target,
      exact              = exact,
      auto_map           = auto_map,
      needs_review       = needs_review,
      is_s1              = is_s1,
      cross              = cross,
      ps_query_asvs_all  = ps_query_asvs_all,
      compare_result     = cr,
      src_x              = src_x,
      src_y              = src_y,
      ps_ref_taxa_hashes   = unname(hash_asvs(ps_ref_asvs)),
      ps_query_taxa_hashes = unname(hash_asvs(ps_query_asvs_all)),
      ps_query_samples     = phyloseq::sample_names(ps_query)
    ),
    class = "projection_plan"
  )
}
