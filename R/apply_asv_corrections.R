# R/apply_asv_corrections.R
#
# Step 3 of ASV harmonization for FoodSeq (trnL / 12SV5) phyloseq data:
# applying user-reviewed merge decisions to a phyloseq object. Steps 1-2
# (pairwise comparison and interactive review) live in the companion file
# plan_harmonization.R; this file's apply_asv_corrections() has no
# dependency on them beyond the data they produce (via a harmonization_plan
# object, or the equivalent raw arguments below) and needs only phyloseq.
#
# Usage
# -----
#   source("plan_harmonization.R")
#   source("apply_asv_corrections.R")
#
#   plan   <- plan_harmonization(ps_list, source_names = names(ps_list))
#   result <- apply_asv_corrections(plan)
#   # result$ps_harmonized, result$decisions, result$asv_summary

#' @title Apply ASV Correction Decisions
#'
#' @description Adds a `glom_name` column to the tax table, applies S2 rank
#'   alignment for user-merged pairs, physically merges ASVs where a
#'   sequence preference was expressed, and applies all manual name updates.
#'
#' @param ps A phyloseq object, or a \code{harmonization_plan} object returned
#'   by \code{plan_harmonization()} (in the companion file
#'   plan_harmonization.R). When a plan is passed, \code{flagged},
#'   \code{tax_rank_cols}, \code{manual_updates}, \code{common_name_updates},
#'   \code{keep_asvs}, \code{s1_skip_keys}, \code{rationale}, and
#'   \code{asv_info} are all taken from the plan and any values passed
#'   explicitly for those arguments are ignored.
#' @param flagged Data frame of flagged pairs, e.g. \code{compare_result$flagged}.
#' @param tax_rank_cols Rank column names; inferred from phyloseq if NULL.
#' @param manual_updates Named character vector c("asv_seq" = "chosen_name")
#'   from \code{resolve_conflicts_interactive()$name_updates}.
#' @param common_name_updates Named character vector (group taxa key ->
#'   common name) from \code{resolve_conflicts_interactive()$common_name_updates}.
#' @param keep_asvs Named character vector c("discard_asv" = "keep_asv") from
#'   \code{resolve_conflicts_interactive()$keep_asvs}. For each pair, the
#'   discard ASV's counts are summed into the keep ASV's counts and the
#'   discard ASV is then removed from the phyloseq object.
#' @param rationale Named character vector c("asv_i|asv_j" = "rationale text")
#'   from \code{resolve_conflicts_interactive()$rationale}. Joined into
#'   \code{$decisions}.
#' @param s1_skip_keys Character vector of "asv_i|asv_j" pair keys from
#'   \code{resolve_conflicts_interactive()$s1_skip_keys}. S1 pairs in this set
#'   are skipped during auto-merge (user chose to keep longer or keep distinct).
#' @param asv_info Per-ASV metadata data frame, e.g. \code{compare_result$asv_info}.
#' @return When \code{ps} is a \code{harmonization_plan}, a list with
#'   \code{ps_harmonized}, \code{decisions}, \code{asv_summary},
#'   \code{unresolved}, and \code{compare_result}. Otherwise a list with:
#'   \itemize{
#'     \item \code{ps} — phyloseq with glom_name in tax table, discarded ASVs pruned
#'     \item \code{decisions} — one row per reviewed pair with scenario, decision, rationale
#'     \item \code{unresolved} — flagged pairs not covered by manual_updates
#'     \item \code{taxtab}, \code{n_s2_corrected}, \code{n_manual}
#'   }
#'   To resolve unresolved pairs manually:
#'   \preformatted{
#'     corrected$taxtab["<asv_seq>", "glom_name"] <- "<name>"
#'     tax_table(corrected$ps) <- as.matrix(corrected$taxtab)
#'   }
#' @export
apply_asv_corrections <- function(ps, flagged = NULL, tax_rank_cols = NULL,
                                  manual_updates      = NULL,
                                  common_name_updates = NULL,
                                  keep_asvs           = NULL,
                                  rationale           = NULL,
                                  s1_skip_keys        = NULL,
                                  asv_info            = NULL) {
  is_plan <- inherits(ps, "harmonization_plan")
  if (is_plan) {
    plan                <- ps
    ps                  <- plan$ps_merged
    flagged             <- plan$flagged
    tax_rank_cols       <- plan$tax_rank_cols
    manual_updates      <- plan$manual$name_updates
    common_name_updates <- plan$manual$common_name_updates
    keep_asvs           <- plan$manual$keep_asvs
    s1_skip_keys        <- plan$manual$s1_skip_keys
    rationale           <- plan$manual$rationale
    asv_info            <- plan$asv_info
  }

  stopifnot(requireNamespace("phyloseq", quietly = TRUE))
  stopifnot(inherits(ps, "phyloseq"))
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

  all_tax_cols    <- colnames(phyloseq::tax_table(ps))
  taxon_set_col   <- if ("taxa"        %in% all_tax_cols) "taxa"
                     else tail(setdiff(all_tax_cols, "glom_name"), 1)
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

  # ── Add `glom_name` column ──────────────────────────────────────────────────────────────────────────
  # Default: each ASV gets its own sequence as a unique glom_name so that
  # no two ASVs can be silently merged unless they were explicitly decided on
  # (S1 auto-merge or user merge decision). Only those ASVs will share a name.
  taxtab$glom_name <- rownames(taxtab)

  # ── Apply S2 rank alignment when a manual merge decision was made ──────
  # When the user merges an S2 pair, the higher-res ASV's rank columns are
  # overwritten with the lower-res ASV's values so the surviving representative
  # row is internally consistent regardless of which ASV has higher abundance.
  s2_rows     <- flagged[!is.na(flagged$scenario) & flagged$scenario == 2L, ]
  n_corrected <- 0L

  if (nrow(s2_rows) > 0 && !is.null(manual_updates) && length(manual_updates) > 0) {
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

  # ── Apply manual updates ───────────────────────────────────────────────────
  n_manual <- 0L
  if (!is.null(manual_updates) && length(manual_updates) > 0) {
    for (asv_seq in names(manual_updates)) {
      if (!asv_seq %in% rownames(taxtab)) {
        warning("manual_updates: ASV not found in tax table: ", asv_seq); next
      }
      taxtab[asv_seq, "glom_name"] <- manual_updates[[asv_seq]]
      n_manual <- n_manual + 1L
    }
  }

  methods::slot(ps, "tax_table", check = FALSE) <- phyloseq::tax_table(as.matrix(taxtab))

  # ── Transfer counts and prune discarded ASVs ────────────────────────────
  # Build a combined discard->keep map from two sources:
  #   (a) S1 auto-merges: the longer ASV is discarded; keep_asv is the shorter.
  #   (b) User decisions from resolve_conflicts_interactive() via keep_asvs.
  # All transfers and pruning happen in one pass.
  combined_keep      <- character(0)  # named: discard_asv -> keep_asv
  combined_keep_type <- character(0)  # named: discard_asv -> decision_type

  s1_rows <- flagged[!is.na(flagged$action) &
                       flagged$action == "auto_merge_keep_shorter", ]
  if (nrow(s1_rows) > 0 && "keep_asv" %in% colnames(s1_rows)) {
    for (k in seq_len(nrow(s1_rows))) {
      row_k      <- s1_rows[k, ]
      pair_key_k <- paste0(row_k$asv_i, "|", row_k$asv_j)
      # Skip pairs where the user overrode the default or chose to keep distinct
      if (!is.null(s1_skip_keys) && pair_key_k %in% s1_skip_keys) next
      keep_k    <- row_k$keep_asv
      if (is.na(keep_k) || !nzchar(keep_k)) next
      discard_k <- if (keep_k == row_k$asv_i) row_k$asv_j else row_k$asv_i
      combined_keep[discard_k]      <- keep_k
      combined_keep_type[discard_k] <- "s1_auto"
    }
  }

  if (!is.null(keep_asvs) && length(keep_asvs) > 0) {
    for (d in names(keep_asvs)) {
      k_asv <- keep_asvs[[d]]
      combined_keep[d] <- k_asv
      is_s2 <- nrow(s2_rows) > 0 &&
                any((s2_rows$asv_i == d & s2_rows$asv_j == k_asv) |
                    (s2_rows$asv_i == k_asv & s2_rows$asv_j == d))
      combined_keep_type[d] <- if (is_s2) "s2_rank_align" else "user_merge"
    }
  }

  taxtab_pre_prune <- taxtab  # capture before any pruning — used for discarded ASV metadata
  n_pruned <- 0L
  if (length(combined_keep) > 0) {
    all_taxa  <- phyloseq::taxa_names(ps)
    rows_mode <- phyloseq::taxa_are_rows(ps)
    otu_mat   <- as(phyloseq::otu_table(ps), "matrix")

    for (discard in names(combined_keep)) {
      keep <- combined_keep[[discard]]
      if (!discard %in% all_taxa) {
        warning("discard ASV not found in phyloseq: ", discard); next
      }
      if (!keep %in% all_taxa) {
        warning("keep ASV not found in phyloseq: ", keep); next
      }
      if (rows_mode) {
        otu_mat[keep, ] <- otu_mat[keep, ] + otu_mat[discard, ]
      } else {
        otu_mat[, keep] <- otu_mat[, keep] + otu_mat[, discard]
      }
    }

    valid_discard <- names(combined_keep)[names(combined_keep) %in% all_taxa]
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

  n_cn_updated <- 0L
  if (!is.null(common_name_updates) && length(common_name_updates) > 0 &&
      !is.null(manual_updates)      && length(manual_updates) > 0) {
    for (taxa_key in names(common_name_updates)) {
      group_asvs <- names(manual_updates)[!is.na(manual_updates) &
                                           nzchar(manual_updates) &
                                           manual_updates == taxa_key]
      survivors  <- intersect(group_asvs, rownames(taxtab))
      if (length(survivors) == 0) {
        warning("common_name_updates: no surviving ASV found for group '",
                taxa_key, "'")
        next
      }
      for (asv_s in survivors) {
        if (taxon_set_col %in% colnames(taxtab))
          taxtab[asv_s, taxon_set_col] <- taxa_key
        if (!is.null(common_name_col) && common_name_col %in% colnames(taxtab))
          taxtab[asv_s, common_name_col] <- common_name_updates[[taxa_key]]
      }
      n_cn_updated <- n_cn_updated + 1L
    }
    methods::slot(ps, "tax_table", check = FALSE) <- phyloseq::tax_table(as.matrix(taxtab))
  }
  if (n_cn_updated > 0)
    message(n_cn_updated, " merge group(s) updated with merged taxa string and custom common name.")

  # ── Unresolved: flagged pairs not covered by manual_updates ──────────────────
  review_actions <- c("flag_lower_resolution", "flag_same_taxon",
                      "flag_split_distinct_names", "flag_taxon_subset",
                      "flag_taxon_overlap", "flag_conflict", "flag_unassigned")
  needs_manual <- flagged$action %in% review_actions
  unresolved   <- flagged[needs_manual, ]

  if (nrow(unresolved) > 0 && length(manual_updates) > 0) {
    resolved_asvs <- names(manual_updates)
    both_covered  <- unresolved$asv_i %in% resolved_asvs &
                     unresolved$asv_j %in% resolved_asvs
    unresolved    <- unresolved[!both_covered, ]
  }

  n_auto <- sum(flagged$action == "auto_merge_keep_shorter", na.rm = TRUE)
  if (n_corrected > 0)
    message(n_corrected, " S2 rank alignment(s) applied (user-merged pairs).")
  if (n_auto > 0)
    message(n_auto, " S1 pair(s) auto-merged: counts transferred, longer ASV pruned.")
  if (n_pruned > 0)
    message(n_pruned, " ASV(s) physically removed (counts transferred to representative).")
  if (n_manual > 0)
    message(n_manual / 2L, " pair(s) resolved manually.")
  if (nrow(unresolved) > 0)
    message(nrow(unresolved), " pair(s) still unresolved (see $unresolved).")
  message("Corrections complete. corrected$ps is ready for downstream analysis.")

  # ── Remove glom_name from tax table ──────────────────────────────────────
  if ("glom_name" %in% colnames(taxtab)) {
    taxtab$glom_name <- NULL
    methods::slot(ps, "tax_table", check = FALSE) <- phyloseq::tax_table(as.matrix(taxtab))
  }

  # ── Build decisions table ─────────────────────────────────────────────────
  discarded_asvs   <- names(combined_keep)
  reps_with_merges <- unique(unname(combined_keep))

  # Rep -> semicolon-delimited list of ASVs merged into it
  if (length(combined_keep) > 0) {
    rep_to_merged <- tapply(names(combined_keep), unname(combined_keep),
                            function(x) paste(sort(x), collapse = "; "))
  } else {
    rep_to_merged <- character(0)
  }

  # Rep -> highest-priority decision_type among all contributing merges
  priority_map      <- c(s2_rank_align = 3L, user_merge = 2L, s1_auto = 1L)
  rep_decision_type <- character(0)
  for (d in names(combined_keep)) {
    k  <- combined_keep[d]
    dt <- combined_keep_type[d]
    if (!k %in% names(rep_decision_type) ||
        priority_map[dt] > priority_map[rep_decision_type[k]])
      rep_decision_type[k] <- dt
  }

  # asv_info row index for O(1) lookup (no closures — avoids JIT issues)
  ai_idx <- if (!is.null(asv_info) && "asv" %in% colnames(asv_info))
               setNames(seq_len(nrow(asv_info)), asv_info$asv)
             else NULL

  # Pre-allocate decision columns as plain vectors
  n_dec             <- length(all_original_asvs)
  dec_asv_status    <- character(n_dec)
  dec_merged_with   <- rep(NA_character_, n_dec)
  dec_merged_into   <- rep(NA_character_, n_dec)
  dec_dtype         <- character(n_dec)
  dec_scenario      <- rep(NA_character_, n_dec)
  dec_common_name   <- rep(NA_character_, n_dec)
  dec_taxa          <- rep(NA_character_, n_dec)
  dec_lowest_rank   <- rep(NA_character_, n_dec)
  dec_lowest_val    <- rep(NA_character_, n_dec)
  dec_batch         <- rep(NA_character_, n_dec)
  dec_reads         <- rep(NA_real_,      n_dec)
  dec_rationale     <- rep(NA_character_, n_dec)

  has_action_col <- "action" %in% colnames(flagged)

  for (idx in seq_along(all_original_asvs)) {
    asv_seq <- all_original_asvs[[idx]]

    # ── Fate / merge columns ───────────────────────────────────────────────
    if (asv_seq %in% discarded_asvs) {
      dec_asv_status[idx]  <- "dropped after merge"
      dec_merged_into[idx] <- unname(combined_keep[[asv_seq]])
      dec_dtype[idx]       <- unname(combined_keep_type[[asv_seq]])
      rep_asv <- unname(combined_keep[[asv_seq]])
      fl <- flagged[((flagged$asv_i == asv_seq & flagged$asv_j == rep_asv) |
                     (flagged$asv_i == rep_asv  & flagged$asv_j == asv_seq)), ]
      if (nrow(fl) > 0) {
        sc <- fl$scenario[[1]]
        dec_scenario[idx] <- if (!is.na(sc)) paste0("S", sc)
                             else if (has_action_col && grepl("conflict",   fl$action[[1]], fixed = TRUE)) "conflict"
                             else if (has_action_col && grepl("unassigned", fl$action[[1]], fixed = TRUE)) "unassigned"
                             else NA_character_
        if (!is.null(rationale) && length(rationale) > 0) {
          k1 <- paste0(fl$asv_i[[1]], "|", fl$asv_j[[1]])
          k2 <- paste0(fl$asv_j[[1]], "|", fl$asv_i[[1]])
          rv <- rationale[k1]; if (is.na(rv)) rv <- rationale[k2]
          if (!is.na(rv)) dec_rationale[idx] <- unname(rv)
        }
      }

    } else if (asv_seq %in% reps_with_merges) {
      dec_asv_status[idx] <- "representative"
      dec_merged_with[idx] <- if (asv_seq %in% names(rep_to_merged))
                                unname(rep_to_merged[[asv_seq]])
                              else NA_character_
      dec_dtype[idx]     <- if (asv_seq %in% names(rep_decision_type))
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
          if (!is.null(rationale) && length(rationale) > 0) {
            k1 <- paste0(fl$asv_i[[1]], "|", fl$asv_j[[1]])
            k2 <- paste0(fl$asv_j[[1]], "|", fl$asv_i[[1]])
            rv <- rationale[k1]; if (is.na(rv)) rv <- rationale[k2]
            if (!is.na(rv)) rat_vals <- c(rat_vals, unname(rv))
          }
        }
      }
      if (length(sc_vals)  > 0) dec_scenario[idx]  <- paste(sc_vals,  collapse = ";")
      if (length(rat_vals) > 0) dec_rationale[idx] <- paste(rat_vals, collapse = " | ")

    } else {
      dec_asv_status[idx] <- "unchanged"
      dec_dtype[idx] <- "unchanged"
    }

    # ── Taxonomy columns ───────────────────────────────────────────────────
    tab <- if (asv_seq %in% discarded_asvs) taxtab_pre_prune else taxtab
    if (asv_seq %in% rownames(tab)) {
      if (!is.null(common_name_col) && common_name_col %in% colnames(tab)) {
        v <- tab[[common_name_col]][[which(rownames(tab) == asv_seq)[[1]]]]
        if (!is.na(v) && nzchar(v)) dec_common_name[idx] <- as.character(v)
      }
      if (taxon_set_col %in% colnames(tab)) {
        v <- tab[[taxon_set_col]][[which(rownames(tab) == asv_seq)[[1]]]]
        if (!is.na(v) && nzchar(v)) dec_taxa[idx] <- as.character(v)
      }
      if (!is.null(ai_idx) && asv_seq %in% names(ai_idx)) {
        ri <- ai_idx[[asv_seq]]
        lr <- asv_info[ri, "deepest_rank"]
        lv <- asv_info[ri, "deepest_name"]
        if (!is.na(lr)) dec_lowest_rank[idx] <- as.character(lr)
        if (!is.na(lv)) dec_lowest_val[idx]  <- as.character(lv)
      } else {
        row_vals <- unlist(tab[asv_seq, tax_rank_cols, drop = FALSE])
        row_vals <- as.character(row_vals)
        non_na   <- which(!is.na(row_vals) & nzchar(row_vals))
        if (length(non_na) > 0) {
          dec_lowest_rank[idx] <- tax_rank_cols[[max(non_na)]]
          dec_lowest_val[idx]  <- row_vals[[max(non_na)]]
        }
      }
    }

    # ── Batch source and reads ─────────────────────────────────────────────
    if (!is.null(ai_idx) && asv_seq %in% names(ai_idx)) {
      ri <- ai_idx[[asv_seq]]
      sv <- asv_info[ri, "sources"]
      if (!is.na(sv)) dec_batch[idx] <- as.character(sv)
    }
    if (asv_seq %in% names(reads_pre)) dec_reads[idx] <- reads_pre[[asv_seq]]
  }

  decisions <- data.frame(
    asv               = all_original_asvs,
    asv_status        = dec_asv_status,
    merged_with       = dec_merged_with,
    merged_into       = dec_merged_into,
    decision_type     = dec_dtype,
    scenario          = dec_scenario,
    common_name       = dec_common_name,
    taxa              = dec_taxa,
    lowest_rank       = dec_lowest_rank,
    lowest_rank_value = dec_lowest_val,
    batch_source      = dec_batch,
    total_reads       = dec_reads,
    rationale         = dec_rationale,
    stringsAsFactors  = FALSE
  )

  if (is_plan) {
    # User-facing return matching plan_harmonization() output structure
    asv_summary <- NULL
    if (!is.null(decisions) && "asv_status" %in% colnames(decisions)) {
      asv_summary <- decisions[decisions$asv_status != "unchanged", , drop = FALSE]
      rownames(asv_summary) <- NULL
    }
    invisible(list(
      ps_harmonized  = ps,
      decisions      = decisions,
      asv_summary    = asv_summary,
      unresolved     = unresolved,
      compare_result = plan$compare_result
    ))
  } else {
    list(
      ps             = ps,
      taxtab         = taxtab,
      unresolved     = unresolved,
      decisions      = decisions,
      n_s2_corrected = n_corrected,
      n_manual       = n_manual / 2L
    )
  }
}


