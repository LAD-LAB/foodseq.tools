# R/plan_harmonization.R
#
# Detects cross-batch (and within-object) ASV redundancy for FoodSeq
# (trnL / 12SV5) phyloseq data and applies user-reviewed merge decisions.
#
# Requires assign_common_names() output (assign_common_names_AS.R), which
# adds two columns to the tax table:
#   common_name -- resolved human-readable food label
#   taxa        -- sorted semicolon-delimited set of matched species
#
# Scenarios detected
# ------------------
#   1  -- substring seq + same species set -> auto-merge; keep shorter ASV
#   1b -- substring seq + same formal LCA + different species sets -> flag
#   2  -- substring seq + different resolution + shared lineage -> flag
#   3  -- non-substring seq + identical species sets -> flag
#   4  -- non-substring seq + species set A subset of B (or B subset of A) -> flag
#   5  -- non-substring seq + species-set overlap >= min_overlap -> flag
#   conflict   -- substring + conflicting formal taxon -> flag
#   unassigned -- substring + missing taxonomy -> flag
#
# Usage
# -----
#   plan   <- plan_harmonization(ps_list, source_names = names(ps_list))
#   result <- apply_harmonization(ps, plan)
#   # result$ps, result$asv_outcomes, result$unresolved, etc.
#
# Returns from compare_asvs()
# -----------------------------------
#   $flagged         -- data frame of all flagged pairs
#   $name_updates    -- named character vector c("asv_seq" = "new_name") for S2
#   $summary         -- frequency table per scenario
#   $asv_info        -- per-ASV metadata including species sets and common names
#   $taxon_set_col   -- name of the column used as species-set key
#   $common_name_col -- name of the common-name column (or NULL)
#   $label_col       -- alias for taxon_set_col (backward compat)

compare_asvs <- function(..., source_names = NULL, tax_rank_cols = NULL,
                               min_overlap = 0.01, flag_common_name = FALSE,
                               verbose = TRUE) {

  # -- 0. Internal helpers ----------------------------------------------------

  # Parse a semicolon-delimited species-set string into a sorted character vector.
  parse_species_set <- function(taxa_str) {
    if (is.na(taxa_str) || !nzchar(taxa_str) || taxa_str == "NA")
      return(character(0))
    sort(trimws(strsplit(taxa_str, ";")[[1]]))
  }

  # Strip infraspecific epithets (var., subsp., ssp., f., cv., agg.) from a
  # character vector of species names, returning base binomials.
  strip_infraspecific <- function(sp_vec) {
    gsub("\\s+(var\\.|subsp\\.|ssp\\.|f\\.|cv\\.|agg\\.).*$", "",
         trimws(sp_vec), perl = TRUE)
  }

  # Classify the relationship between two species sets.
  # Returns list(type, overlap_coef).
  #   type: "equal" | "a_subset_b" | "b_subset_a" | "overlap" | "disjoint" | "unassigned"
  # After exact-string checks, falls back to a normalised-infraspecific check so
  # that e.g. {"Brassica oleracea var. italica"} vs {"Brassica oleracea"} is
  # classified as "a_subset_b" rather than "disjoint".
  taxa_set_relation <- function(set_a, set_b) {
    if (length(set_a) == 0 || length(set_b) == 0)
      return(list(type = "unassigned", overlap_coef = NA_real_))
    inter   <- length(intersect(set_a, set_b))
    oc      <- inter / min(length(set_a), length(set_b))
    if (inter == 0) {
      # Exact strings share nothing -- try normalised infraspecific comparison.
      norm_a  <- unique(strip_infraspecific(set_a))
      norm_b  <- unique(strip_infraspecific(set_b))
      ni      <- length(intersect(norm_a, norm_b))
      if (ni == 0)
        return(list(type = "disjoint", overlap_coef = 0))
      n_oc <- ni / min(length(norm_a), length(norm_b))
      # Deliberately never return "equal" here, even when norm_a and norm_b
      # are setequal: the two sides differ in their exact (infraspecific)
      # names, so treating them as fully identical would let this pair
      # qualify for scenario 1's *unreviewed* auto-merge (see same_comparison
      # below). Falling into a_subset_b/b_subset_a/overlap instead keeps the
      # pair on a human-reviewed path (S4/S5).
      if (all(norm_a %in% norm_b))
        return(list(type = "a_subset_b", overlap_coef = n_oc))
      if (all(norm_b %in% norm_a))
        return(list(type = "b_subset_a", overlap_coef = n_oc))
      return(list(type = "overlap", overlap_coef = n_oc))
    }
    if (setequal(set_a, set_b))
      return(list(type = "equal",       overlap_coef = 1))
    if (all(set_a %in% set_b))
      return(list(type = "a_subset_b",  overlap_coef = oc))
    if (all(set_b %in% set_a))
      return(list(type = "b_subset_a",  overlap_coef = oc))
    list(type = "overlap", overlap_coef = oc)
  }

  # -- 1. Collect inputs ------------------------------------------------------
  stopifnot(requireNamespace("phyloseq", quietly = TRUE))

  inputs <- list(...)
  if (length(inputs) == 1 && is.list(inputs[[1]]) &&
      !inherits(inputs[[1]], "phyloseq")) {
    inputs <- inputs[[1]]
  }
  stopifnot(length(inputs) >= 1)
  for (i in seq_along(inputs)) stopifnot(inherits(inputs[[i]], "phyloseq"))

  if (is.null(source_names)) {
    source_names <- as.character(seq_along(inputs))
  } else {
    stopifnot(length(source_names) == length(inputs))
  }

  # -- 2. Identify columns ----------------------------------------------------
  # Prefer explicit assign_common_names_AS.R column names; fall back to
  # last-column convention for compatibility with other naming functions.
  all_tax_cols <- colnames(phyloseq::tax_table(inputs[[1]]))

  taxon_set_col   <- if ("taxa"        %in% all_tax_cols) "taxa"
                     else utils::tail(setdiff(all_tax_cols, "glom_name"), 1)
  common_name_col <- if ("common_name" %in% all_tax_cols) "common_name" else NULL
  label_col       <- taxon_set_col   # backward-compat alias

  if (is.null(tax_rank_cols)) {
    exclude       <- c(taxon_set_col, common_name_col, "glom_name")
    tax_rank_cols <- setdiff(all_tax_cols, exclude)
  }

  # -- 3. Pool unique ASVs across all objects ---------------------------------
  # taxon_set_col/common_name_col are derived from inputs[[1]] alone (above),
  # so a later input missing one of those columns entirely -- e.g. a raw new
  # batch that hasn't been through assign_common_names() yet, exactly the
  # case apply_projection() hands in -- must still end up with the SAME column
  # set as every other input's rows, or the final rbind() below breaks with
  # a bare "numbers of columns do not match". Pad any column absent from a
  # given input's own tax_table with NA rather than silently omitting it.
  extract_cols <- c(tax_rank_cols, taxon_set_col,
                    if (!is.null(common_name_col)) common_name_col)
  asv_rows <- list()
  for (i in seq_along(inputs)) {
    ps   <- inputs[[i]]
    otu  <- phyloseq::otu_table(ps)
    seqs <- if (phyloseq::taxa_are_rows(ps)) rownames(otu) else colnames(otu)
    tax  <- as.data.frame(phyloseq::tax_table(ps)@.Data, stringsAsFactors = FALSE)

    present_cols <- intersect(extract_cols, colnames(tax))
    tax  <- tax[seqs, present_cols, drop = FALSE]
    df   <- data.frame(asv = seqs, stringsAsFactors = FALSE)
    df[present_cols] <- tax
    for (mc in setdiff(extract_cols, present_cols)) df[[mc]] <- NA_character_
    df   <- df[, c("asv", extract_cols), drop = FALSE]
    df$source <- source_names[i]
    asv_rows[[i]] <- df
  }
  all_asvs <- do.call(rbind, asv_rows)

  # Deduplicate: prefer rows with non-all-NA formal taxonomy per unique ASV
  dedup <- function(group) {
    all_na   <- apply(group[, tax_rank_cols, drop = FALSE], 1,
                      function(r) all(is.na(r)))
    assigned <- group[!all_na, , drop = FALSE]
    if (nrow(assigned) > 0) assigned[1, ] else group[1, ]
  }
  asv_info <- do.call(rbind, lapply(split(all_asvs, all_asvs$asv), dedup))

  # sources: comma-separated batch names each ASV appeared in
  source_map <- tapply(all_asvs$source, all_asvs$asv,
                       function(s) paste(sort(unique(s)), collapse = ","))
  asv_info$sources <- source_map[asv_info$asv]
  rownames(asv_info) <- NULL

  # -- 4. Per-ASV deepest resolved rank, species set, and display name --------
  get_deepest <- function(tax_row) {
    vals     <- as.character(tax_row[tax_rank_cols])
    assigned <- which(!is.na(vals) & vals != "" & vals != "NA")
    if (length(assigned) == 0)
      return(list(level = NA_integer_, rank = NA_character_,
                  name  = NA_character_, lineage = NA_character_))
    deepest <- max(assigned)
    list(
      level   = deepest,
      rank    = tax_rank_cols[deepest],
      name    = vals[deepest],
      lineage = paste(vals[assigned], collapse = "|")
    )
  }

  deepest_list <- lapply(seq_len(nrow(asv_info)), function(i)
    get_deepest(asv_info[i, tax_rank_cols]))
  asv_info$deepest_level <- sapply(deepest_list, `[[`, "level")
  asv_info$deepest_rank  <- sapply(deepest_list, `[[`, "rank")
  asv_info$deepest_name  <- sapply(deepest_list, `[[`, "name")
  asv_info$lineage       <- sapply(deepest_list, `[[`, "lineage")
  asv_info$len           <- nchar(asv_info$asv)

  # -- Per-ASV read count stats across all inputs ----------------------------
  # Build a list of transposed OTU matrices (ASV x sample) for each input.
  reads_mats <- lapply(inputs, function(ps_i) {
    m <- as.matrix(phyloseq::otu_table(ps_i))
    if (!phyloseq::taxa_are_rows(ps_i)) t(m) else m   # ensure rows = ASVs
  })
  names(reads_mats) <- source_names

  # Total reads and sample counts per source object (denominators for % and prevalence).
  source_totals  <- vapply(reads_mats, function(m) sum(m),     numeric(1))
  source_n_samps <- vapply(reads_mats, function(m) ncol(m),    integer(1))

  reads_stats <- lapply(asv_info$asv, function(asv) {
    src_names <- trimws(strsplit(source_map[asv], ",")[[1]])
    src_mats  <- reads_mats[intersect(src_names, names(reads_mats))]

    # Total reads across all inputs.
    vals_all  <- unlist(lapply(reads_mats, function(m) {
      if (asv %in% rownames(m)) as.numeric(m[asv, ]) else numeric(0)
    }))
    total_all <- if (length(vals_all) == 0) 0L else sum(vals_all)
    max_all   <- if (length(vals_all) == 0) 0L else max(vals_all)

    # % of reads: denominator = total reads in source object(s).
    denom_reads <- sum(source_totals[names(src_mats)])
    pct         <- if (denom_reads > 0) round(100 * total_all / denom_reads, 2) else NA_real_

    # Prevalence: % of samples with >0 reads, within source object(s).
    n_detected  <- sum(unlist(lapply(src_mats, function(m) {
      if (asv %in% rownames(m)) sum(m[asv, ] > 0L) else 0L
    })))
    n_samps     <- sum(source_n_samps[names(src_mats)])
    prev_pct    <- if (n_samps > 0) round(100 * n_detected / n_samps, 1) else NA_real_

    list(reads_total = total_all, reads_max = max_all,
         reads_pct = pct, reads_prev = prev_pct)
  })
  asv_info$reads_total <- vapply(reads_stats, `[[`, numeric(1), "reads_total")
  asv_info$reads_max   <- vapply(reads_stats, `[[`, numeric(1), "reads_max")
  asv_info$reads_pct   <- vapply(reads_stats, `[[`, numeric(1), "reads_pct")
  asv_info$reads_prev  <- vapply(reads_stats, `[[`, numeric(1), "reads_prev")

  # Species-set string (taxa column from assign_common_names_AS.R)
  if (taxon_set_col %in% colnames(asv_info)) {
    taxon_strs <- as.character(asv_info[[taxon_set_col]])
    taxon_strs[is.na(taxon_strs) | taxon_strs == "" | taxon_strs == "NA"] <- NA_character_
  } else {
    taxon_strs <- rep(NA_character_, nrow(asv_info))
  }
  asv_info$taxon_str <- taxon_strs

  # Common name (common_name column from assign_common_names_AS.R)
  if (!is.null(common_name_col) && common_name_col %in% colnames(asv_info)) {
    cn_vals <- as.character(asv_info[[common_name_col]])
    cn_vals[is.na(cn_vals) | cn_vals == "" | cn_vals == "NA"] <- NA_character_
  } else {
    cn_vals <- rep(NA_character_, nrow(asv_info))
  }
  asv_info$common_name_val <- cn_vals

  # label_val: taxon_str (backward-compat alias)
  asv_info$label_val <- taxon_strs

  # comparison_name: human-readable display label (common_name > taxon_str > deepest_name)
  asv_info$comparison_name <- ifelse(!is.na(cn_vals), cn_vals,
                               ifelse(!is.na(taxon_strs), taxon_strs,
                                      asv_info$deepest_name))

  n_asvs <- nrow(asv_info)
  if (verbose) {
    message("Comparing ", n_asvs, " unique ASVs across ",
            length(inputs), " phyloseq object(s).")
    message("Taxon-set column : '", taxon_set_col, "'")
    if (!is.null(common_name_col))
      message("Common-name column: '", common_name_col, "'")
    message("min_overlap threshold (S5): ", min_overlap,
            " | flag_common_name: ", flag_common_name)
  }

  # -- 5. Vectorised substring scan ------------------------------------------
  ord        <- order(asv_info$len)
  asv_sorted <- asv_info[ord, ]
  seqs       <- asv_sorted$asv

  sub_pairs   <- vector("list", n_asvs)
  sub_key_set <- character(0)

  for (i in seq_len(n_asvs - 1)) {
    candidates <- (i + 1):n_asvs
    candidates <- candidates[asv_sorted$len[candidates] > asv_sorted$len[i]]
    if (length(candidates) == 0) next
    hits <- candidates[grepl(seqs[i], seqs[candidates], fixed = TRUE)]
    if (length(hits) > 0) {
      sub_pairs[[i]] <- hits
      sub_key_set    <- c(sub_key_set, seqs[i], seqs[hits])
    }
  }
  sub_key_set <- unique(sub_key_set)

  # -- 6. Categorise each substring pair -------------------------------------
  make_row <- function(i, j, scenario, action, keep_asv, keep_name, note) {
    ai <- asv_sorted[i, ]; aj <- asv_sorted[j, ]
    set_ai  <- parse_species_set(ai$taxon_str)
    set_aj  <- parse_species_set(aj$taxon_str)
    rel     <- taxa_set_relation(set_ai, set_aj)
    same_cn <- !is.na(ai$common_name_val) && !is.na(aj$common_name_val) &&
               ai$common_name_val == aj$common_name_val
    ka_rank <- if (!is.na(keep_asv) && keep_asv == ai$asv) ai$deepest_rank
               else if (!is.na(keep_asv))                   aj$deepest_rank
               else NA_character_
    data.frame(
      asv_i             = ai$asv,
      asv_j             = aj$asv,
      len_i             = ai$len,
      len_j             = aj$len,
      source_i          = ai$sources,
      source_j          = aj$sources,
      reads_total_i     = ai$reads_total,
      reads_total_j     = aj$reads_total,
      reads_max_i       = ai$reads_max,
      reads_max_j       = aj$reads_max,
      reads_pct_i       = ai$reads_pct,
      reads_pct_j       = aj$reads_pct,
      reads_prev_i      = ai$reads_prev,
      reads_prev_j      = aj$reads_prev,
      is_substring      = TRUE,
      rank_i            = ai$deepest_rank,
      rank_j            = aj$deepest_rank,
      deepest_name_i    = ai$deepest_name,
      deepest_name_j    = aj$deepest_name,
      taxon_i           = ai$taxon_str,
      taxon_j           = aj$taxon_str,
      common_name_i     = ai$common_name_val,
      common_name_j     = aj$common_name_val,
      same_common_name  = same_cn,
      overlap_coef      = if (!is.null(rel$overlap_coef)) rel$overlap_coef else NA_real_,
      same_resolution   = !is.na(ai$deepest_level) && !is.na(aj$deepest_level) &&
                          ai$deepest_level == aj$deepest_level,
      same_identity     = rel$type == "equal",
      scenario          = scenario,
      action            = action,
      keep_asv          = keep_asv,
      keep_name         = keep_name,
      keep_taxonomy     = if (!is.na(keep_asv) && !is.na(ka_rank) && !is.na(keep_name))
                            paste0(ka_rank, ": ", keep_name) else NA_character_,
      note              = note,
      stringsAsFactors  = FALSE
    )
  }

  flagged_rows <- list()

  for (i in seq_len(n_asvs - 1)) {
    if (is.null(sub_pairs[[i]])) next
    ai <- asv_sorted[i, ]
    for (j in sub_pairs[[i]]) {
      aj <- asv_sorted[j, ]

      both_assigned <- !is.na(ai$deepest_level) && !is.na(aj$deepest_level)
      if (!both_assigned) {
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = NA_integer_, action = "flag_unassigned",
          keep_asv = NA_character_, keep_name = NA_character_,
          note = paste0("Substring pair; one or both ASVs lack taxonomy. ",
                        "i assigned=", !is.na(ai$deepest_level),
                        ", j assigned=", !is.na(aj$deepest_level))
        )
        next
      }

      set_ai <- parse_species_set(ai$taxon_str)
      set_aj <- parse_species_set(aj$taxon_str)
      rel    <- taxa_set_relation(set_ai, set_aj)

      # same_comparison: species-set equality when sets available;
      # fall back to comparison_name string equality when both lack taxon_str.
      same_comparison <- if (length(set_ai) > 0 || length(set_aj) > 0) {
        rel$type == "equal"
      } else {
        !is.na(ai$comparison_name) && !is.na(aj$comparison_name) &&
        ai$comparison_name == aj$comparison_name
      }

      same_formal_name <- !is.na(ai$deepest_name) && !is.na(aj$deepest_name) &&
                          ai$deepest_name == aj$deepest_name
      same_res         <- ai$deepest_level == aj$deepest_level
      shared_branch    <- !is.na(ai$lineage) && !is.na(aj$lineage) &&
                          (startsWith(aj$lineage, ai$lineage) ||
                           startsWith(ai$lineage, aj$lineage))

      if (same_res && same_comparison) {
        taxa_label <- if (!is.na(ai$taxon_str)) ai$taxon_str else ai$comparison_name
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = 1L, action = "auto_merge_keep_shorter",
          keep_asv = ai$asv, keep_name = ai$comparison_name,
          note = paste0("Substring + same species set. taxa: '", taxa_label,
                        "'. apply_harmonization will transfer counts and prune the longer ASV.")
        )
      } else if (same_res && same_formal_name && !same_comparison) {
        taxa_i <- if (!is.na(ai$taxon_str)) ai$taxon_str else ai$comparison_name
        taxa_j <- if (!is.na(aj$taxon_str)) aj$taxon_str else aj$comparison_name
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = NA_integer_, action = "flag_split_distinct_names",
          keep_asv = NA_character_, keep_name = NA_character_,
          note = paste0("Substring + same formal taxonomy ('", ai$deepest_name,
                        "') but different species sets. taxa_i: '", taxa_i,
                        "' | taxa_j: '", taxa_j, "'.")
        )
      } else if (!same_res && shared_branch) {
        if (ai$deepest_level <= aj$deepest_level) {
          keep <- ai; other <- aj
        } else {
          keep <- aj; other <- ai
        }
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = 2L, action = "flag_lower_resolution",
          keep_asv = keep$asv, keep_name = keep$comparison_name,
          note = paste0("Substring + shared lineage, different resolution. ",
                        "Lower-res: '", keep$comparison_name, "' (",
                        keep$deepest_rank, ") | Higher-res: '",
                        other$comparison_name, "' (", other$deepest_rank, ").")
        )
      } else if (same_res && !same_formal_name) {
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = NA_integer_, action = "flag_conflict",
          keep_asv = NA_character_, keep_name = NA_character_,
          note = paste0("Substring pair with conflicting taxonomy at same rank (",
                        ai$deepest_rank, "): '", ai$deepest_name,
                        "' vs '", aj$deepest_name, "'. Manual inspection needed.")
        )
      } else {
        flagged_rows[[length(flagged_rows) + 1]] <- make_row(
          i, j, scenario = NA_integer_, action = "flag_conflict",
          keep_asv = NA_character_, keep_name = NA_character_,
          note = "Substring pair; different taxonomic branches. Manual inspection needed."
        )
      }
    }
  }

  # -- 7. Non-substring pairs: S3 (exact), S4 (subset), S5 (overlap) ---------
  # Iterates all pairs. For pairs with taxon_str, uses set-based comparison.
  # For pairs without taxon_str, falls back to comparison_name string equality
  # for S3 only (cannot compute S4/S5 without species sets).

  ns_row <- function(a, b, rel, scenario, action, keep_name, note) {
    same_cn <- !is.na(a$common_name_val) && !is.na(b$common_name_val) &&
               a$common_name_val == b$common_name_val
    data.frame(
      asv_i             = a$asv,
      asv_j             = b$asv,
      len_i             = a$len,
      len_j             = b$len,
      source_i          = a$sources,
      source_j          = b$sources,
      reads_total_i     = a$reads_total,
      reads_total_j     = b$reads_total,
      reads_max_i       = a$reads_max,
      reads_max_j       = b$reads_max,
      reads_pct_i       = a$reads_pct,
      reads_pct_j       = b$reads_pct,
      reads_prev_i      = a$reads_prev,
      reads_prev_j      = b$reads_prev,
      is_substring      = FALSE,
      rank_i            = a$deepest_rank,
      rank_j            = b$deepest_rank,
      deepest_name_i    = a$deepest_name,
      deepest_name_j    = b$deepest_name,
      taxon_i           = a$taxon_str,
      taxon_j           = b$taxon_str,
      common_name_i     = a$common_name_val,
      common_name_j     = b$common_name_val,
      same_common_name  = same_cn,
      overlap_coef      = if (!is.null(rel$overlap_coef)) rel$overlap_coef else NA_real_,
      same_resolution   = !is.na(a$deepest_level) && !is.na(b$deepest_level) &&
                          a$deepest_level == b$deepest_level,
      same_identity     = rel$type == "equal",
      scenario          = scenario,
      action            = action,
      keep_asv          = NA_character_,
      keep_name         = keep_name,
      keep_taxonomy     = NA_character_,
      note              = note,
      stringsAsFactors  = FALSE
    )
  }

  asv_assigned <- asv_info[!is.na(asv_info$deepest_level), ]

  if (nrow(asv_assigned) >= 2) {
    for (pi in seq_len(nrow(asv_assigned) - 1)) {
      for (pj in (pi + 1):nrow(asv_assigned)) {
        a <- asv_assigned[pi, ]; b <- asv_assigned[pj, ]

        already_substring <-
          a$asv %in% sub_key_set && b$asv %in% sub_key_set &&
          (grepl(a$asv, b$asv, fixed = TRUE) || grepl(b$asv, a$asv, fixed = TRUE))
        if (already_substring) next

        has_sets <- !is.na(a$taxon_str) && !is.na(b$taxon_str)
        same_cn  <- !is.na(a$common_name_val) && !is.na(b$common_name_val) &&
                    a$common_name_val == b$common_name_val

        if (has_sets) {
          set_a <- parse_species_set(a$taxon_str)
          set_b <- parse_species_set(b$taxon_str)
          rel   <- taxa_set_relation(set_a, set_b)

          if (rel$type == "equal") {
            # S3: identical species sets -> flag for user decision
            taxa_label <- if (!is.na(a$taxon_str)) a$taxon_str else a$comparison_name
            flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
              a, b, rel, scenario = 3L, action = "flag_same_taxon",
              keep_name = NA_character_,
              note = paste0("Non-substring ASVs with identical species set. taxa: '",
                            taxa_label, "'.")
            )

          } else if (rel$type %in% c("a_subset_b", "b_subset_a")) {
            # S4: strict subset -> flag for user decision
            if (rel$type == "a_subset_b") {
              sub_lbl <- "ASV i"; sup_lbl <- "ASV j"
              unique_sp <- paste(setdiff(set_b, set_a), collapse = "; ")
            } else {
              sub_lbl <- "ASV j"; sup_lbl <- "ASV i"
              unique_sp <- paste(setdiff(set_a, set_b), collapse = "; ")
            }
            cn_note <- if (same_cn)
              paste0(" Same common name: '", a$common_name_val, "'.")
            else ""
            flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
              a, b, rel, scenario = 4L, action = "flag_taxon_subset",
              keep_name = NA_character_,
              note = paste0(sub_lbl, " species set is a strict subset of ",
                            sup_lbl, ". Species unique to superset: ",
                            unique_sp, ".", cn_note,
                            " Review whether these represent the same food item.")
            )

          } else if (rel$type == "overlap" &&
                     (rel$overlap_coef >= min_overlap ||
                      (flag_common_name && same_cn))) {
            # S5: partial overlap >= threshold, or same common name if flag_common_name -> flag
            shared_sp <- paste(intersect(set_a, set_b), collapse = "; ")
            cn_note   <- if (same_cn)
              paste0(" Same common name: '", a$common_name_val, "'.")
            else ""
            flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
              a, b, rel, scenario = 5L, action = "flag_taxon_overlap",
              keep_name = NA_character_,
              note = paste0("Non-substring ASVs with overlapping species sets ",
                            sprintf("(overlap=%.2f). ", rel$overlap_coef),
                            "Shared species: ", shared_sp, ".", cn_note)
            )

          } else if (flag_common_name && same_cn) {
            # Same common name but low/no species overlap -> flag conflict
            flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
              a, b, rel, scenario = NA_integer_, action = "flag_conflict",
              keep_name = NA_character_,
              note = paste0("Same common name ('", a$common_name_val,
                            "') but low species overlap ",
                            sprintf("(overlap=%.2f). Manual inspection needed.",
                                    rel$overlap_coef))
            )
          }

        } else if (!has_sets && flag_common_name && same_cn) {
          # No species sets; same common name -> flag for user decision
          rel_na <- list(type = NA_character_, overlap_coef = NA_real_)
          flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
            a, b, rel_na, scenario = 3L, action = "flag_same_taxon",
            keep_name = NA_character_,
            note = paste0("Non-substring ASVs with same common name ('",
                          a$common_name_val, "') and no species-set data.")
          )

        } else if (!has_sets &&
                   !is.na(a$comparison_name) && !is.na(b$comparison_name) &&
                   a$comparison_name == b$comparison_name) {
          # No species sets; same comparison_name (deepest formal rank) -> flag
          rel_na <- list(type = NA_character_, overlap_coef = NA_real_)
          flagged_rows[[length(flagged_rows) + 1]] <- ns_row(
            a, b, rel_na, scenario = 3L, action = "flag_same_taxon",
            keep_name = NA_character_,
            note = paste0("Non-substring ASVs with same comparison_name ('",
                          a$comparison_name, "') and no species-set data.")
          )
        }
      }
    }
  }

  # -- 8. Compile flagged data frame -----------------------------------------
  if (length(flagged_rows) == 0) {
    flagged_df <- data.frame(
      asv_i = character(), asv_j = character(),
      len_i = integer(), len_j = integer(),
      source_i = character(), source_j = character(),
      reads_total_i = integer(), reads_total_j = integer(),
      reads_max_i = integer(), reads_max_j = integer(),
      reads_pct_i = numeric(), reads_pct_j = numeric(),
      reads_prev_i = numeric(), reads_prev_j = numeric(),
      is_substring = logical(),
      rank_i = character(), rank_j = character(),
      deepest_name_i = character(), deepest_name_j = character(),
      taxon_i = character(), taxon_j = character(),
      common_name_i = character(), common_name_j = character(),
      same_common_name = logical(),
      overlap_coef = numeric(),
      same_resolution = logical(), same_identity = logical(),
      scenario = integer(), action = character(),
      keep_asv = character(), keep_name = character(),
      keep_taxonomy = character(), note = character(),
      stringsAsFactors = FALSE
    )
  } else {
    flagged_df <- do.call(rbind, flagged_rows)
    rownames(flagged_df) <- NULL
  }

  # -- 9. name_updates -------------------------------------------------------
  # S2 is now flagged (not auto); all auto-corrections come from S1 only.
  # name_updates is kept as an empty vector for API compatibility.
  name_updates <- character(0)

  # -- 10. Summary table -----------------------------------------------------
  action_levels <- c(
    "auto_merge_keep_shorter",
    "flag_lower_resolution",
    "flag_same_taxon",
    "flag_split_distinct_names",
    "flag_taxon_subset",
    "flag_taxon_overlap",
    "flag_conflict",
    "flag_unassigned"
  )
  scenario_labels <- c(
    "auto_merge_keep_shorter"  = "S1  substring, same species set -- auto-merge, keep shorter",
    "flag_lower_resolution"    = "S2  substring, diff resolution, shared branch -- review",
    "flag_same_taxon"          = "S3  non-substring, same species set -- review",
    "flag_split_distinct_names"= "S1b substring, same LCA, diff species sets -- review",
    "flag_taxon_subset"        = "S4  non-substring, species set subset -- review",
    "flag_taxon_overlap"       = "S5  non-substring, species set overlap -- review",
    "flag_conflict"            = "Conflict -- taxonomy mismatch or name conflict",
    "flag_unassigned"          = "Unassigned -- missing taxonomy"
  )
  n_total    <- nrow(flagged_df)
  summary_df <- if (n_total == 0) {
    data.frame(action = action_levels, description = scenario_labels[action_levels],
               n_pairs = 0L, pct_of_total = NA_real_, stringsAsFactors = FALSE)
  } else {
    counts <- table(factor(flagged_df$action, levels = action_levels))
    data.frame(
      action       = action_levels,
      description  = scenario_labels[action_levels],
      n_pairs      = as.integer(counts),
      pct_of_total = round(100 * as.integer(counts) / n_total, 1),
      stringsAsFactors = FALSE
    )
  }

  review_actions_all <- c("flag_lower_resolution", "flag_same_taxon",
                          "flag_split_distinct_names", "flag_taxon_subset",
                          "flag_taxon_overlap", "flag_conflict", "flag_unassigned")

  if (verbose) {
    message("\nASV cross-batch scenario summary:")
    for (i in seq_len(nrow(summary_df))) {
      message(sprintf("  %-62s  n = %d  (%.1f%%)",
                      summary_df$description[i],
                      summary_df$n_pairs[i],
                      if (is.na(summary_df$pct_of_total[i])) 0
                      else summary_df$pct_of_total[i]))
    }
    n_auto   <- sum(flagged_df$action == "auto_merge_keep_shorter", na.rm = TRUE)
    n_review <- sum(flagged_df$action %in% review_actions_all,      na.rm = TRUE)
    message(sprintf("\n  %d pair(s) handled automatically, %d require review.",
                    n_auto, n_review))
  }

  asv_info_cols <- intersect(
    c("asv", "len", "sources", tax_rank_cols,
      "deepest_level", "deepest_rank", "deepest_name", "lineage",
      "reads_total", "reads_max", "reads_pct", "reads_prev",
      "taxon_str", "common_name_val", "label_val", "comparison_name"),
    colnames(asv_info)
  )

  list(
    flagged         = flagged_df,
    name_updates    = name_updates,
    summary         = summary_df,
    asv_info        = asv_info[, asv_info_cols],
    taxon_set_col   = taxon_set_col,
    common_name_col = common_name_col,
    label_col       = label_col    # backward-compat alias
  )
}


# resolve_conflicts_interactive -----------------------------------------------
#
# Reviews flagged pairs via a Shiny gadget (RStudio Viewer pane).
# Falls back to a warning if shiny/miniUI are not installed.
#
# Five cases per pair:
#   P -- touches a prior_decisions group (see below) -> join group or keep distinct
#   A -- both frozen, same group  -> auto-advance
#   B -- both frozen, diff groups -> merge or keep separate
#   C -- one frozen, one not      -> join group or keep separate
#   D -- neither frozen           -> standard options
#
# prior_asv_final (from plan_harmonization()'s `prior_decisions` argument):
# a named character vector, ASV -> the final representative ASV of the
# established (>=2-member) prior-decisions group it belongs to (covers both
# discards and the representative itself, which maps to itself). A row where
# exactly one of asv_i/asv_j is a name in this vector is case P -- the
# review only ever decides whether the OTHER (new) ASV joins that group, via
# keep_asvs (the same discard->keep mechanism ordinary merges use), never by
# touching name_updates/the group's own taxa-string state. A row where BOTH
# sides are prior_asv_final members resolving to DIFFERENT final reps is
# also case P, but only offers an "acknowledge" action -- reconciling two
# already-decided groups isn't supported incrementally (re-run
# plan_harmonization() on the full dataset instead). Rows that are an exact
# repeat of an already-resolved prior_decisions pair never reach here at all
# -- plan_harmonization() carries those forward directly into the final
# decisions table without any review step.
#
# Returns:
#   $name_updates        -- named character vector -> apply_harmonization()
#   $common_name_updates -- named character vector (group taxa key -> common name)
#   $keep_asvs           -- named character vector c("discard_asv" = "keep_asv")
#   $s1_skip_keys        -- character vector of overridden S1 pair keys
#   $rationale           -- named character vector c("asv_i|asv_j" = "rationale")
#   $skipped             -- data frame of skipped pairs


resolve_conflicts_interactive <- function(flagged, compare_result = NULL,
                                          startup_warning = NULL,
                                          prior_asv_final = character(0),
                                          exact_repeat_keys = character(0)) {
  stopifnot(is.data.frame(flagged))

  # A row touches a prior-decisions group when either side is a member of
  # one (see header comment). Computed once, up front, and consulted both
  # for to_resolve's own "P" bucket and by pair_info()/option_buttons()
  # per-row below.
  touches_prior <- if (length(prior_asv_final) > 0) {
    flagged$asv_i %in% names(prior_asv_final) | flagged$asv_j %in% names(prior_asv_final)
  } else {
    rep(FALSE, nrow(flagged))
  }
  flagged$.touches_prior <- touches_prior

  # Rows identical to an already-resolved prior_decisions pair never need
  # review here -- plan_harmonization() carries their decision forward
  # directly into the final decisions table (see .build_decisions_table()'s
  # prior_lookup argument). Excluded from `flagged` for the review UI's own
  # purposes, but plan_harmonization() re-attaches them to the full table
  # afterward.
  is_exact_repeat <- if (length(exact_repeat_keys) > 0) {
    paste0(flagged$asv_i, "|", flagged$asv_j) %in% exact_repeat_keys |
      paste0(flagged$asv_j, "|", flagged$asv_i) %in% exact_repeat_keys
  } else {
    rep(FALSE, nrow(flagged))
  }
  flagged <- flagged[!is_exact_repeat, ]

  review_actions <- c("flag_lower_resolution", "flag_same_taxon",
                      "flag_split_distinct_names", "flag_taxon_subset",
                      "flag_taxon_overlap", "flag_conflict", "flag_unassigned")
  to_resolve <- flagged[flagged$action %in% review_actions, ]
  s1_all     <- flagged[!is.na(flagged$action) &
                          flagged$action == "auto_merge_keep_shorter", ]

  if (nrow(to_resolve) == 0) {
    message("No pairs require review.")
    return(invisible(list(name_updates        = character(0),
                          common_name_updates = character(0),
                          keep_asvs           = character(0),
                          s1_skip_keys        = character(0),
                          prior_keep_distinct_keys = character(0),
                          rationale           = character(0),
                          skipped             = to_resolve)))
  }

  if (!interactive()) {
    warning("resolve_conflicts_interactive() called in a non-interactive session. ",
            "All pairs returned as skipped.")
    return(invisible(list(name_updates        = character(0),
                          common_name_updates = character(0),
                          keep_asvs           = character(0),
                          s1_skip_keys        = character(0),
                          prior_keep_distinct_keys = character(0),
                          rationale           = character(0),
                          skipped             = to_resolve)))
  }


  # name_updates value of "" means explicitly left unassigned (frozen but unnamed).
  # NA (absent from the vector) means not yet reviewed.

  # \u2500\u2500 Build ASV metadata lookup \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  has_reads <- "reads_total_i" %in% colnames(flagged)
  meta_cols_i <- c("asv_i", "len_i", "source_i",
                   if (has_reads) c("reads_total_i", "reads_max_i",
                                    "reads_pct_i",   "reads_prev_i"),
                   "taxon_i", "common_name_i", "rank_i", "deepest_name_i")
  meta_cols_j <- c("asv_j", "len_j", "source_j",
                   if (has_reads) c("reads_total_j", "reads_max_j",
                                    "reads_pct_j",   "reads_prev_j"),
                   "taxon_j", "common_name_j", "rank_j", "deepest_name_j")
  base_names <- c("asv", "len", "source",
                  if (has_reads) c("reads_total", "reads_max", "reads_pct", "reads_prev"),
                  "taxon", "common_name", "rank", "deepest_name")

  mi <- flagged[, intersect(meta_cols_i, colnames(flagged))]
  mj <- flagged[, intersect(meta_cols_j, colnames(flagged))]
  colnames(mi) <- base_names[seq_len(ncol(mi))]
  colnames(mj) <- base_names[seq_len(ncol(mj))]
  asv_meta <- unique(rbind(mi, mj))
  asv_meta  <- asv_meta[!duplicated(asv_meta$asv), ]
  rownames(asv_meta) <- asv_meta$asv

  # \u2500\u2500 Helpers \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  default_name <- function(asv) {
    a <- asv_meta[asv, ]
    if ("taxon"        %in% colnames(asv_meta) && !is.na(a$taxon)        && nzchar(a$taxon))
      return(a$taxon)
    if ("deepest_name" %in% colnames(asv_meta) && !is.na(a$deepest_name) && nzchar(a$deepest_name))
      return(a$deepest_name)
    NA_character_
  }

  # Pure version of merge_groups_to_name: returns updated nu.
  do_merge_groups <- function(nu, old_name, new_name) {
    targets <- names(nu)[!is.na(nu) & nzchar(nu) & nu == old_name]
    nu[targets] <- new_name
    nu
  }

  # Merges semicolon-delimited taxa strings into one alphabetically sorted union.
  # NA and empty strings are silently dropped; returns NA if nothing is left.
  merge_taxa_strings <- function(taxa_vec) {
    valid  <- taxa_vec[!is.na(taxa_vec) & nzchar(taxa_vec)]
    if (length(valid) == 0) return(NA_character_)
    all_sp <- unique(trimws(unlist(strsplit(valid, ";"))))
    all_sp <- all_sp[nzchar(all_sp)]
    if (length(all_sp) == 0) return(NA_character_)
    paste(sort(all_sp), collapse = "; ")
  }

  # Pure version: takes nu (name_updates vector) as argument.
  get_group_taxa_strs <- function(group_name, nu) {
    asvs <- names(nu)[!is.na(nu) & nzchar(nu) & nu == group_name]
    taxa_strs <- if ("taxon" %in% colnames(asv_meta) && length(asvs) > 0)
      asv_meta[intersect(asvs, rownames(asv_meta)), "taxon"]
    else character(0)
    c(group_name, taxa_strs)
  }

  # Determines case (P/P_conflict/A/B/B_mixed/C/C_unassigned/D) from prior-
  # decisions membership first, then frozen (name_updates) state. A row
  # touching a prior_decisions group is ALWAYS case P/P_conflict, checked
  # before anything else -- prior-group membership is a fact about the ASV
  # (from an earlier round), not something this round's naming-phase state
  # (nu) should ever override or reclassify.
  pair_info <- function(row, nu) {
    if (isTRUE(row$.touches_prior)) {
      pside <- row$asv_i %in% names(prior_asv_final)
      qside <- row$asv_j %in% names(prior_asv_final)
      if (pside && qside) {
        return(list(case = "P_conflict",
                    prior_i = unname(prior_asv_final[[row$asv_i]]),
                    prior_j = unname(prior_asv_final[[row$asv_j]])))
      }
      new_asv   <- if (pside) row$asv_j else row$asv_i
      prior_asv <- if (pside) row$asv_i else row$asv_j
      return(list(case = "P_join", new_side = if (pside) "j" else "i",
                  new_asv = new_asv, prior_asv = prior_asv,
                  prior_rep = unname(prior_asv_final[[prior_asv]])))
    }
    fi <- !is.na(nu[row$asv_i]); fj <- !is.na(nu[row$asv_j])
    gi <- if (fi) nu[row$asv_i] else NA_character_
    gj <- if (fj) nu[row$asv_j] else NA_character_
    ni <- fi && !is.na(gi) && nzchar(gi)
    nj <- fj && !is.na(gj) && nzchar(gj)
    ui <- fi && !ni; uj <- fj && !nj
    list(fi=fi,fj=fj,gi=gi,gj=gj,ni=ni,nj=nj,ui=ui,uj=uj,
         case = if (!fi && !fj) "D"
                else if (fi && fj && ni && nj && gi == gj) "A"
                else if (fi && fj && !ni && !nj) "A"
                else if (fi && fj && ((ni && uj) || (nj && ui))) "B_mixed"
                else if (fi && fj) "B"
                else if (fi) (if (ui) "C_unassigned" else "C")
                else          (if (uj) "C_unassigned" else "C"))
  }



  # -- Gadget dependency check ------------------------------------------------
  if (!requireNamespace("shiny",  quietly = TRUE) ||
      !requireNamespace("miniUI", quietly = TRUE) ||
      !requireNamespace("DT",     quietly = TRUE)) {
    message("Install required packages for the interactive gadget:\n",
            "  install.packages(c('shiny','miniUI','DT'))")
    return(invisible(list(name_updates        = character(0),
                          common_name_updates = character(0),
                          keep_asvs           = character(0),
                          s1_skip_keys        = character(0),
                          prior_keep_distinct_keys = character(0),
                          rationale           = character(0),
                          skipped             = to_resolve)))
  }

  # -- Scenario metadata ------------------------------------------------------
  action_to_scn <- c(
    auto_merge_keep_shorter   = "S1",  flag_lower_resolution     = "S2",
    flag_same_taxon           = "S3",  flag_split_distinct_names = "S1b",
    flag_taxon_subset         = "S4",  flag_taxon_overlap        = "S5",
    flag_conflict             = "conflict", flag_unassigned       = "unassigned"
  )
  scn_desc <- c(
    S1 = "Substring -- same species set (auto-merge)",
    S1b = "Substring -- same LCA, different species sets",
    S2 = "Substring -- different resolution, shared branch",
    S3 = "Non-substring -- identical species sets",
    S4 = "Non-substring -- species set subset",
    S5 = "Non-substring -- overlapping species sets",
    conflict = "Substring -- conflicting taxonomy at same rank",
    unassigned = "Substring -- missing taxonomy"
  )
  scn_colour <- c(
    S1 = "#dff0d8", S1b = "#fcf8e3", S2 = "#fcf8e3",
    S3 = "#d9edf7", S4 = "#f2dede", S5 = "#f5f5f5",
    conflict = "#f2dede", unassigned = "#f5f5f5"
  )

  # -- CSS --------------------------------------------------------------------
  gadget_css <- shiny::tags$style(shiny::HTML(
    'body,.container-fluid{padding:0 8px;font-size:13px;}
     .hdr{border-bottom:2px solid #337ab7;padding-bottom:5px;margin-bottom:10px;padding-right:115px;}
     .hdr h4{margin:0;display:inline;} .prog{float:right;font-size:11px;color:#888;margin-top:4px;}
     .apanel{border:1px solid #ddd;border-radius:4px;padding:8px;background:#fafafa;
              height:260px;overflow-y:auto;}
     .apanel.frozen{border-color:#337ab7;background:#f0f7ff;}
     .apanel.unass{border-color:#aaa;background:#fefefe;}
     .seq{font-family:monospace;font-size:10px;color:#666;word-break:break-all;margin:2px 0 4px;}
     .splist{font-size:11px;padding-left:16px;margin:3px 0;}
     .fbadge{background:#d9edf7;border:1px solid #bce8f1;border-radius:3px;
              padding:2px 5px;font-size:10px;margin-top:4px;word-break:break-word;}
     .ubadge{background:#fcf8e3;border:1px solid #faebcc;border-radius:3px;
              padding:2px 5px;font-size:10px;margin-top:4px;}
     .scnbar{border-radius:3px;padding:5px 8px;margin-bottom:6px;font-size:12px;font-weight:bold;}
     .note{font-size:11px;color:#666;font-style:italic;margin-bottom:5px;}
     .pair-subdesc{font-size:12px;color:#333;margin-bottom:6px;padding:2px 0;}
     .btn-choice{display:block; width:fit-content; margin:2px 0;} hr.thin{margin:6px 0;}
     .pair-wrap{display:flex; flex-direction:column;}
     .pair-context{flex:0 0 auto;}
     .pair-actions{flex:0 0 auto; border-top:2px solid #e0e0e0; padding-top:8px; background:#fff;}
     .global-btn{position:fixed; top:8px; right:12px; z-index:9999;}
     #flagged_dt table.dataTable td{max-width:220px !important; overflow:hidden !important;
       text-overflow:ellipsis !important; white-space:nowrap !important; cursor:text;}'
  ))

  # -- ASV panel --------------------------------------------------------------
  asv_panel <- function(asv_seq, label, nu, cn_grp) {
    if (!asv_seq %in% rownames(asv_meta))
      return(shiny::div(class = "apanel", label, ": not found"))
    a  <- asv_meta[asv_seq, ]
    fi <- !is.na(nu[asv_seq])
    gi <- if (fi) nu[asv_seq] else NA_character_
    ni <- fi && !is.na(gi) && nzchar(gi)
    cls <- paste("apanel", if (ni) "frozen" else if (fi) "unass" else "")
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
    badge <- NULL
    if (ni) {
      cn  <- if (!is.null(cn_grp) && !is.na(cn_grp[gi])) cn_grp[[gi]] else NA_character_
      cntxt <- if (!is.na(cn) && nzchar(cn)) sprintf(' ("%s")', cn) else ""
      badge <- shiny::div(class = "fbadge",
                          shiny::tags$strong("MERGED: "), paste0(gi, cntxt))
    } else if (fi) {
      badge <- shiny::div(class = "ubadge", "Left unassigned")
    }
    shiny::div(class = cls,
      shiny::tags$h6(style = "margin-top:0;", label),
      shiny::div(class = "seq", substr(asv_seq, 1, 55),
                 if (nchar(asv_seq) > 55) "..." else ""),
      shiny::div(shiny::tags$strong(sprintf("%d bp", a$len)), " | Source: ", a$source),
      reads_div, shiny::tags$hr(class = "thin"),
      shiny::div(sprintf("%s: %s",
                         if (!is.na(a$rank)) a$rank else "unassigned",
                         if (!is.na(a$deepest_name)) a$deepest_name else "(unassigned)")),
      shiny::div(shiny::tags$em("Common name: "),
                 if (!is.na(a$common_name)) a$common_name else "(none)"),
      sp_ul, badge)
  }

  # -- Flagged-pairs column defaults (for DT modal) ---------------------------
  all_flagged_cols <- if (!is.null(compare_result)) colnames(compare_result$flagged) else character(0)
  default_flagged_cols <- intersect(
    c("asv_i", "asv_j", "len_i", "len_j", "scenario",
      "source_i", "source_j",
      "reads_total_i", "reads_total_j", "reads_prev_i", "reads_prev_j",
      "taxon_i", "taxon_j", "common_name_i", "common_name_j",
      "overlap_coef", "action"),
    all_flagged_cols)

  # -- UI ---------------------------------------------------------------------
  ui <- shiny::fluidPage(
    gadget_css,
    if (!is.null(compare_result))
      shiny::div(class = "global-btn",
                 shiny::actionButton("show_flagged", "Flagged pairs",
                                     class = "btn-default btn-sm",
                                     icon  = shiny::icon("table"))),
    shiny::uiOutput("main_ui")
  )

  # -- Server -----------------------------------------------------------------
  server <- function(input, output, session) {

    rv <- shiny::reactiveValues(
      phase                 = if (!is.null(compare_result)) "summary"
                              else if (nrow(s1_all) > 0) "s1_overview"
                              else "scn_filter",
      s1_idx = 1L, s1_snaps = list(),
      keep_asvs    = character(0), s1_skip_keys = character(0),
      prior_keep_distinct_keys = character(0),
      name_updates          = character(0),
      common_name_for_group = character(0),
      rationale             = character(0),
      k = 1L, snapshots = list(), skipped_rows = list(),
      to_res  = to_resolve, n_total = nrow(to_resolve),
      merge_groups = list(), gk = 1L, seq_snaps = list(),
      auto_keep    = character(0), auto_reasons = list(),
      pending_taxa = NULL,         pending_pair = NULL
    )

    # -- Phase renderers ----------------------------------------------------
    output$main_ui <- shiny::renderUI({
      switch(rv$phase,
        summary = {
          shiny::req(!is.null(compare_result))
          summ    <- compare_result$summary
          n_asvs  <- if (!is.null(compare_result$asv_info)) nrow(compare_result$asv_info) else NA
          n_pairs <- nrow(compare_result$flagged)
          shiny::tagList(
            if (!is.null(startup_warning))
              shiny::div(class = "alert alert-warning",
                         style = "white-space:pre-wrap; font-size:0.9em;",
                         shiny::tags$strong("Warning: "), startup_warning),
            shiny::div(class = "hdr", shiny::tags$h4("compare_asvs -- Summary")),
            if (!is.na(n_asvs))
              shiny::p(class = "text-muted small",
                       sprintf("%d unique ASVs compared | %d total pairs classified",
                               n_asvs, n_pairs)),
            DT::DTOutput("summary_scenario_dt"),
            shiny::tags$hr(class = "thin"),
            shiny::p("Would you like to inspect the full flagged pairs table before proceeding?"),
            shiny::p(class = "text-muted small",
                     "(The table is also accessible throughout the session via the",
                     shiny::tags$strong("'Flagged pairs'"), "button.)"),
            shiny::actionButton("summary_view_table", "View full table",
                                class = "btn-default"),
            shiny::actionButton("summary_proceed", "Proceed to decision phase ->",
                                class = "btn-primary", style = "margin-left:8px;"))
        },
        s1_overview = {
          n <- nrow(s1_all)
          shiny::tagList(
            if (!is.null(startup_warning))
              shiny::div(class = "alert alert-warning",
                         style = "white-space:pre-wrap; font-size:0.9em;",
                         shiny::tags$strong("Warning: "), startup_warning),
            shiny::div(class = "hdr", shiny::tags$h4("S1 Auto-merge"),
                       shiny::div(class = "prog", "Phase 1 / 4")),
            shiny::p(sprintf("%d pair(s) will be auto-merged (shorter kept, longer pruned).", n)),
            shiny::p("Review any of these decisions before proceeding?"),
            shiny::actionButton("s1_proceed",    "Proceed with all auto-merges", class = "btn-success"),
            shiny::actionButton("s1_review_btn", "Review S1 pairs",
                                class = "btn-default", style = "margin-left:6px;"))
        },
        s1_review = {
          idx <- rv$s1_idx; n <- nrow(s1_all); row <- s1_all[idx, ]
          shorter <- row$keep_asv
          longer  <- if (!is.na(shorter) && shorter == row$asv_i) row$asv_j else row$asv_i
          shiny::tagList(
            shiny::div(class = "hdr",
                       shiny::tags$h4(sprintf("S1 Review -- %d / %d", idx, n)),
                       shiny::div(class = "prog", "Phase 1 / 4")),
            shiny::div(class = "pair-wrap",
              shiny::div(class = "pair-context",
                shiny::div(class = "scnbar", style = "background:#dff0d8;",
                           "S1: same species set -- auto-merge (shorter kept)"),
                if (!is.na(row$note) && nzchar(row$note)) shiny::div(class = "note", row$note),
                shiny::fluidRow(
                  shiny::column(6, asv_panel(shorter, "Shorter ASV (DEFAULT keep)",
                                             rv$name_updates, rv$common_name_for_group)),
                  shiny::column(6, asv_panel(longer,  "Longer ASV (default discard)",
                                             rv$name_updates, rv$common_name_for_group)))),
              shiny::div(class = "pair-actions",
                if (idx > 1) shiny::actionButton("s1_back", "<- Back",
                                                  class = "btn-default btn-sm"),
                shiny::actionButton("s1_shorter", "Keep shorter (DEFAULT)", class = "btn-success btn-choice"),
                shiny::actionButton("s1_longer",  "Keep longer (override)", class = "btn-warning btn-choice"),
                shiny::actionButton("s1_both",    "Keep both distinct",     class = "btn-default btn-choice"),
                shiny::textInput("s1_rat", "Rationale (optional):", value = "", width = "100%"),
                shiny::actionButton("s1_done_all", "Skip remaining S1 ->", class = "btn-link btn-sm"))))
        },
        scn_filter = {
          acts <- unique(to_resolve$action)
          scns <- action_to_scn[acts]; scns <- scns[!is.na(scns)]
          rows <- lapply(seq_along(scns), function(i) {
            lbl <- scns[i]; act <- names(scns)[i]
            cnt <- sum(to_resolve$action == act, na.rm = TRUE)
            desc <- if (!is.na(scn_desc[lbl])) scn_desc[lbl] else lbl
            shiny::fluidRow(style = "margin-bottom:4px;",
              shiny::column(1, shiny::checkboxInput(paste0("skip_", act), NULL, FALSE)),
              shiny::column(2, shiny::tags$strong(lbl)),
              shiny::column(1, sprintf("(%d)", cnt)),
              shiny::column(8, shiny::tags$span(class = "text-muted", desc)))
          })
          n_prior <- sum(to_resolve$.touches_prior, na.rm = TRUE)
          prior_row <- if (n_prior > 0) shiny::fluidRow(
            style = "margin-bottom:4px;",
            shiny::column(1, shiny::checkboxInput("skip_vs_prior", NULL, FALSE)),
            shiny::column(2, shiny::tags$strong("vs. prior decision")),
            shiny::column(1, sprintf("(%d)", n_prior)),
            shiny::column(8, shiny::tags$span(class = "text-muted",
              "Touches an ASV already resolved by prior_decisions -- skipping defaults these to keep distinct (not left unresolved).")))
          else NULL
          shiny::tagList(
            shiny::div(class = "hdr", shiny::tags$h4("Scenario Filter"),
                       shiny::div(class = "prog", "Phase 2 / 4")),
            shiny::p("Check scenario types to skip entirely:"),
            do.call(shiny::tagList, rows), prior_row, shiny::tags$hr(class = "thin"),
            shiny::actionButton("filter_go", "Proceed to decision phase ->",
                                class = "btn-primary btn-lg"))
        },
        name_assign = {
          k <- rv$k; n <- rv$n_total
          if (k > n)
            return(shiny::tagList(
              shiny::div(class = "alert alert-success",
                         shiny::tags$strong("Name assignment complete."),
                         sprintf(" %d pair(s) reviewed.", n)),
              shiny::actionButton("to_seq", "Continue to sequence selection ->",
                                  class = "btn-primary btn-lg")))
          row  <- rv$to_res[k, ]
          pi   <- pair_info(row, rv$name_updates)
          scn  <- action_to_scn[row$action]; if (is.na(scn)) scn <- row$action
          bg   <- if (!is.na(scn_colour[scn])) scn_colour[scn] else "#f5f5f5"
          can_back <- k > 1 && !is.null(tryCatch(rv$snapshots[[k - 1]], error = function(e) NULL))
          shiny::tagList(
            shiny::div(class = "hdr",
                       shiny::tags$h4(sprintf("Name Assignment -- %d / %d", k, n)),
                       shiny::div(class = "prog", "Phase 3 / 4")),
            shiny::div(class = "pair-wrap",
              shiny::div(class = "pair-context",
                shiny::div(class = "scnbar", style = sprintf("background:%s;", bg), scn),
                shiny::div(class = "pair-subdesc", pair_short_desc(row)),
                shiny::fluidRow(
                  shiny::column(6, asv_panel(row$asv_i, "ASV i",
                                             rv$name_updates, rv$common_name_for_group)),
                  shiny::column(6, asv_panel(row$asv_j, "ASV j",
                                             rv$name_updates, rv$common_name_for_group)))),
              shiny::div(class = "pair-actions",
                option_buttons(row, pi, rv$name_updates),
                shiny::textInput("na_rat", "Rationale (optional):", value = "", width = "100%"),
                if (can_back)
                  shiny::actionButton("na_back", "<- Back", class = "btn-default btn-sm"))))
        },
        seq_sel_overview = {
          mg <- rv$merge_groups; n <- length(mg)
          rows <- lapply(seq_along(mg), function(i) {
            gn   <- names(mg)[i]; asvs <- mg[[gn]]
            ka   <- rv$auto_keep
            disc <- names(ka)[names(ka) %in% asvs]
            rep_seq <- if (length(disc) > 0) unique(ka[disc])[1] else asvs[1]
            reason  <- if (!is.null(rv$auto_reasons[[gn]])) rv$auto_reasons[[gn]] else "--"
            cn      <- rv$common_name_for_group[gn]
            cn_txt  <- if (!is.null(cn) && !is.na(cn) && nzchar(cn))
                         sprintf(' ("%s")', cn) else ""
            shiny::fluidRow(style = "margin-bottom:6px; border-bottom:1px solid #eee; padding-bottom:4px;",
              shiny::column(4, shiny::tags$strong(
                substr(gn, 1, 35), if (nchar(gn) > 35) "..." else "",
                shiny::tags$em(cn_txt))),
              shiny::column(2, sprintf("%d ASVs", length(asvs))),
              shiny::column(6, shiny::tags$span(class = "text-muted small",
                shiny::tags$em("Auto: "), reason)))
          })
          shiny::tagList(
            shiny::div(class = "hdr",
              shiny::tags$h4("Sequence Selection"),
              shiny::div(class = "prog", "Phase 4 / 4")),
            shiny::p(sprintf(
              "%d merge group(s). Representatives auto-selected by decision tree:", n)),
            shiny::p(class = "text-muted small",
              "(1) Lowest taxonomic resolution -> (2) shortest sequence -> (3) highest total reads"),
            do.call(shiny::tagList, rows),
            shiny::tags$hr(class = "thin"),
            shiny::actionButton("seq_accept_all", "Accept all auto-selections",
                                class = "btn-success btn-lg"),
            shiny::actionButton("seq_review_all", "Review group by group",
                                class = "btn-default btn-lg", style = "margin-left:8px;"))
        },
        seq_sel = {
          gk <- rv$gk; mg <- rv$merge_groups; n <- length(mg)
          if (gk > n)
            return(shiny::tagList(
              shiny::div(class = "alert alert-success", "Sequence selection complete."),
              shiny::actionButton("finish", "Done -- return results",
                                  class = "btn-primary btn-lg")))
          grp_name <- names(mg)[gk]; asvs <- mg[[grp_name]]
          cn <- rv$common_name_for_group[grp_name]
          cn_txt <- if (!is.null(cn) && !is.na(cn) && nzchar(cn))
                      sprintf(' ("%s")', cn) else ""
          can_back <- gk > 1 && !is.null(tryCatch(rv$seq_snaps[[gk - 1]], error = function(e) NULL))
          # Identify auto-selected representative for this group
          ka      <- rv$auto_keep
          disc    <- names(ka)[names(ka) %in% asvs]
          auto_rep <- if (length(disc) > 0) unique(ka[disc])[1] else asvs[1]
          reason   <- if (!is.null(rv$auto_reasons[[grp_name]])) rv$auto_reasons[[grp_name]] else "auto"
          asv_rows <- lapply(seq_along(asvs), function(i) {
            a_seq   <- asvs[i]
            a       <- if (a_seq %in% rownames(asv_meta)) asv_meta[a_seq, ] else NULL
            rds     <- if (!is.null(a) && has_reads && !is.na(a$reads_total))
              sprintf(" | %s reads", format(a$reads_total, big.mark = ",")) else ""
            is_auto <- identical(a_seq, auto_rep)

            tax_str <- if (!is.null(a) && !is.na(a$rank) && nzchar(a$rank))
              sprintf("%s: %s", a$rank,
                      if (!is.na(a$deepest_name)) a$deepest_name else "?")
            else "(unassigned)"

            taxa_str <- if (!is.null(a) && "taxon" %in% colnames(asv_meta) &&
                            !is.na(a$taxon) && nzchar(a$taxon)) {
              sp <- trimws(strsplit(a$taxon, ";")[[1]])
              if (length(sp) <= 4)
                paste(sp, collapse = "; ")
              else
                paste0(paste(sp[1:3], collapse = "; "), sprintf(" ... (+%d more)", length(sp) - 3))
            } else NULL

            cn_str <- if (!is.null(a) && "common_name" %in% colnames(asv_meta) &&
                          !is.na(a$common_name) && nzchar(a$common_name))
              a$common_name else NULL

            shiny::fluidRow(
              style = sprintf("margin-bottom:8px; padding-bottom:6px; border-bottom:1px solid %s;",
                              if (is_auto) "#c3e6cb" else "#eee"),
              shiny::column(9,
                shiny::div(style = "font-family:monospace;font-size:10px; color:#555;",
                           sprintf("[%d] ", i), substr(a_seq, 1, 45), "..."),
                if (!is.null(a))
                  shiny::div(class = "text-muted small",
                             sprintf("%d bp | Source: %s%s", a$len, a$source, rds)),
                shiny::div(class = "small",
                           shiny::tags$strong("Resolution: "), tax_str,
                           if (!is.null(cn_str))
                             shiny::tagList(" | ", shiny::tags$em(cn_str))),
                if (!is.null(taxa_str))
                  shiny::div(class = "text-muted small",
                             shiny::tags$strong("Taxa: "), taxa_str),
                if (is_auto) shiny::div(class = "text-success small",
                                        shiny::tags$strong("<- AUTO-SELECTED: "), reason)),
              shiny::column(3,
                shiny::actionButton(sprintf("seq_pick_%d", i),
                  if (is_auto) sprintf("Select [%d] (default)", i)
                  else sprintf("Select [%d]", i),
                  class = if (is_auto) "btn-success btn-xs" else "btn-primary btn-xs")))
          })
          shiny::tagList(
            shiny::div(class = "hdr",
                       shiny::tags$h4(sprintf("Sequence Selection -- %d / %d", gk, n)),
                       shiny::div(class = "prog", "Phase 4 / 4")),
            shiny::p(shiny::tags$strong("Group: "), grp_name, shiny::tags$em(cn_txt)),
            shiny::p(sprintf("%d ASVs -- choose representative (green = auto-selected):",
                             length(asvs))),
            do.call(shiny::tagList, asv_rows),
            shiny::tags$hr(class = "thin"),
            if (can_back) shiny::actionButton("seq_back", "<- Back",
                                               class = "btn-default btn-sm"),
            shiny::actionButton("seq_exit", "Exit (use auto for remaining)",
                                class = "btn-link btn-sm", style = "margin-left:8px;"))
        },
        done = shiny::tagList(
          shiny::div(class = "alert alert-success", shiny::tags$h4("Review complete.")),
          shiny::actionButton("finish", "Done -- return results",
                              class = "btn-primary btn-lg"))
      )
    })

    # -- Option buttons -----------------------------------------------------
    option_buttons <- function(row, pi, nu) {
      skip_exit <- shiny::tagList(shiny::tags$br(),
        shiny::actionButton("na_skip", "Skip", class = "btn-link btn-sm"),
        shiny::actionButton("na_exit", "Exit", class = "btn-link btn-sm",
                            style = "margin-left:6px;"))

      if (pi$case == "P_join") {
        rep_label <- default_name(pi$prior_rep)
        rep_txt   <- if (!is.na(rep_label)) sprintf('"%s"', rep_label)
                     else sprintf("ASV %s", substr(pi$prior_rep, 1, 20))
        return(shiny::tagList(
          shiny::div(class = "alert alert-info",
                     style = "padding:5px 8px; font-size:11.5px;",
                     sprintf("ASV %s was already resolved by a prior decision (part of the group represented by %s). This review only decides whether the NEW ASV joins that group -- to change anything about the group itself, re-run plan_harmonization() on the full dataset instead.",
                             pi$prior_asv, rep_txt)),
          shiny::actionButton("na_prior_join",
            sprintf("Join ASV %s's group -> %s", pi$new_side, rep_txt),
            class = "btn-primary btn-choice"),
          shiny::actionButton("na_prior_distinct",
            sprintf("Keep ASV %s distinct", pi$new_side),
            class = "btn-default btn-choice"),
          skip_exit))
      }
      if (pi$case == "P_conflict") {
        ri_label <- default_name(pi$prior_i); rj_label <- default_name(pi$prior_j)
        ri_txt <- if (!is.na(ri_label)) sprintf('"%s"', ri_label) else substr(pi$prior_i, 1, 20)
        rj_txt <- if (!is.na(rj_label)) sprintf('"%s"', rj_label) else substr(pi$prior_j, 1, 20)
        return(shiny::tagList(
          shiny::div(class = "alert alert-warning",
                     style = "padding:5px 8px; font-size:11.5px;",
                     sprintf("Both ASVs here already belong to separate, already-decided prior groups (%s and %s). Reconciling two existing groups isn't supported in this incremental review -- re-run plan_harmonization() on the full dataset if they should actually be merged.", ri_txt, rj_txt)),
          shiny::actionButton("na_prior_ack", "Acknowledge (keep distinct)",
                              class = "btn-default btn-choice"),
          skip_exit))
      }

      if (pi$case == "A")
        return(shiny::div(class = "alert alert-info",
          if (pi$ni && pi$nj) sprintf('Both in group "%s". Click Continue.', pi$gi)
          else "Both unassigned. Click Continue.",
          shiny::actionButton("na_advance", "Continue ->", class = "btn-primary btn-sm",
                              style = "margin-left:8px;")))
      if (pi$case == "B") {
        sz_i <- sum(!is.na(nu) & nzchar(nu) & nu == pi$gi)
        sz_j <- sum(!is.na(nu) & nzchar(nu) & nu == pi$gj)
        return(shiny::tagList(
          shiny::actionButton("na_use_i",
            sprintf('Merge under taxa string i: "%s" (%d ASVs)', pi$gi, sz_i),
            class = "btn-primary btn-choice"),
          shiny::actionButton("na_use_j",
            sprintf('Merge under taxa string j: "%s" (%d ASVs)', pi$gj, sz_j),
            class = "btn-primary btn-choice"),
          shiny::actionButton("na_custom", "Merge under custom name",
                              class = "btn-info btn-choice"),
          shiny::actionButton("na_keep", "Keep groups separate",
                              class = "btn-default btn-choice"),
          skip_exit))
      }
      if (pi$case == "B_mixed") {
        grp   <- if (pi$ni) pi$gi else pi$gj
        uside <- if (pi$ni) "j" else "i"
        sz    <- sum(!is.na(nu) & nzchar(nu) & nu == grp)
        return(shiny::tagList(
          shiny::actionButton("na_assign_grp",
            sprintf('Add ASV %s -> "%s" (%d)', uside, grp, sz),
            class = "btn-primary btn-choice"),
          shiny::actionButton("na_keep_unass",
            sprintf("Leave ASV %s unassigned", uside), class = "btn-default btn-choice"),
          skip_exit))
      }
      if (pi$case == "C_unassigned") {
        fside <- if (pi$fi) "i" else "j"; nside <- if (pi$fi) "j" else "i"
        new_asv <- if (pi$fi) row$asv_j else row$asv_i
        nd <- default_name(new_asv)
        return(shiny::tagList(
          shiny::actionButton("na_leave_unass",
            sprintf("Leave ASV %s unassigned (same as ASV %s)", nside, fside),
            class = "btn-warning btn-choice"),
          shiny::actionButton("na_give_name",
            if (!is.na(nd)) sprintf('Freeze ASV %s: "%s"', nside, nd)
            else sprintf("Give ASV %s a name", nside),
            class = "btn-primary btn-choice"),
          skip_exit))
      }
      if (pi$case == "C") {
        fside <- if (pi$fi) "i" else "j"; nside <- if (pi$fi) "j" else "i"
        frz   <- if (pi$fi) pi$gi else pi$gj
        new_asv <- if (pi$fi) row$asv_j else row$asv_i
        nd    <- default_name(new_asv); sz <- sum(!is.na(nu) & nzchar(nu) & nu == frz)
        new_t <- if ("taxon" %in% colnames(asv_meta) && new_asv %in% rownames(asv_meta))
                   asv_meta[new_asv, "taxon"] else NA_character_
        expand_note <- if (!is.na(new_t) && nzchar(new_t) &&
                           frz %in% names(rv$common_name_for_group)) {
          mk <- merge_taxa_strings(c(frz, new_t))
          if (!identical(mk, frz))
            shiny::div(class = "alert alert-info",
                       style = "padding:3px 6px;font-size:11px;margin-bottom:4px;",
                       "Adding this ASV will expand the group's taxa string.") else NULL
        } else NULL
        sep_label <- if (!is.na(nd) && nzchar(nd) && nd != frz)
          sprintf('Keep ASV %s separate ("%s")', nside, nd)
        else sprintf("Keep ASV %s separate", nside)
        return(shiny::tagList(expand_note,
          shiny::actionButton("na_add_grp",
            sprintf('Add ASV %s -> "%s" (%d)', nside, frz, sz),
            class = "btn-primary btn-choice"),
          shiny::actionButton("na_keep_sep", sep_label, class = "btn-default btn-choice"),
          skip_exit))
      }
      # Case D
      ni <- default_name(row$asv_i); nj <- default_name(row$asv_j)
      shiny::tagList(
        shiny::actionButton("na_use_i",
          if (!is.na(ni)) sprintf('Merge under taxa string i: "%s"', ni)
          else "Merge under ASV i's taxa string",
          class = "btn-primary btn-choice"),
        shiny::actionButton("na_use_j",
          if (!is.na(nj)) sprintf('Merge under taxa string j: "%s"', nj)
          else "Merge under ASV j's taxa string",
          class = "btn-primary btn-choice"),
        shiny::actionButton("na_custom", "Merge under custom name",
                            class = "btn-info btn-choice"),
        shiny::actionButton("na_keep", "Keep as distinct", class = "btn-default btn-choice"),
        skip_exit)
    }

    # One-line pair description for the scenario banner.
    pair_short_desc <- function(row) {
      sub_txt <- if (isTRUE(row$is_substring)) "ASV i is a substring of ASV j"
                 else "non-nested sequences"
      switch(row$action,
        auto_merge_keep_shorter   = paste0(sub_txt, " -- identical species sets"),
        flag_split_distinct_names = paste0(sub_txt, " -- same formal taxon, different species sets"),
        flag_lower_resolution     = {
          ri <- if (!is.na(row$rank_i) && nzchar(row$rank_i)) row$rank_i else "?"
          rj <- if (!is.na(row$rank_j) && nzchar(row$rank_j)) row$rank_j else "?"
          paste0(sub_txt, sprintf(" -- ASV i at %s level, ASV j at %s level", ri, rj))
        },
        flag_same_taxon           = paste0(sub_txt, " -- identical species sets"),
        flag_taxon_subset         = {
          si <- trimws(strsplit(if (!is.na(row$taxon_i) && nzchar(row$taxon_i))
                                  row$taxon_i else "", ";")[[1]])
          sj <- trimws(strsplit(if (!is.na(row$taxon_j) && nzchar(row$taxon_j))
                                  row$taxon_j else "", ";")[[1]])
          if (length(si) > 0 && length(sj) > 0 && all(si %in% sj))
            "Non-nested -- ASV i species set is a subset of ASV j"
          else
            "Non-nested -- ASV j species set is a subset of ASV i"
        },
        flag_taxon_overlap        = {
          oc <- if ("overlap_coef" %in% colnames(row) && !is.na(row$overlap_coef))
                  sprintf(" (overlap = %.2f)", row$overlap_coef) else ""
          paste0(sub_txt, " -- partially overlapping species sets", oc)
        },
        flag_conflict             = {
          rk <- if (!is.na(row$rank_i) && nzchar(row$rank_i)) row$rank_i else "same rank"
          paste0(sub_txt, sprintf(" -- conflicting %s assignments", rk))
        },
        flag_unassigned           = paste0(sub_txt, " -- one or both ASVs lack taxonomy"),
        sub_txt
      )
    }

    # -- Snapshot helpers ---------------------------------------------------
    save_snap <- function() {
      rv$snapshots[[rv$k]] <- list(
        name_updates          = rv$name_updates,
        common_name_for_group = rv$common_name_for_group,
        rationale             = rv$rationale,
        n_skipped             = length(rv$skipped_rows),
        keep_asvs             = rv$keep_asvs,
        prior_keep_distinct_keys = rv$prior_keep_distinct_keys)
    }
    record_rat <- function(row) {
      r <- if (!is.null(input$na_rat)) trimws(input$na_rat) else ""
      if (nzchar(r)) rv$rationale[paste0(row$asv_i, "|", row$asv_j)] <- r
    }
    # Advance k, auto-skipping Case A pairs.
    advance <- function() {
      new_k <- rv$k + 1L; nu <- rv$name_updates
      to_r  <- rv$to_res;  n  <- rv$n_total
      while (new_k <= n) {
        pi <- pair_info(to_r[new_k, ], nu)
        if (pi$case == "A") {
          rv$snapshots[[new_k]] <- list(
            name_updates          = nu,
            common_name_for_group = rv$common_name_for_group,
            rationale             = rv$rationale,
            n_skipped             = length(rv$skipped_rows),
            keep_asvs             = rv$keep_asvs,
            prior_keep_distinct_keys = rv$prior_keep_distinct_keys)
          new_k <- new_k + 1L
        } else break
      }
      rv$k <- new_k
      shiny::updateTextInput(session, "na_rat", value = "")
    }
    cur_row <- function() rv$to_res[rv$k, ]
    cur_pi  <- function() pair_info(cur_row(), rv$name_updates)

    # -- Phase 0: Summary + on-demand flagged pairs table ------------------
    next_phase_after_summary <- function()
      if (nrow(s1_all) > 0) "s1_overview" else "scn_filter"

    shiny::observeEvent(input$summary_proceed, {
      rv$phase <- next_phase_after_summary()
    })

    open_flagged_modal <- function() {
      shiny::showModal(shiny::modalDialog(
        title = "Flagged pairs",
        size  = "l",
        shiny::checkboxGroupInput(
          "flagged_cols", "Columns to display:",
          choices  = all_flagged_cols,
          selected = default_flagged_cols,
          inline   = TRUE),
        shiny::tags$hr(class = "thin"),
        DT::DTOutput("flagged_dt"),
        shiny::tags$hr(class = "thin"),
        shiny::fluidRow(
          shiny::column(8, shiny::textInput(
            "csv_path", NULL, value = "",
            placeholder = "path/to/flagged_pairs.csv", width = "100%")),
          shiny::column(4, shiny::actionButton(
            "export_csv", "Export CSV", class = "btn-default btn-sm",
            style = "margin-top:1px;"))),
        shiny::p(class = "text-muted small", style = "margin-top:4px;",
                 "The full flagged-pairs table is returned as ",
                 shiny::tags$code("result$compare_result$flagged"),
                 " -- export here if you need it before the session ends."),
        footer    = shiny::modalButton("Close"),
        easyClose = TRUE
      ))
    }

    shiny::observeEvent(input$summary_view_table, { open_flagged_modal() })
    shiny::observeEvent(input$show_flagged,        { open_flagged_modal() })

    output$flagged_dt <- DT::renderDT({
      shiny::req(!is.null(compare_result))
      cols <- if (!is.null(input$flagged_cols) && length(input$flagged_cols) > 0)
                intersect(input$flagged_cols, colnames(compare_result$flagged))
              else default_flagged_cols
      df <- compare_result$flagged[, cols, drop = FALSE]
      for (col in intersect(c("asv_i", "asv_j"), names(df)))
        df[[col]] <- paste0(substr(df[[col]], 1, 25), "...")
      DT::datatable(df,
                    filter   = "top",
                    rownames = FALSE,
                    class    = "nowrap",
                    options  = list(
                      scrollX    = TRUE,
                      pageLength = 15,
                      autoWidth  = TRUE,
                      createdRow = DT::JS(
                        "function(row, data, index) {
                           $('td', row).each(function() {
                             var txt = $.trim($(this).text());
                             if (txt.length > 0) $(this).attr('title', txt);
                           });
                         }"
                      )
                    ))
    }, server = TRUE)

    output$summary_scenario_dt <- DT::renderDT({
      shiny::req(!is.null(compare_result))
      summ <- compare_result$summary
      DT::datatable(
        summ[, c("description", "n_pairs", "pct_of_total")],
        colnames  = c("Scenario", "N pairs", "% of total"),
        rownames  = FALSE,
        options   = list(dom = "t", pageLength = nrow(summ), ordering = FALSE),
        selection = "none")
    }, server = FALSE)


    shiny::observeEvent(input$export_csv, {
      shiny::req(!is.null(compare_result))
      path <- trimws(input$csv_path)
      if (!nzchar(path))
        path <- sprintf("flagged_pairs_%s.csv",
                        format(Sys.time(), "%Y%m%d_%H%M%S"))
      cols <- if (!is.null(input$flagged_cols) && length(input$flagged_cols) > 0)
                intersect(input$flagged_cols, colnames(compare_result$flagged))
              else colnames(compare_result$flagged)
      tryCatch({
        utils::write.csv(compare_result$flagged[, cols, drop = FALSE],
                  file = path, row.names = FALSE)
        shiny::showNotification(paste("Exported:", path), type = "message",
                                duration = 5)
      }, error = function(e)
        shiny::showNotification(paste("Export failed:", e$message),
                                type = "error", duration = 8))
    })

    # -- Phase 1: S1 overview -----------------------------------------------
    shiny::observeEvent(input$s1_proceed,    { rv$phase <- "scn_filter" })
    shiny::observeEvent(input$s1_done_all,   { rv$phase <- "scn_filter" })
    shiny::observeEvent(input$s1_review_btn, {
      rv$s1_idx <- 1L; rv$s1_snaps <- list(); rv$phase <- "s1_review" })

    apply_s1 <- function(choice) {
      idx <- rv$s1_idx; row <- s1_all[idx, ]
      pk  <- paste0(row$asv_i, "|", row$asv_j)
      asv_s <- row$keep_asv
      asv_l <- if (!is.na(asv_s) && asv_s == row$asv_i) row$asv_j else row$asv_i
      rv$s1_snaps[[idx]] <- list(keep_asvs    = rv$keep_asvs,
                                  s1_skip_keys = rv$s1_skip_keys,
                                  rationale    = rv$rationale)
      rat <- if (!is.null(input$s1_rat)) trimws(input$s1_rat) else ""
      if (nzchar(rat)) rv$rationale[pk] <- rat
      if (choice == "longer") {
        rv$keep_asvs[asv_s] <- asv_l; rv$s1_skip_keys <- c(rv$s1_skip_keys, pk)
      } else if (choice == "both") {
        rv$s1_skip_keys <- c(rv$s1_skip_keys, pk)
      }
      shiny::updateTextInput(session, "s1_rat", value = "")
      if (rv$s1_idx >= nrow(s1_all)) rv$phase <- "scn_filter"
      else rv$s1_idx <- rv$s1_idx + 1L
    }
    shiny::observeEvent(input$s1_shorter, { apply_s1("shorter") })
    shiny::observeEvent(input$s1_longer,  { apply_s1("longer")  })
    shiny::observeEvent(input$s1_both,    { apply_s1("both")    })
    shiny::observeEvent(input$s1_back, {
      idx <- rv$s1_idx; if (idx <= 1) return()
      prev <- tryCatch(rv$s1_snaps[[idx - 1]], error = function(e) NULL)
      if (is.null(prev)) return()
      rv$keep_asvs    <- prev$keep_asvs
      rv$s1_skip_keys <- prev$s1_skip_keys
      rv$rationale    <- prev$rationale
      rv$s1_snaps[[idx - 1]] <- NULL
      rv$s1_idx <- idx - 1L
      shiny::updateTextInput(session, "s1_rat", value = "")
    })

    # -- Phase 2: Scenario filter -------------------------------------------
    shiny::observeEvent(input$filter_go, {
      skip_acts <- character(0)
      for (act in unique(to_resolve$action))
        if (isTRUE(input[[paste0("skip_", act)]])) skip_acts <- c(skip_acts, act)
      remaining <- to_resolve
      if (length(skip_acts) > 0) {
        rv$skipped_rows <- c(rv$skipped_rows,
                             list(to_resolve[to_resolve$action %in% skip_acts, ]))
        remaining <- remaining[!remaining$action %in% skip_acts, ]
      }
      # "vs. prior decision" skip is NOT a plain shelve-as-unresolved like
      # the scenario-type checkboxes above -- per plan_harmonization()'s
      # design, skipping these defaults every one to keep_distinct (the new
      # ASV stays out of the prior group) rather than leaving them blank.
      if (isTRUE(input$skip_vs_prior)) {
        prior_rows <- remaining[remaining$.touches_prior, ]
        if (nrow(prior_rows) > 0) {
          rv$prior_keep_distinct_keys <- c(rv$prior_keep_distinct_keys,
            paste0(prior_rows$asv_i, "|", prior_rows$asv_j))
        }
        remaining <- remaining[!remaining$.touches_prior, ]
      }
      rv$to_res <- remaining
      rv$n_total <- nrow(rv$to_res); rv$k <- 1L
      # Auto-skip initial Case A pairs
      nu <- rv$name_updates
      while (rv$k <= rv$n_total) {
        pi <- pair_info(rv$to_res[rv$k, ], nu)
        if (pi$case != "A") break
        rv$snapshots[[rv$k]] <- list(name_updates = nu,
          common_name_for_group = rv$common_name_for_group,
          rationale = rv$rationale, n_skipped = length(rv$skipped_rows),
          keep_asvs = rv$keep_asvs,
          prior_keep_distinct_keys = rv$prior_keep_distinct_keys)
        rv$k <- rv$k + 1L
      }
      rv$phase <- "name_assign"
    })

    # -- Phase 3: Name assignment -------------------------------------------
    shiny::observeEvent(input$na_back, {
      k    <- rv$k
      snap <- tryCatch(rv$snapshots[[k - 1]], error = function(e) NULL)
      if (is.null(snap)) return()
      rv$name_updates          <- snap$name_updates
      rv$common_name_for_group <- if (!is.null(snap$common_name_for_group))
                                     snap$common_name_for_group else character(0)
      rv$rationale             <- snap$rationale
      rv$skipped_rows          <- rv$skipped_rows[seq_len(snap$n_skipped)]
      rv$keep_asvs             <- if (!is.null(snap$keep_asvs))
                                     snap$keep_asvs else character(0)
      rv$prior_keep_distinct_keys <- if (!is.null(snap$prior_keep_distinct_keys))
                                        snap$prior_keep_distinct_keys else character(0)
      rv$snapshots[[k - 1]]    <- NULL
      rv$k <- k - 1L
    })

    shiny::observeEvent(input$na_advance, { save_snap(); advance() })

    shiny::observeEvent(input$na_use_i, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      if (pi$case == "B") {
        rv$name_updates <- do_merge_groups(
          do_merge_groups(rv$name_updates, pi$gj, pi$gi), pi$gi, pi$gi)
      } else {
        chosen <- default_name(row$asv_i)
        if (!is.na(chosen)) {
          nu <- rv$name_updates; nu[row$asv_i] <- chosen; nu[row$asv_j] <- chosen
          rv$name_updates <- nu
        }
      }
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_use_j, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      if (pi$case == "B") {
        rv$name_updates <- do_merge_groups(
          do_merge_groups(rv$name_updates, pi$gi, pi$gj), pi$gj, pi$gj)
      } else {
        chosen <- default_name(row$asv_j)
        if (!is.na(chosen)) {
          nu <- rv$name_updates; nu[row$asv_i] <- chosen; nu[row$asv_j] <- chosen
          rv$name_updates <- nu
        }
      }
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_custom, {
      row <- cur_row(); pi <- cur_pi()
      if (pi$case == "B") {
        all_t <- c(get_group_taxa_strs(pi$gi, rv$name_updates),
                   get_group_taxa_strs(pi$gj, rv$name_updates))
      } else {
        ti <- if ("taxon" %in% colnames(asv_meta) && row$asv_i %in% rownames(asv_meta))
               asv_meta[row$asv_i, "taxon"] else NA_character_
        tj <- if ("taxon" %in% colnames(asv_meta) && row$asv_j %in% rownames(asv_meta))
               asv_meta[row$asv_j, "taxon"] else NA_character_
        all_t <- c(ti, tj)
      }
      rv$pending_taxa <- merge_taxa_strings(all_t)
      rv$pending_pair <- list(row = row, pi = pi)
      shiny::showModal(shiny::modalDialog(title = "Custom merge",
        shiny::p(shiny::tags$strong("Auto-computed taxa string:")),
        shiny::div(style = "background:#f5f5f5;padding:6px;border-radius:3px;font-size:11px;word-break:break-word;",
                   if (!is.na(rv$pending_taxa)) rv$pending_taxa
                   else "(no taxa data -- enter below)"),
        if (is.na(rv$pending_taxa))
          shiny::textInput("custom_taxa_manual", "Taxa string (semicolon-separated):",
                           width = "100%"),
        shiny::textInput("custom_cn", "Common name for merged group:", width = "100%"),
        footer = shiny::tagList(shiny::modalButton("Cancel"),
          shiny::actionButton("confirm_custom", "Confirm merge", class = "btn-primary"))))
    })

    shiny::observeEvent(input$confirm_custom, {
      row  <- rv$pending_pair$row; pi <- rv$pending_pair$pi
      taxa <- rv$pending_taxa
      if (is.na(taxa)) {
        manual <- trimws(input$custom_taxa_manual)
        taxa   <- if (nzchar(manual)) manual else "(no taxa)"
      }
      cn <- trimws(input$custom_cn)
      shiny::removeModal(); save_snap()
      nu <- rv$name_updates; cn_grp <- rv$common_name_for_group
      if (pi$case == "B") {
        nu <- do_merge_groups(nu, pi$gi, taxa)
        nu <- do_merge_groups(nu, pi$gj, taxa)
      } else {
        nu[row$asv_i] <- taxa; nu[row$asv_j] <- taxa
      }
      if (nzchar(cn)) cn_grp[taxa] <- cn
      rv$name_updates <- nu; rv$common_name_for_group <- cn_grp
      record_rat(row); advance()
    })

    # -- Prior-decisions review (case P_join/P_conflict) --------------------
    # Feeds keep_asvs directly (the same discard->keep mechanism ordinary
    # merges use) rather than name_updates -- joining an existing prior
    # group must never touch that group's own taxa-string/name state, only
    # attach the new ASV to it.
    shiny::observeEvent(input$na_prior_join, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      rv$keep_asvs[pi$new_asv] <- pi$prior_rep
      record_rat(row); advance()
    })
    shiny::observeEvent(input$na_prior_distinct, {
      row <- cur_row(); save_snap()
      rv$prior_keep_distinct_keys <- c(rv$prior_keep_distinct_keys,
                                       paste0(row$asv_i, "|", row$asv_j))
      record_rat(row); advance()
    })
    shiny::observeEvent(input$na_prior_ack, {
      row <- cur_row(); save_snap()
      rv$prior_keep_distinct_keys <- c(rv$prior_keep_distinct_keys,
                                       paste0(row$asv_i, "|", row$asv_j))
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_keep, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      if (pi$case != "B") {
        ni <- default_name(row$asv_i); nj <- default_name(row$asv_j)
        if (is.na(ni)) ni <- ""; if (is.na(nj)) nj <- ""
        if (nzchar(ni) && nzchar(nj) && ni == nj) {
          shiny::showModal(shiny::modalDialog(title = "Name collision",
            shiny::p(sprintf('Both default to "%s". Enter different name for ASV j:', ni)),
            shiny::textInput("distinct_nj", "Name for ASV j (empty = unassigned):",
                             width = "100%"),
            footer = shiny::tagList(shiny::modalButton("Cancel"),
              shiny::actionButton("confirm_distinct", "Confirm", class = "btn-primary"))))
          return()
        }
        nu <- rv$name_updates; nu[row$asv_i] <- ni; nu[row$asv_j] <- nj
        rv$name_updates <- nu
      }
      record_rat(row); advance()
    })
    shiny::observeEvent(input$confirm_distinct, {
      row <- cur_row()
      ni  <- default_name(row$asv_i); if (is.na(ni)) ni <- ""
      nj  <- trimws(input$distinct_nj)
      shiny::removeModal()
      nu  <- rv$name_updates; nu[row$asv_i] <- ni; nu[row$asv_j] <- nj
      rv$name_updates <- nu; record_rat(row); advance()
    })

    shiny::observeEvent(input$na_assign_grp, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      u_asv <- if (pi$ui) row$asv_i else row$asv_j
      grp   <- if (pi$ni) pi$gi else pi$gj
      new_t <- if ("taxon" %in% colnames(asv_meta) && u_asv %in% rownames(asv_meta))
                 asv_meta[u_asv, "taxon"] else NA_character_
      nu <- rv$name_updates; cn_grp <- rv$common_name_for_group
      if (!is.na(new_t) && nzchar(new_t) && grp %in% names(cn_grp)) {
        mk <- merge_taxa_strings(c(grp, new_t))
        if (!identical(mk, grp)) {
          nu <- do_merge_groups(nu, grp, mk)
          cn_grp[mk] <- cn_grp[grp]; cn_grp <- cn_grp[names(cn_grp) != grp]; grp <- mk
        }
      }
      nu[u_asv] <- grp; rv$name_updates <- nu; rv$common_name_for_group <- cn_grp
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_keep_unass, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      u_asv <- if (pi$ui) row$asv_i else row$asv_j
      nu <- rv$name_updates; nu[u_asv] <- ""; rv$name_updates <- nu
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_leave_unass, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      new_asv <- if (pi$fi) row$asv_j else row$asv_i
      nu <- rv$name_updates; nu[new_asv] <- ""; rv$name_updates <- nu
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_give_name, {
      row <- cur_row(); pi <- cur_pi()
      new_asv <- if (pi$fi) row$asv_j else row$asv_i
      nd      <- default_name(new_asv)
      if (!is.na(nd) && nzchar(nd)) {
        save_snap()
        nu <- rv$name_updates; nu[new_asv] <- nd; rv$name_updates <- nu
        record_rat(row); advance()
      } else {
        shiny::showModal(shiny::modalDialog(title = "Name for ASV",
          shiny::textInput("give_name_in", "Name (empty = leave unassigned):", width = "100%"),
          footer = shiny::tagList(shiny::modalButton("Cancel"),
            shiny::actionButton("confirm_give", "Confirm", class = "btn-primary"))))
      }
    })
    shiny::observeEvent(input$confirm_give, {
      row <- cur_row(); pi <- cur_pi()
      new_asv <- if (pi$fi) row$asv_j else row$asv_i
      nm <- trimws(input$give_name_in)
      shiny::removeModal(); save_snap()
      nu <- rv$name_updates; nu[new_asv] <- nm; rv$name_updates <- nu
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_add_grp, {
      row <- cur_row(); pi <- cur_pi(); save_snap()
      new_asv  <- if (pi$fi) row$asv_j else row$asv_i
      frz_name <- if (pi$fi) pi$gi else pi$gj
      new_t    <- if ("taxon" %in% colnames(asv_meta) && new_asv %in% rownames(asv_meta))
                    asv_meta[new_asv, "taxon"] else NA_character_
      nu <- rv$name_updates; cn_grp <- rv$common_name_for_group; eff <- frz_name
      if (!is.na(new_t) && nzchar(new_t) && frz_name %in% names(cn_grp)) {
        mk <- merge_taxa_strings(c(frz_name, new_t))
        if (!identical(mk, frz_name)) {
          nu <- do_merge_groups(nu, frz_name, mk)
          cn_grp[mk] <- cn_grp[frz_name]; cn_grp <- cn_grp[names(cn_grp) != frz_name]
          eff <- mk
        }
      }
      nu[new_asv] <- eff; rv$name_updates <- nu; rv$common_name_for_group <- cn_grp
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_keep_sep, {
      row <- cur_row(); pi <- cur_pi()
      new_asv  <- if (pi$fi) row$asv_j else row$asv_i
      frz_name <- if (pi$fi) pi$gi else pi$gj
      nd       <- default_name(new_asv)
      if (!is.na(nd) && nzchar(nd) && nd == frz_name) {
        shiny::showModal(shiny::modalDialog(title = "Name collision",
          shiny::p(sprintf('Default "%s" matches frozen group. Enter different name or leave empty to unassign.', frz_name)),
          shiny::textInput("sep_name_in", "Name:", width = "100%"),
          footer = shiny::tagList(shiny::modalButton("Cancel"),
            shiny::actionButton("confirm_sep", "Confirm", class = "btn-primary"))))
      } else {
        save_snap()
        resolved <- if (!is.na(nd) && nzchar(nd)) nd else ""
        nu <- rv$name_updates; nu[new_asv] <- resolved; rv$name_updates <- nu
        record_rat(row); advance()
      }
    })
    shiny::observeEvent(input$confirm_sep, {
      row <- cur_row(); pi <- cur_pi()
      new_asv <- if (pi$fi) row$asv_j else row$asv_i
      nm <- trimws(input$sep_name_in)
      shiny::removeModal(); save_snap()
      nu <- rv$name_updates; nu[new_asv] <- nm; rv$name_updates <- nu
      record_rat(row); advance()
    })

    shiny::observeEvent(input$na_skip, {
      row <- cur_row(); save_snap()
      rv$skipped_rows[[length(rv$skipped_rows) + 1]] <- row; advance()
    })
    shiny::observeEvent(input$na_exit, {
      k <- rv$k; n <- rv$n_total; to_r <- rv$to_res
      if (k <= n) rv$skipped_rows <- c(rv$skipped_rows, list(to_r[k:n, , drop = FALSE]))
      rv$k <- n + 1L
    })

    shiny::observeEvent(input$to_seq, {
      valid_u <- rv$name_updates[!is.na(rv$name_updates) & nzchar(rv$name_updates)]
      grps <- if (length(valid_u) > 0) {
        g <- split(names(valid_u), valid_u); g <- lapply(g, unique)
        g[vapply(g, length, integer(1)) >= 2L]
      } else list()

      # -- Auto-selection decision tree -------------------------------------
      rank_spec_order <- c("superkingdom","phylum","class","order","family",
                           "genus","species","subspecies","varietas","forma")
      get_rank_spec <- function(rank_name) {
        if (is.null(rank_name) || is.na(rank_name) || !nzchar(rank_name)) return(0L)
        idx <- match(tolower(trimws(rank_name)), rank_spec_order)
        if (is.na(idx)) 0L else as.integer(idx)
      }
      auto_pick_rep <- function(asvs) {
        stats <- lapply(asvs, function(a) {
          d <- if (a %in% rownames(asv_meta)) asv_meta[a, ] else NULL
          list(
            asv       = a,
            rank_spec = if (!is.null(d) && !is.na(d$rank)) get_rank_spec(d$rank) else 0L,
            rank_name = if (!is.null(d) && !is.na(d$rank)) as.character(d$rank) else NA_character_,
            len       = if (!is.null(d)) as.integer(d$len) else nchar(a),
            reads     = if (!is.null(d) && has_reads && !is.na(d$reads_total))
                          as.numeric(d$reads_total) else 0
          )
        })
        min_rs <- min(sapply(stats, `[[`, "rank_spec"))
        cands  <- stats[sapply(stats, `[[`, "rank_spec") == min_rs]
        if (length(cands) == 1) {
          rn <- cands[[1]]$rank_name
          return(list(keep   = cands[[1]]$asv,
                      reason = sprintf("lowest resolution (%s)",
                                       if (!is.na(rn)) rn else "unassigned")))
        }
        min_len <- min(sapply(cands, `[[`, "len"))
        cands   <- cands[sapply(cands, `[[`, "len") == min_len]
        if (length(cands) == 1) {
          rn <- cands[[1]]$rank_name
          return(list(keep   = cands[[1]]$asv,
                      reason = sprintf("tied resolution (%s) -> shortest (%d bp)",
                                       if (!is.na(rn)) rn else "unassigned", min_len)))
        }
        max_rd <- max(sapply(cands, `[[`, "reads"))
        cands  <- cands[sapply(cands, `[[`, "reads") == max_rd]
        rn <- cands[[1]]$rank_name
        list(keep   = cands[[1]]$asv,
             reason = sprintf("tied resolution/length -> most reads (%s)",
                              format(max_rd, big.mark = ",")))
      }

      # Seed from rv$keep_asvs, not character(0) -- na_prior_join() (case
      # P_join, above) writes directly into rv$keep_asvs during the
      # name_assign phase, and rv$auto_keep is the baseline every later
      # "reset to auto defaults" path (seq_accept_all/seq_review_all/
      # seq_exit, all below) resets rv$keep_asvs back to -- losing this
      # seed here would silently drop every prior-group-join decision the
      # moment sequence selection starts.
      auto_ka <- rv$keep_asvs; auto_rs <- list()
      for (gn in names(grps)) {
        res      <- auto_pick_rep(grps[[gn]])
        keep_seq <- res$keep
        for (d in grps[[gn]][grps[[gn]] != keep_seq]) auto_ka[d] <- keep_seq
        auto_rs[[gn]] <- res$reason
      }

      rv$merge_groups  <- grps
      rv$auto_keep     <- auto_ka
      rv$auto_reasons  <- auto_rs
      rv$keep_asvs     <- auto_ka   # default: accept auto selections
      rv$gk            <- 1L
      rv$seq_snaps     <- list()
      rv$phase <- if (length(grps) > 0) "seq_sel_overview" else "done"
    })

    # -- Phase 4: Sequence selection ----------------------------------------
    # Overview: accept all auto-selections or enter group-by-group review
    shiny::observeEvent(input$seq_accept_all, {
      rv$keep_asvs <- rv$auto_keep
      rv$phase     <- "done"
    })
    shiny::observeEvent(input$seq_review_all, {
      rv$keep_asvs <- rv$auto_keep   # start review from auto defaults
      rv$gk        <- 1L
      rv$seq_snaps <- list()
      rv$phase     <- "seq_sel"
    })

    # Group-by-group review
    for (.i in 1:20) local({
      idx <- .i
      shiny::observeEvent(input[[sprintf("seq_pick_%d", idx)]], {
        if (rv$phase != "seq_sel") return()
        gk <- rv$gk; mg <- rv$merge_groups
        if (gk > length(mg)) return()
        asvs <- mg[[gk]]; if (idx > length(asvs)) return()
        rv$seq_snaps[[gk]] <- rv$keep_asvs
        keep_seq <- asvs[idx]
        for (d in asvs[asvs != keep_seq]) rv$keep_asvs[d] <- keep_seq
        rv$gk <- gk + 1L
      }, ignoreInit = TRUE)
    })

    shiny::observeEvent(input$seq_back, {
      gk <- rv$gk
      snap_ka <- tryCatch(rv$seq_snaps[[gk - 1]], error = function(e) NULL)
      if (is.null(snap_ka)) return()
      rv$keep_asvs <- snap_ka; rv$seq_snaps[[gk - 1]] <- NULL; rv$gk <- gk - 1L
    })
    shiny::observeEvent(input$seq_exit, {
      # Apply auto defaults for remaining unreviewed groups then finish
      rv$keep_asvs <- rv$auto_keep
      rv$phase     <- "done"
    })

    # -- Finish -------------------------------------------------------------
    shiny::observeEvent(input$finish, {
      skipped_df <- if (length(rv$skipped_rows) > 0)
        do.call(rbind, lapply(rv$skipped_rows,
                              function(x) if (is.data.frame(x)) x else as.data.frame(x)))
      else flagged[0L, ]
      shiny::stopApp(invisible(list(
        name_updates             = rv$name_updates,
        common_name_updates      = rv$common_name_for_group,
        keep_asvs                = rv$keep_asvs,
        s1_skip_keys             = rv$s1_skip_keys,
        prior_keep_distinct_keys = rv$prior_keep_distinct_keys,
        rationale                = rv$rationale,
        skipped                  = skipped_df)))
    })
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
  shiny::runGadget(ui, server, viewer = shiny::paneViewer(minHeight = 600))
}




# Builds the per-pair decisions table: `flagged` plus decision/chosen_name/
# rationale columns. S1 (auto-merge) rows are pre-filled with the
# algorithmic default (overridable by editing `decision`); every other row
# starts blank unless `manual` (from resolve_conflicts_interactive(), if that
# interactive step ran) already carried a decision for it. Keyed directly by
# full ASV sequence (flagged$asv_i/asv_j) -- no separate ID-lookup table
# needed, unlike the CSV mechanism this replaced.
#
# `prior_lookup` (optional, from plan_harmonization()'s `prior_decisions`):
# a named list, keyed by "asv_i|asv_j" in the PRIOR round's own orientation,
# each element list(decision=, chosen_name=, rationale=). A row here that
# matches one of these keys (in either orientation) is an exact repeat of an
# already-resolved prior pair -- its decision is carried forward directly,
# short-circuiting every other fill rule below (highest priority; these
# rows are never shown to the reviewer at all -- see
# resolve_conflicts_interactive()'s exact_repeat_keys filtering). Carrying
# forward is orientation-aware: compare_asvs() doesn't guarantee asv_i/asv_j
# order is stable across separate runs (e.g. the same two ASVs could sort
# into opposite i/j slots depending on what else is being compared this
# round), so a bare "merge_keep_i"/"merge_keep_j" from the prior round is
# reoriented to match THIS round's i/j assignment before being applied --
# otherwise the merge direction could silently flip onto the wrong ASV.
.build_decisions_table <- function(flagged, manual, prior_lookup = NULL,
                                    prior_keep_asvs = character(0)) {
  n <- nrow(flagged)
  decision      <- rep(NA_character_, n)
  chosen_name   <- rep(NA_character_, n)
  rationale_out <- rep(NA_character_, n)

  pair_key     <- paste0(flagged$asv_i, "|", flagged$asv_j)
  pair_key_rev <- paste0(flagged$asv_j, "|", flagged$asv_i)

  reorient_decision <- function(dec, reversed) {
    if (!reversed || is.na(dec)) return(dec)
    if (identical(dec, "merge_keep_i")) return("merge_keep_j")
    if (identical(dec, "merge_keep_j")) return("merge_keep_i")
    dec   # "keep_distinct" or a literal external ASV sequence -- orientation-independent
  }

  default_name_for <- function(i, side) {
    for (col in paste0(c("common_name_", "taxon_", "deepest_name_"), side)) {
      v <- flagged[[col]][i]
      if (!is.null(v) && !is.na(v) && nzchar(v)) return(v)
    }
    NA_character_
  }

  # Resolve an ASV through a (possibly chained/hub) discard->keep map to its
  # final representative -- the same walk apply_harmonization() itself
  # uses for combined_keep. resolve_conflicts_interactive() builds
  # multi-ASV merge groups as a hub: every non-representative member maps
  # directly to one shared representative, not pairwise to each other. A
  # plain "does asv_i map straight to asv_j (or vice versa)" check misses
  # that, so both sides must be resolved to their final target before
  # comparing.
  resolve_final <- function(asv, map) {
    seen <- character(0)
    cur  <- asv
    while (cur %in% names(map)) {
      if (cur %in% seen) return(NA_character_)  # cycle guard
      seen <- c(seen, cur)
      cur  <- unname(map[[cur]])
    }
    cur
  }

  # Merged with the PRIOR round's own discard->keep chain (from
  # plan_harmonization()'s prior_decisions, if any) -- not just this round's
  # fresh manual$keep_asvs. Needed so that a P_join decision (new_asv ->
  # prior_rep, in manual$keep_asvs) resolves correctly even when the prior
  # ASV touched by this round's pair was itself a non-representative member
  # of its prior group (not the group's own representative) -- otherwise
  # resolve_final() below has no way to walk from that member on to its
  # true final representative, and this row would be wrongly left
  # unresolved despite the join already having been decided. This round's
  # own entries take precedence on the rare chance a key appears in both.
  ka <- c(manual$keep_asvs, prior_keep_asvs)
  ka <- ka[!duplicated(names(ka))]

  for (i in seq_len(n)) {
    # Exact repeat of an already-resolved prior_decisions pair -- highest
    # priority, short-circuits every other fill rule for this row.
    if (!is.null(prior_lookup)) {
      hit <- prior_lookup[[pair_key[i]]]
      reversed <- FALSE
      if (is.null(hit)) {
        hit      <- prior_lookup[[pair_key_rev[i]]]
        reversed <- TRUE
      }
      if (!is.null(hit)) {
        decision[i]      <- reorient_decision(hit$decision, reversed)
        chosen_name[i]   <- hit$chosen_name
        rationale_out[i] <- hit$rationale
        next
      }
    }

    if (identical(flagged$action[i], "auto_merge_keep_shorter")) {
      keep <- flagged$keep_asv[i]
      if (!is.na(keep) && identical(keep, flagged$asv_i[i])) {
        decision[i]    <- "merge_keep_i"
        chosen_name[i] <- default_name_for(i, "i")
      } else if (!is.na(keep) && identical(keep, flagged$asv_j[i])) {
        decision[i]    <- "merge_keep_j"
        chosen_name[i] <- default_name_for(i, "j")
      }
    }

    # Carry over decisions already made via resolve_conflicts_interactive().
    # Resolve both sides of the pair through the full keep_asvs map (not
    # just a direct i<->j lookup) so multi-ASV merge groups are recognized
    # as resolved even when neither asv_i nor asv_j is the group's chosen
    # representative -- in that case the literal representative ASV
    # sequence is recorded rather than forcing it into the i/j vocabulary.
    if (length(ka) > 0) {
      final_i <- resolve_final(flagged$asv_i[i], ka)
      final_j <- resolve_final(flagged$asv_j[i], ka)
      if (!is.na(final_i) && !is.na(final_j) && identical(final_i, final_j)) {
        if (identical(final_i, flagged$asv_i[i])) {
          decision[i] <- "merge_keep_i"
        } else if (identical(final_i, flagged$asv_j[i])) {
          decision[i] <- "merge_keep_j"
        } else {
          decision[i] <- final_i
        }
      }
    }
    if (length(manual$name_updates) > 0) {
      nm_i <- manual$name_updates[flagged$asv_i[i]]
      nm_j <- manual$name_updates[flagged$asv_j[i]]
      # Both sides having a (non-NA) name_updates entry means this pair was
      # actually visited during resolve_conflicts_interactive()'s naming
      # phase -- name_updates starts as an empty vector, so an unvisited ASV
      # has no entry at all (NA on lookup), unlike a visited-but-left-
      # unassigned one (explicitly set to ""). Only a shared, non-empty name
      # actually merges the pair (mirrors the exact filter
      # resolve_conflicts_interactive() itself uses to build merge groups:
      # nzchar() on both, dropping empty/unassigned entries first) -- that
      # case is already handled by decision[i] via keep_asvs above. Anything
      # else where both sides *were* visited (matching names but empty/
      # unassigned, or genuinely different names) means the user reviewed
      # this pair and explicitly chose not to merge it, so it's resolved as
      # "keep_distinct" rather than left blank/unresolved.
      if (length(nm_i) && length(nm_j) && !is.na(nm_i) && !is.na(nm_j)) {
        if (nzchar(nm_i) && nzchar(nm_j) && identical(unname(nm_i), unname(nm_j))) {
          chosen_name[i] <- unname(nm_i)
        } else if (is.na(decision[i])) {
          decision[i] <- "keep_distinct"
        }
      }
    }
    if (length(manual$rationale) > 0) {
      rv <- manual$rationale[pair_key[i]]
      if (is.na(rv)) rv <- manual$rationale[pair_key_rev[i]]
      if (length(rv) && !is.na(rv)) rationale_out[i] <- unname(rv)
    }
    if (!is.null(manual$s1_skip_keys) && pair_key[i] %in% manual$s1_skip_keys) {
      decision[i] <- "keep_distinct"
    }
    # "vs. prior decision" pairs explicitly kept distinct (na_prior_distinct/
    # na_prior_ack) or bulk-skipped via scn_filter's own checkbox -- both
    # populate this the same way, always resolving to keep_distinct rather
    # than left blank/unresolved.
    if (!is.null(manual$prior_keep_distinct_keys) &&
        pair_key[i] %in% manual$prior_keep_distinct_keys) {
      decision[i] <- "keep_distinct"
    }
  }

  cbind(flagged, data.frame(decision = decision, chosen_name = chosen_name,
                            rationale = rationale_out, stringsAsFactors = FALSE))
}

#' @title Plan ASV Harmonization Across Batches
#'
#' @description Runs \code{compare_asvs()} to detect cross-batch (and
#'   within-object) ASV redundancy, then \code{resolve_conflicts_interactive()}
#'   for interactive review, and packages the result as a
#'   \code{harmonization_plan} object: a per-pair \code{decisions} table
#'   (one row per flagged pair, with \code{decision}/\code{chosen_name}/
#'   \code{rationale} columns -- pre-filled for algorithmic S1 auto-merges
#'   and anything already decided interactively, blank otherwise) plus a
#'   \code{summary} table of pair counts by scenario. Pass the plan (or its
#'   \code{$decisions} table directly, or a path to it saved as a CSV) to
#'   \code{\link{apply_harmonization}} to apply the decisions.
#'
#'   The returned plan intentionally does **not** carry a phyloseq object --
#'   \code{\link{apply_harmonization}} takes the same \code{ps} (or the
#'   merged result of your own \code{ps_list}) directly, so the plan itself
#'   stays a small, self-contained, easily-saved (e.g. \code{write.csv()} on
#'   \code{$decisions}) artifact.
#'
#' @param ps_list Named list of phyloseq objects, or a single phyloseq
#'   object.
#' @param source_names Optional batch labels; length must equal
#'   \code{ps_list}. Defaults to \code{names(ps_list)} or sequential
#'   integers.
#' @param min_overlap Numeric; S5 species-set overlap threshold. Default
#'   \code{0.01}.
#' @param flag_common_name If \code{TRUE}, also flag S5 pairs that share a
#'   common name. Default \code{FALSE}.
#' @param verbose If \code{TRUE}, print a scenario summary. Default
#'   \code{TRUE}.
#' @param tax_rank_cols Character vector of formal taxonomic rank column
#'   names; inferred from the phyloseq tax table when \code{NULL}.
#' @param save_path Optional file path to write \code{$decisions} to as a
#'   CSV immediately (\code{write.csv(..., row.names = FALSE)}), for a
#'   ready-to-hand-edit file on disk without a separate step. If \code{NULL}
#'   (default), nothing is written -- \code{$decisions} is only returned
#'   in-memory.
#' @param prior_decisions Optional decisions from an earlier
#'   \code{plan_harmonization()} round on a related (not necessarily
#'   identical) set of objects -- the list it returned, its
#'   \code{$decisions} data frame, or a path to one saved as a CSV. Reduces
#'   review work and keeps a continuous project's decisions consistent:
#'   a pair identical to an already-resolved \code{prior_decisions} row is
#'   carried forward automatically (never re-shown for review); a
#'   genuinely new pair that touches an ASV \code{prior_decisions} already
#'   resolved is flagged with \code{touches_prior_decision = TRUE} and
#'   surfaced as its own skippable category in the interactive review --
#'   reviewing one only ever decides whether the \emph{new} ASV joins that
#'   already-established group or stays distinct, never anything about the
#'   existing group itself (a note in the review UI recommends a full
#'   re-run for that); skipping the category defaults every such pair to
#'   \code{"keep_distinct"} rather than leaving it unresolved. If
#'   \code{NULL} (default), every pair is reviewed fresh, as before.
#'
#' @return A \code{harmonization_plan} object (a list with class
#'   \code{"harmonization_plan"}) with two elements: \code{decisions} (one
#'   row per flagged pair -- including \code{source_i}/\code{source_j}
#'   source-batch columns; \code{touches_prior_decision}, set when
#'   \code{prior_decisions} was supplied and this pair touches an ASV it
#'   already resolved; and two manifest columns, \code{plan_samples}
#'   and \code{plan_taxa_hashes}, recording \code{ps_list}'s full original
#'   sample names and a compact hash of every original ASV, used by
#'   \code{\link{apply_harmonization}}'s sample/taxa checkpoints -- ready
#'   to hand-edit and pass to \code{\link{apply_harmonization}}) and
#'   \code{summary} (pair counts by scenario). When \code{prior_decisions}
#'   is supplied, \code{decisions} is always the full cumulative table
#'   (everything from \code{prior_decisions} plus whatever this round
#'   resolved) -- ready to pass straight back in as \code{prior_decisions}
#'   for the next round, with no manual stitching required.
#'
#' @examples
#' \dontrun{
#' plan   <- plan_harmonization(ps_list, source_names = names(ps_list))
#' # review/edit plan$decisions (in R, or via write.csv() + a spreadsheet editor)
#' result <- apply_harmonization(ps, plan)
#' # result$ps, result$asv_outcomes, result$unresolved, etc.
#'
#' # Later, a new batch arrives -- reuse plan's decisions, only review what's new:
#' plan2 <- plan_harmonization(ps_list2, prior_decisions = plan)
#' }
#'
#' @export
plan_harmonization <- function(ps_list,
                                source_names     = NULL,
                                min_overlap      = 0.01,
                                flag_common_name = FALSE,
                                verbose          = TRUE,
                                tax_rank_cols    = NULL,
                                save_path        = NULL,
                                prior_decisions  = NULL) {

  if (inherits(ps_list, "phyloseq")) ps_list <- list(ps_list)
  stopifnot(is.list(ps_list), length(ps_list) >= 1L)
  for (i in seq_along(ps_list))
    stopifnot(inherits(ps_list[[i]], "phyloseq"))

  .strip_df <- function(df) {
    as.data.frame(
      lapply(df, function(x) {
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
  }

  # -- Normalize prior_decisions (same 3 input types apply_harmonization()
  # accepts for `decisions`) --------------------------------------------------
  prior_df <- NULL
  if (!is.null(prior_decisions)) {
    prior_df <- if (inherits(prior_decisions, "harmonization_plan")) {
      prior_decisions$decisions
    } else if (is.data.frame(prior_decisions)) {
      prior_decisions
    } else if (is.character(prior_decisions) && length(prior_decisions) == 1) {
      utils::read.csv(prior_decisions, stringsAsFactors = FALSE, na.strings = "")
    } else {
      stop("`prior_decisions` must be the list returned by plan_harmonization(), a ",
           "decisions data.frame, or a file path to a decisions CSV.")
    }
    stopifnot(is.data.frame(prior_df))
    prior_df <- .strip_df(prior_df)
  }

  # Step 1: detect redundant pairs
  compare_result <- compare_asvs(
    ps_list,
    source_names     = source_names,
    min_overlap      = min_overlap,
    flag_common_name = flag_common_name,
    verbose          = verbose,
    tax_rank_cols    = tax_rank_cols
  )
  compare_result$flagged <- .strip_df(compare_result$flagged)

  # -- prior_decisions: exact-repeat lookup + prior-group membership --------
  # Reuses apply_harmonization.R's .decisions_to_vectors() to turn
  # prior_df's own `decision` column into a discard->keep map -- exactly the
  # same construction/validation already used to apply a decisions table,
  # just borrowed here to characterize what a PRIOR round already resolved.
  prior_lookup    <- NULL   # "asv_i|asv_j" (prior orientation) -> decision/chosen_name/rationale
  prior_asv_final <- character(0)   # ASV -> final rep, for established (>=2-member) prior groups
  prior_keep_asvs <- character(0)   # prior round's own raw discard->keep map (for resolve_final chains)

  if (!is.null(prior_df) && nrow(prior_df) > 0) {
    decided <- prior_df[!is.na(prior_df$decision), ]
    if (nrow(decided) > 0) {
      prior_lookup <- stats::setNames(
        lapply(seq_len(nrow(decided)), function(i) list(
          decision    = decided$decision[i],
          chosen_name = if ("chosen_name" %in% colnames(decided)) decided$chosen_name[i] else NA_character_,
          rationale   = if ("rationale"   %in% colnames(decided)) decided$rationale[i]   else NA_character_
        )),
        paste0(decided$asv_i, "|", decided$asv_j)
      )
    }

    prior_vecs      <- .decisions_to_vectors(prior_df)
    prior_keep_asvs <- prior_vecs$keep_asvs

    if (length(prior_keep_asvs) > 0) {
      resolve_final_prior <- function(asv, map) {
        seen <- character(0); cur <- asv
        while (cur %in% names(map)) {
          if (cur %in% seen) return(NA_character_)
          seen <- c(seen, cur); cur <- unname(map[[cur]])
        }
        cur
      }
      finals <- stats::setNames(
        vapply(names(prior_keep_asvs), resolve_final_prior, character(1), map = prior_keep_asvs),
        names(prior_keep_asvs)
      )
      reps <- unique(unname(finals))
      prior_asv_final <- c(finals, stats::setNames(reps, reps))
      prior_asv_final <- prior_asv_final[!duplicated(names(prior_asv_final))]
    }
  }

  exact_repeat_keys <- if (!is.null(prior_lookup))
    unique(c(names(prior_lookup),
             vapply(strsplit(names(prior_lookup), "|", fixed = TRUE),
                    function(x) paste0(x[2], "|", x[1]), character(1))))
  else character(0)

  # touches_prior_decision: kept as its own column, orthogonal to `scenario`
  # -- a pair's underlying taxonomic relationship (S1-S5/conflict/
  # unassigned) never changes based on prior_decisions; this only tags
  # whether it also happens to connect to already-resolved prior work.
  # Exact repeats are NOT tagged here -- they're not "new" information for
  # the reviewer, just a straight carry-forward (see prior_lookup above).
  touches_prior <- if (length(prior_asv_final) > 0) {
    (compare_result$flagged$asv_i %in% names(prior_asv_final) |
       compare_result$flagged$asv_j %in% names(prior_asv_final)) &
      !(paste0(compare_result$flagged$asv_i, "|", compare_result$flagged$asv_j) %in% exact_repeat_keys)
  } else {
    rep(FALSE, nrow(compare_result$flagged))
  }
  compare_result$flagged$touches_prior_decision <- touches_prior

  # Step 2: interactive review (Shiny gadget)
  manual <- resolve_conflicts_interactive(
    flagged           = compare_result$flagged,
    compare_result    = compare_result,
    prior_asv_final   = prior_asv_final,
    exact_repeat_keys = exact_repeat_keys
  )

  decisions <- .build_decisions_table(compare_result$flagged, manual,
                                      prior_lookup    = prior_lookup,
                                      prior_keep_asvs = prior_keep_asvs)

  # -- Sample/taxa manifest, for apply_harmonization()'s checkpoints ------
  # Captures the FULL original sample and taxa universe from ps_list (not
  # just the flagged subset) -- a taxon that was never flagged is still
  # part of the plan and shouldn't later look "foreign" to apply_asv_
  # corrections(). Stored directly on `decisions` (one delimited string per
  # column, repeated identically on every row) rather than a companion
  # file, so the checkpoint travels with the CSV itself and can't drift out
  # of sync with a separately-stored manifest. Sample names are stored
  # directly (small, human-readable, useful for debugging a mismatch); ASV
  # sequences are stored as short digest::digest() hashes instead of the
  # raw sequences themselves -- apply_harmonization() only ever needs a
  # membership test (is *this* ASV, which it already has in hand, a member
  # of the plan's original set?), not the ability to reconstruct the list,
  # and raw sequences for a few thousand ASVs would risk exceeding Excel's
  # per-cell character limit where the hashes comfortably don't.
  all_samples <- sort(unique(unlist(lapply(ps_list, phyloseq::sample_names))))
  all_taxa    <- sort(unique(unlist(lapply(ps_list, phyloseq::taxa_names))))
  taxa_hashes <- vapply(all_taxa, function(x)
    digest::digest(x, algo = "xxhash32", serialize = FALSE), character(1))

  n_dec <- nrow(decisions)
  decisions$plan_samples      <- rep(paste(all_samples, collapse = ";"), n_dec)
  decisions$plan_taxa_hashes  <- rep(paste(unname(taxa_hashes), collapse = ";"), n_dec)

  if (!is.null(save_path)) {
    utils::write.csv(decisions, save_path, row.names = FALSE, na = "")
    if (verbose) message("Decisions written to '", save_path, "'.")
  }

  structure(
    list(
      decisions = decisions,
      summary   = compare_result$summary
    ),
    class = "harmonization_plan"
  )
}

