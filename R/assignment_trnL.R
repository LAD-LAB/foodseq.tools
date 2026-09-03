# This file contains modified code originally from the dada2 project.
# Original copyright: (c) 2015-2025 Benjamin J. Callahan and contributors
# Original license: GNU Lesser General Public License (LGPL-2)
#
# Modifications made by Ashish Subramanian, 2025.
#
# This file is therefore licensed under the LGPL-2,
# in compliance with the original dada2 license.
#
# See the LICENSE file in this repository for full details.

assignSpecies_Taxonomy <- function(
    seqs,
    refFasta,
    tryRC = FALSE,
    n = 2000,
    ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "varietas", "forma"),
    sep = ";",
    drop_nohit = FALSE
) {
  na_tokens <- c("NA", "N/A", "", "Unassigned", "unassigned", "Unclassified", "unclassified", "Unknown", "unknown", "nan", "NaN", "NULL", "None")
  # deps: ShortRead, Biostrings, dplyr, tidyr, tibble, purrr
  seqs <- dada2::getSequences(seqs)

  # Read references
  refsr   <- ShortRead::readFasta(refFasta)
  headers <- methods::as(ShortRead::id(refsr), "character")
  refs    <- ShortRead::sread(refsr)

  # Detect if headers look like assignTaxonomy (semicolon-delimited)
  header_parts_1 <- strsplit(headers[[1]], sep, fixed = FALSE)[[1]]
  looks_taxonomy <- (length(header_parts_1) >= 3) || grepl(";", headers[[1]], fixed = TRUE)

  # Build reference metadata
  # - If taxonomy-style, keep the whole semicolon string
  # - Also derive a fallback "Species" from first two tokens for non-taxonomy refs
  tokens_12 <- strsplit(headers, "\\s+")
  species_fallback <- vapply(tokens_12, function(x) paste(utils::head(x, 2), collapse = " "), character(1))

  ref_meta <- tibble::tibble(
    ref_idx = seq_along(headers),
    header  = headers,
    taxonomy_string = if (looks_taxonomy) headers else NA_character_,
    Species = if (!looks_taxonomy) species_fallback else NA_character_
  )

  # Prepare hits container and group queries by length (performance)
  hits <- vector("list", length(seqs))
  lens <- nchar(seqs)

  counts <- as.integer(table(lens))
  total_batches <- sum(ceiling(counts / n))

  pb <- utils::txtProgressBar(min = 0, max = total_batches, initial = 0, style = 3)
  on.exit(close(pb), add = TRUE)
  done <- 0

  for (len in unique(lens)) {
    i.len <- which(lens == len)
    n.len <- length(i.len)
    j.lo  <- 1
    j.hi  <- min(n, n.len)

    while (j.lo <= n.len) {
      i.loop <- i.len[j.lo:j.hi]

      seqdict <- Biostrings::PDict(seqs[i.loop])
      vhit <- (Biostrings::vcountPDict(seqdict, refs) > 0)
      if (tryRC) {
        vhit_rc <- (Biostrings::vcountPDict(seqdict, Biostrings::reverseComplement(refs)) > 0)
        vhit <- vhit | vhit_rc
      }

      # Store a logical vector of hits per query (by original index)
      hits[i.loop] <- lapply(seq_len(nrow(vhit)), function(x) vhit[x, ])

      done <- done + 1
      utils::setTxtProgressBar(pb, done)

      j.lo <- j.hi + 1
      j.hi <- min(j.hi + n, n.len)
      rm(seqdict); gc()
    }
  }

  close(pb)

  # Convert hits list -> long data.frame of asv x ref_idx matches
  rows <- lapply(seq_along(seqs), function(i) {
    lh <- hits[[i]]
    if (is.null(lh)) lh <- rep(FALSE, length(headers))
    idx <- which(lh)
    if (length(idx) == 0L) {
      if (drop_nohit) {
        return(NULL)
      } else {
        return(data.frame(asv = seqs[i], ref_idx = NA_integer_, stringsAsFactors = FALSE))
      }
    } else {
      data.frame(asv = seqs[i], ref_idx = idx, stringsAsFactors = FALSE)
    }
  })
  m <- dplyr::bind_rows(rows)

  if (is.null(m) || nrow(m) == 0L) {
    # No hits at all
    out_empty <- tibble::tibble(asv = character())
    if (looks_taxonomy) {
      for (r in ranks) out_empty[[r]] <- character()
    } else {
      out_empty[["Species"]] <- character()
    }
    return(out_empty)
  }

  # Join with reference metadata to get taxonomy string or species fallback
  m <- dplyr::left_join(m, ref_meta, by = "ref_idx")

  if (looks_taxonomy) {
    # Split taxonomy_string into rank columns
    # Determine max rank length across refs to avoid truncation
    split_vec <- strsplit(m$taxonomy_string, sep, fixed = FALSE)
    max_len   <- max(vapply(split_vec, length, integer(1), USE.NAMES = FALSE))

    # If provided ranks < max_len, pad with generic names; if > max_len, extra ranks become NA
    if (length(ranks) < max_len) {
      ranks <- c(ranks, paste0("Rank", seq_len(max_len - length(ranks))))
    } else if (length(ranks) > max_len) {
      ranks <- ranks[seq_len(max_len)]
    }

    # Build a data.frame of split ranks
    pad_to <- function(x, n) { length(x) <- n; x }
    rank_mat <- do.call(rbind, lapply(split_vec, pad_to, n = length(ranks)))

    rank_mat <- apply(rank_mat, 2, trimws)
    na_upper <- toupper(na_tokens)
    is_na_token <- is.na(rank_mat) | (toupper(rank_mat) %in% na_upper)
    rank_mat[is_na_token] <- NA_character_

    rank_df  <- as.data.frame(rank_mat, stringsAsFactors = FALSE)
    colnames(rank_df) <- ranks

    out <- dplyr::bind_cols(
      m[, c("asv"), drop = FALSE],
      rank_df
    ) %>%
      dplyr::distinct()  # in case of duplicate identical matches

    return(out)
  } else {
    # Non-taxonomy headers: return ASV + Species (first two tokens)
    out <- m %>%
      dplyr::transmute(asv, Species) %>%
      dplyr::distinct()
    return(out)
  }
}

#' @title trnL Assignment
#'
#' @description This function creates a taxonomy table for trnL ASVs.
#'
#' @param qiime_asvtab The qiime-formatted ASV table
#' @param ref_fasta The path to the trnL reference FASTA
#' @param out_dir The path to the output directory
#' @param tryRC TRUE or FALSE; whether to check reverse complements
#' @param manual_taxids NULL; manual additions to be made
#' @param verbose TRUE or FALSE; whether to output messages
#' @return A taxonomy table, formatted as a character matrix; a list of statistics regarding assignment,
#' and a CSV of NAs saved to the output directory.
#' @export
assignment_trnL <- function(qiime_asvtab = qiime.asvtab.trnL,
                             ref_fasta = ref.trnL,
                             out_dir = qiime.dir.trnL,
                             tryRC = TRUE,
                             manual_taxids = NULL,
                             verbose = TRUE) {

  if (is.null(qiime_asvtab)) {
    qiime_asvtab <- get("qiime.asvtab.trnL", envir = parent.frame())
  }

  if (is.null(ref_fasta)) {
    ref_fasta <- get("ref.trnL", envir = parent.frame())
  }

  ranks <- c("superkingdom","phylum","class","order","family",
             "genus","species","subspecies","varietas","forma")

  # Adaptation of assignSpecies_Taxonomy for an assignTaxonomy-formatted reference
  table <- assignSpecies_Taxonomy(qiime_asvtab, ref_fasta, tryRC = tryRC)

  # Clean up subspecific ranks (subspecies, varietas, forma) BEFORE condenseTaxa
  # These should only be kept if ALL matches for an ASV agree on them
  subspecific_ranks <- c("subspecies", "varietas", "forma")
  table_clean <- table %>%
    dplyr::group_by(asv) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(subspecific_ranks),
      ~ if (dplyr::n() > 1) {
        # If multiple matches exist, only keep subspecific rank if all agree (ignoring NAs)
        non_na_vals <- .[!is.na(.)]
        if (length(non_na_vals) > 0 && length(unique(non_na_vals)) > 1) {
          NA_character_  # Disagreement: set to NA
        } else {
          .  # Agreement or all NA: keep as is
        }
      } else {
        .  # Single match: keep as is
      }
    )) %>%
    dplyr::ungroup()

  taxaTable <- as.matrix(table_clean[, ranks, drop = FALSE])
  groupings <- table_clean$asv
  assignments <- taxonomizr::condenseTaxa(taxaTable, groupings = groupings)


  # Normalize orientation: ensure rows = ASVs, cols = ranks
  if (is.null(colnames(assignments)) || !all(ranks %in% colnames(assignments))) {
    assignments <- t(assignments)
  }
  # Ensure all target ranks exist and ordered
  missing_cols <- setdiff(ranks, colnames(assignments))
  if (length(missing_cols)) {
    for (mc in missing_cols) {
      assignments <- cbind(assignments, stats::setNames(rep(NA_character_, nrow(assignments)), mc))
    }
  }
  assignments <- assignments[, ranks, drop = FALSE]

  # ---- 2) Reorder rows to match ASV order in qiime_asvtab (if present) ----
  asv_order <- colnames(qiime_asvtab)
  if (!is.null(asv_order)) {
    keep <- intersect(asv_order, rownames(assignments))
    rest <- setdiff(rownames(assignments), keep)
    assignments <- rbind(assignments[keep, , drop = FALSE],
                         assignments[rest, , drop = FALSE])
  }

  # ---- 3) Coerce to character matrix ----
  assignments_mat <- as.matrix(assignments)
  storage.mode(assignments_mat) <- "character"

  # ---- 4) ASVs with NA at *superkingdom* + counts; write CSV ----
  asvs_na <- rownames(assignments_mat)[is.na(assignments_mat[, "superkingdom"])]

  if (length(asvs_na)) {
    # Number of samples each ASV appears in (reads > 0)
    samples_with_reads <- colSums(qiime_asvtab[, asvs_na, drop = FALSE] > 0, na.rm = TRUE)
    # Total reads per ASV across all samples
    total_reads_by_asv <- colSums(qiime_asvtab[, asvs_na, drop = FALSE], na.rm = TRUE)

    # align orders
    samples_with_reads <- samples_with_reads[asvs_na]
    total_reads_by_asv <- total_reads_by_asv[asvs_na]

    unassigned_df <- data.frame(
      ASV                = asvs_na,
      samples_with_reads = as.integer(samples_with_reads),
      total_reads        = as.numeric(total_reads_by_asv),
      stringsAsFactors   = FALSE
    )
  } else {
    unassigned_df <- data.frame(
      ASV = character(0),
      samples_with_reads = integer(0),
      total_reads = numeric(0)
    )
  }

  if (!is.null(out_dir) && !is.na(out_dir) && dir.exists(out_dir)) {
    readr::write_csv(unassigned_df, file.path(out_dir, "trnL_NAs.csv"))
  } else if (verbose) {
    message("Skipping CSV export: Output directory is not set.")
  }

  # ---- 5) Metrics ----
  # Unassigned ASVs (Species NA before any LCA condensing)
  unassigned_asvs <- table$asv[is.na(table$species)]

  n_asvs <- ncol(qiime_asvtab)
  pct_asvs_assigned <- 100 * (1 - (length(unassigned_asvs)) / n_asvs)

  total_reads <- sum(qiime_asvtab)
  reads_unassigned <- if (length(unassigned_asvs)) {
    sum(qiime_asvtab[, unassigned_asvs, drop = FALSE])
  } else 0
  pct_reads_unassigned <- 100 * reads_unassigned / total_reads
  pct_reads_assigned   <- 100 - pct_reads_unassigned

  if (verbose) cat(sprintf(
    "%.2f%% of ASVs have been assigned, covering %.2f%% of reads in the dataset.\n",
    pct_asvs_assigned, pct_reads_assigned
  ))

  per_rank_percent_assigned <- colSums(!is.na(assignments_mat)) / nrow(assignments_mat) * 100
  per_rank_percent_assigned <- per_rank_percent_assigned[ranks]

  # print names + values, one per line
  if (verbose) cat("\nBy taxonomic rank: \n")
  if (verbose) cat(paste0(
    names(per_rank_percent_assigned), ": ",
    sprintf("%.2f", per_rank_percent_assigned), "%"
  ), sep = "\n")
  cat("\n")

  # ---- 6) Return: character matrix + metrics ----
  list(
    assignments = assignments_mat,                  # character matrix (ASV x ranks)
    metrics = list(
      percent_asvs_assigned       = pct_asvs_assigned,  # scalar %
      percent_reads_assigned      = pct_reads_assigned, # scalar %
      per_rank_percent_assigned   = per_rank_percent_assigned # named vector (% per rank)
    ),
    NAs_trnL = unassigned_df
  )
}
