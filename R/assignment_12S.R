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

assignTaxonomy_mod <- function (seqs, refFasta, minBoot = 50, tryRC = FALSE, outputBootstraps = FALSE,
                                taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family",
                                              "Genus", "Species"), multithread = FALSE, verbose = FALSE)
{
  seqs <- dada2::getSequences(seqs)
  refsr <- ShortRead::readFasta(refFasta)
  lens <- Biostrings::width(ShortRead::sread(refsr))
  refs <- as.character(ShortRead::sread(refsr))
  tax <- as.character(ShortRead::id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x))
  if (!grepl(";", tax[[1]])) {
    if (length(unlist(strsplit(tax[[1]], "\\s"))) == 3) {
      stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
    }
    else {
      stop("Incorrect reference file format for assignTaxonomy.")
    }
  }
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- max(tax.depth)
  for (i in seq(length(tax))) {
    if (tax.depth[[i]] < td) {
      for (j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol = td,
                    byrow = TRUE)
  tax.df <- as.data.frame(tax.mat)
  for (i in seq(ncol(tax.df))) {
    tax.df[, i] <- factor(tax.df[, i])
    tax.df[, i] <- as.integer(tax.df[, i])
  }
  tax.mat.int <- as.matrix(tax.df)
  if (is.logical(multithread)) {
    if (multithread == TRUE) {
      RcppParallel::setThreadOptions(numThreads = "auto")
    }
    else {
      RcppParallel::setThreadOptions(numThreads = 1)
    }
  }
  else if (is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  assignment <- dada2:::C_assign_taxonomy2(seqs, dada2::rc(seqs), refs, ref.to.genus,
                                   tax.mat.int, tryRC, verbose)
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,
  ] >= minBoot])
  tax.out <- matrix(NA_character_, nrow = length(seqs), ncol = td)
  for (i in seq(length(seqs))) {
    if (length(taxes[[i]]) > 0) {
      tax.out[i, 1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
  tax.out[tax.out == "_DADA2_UNSPECIFIED"] <- NA_character_

  ## --- Added NA-token cleanup (only new code) ---
  na_tokens <- c("NA","N/A","","Unassigned","unassigned","Unclassified","unclassified",
                 "Unknown","unknown","nan","NaN","NULL","None","_DADA2_UNSPECIFIED")
  # trim whitespace, then replace any token (case-insensitive) with NA
  tax.out <- apply(tax.out, 2, function(col) trimws(col))
  mask <- tolower(tax.out) %in% tolower(na_tokens)
  tax.out[mask] <- NA_character_
  ## ----------------------------------------------

  if (outputBootstraps) {
    boots.out <- matrix(boots, nrow = length(seqs), ncol = td)
    rownames(boots.out) <- seqs
    colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
    list(tax = tax.out, boot = boots.out)
  }
  else {
    tax.out
  }
}

#' @title 12S Assignment
#'
#' @description This function creates a taxonomy table for 12S ASVs.
#'
#' @param qiime_asvtab The qiime-formatted ASV table
#' @param ref_fasta The path to the 12S reference FASTA
#' @param out_dir The path to the output directory
#' @param tryRC TRUE or FALSE; whether to check reverse complements
#' @param verbose TRUE or FALSE; whether to output messages
#' @return A taxonomy table, formatted as a character matrix; a list of statistics regarding assignment,
#' and a CSV of NAs saved to the output directory.
#' @export
assignment_12S <- function(qiime_asvtab = qiime.asvtab.12S,
                                ref_fasta = ref.12S,
                                out_dir = qiime.dir.12S,
                                tryRC = TRUE,
                                verbose = TRUE) {

  taxLevels <- c("kingdom","phylum","class","order",
                "family","genus","species","subspecies")

  # 1) Taxonomy
  taxtab <- assignTaxonomy_mod(qiime_asvtab,
                                  refFasta = ref_fasta,
                                  tryRC = tryRC,
                                  taxLevels = taxLevels)

  # Coerce to data.frame and ensure rank columns exist (in case some are missing)
  taxonomy <- as.data.frame(taxtab, stringsAsFactors = FALSE)
  for (r in taxLevels) if (!r %in% names(taxonomy)) taxonomy[[r]] <- NA_character_

  # Ensure ASV IDs are available as rownames (DADA2 uses rownames = ASV IDs)
  asv_ids <- rownames(taxonomy)
  if (is.null(asv_ids)) stop("Row names (ASV IDs) are missing in the DADA2 taxonomy result.")
  taxonomy$asv <- asv_ids

  # 2) Compute metrics
  unassigned_asvs <- taxonomy$asv[is.na(taxonomy$kingdom)]

  n_asvs <- ncol(qiime_asvtab)
  pct_asvs_assigned <- 100 * (1 - length(unassigned_asvs) / n_asvs)

  total_reads <- sum(qiime_asvtab)
  reads_unassigned <- if (length(unassigned_asvs)) {
    # assume ASV column names of the ASV table match rownames(taxonomy)
    present <- intersect(unassigned_asvs, colnames(qiime_asvtab))
    if (length(present)) sum(qiime_asvtab[, present, drop = FALSE]) else 0
  } else 0
  pct_reads_unassigned <- 100 * reads_unassigned / total_reads
  pct_reads_assigned   <- 100 - pct_reads_unassigned

  # 3) Clean up subspecific ranks BEFORE condenseTaxa
  # Subspecies should only be kept if there's no ambiguity
  # If multiple matches exist for an ASV, subspecies is NA'd out if they disagree
  subspecific_ranks <- c("subspecies")
  taxonomy_clean <- taxonomy %>%
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

  # Condense to last common ancestor by ASV
  # taxonomizr::condenseTaxa expects a data.frame of ranks; group by ASV ids
  rank_df <- taxonomy_clean[, taxLevels, drop = FALSE]
  assignments <- taxonomizr::condenseTaxa(rank_df, groupings = taxonomy_clean$asv)

  # Normalize orientation to rows = ASVs, cols = ranks
  if (is.null(colnames(assignments)) || !all(taxLevels %in% colnames(assignments))) {
    assignments <- t(assignments)
  }
  # Add any missing rank columns and reorder
  missing_cols <- setdiff(taxLevels, colnames(assignments))
  if (length(missing_cols)) {
    for (mc in missing_cols) {
      assignments <- cbind(assignments, stats::setNames(rep(NA_character_, nrow(assignments)), mc))
    }
  }
  assignments <- assignments[, taxLevels, drop = FALSE]

  # target_asv <- "ATCGCTCGGCCTGCCCCATACAACAAACCTTCACCTCGTCGTGAGGTGAGTCTCAGTCATACGTAAAAATCTAGCGTCTAGCCTACCTAATATCGGTTTG"
  #
  # manual_row <- c(
  #   "synthetic 12S ASV",  # kingdom
  #   "synthetic 12S ASV",  # phylum
  #   "synthetic 12S ASV",  # class
  #   "synthetic 12S ASV",  # order
  #   "synthetic 12S ASV",  # family
  #   "synthetic 12S ASV",  # genus
  #   "synthetic 12S ASV",  # species
  #   NA                       # subspecies
  # )
  #
  # names(manual_row) <- taxLevels  # make sure names match columns
  #
  # if (target_asv %in% rownames(assignments)) {
  #   assignments[target_asv, names(manual_row)] <-
  #     manual_row[names(manual_row) %in% colnames(assignments)]
  # } else if (verbose) {
  #   message(sprintf("Manual override skipped: ASV '%s', the synthetic 12S ASV, is not found in assignments.", target_asv))
  # }

  # Reorder rows to match ASV order in qiime_asvtab (nice for downstream joins)
  asv_order <- colnames(qiime_asvtab)
  if (!is.null(asv_order)) {
    keep <- intersect(asv_order, rownames(assignments))
    rest <- setdiff(rownames(assignments), keep)
    assignments <- rbind(assignments[keep, , drop = FALSE],
                         assignments[rest, , drop = FALSE])
  }

  # Coerce to character matrix
  assignments_mat <- as.matrix(assignments)
  storage.mode(assignments_mat) <- "character"

  # Which ASVs have NA at both 'order' and 'family'?
  asvs_na <- rownames(assignments_mat)[
    is.na(assignments_mat[, "order"]) & is.na(assignments_mat[, "family"])
  ]

  # Build the report (samples with reads > 0, and total reads)
  if (length(asvs_na)) {
    samples_with_reads <- colSums(qiime_asvtab[, asvs_na, drop = FALSE] > 0, na.rm = TRUE)
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
    readr::write_csv(unassigned_df, file.path(out_dir, "12S_NAs.csv"))
  } else if (verbose) {
    message("Skipping CSV export: Output directory is not set.")
  }

  # 4) Per-rank % assigned (non-NA across ASVs)
  per_rank_percent_assigned <- colSums(!is.na(assignments_mat)) / nrow(assignments_mat) * 100
  per_rank_percent_assigned <- per_rank_percent_assigned[taxLevels]

  # 5) Optional printing
  if (verbose) {
    cat(sprintf(
      "%.2f%% of ASVs have been assigned, covering %.2f%% of reads in the dataset.\n",
      pct_asvs_assigned, pct_reads_assigned
    ))
    cat("\nBy taxonomic rank: \n")
    cat(paste0(
      names(per_rank_percent_assigned), ": ",
      sprintf("%.2f", per_rank_percent_assigned), "%"
    ), sep = "\n")
    cat("\n")
  }

  # 6) Return
  list(
    assignments = assignments_mat,                  # character matrix (ASV x ranks)
    metrics = list(
      percent_asvs_assigned     = pct_asvs_assigned,
      percent_reads_assigned    = pct_reads_assigned,
      per_rank_percent_assigned = per_rank_percent_assigned
    ),
    NAs_12S = unassigned_df
  )
}
