#' @title Assign Common Names to a Phyloseq Object
#'
#' @description Matches ASV sequences in a phyloseq's tax_table to a
#'   reference CSV of common names, adding `common_name` and `taxa` columns.
#'   Handles multi-match conflicts via superset logic, genus-level
#'   resolution, and smart name merging. For ASVs with no CSV match but a
#'   species-level (or below) taxonomy assignment, common names are
#'   propagated from sibling ASVs that share the same taxonomic
#'   classification, provided all siblings agree on a single unambiguous
#'   name.
#'
#' @param physeq A phyloseq object
#' @param common_names_csv Path to a CSV with columns: `asv`, `conventional_name`,
#'   `taxon`. Optional columns: `genus`, `genus_conventional_name` (for
#'   genus-level conflict resolution).
#' @param report_conflicts If TRUE, print a summary of multi-match conflicts
#' @param report_all_conflicts If TRUE, report all conflicts; if FALSE, only
#'   report unresolved ones
#' @param concatenate_conflicts If TRUE, concatenate unresolved conventional
#'   names using smart merge; if FALSE, use the first match
#' @return The phyloseq object with `common_name` and `taxa` columns added to
#'   its tax_table. Multi-match conflict details are stored as an attribute
#'   `"common_name_conflicts"`.
#' @export
assign_common_names <- function(physeq, common_names_csv,
                                report_conflicts = TRUE,
                                report_all_conflicts = TRUE,
                                concatenate_conflicts = TRUE) {
  # Read the common names CSV
  common_names <- utils::read.csv(common_names_csv, stringsAsFactors = FALSE)

  # Build genus lookup from ALL rows in CSV (only non-NA values)
  genus_key <- NULL
  if ("genus" %in% colnames(common_names) && "genus_conventional_name" %in% colnames(common_names)) {
    valid_genus <- !is.na(common_names$genus) &
      !is.na(common_names$genus_conventional_name) &
      common_names$genus != "" &
      common_names$genus_conventional_name != ""

    # For each genus entry, split by ";" in case multiple genera are listed
    genus_list <- list()
    for (idx in which(valid_genus)) {
      genera <- strsplit(common_names$genus[idx], ";")[[1]]
      conv_names <- strsplit(common_names$genus_conventional_name[idx], ";")[[1]]

      genera <- trimws(genera)
      conv_names <- trimws(conv_names)

      for (j in seq_along(genera)) {
        if (j <= length(conv_names)) {
          genus_list[[genera[j]]] <- conv_names[j]
        }
      }
    }

    if (length(genus_list) > 0) {
      genus_key <- unlist(genus_list)
    }
  }

  # Get the tax_table from phyloseq
  tax_tab <- as.data.frame(phyloseq::tax_table(physeq))

  # Get ASV sequences
  if ("ASV" %in% colnames(tax_tab)) {
    asv_seqs <- tax_tab$ASV
  } else {
    asv_seqs <- rownames(tax_tab)
  }

  # Initialize taxa and common_name columns (in this order so common_name
  # ends up attached right after taxa in the output tax table)
  tax_tab$taxa <- NA_character_
  tax_tab$common_name <- NA_character_

  # Store all conflicts for reporting
  all_conflicts <- list()

  # --- Internal helpers ---

  extract_genera <- function(taxon_string) {
    species_list <- strsplit(taxon_string, "; ")[[1]]
    genera <- sapply(strsplit(species_list, " "), `[`, 1)
    return(unique(genera))
  }

  merge_taxa_strings <- function(taxa_vec) {
    all_taxa <- c()
    for (taxa_string in taxa_vec) {
      if (!is.na(taxa_string) && taxa_string != "") {
        split_taxa <- strsplit(taxa_string, "; ")[[1]]
        split_taxa <- trimws(split_taxa)
        all_taxa <- c(all_taxa, split_taxa)
      }
    }
    unique_taxa <- sort(unique(all_taxa))
    if (length(unique_taxa) > 0) paste(unique_taxa, collapse = "; ") else NA_character_
  }

  singularize <- function(word) {
    word <- tolower(trimws(word))
    special_cases <- c(
      "berries" = "berry", "cherries" = "cherry", "strawberries" = "strawberry",
      "raspberries" = "raspberry", "blackberries" = "blackberry",
      "blueberries" = "blueberry", "cranberries" = "cranberry",
      "mulberries" = "mulberry", "gooseberries" = "gooseberry"
    )
    if (word %in% names(special_cases)) return(special_cases[word])
    if (grepl("sses$", word)) return(sub("sses$", "ss", word))
    if (grepl("shes$", word)) return(sub("shes$", "sh", word))
    if (grepl("ches$", word)) return(sub("ches$", "ch", word))
    if (grepl("xes$", word)) return(sub("xes$", "x", word))
    if (grepl("ies$", word) && !grepl("[aeiou]ies$", word)) return(sub("ies$", "y", word))
    if (grepl("ves$", word)) return(sub("ves$", "f", word))
    if (grepl("s$", word) && !grepl("ss$", word) && !grepl("us$", word)) return(sub("s$", "", word))
    word
  }

  pluralize <- function(word) {
    word_lower <- tolower(trimws(word))
    special_cases <- c("potato" = "potatoes", "tomato" = "tomatoes", "cactus" = "cacti")
    if (word_lower %in% names(special_cases)) {
      if (substring(word, 1, 1) == toupper(substring(word, 1, 1))) {
        return(paste0(toupper(substring(special_cases[word_lower], 1, 1)),
                      substring(special_cases[word_lower], 2)))
      }
      return(special_cases[word_lower])
    }
    if (grepl("s$", word_lower)) return(word)
    if (grepl("y$", word_lower) && !grepl("[aeiou]y$", word_lower)) return(sub("y$", "ies", word))
    if (grepl("(s|ss|sh|ch|x|z)$", word_lower)) return(paste0(word, "es"))
    if (grepl("f$", word_lower)) return(sub("f$", "ves", word))
    if (grepl("fe$", word_lower)) return(sub("fe$", "ves", word))
    paste0(word, "s")
  }

  is_subtype <- function(specific, general) {
    specific <- tolower(trimws(specific))
    general <- tolower(trimws(general))
    specific_sing <- singularize(specific)
    general_sing <- singularize(general)
    if (grepl(paste0("\\s", general_sing, "$"), specific_sing)) return(TRUE)
    if (grepl(paste0("\\s", general, "$"), specific)) return(TRUE)
    FALSE
  }

  consolidate_subtypes <- function(names_vec) {
    if (length(names_vec) <= 1) return(names_vec)
    is_subtype_of <- rep(FALSE, length(names_vec))
    general_terms <- character(length(names_vec))
    for (i in seq_along(names_vec)) {
      for (j in seq_along(names_vec)) {
        if (i != j && is_subtype(names_vec[i], names_vec[j])) {
          is_subtype_of[i] <- TRUE
          general_terms[i] <- names_vec[j]
          break
        }
      }
    }
    if (any(is_subtype_of)) {
      unique_generals <- unique(general_terms[is_subtype_of])
      result <- character(0)
      for (general in unique_generals) result <- c(result, pluralize(general))
      for (i in seq_along(names_vec)) {
        if (!is_subtype_of[i] && !(names_vec[i] %in% general_terms)) {
          result <- c(result, names_vec[i])
        }
      }
      return(unique(result))
    }
    names_vec
  }

  deduplicate_with_plurals <- function(names_vec) {
    if (length(names_vec) == 0) return(character(0))
    singular_map <- sapply(names_vec, singularize)
    unique_singulars <- unique(singular_map)
    result <- character(length(unique_singulars))
    for (i in seq_along(unique_singulars)) {
      first_idx <- which(singular_map == unique_singulars[i])[1]
      result[i] <- names_vec[first_idx]
    }
    result
  }

  merge_like_patterns <- function(names_vec) {
    like_pattern <- "^(.+?)\\s+like\\s+(.+)$"
    like_indices <- grep(like_pattern, names_vec, ignore.case = TRUE)
    if (length(like_indices) == 0) return(names_vec)
    like_groups <- list()
    for (idx in like_indices) {
      m <- regmatches(names_vec[idx], regexec(like_pattern, names_vec[idx], ignore.case = TRUE))[[1]]
      x_part <- trimws(m[2]); y_part <- trimws(m[3])
      if (x_part %in% names(like_groups)) {
        like_groups[[x_part]] <- c(like_groups[[x_part]], y_part)
      } else {
        like_groups[[x_part]] <- y_part
      }
    }
    merged_likes <- character(0)
    for (x in names(like_groups)) {
      y_parts <- unique(like_groups[[x]])
      if (length(y_parts) == 1) {
        merged_likes <- c(merged_likes, paste(x, "like", y_parts))
      } else {
        merged_likes <- c(merged_likes, paste(x, "like",
          paste(y_parts[-length(y_parts)], collapse = ", "), "and", y_parts[length(y_parts)]))
      }
    }
    c(names_vec[-like_indices], merged_likes)
  }

  process_suffix_patterns <- function(names_vec) {
    suffix_patterns <- list(
      other = ",?\\s*and other ([^,]+)$",
      relatives = ",?\\s*and relatives(?: in the ([^,]+))?$"
    )
    suffixes <- list()
    core_names <- character(0)
    for (name in names_vec) {
      found_suffix <- FALSE
      for (pattern_type in names(suffix_patterns)) {
        pattern <- suffix_patterns[[pattern_type]]
        if (grepl(pattern, name, ignore.case = TRUE)) {
          m <- regmatches(name, regexec(pattern, name, ignore.case = TRUE))[[1]]
          suffix_text <- sub("^,?\\s*", "", m[1])
          if (pattern_type %in% names(suffixes)) {
            suffixes[[pattern_type]] <- c(suffixes[[pattern_type]], suffix_text)
          } else {
            suffixes[[pattern_type]] <- suffix_text
          }
          core <- sub(pattern, "", name, ignore.case = TRUE)
          core <- sub(",?\\s*$", "", core)
          core_parts <- unlist(strsplit(core, ";"))
          for (part in core_parts) {
            comma_parts <- unlist(strsplit(part, ","))
            for (comma_part in comma_parts) {
              and_parts <- trimws(unlist(strsplit(comma_part, "\\s+and\\s+")))
              core_names <- c(core_names, and_parts[and_parts != ""])
            }
          }
          found_suffix <- TRUE
          break
        }
      }
      if (!found_suffix) {
        parts <- trimws(unlist(strsplit(name, ",|;| and ")))
        core_names <- c(core_names, parts[parts != ""])
      }
    }
    list(cores = unique(core_names), suffixes = suffixes)
  }

  format_with_and <- function(items, suffix = NULL) {
    if (length(items) == 0) return("")
    has_suffix <- !is.null(suffix) && length(suffix) > 0
    if (length(items) == 1) {
      result <- items[1]
    } else if (length(items) == 2) {
      result <- if (has_suffix) paste(items, collapse = ", ") else paste(items, collapse = " and ")
    } else {
      if (has_suffix) {
        result <- paste(items, collapse = ", ")
      } else {
        result <- paste(paste(items[-length(items)], collapse = ", "),
                        items[length(items)], sep = ", and ")
      }
    }
    if (has_suffix) result <- paste0(result, ", ", suffix[1])
    result
  }

  smart_merge_names <- function(names_vec) {
    if (length(names_vec) == 0) return("")
    names_vec <- merge_like_patterns(names_vec)
    processed <- process_suffix_patterns(names_vec)
    cores <- consolidate_subtypes(processed$cores)
    cores <- deduplicate_with_plurals(cores)
    final_suffix <- NULL
    if (length(processed$suffixes) > 0) {
      if ("other" %in% names(processed$suffixes)) {
        final_suffix <- processed$suffixes$other[1]
      } else if ("relatives" %in% names(processed$suffixes)) {
        final_suffix <- processed$suffixes$relatives[1]
      }
    }
    format_with_and(cores, final_suffix)
  }

  # --- Main matching loop ---

  for (i in seq_along(asv_seqs)) {
    query_seq <- asv_seqs[i]
    matches <- grepl(query_seq, common_names$asv, fixed = TRUE)

    if (sum(matches) == 0) {
      next
    } else if (sum(matches) == 1) {
      tax_tab$common_name[i] <- common_names$conventional_name[matches]
      tax_tab$taxa[i] <- common_names$taxon[matches]
    } else {
      matched_indices <- which(matches)
      matched_taxa <- common_names$taxon[matched_indices]
      matched_conv_names <- common_names$conventional_name[matched_indices]

      resolution_method <- "unresolved"
      assigned_name <- NA
      resolved <- FALSE

      tax_tab$taxa[i] <- merge_taxa_strings(matched_taxa)

      # Superset check
      for (j in seq_along(matched_taxa)) {
        species_j <- strsplit(matched_taxa[j], "; ")[[1]]
        is_superset_of_all <- TRUE
        for (k in seq_along(matched_taxa)) {
          if (j == k) next
          species_k <- strsplit(matched_taxa[k], "; ")[[1]]
          if (!all(species_k %in% species_j)) { is_superset_of_all <- FALSE; break }
        }
        if (is_superset_of_all) {
          tax_tab$common_name[i] <- matched_conv_names[j]
          assigned_name <- matched_conv_names[j]
          resolution_method <- "superset"
          resolved <- TRUE
          break
        }
      }

      # Genus-level resolution
      if (!resolved && !is.null(genus_key)) {
        all_genera <- unique(unlist(lapply(matched_taxa, extract_genera)))
        genus_common_names <- genus_key[all_genera]
        genus_common_names <- genus_common_names[!is.na(genus_common_names)]

        if (length(genus_common_names) == length(all_genera) && length(all_genera) > 0) {
          if (length(unique(genus_common_names)) == 1) {
            tax_tab$common_name[i] <- unique(genus_common_names)[1]
            assigned_name <- unique(genus_common_names)[1]
            resolution_method <- "genus_single"
            resolved <- TRUE
          } else {
            formatted_name <- smart_merge_names(genus_common_names)
            tax_tab$common_name[i] <- formatted_name
            assigned_name <- formatted_name
            resolution_method <- "genus_multiple"
            resolved <- TRUE
          }
        } else if (length(genus_common_names) > 0) {
          genera_without_mapping <- all_genera[!all_genera %in% names(genus_key)]
          conv_names_for_unmapped <- character(0)
          for (idx in matched_indices) {
            row_genera <- extract_genera(common_names$taxon[idx])
            if (any(row_genera %in% genera_without_mapping)) {
              conv_name <- common_names$conventional_name[idx]
              if (!is.na(conv_name) && conv_name != "") {
                conv_names_for_unmapped <- c(conv_names_for_unmapped, conv_name)
              }
            }
          }
          all_names <- c(genus_common_names, conv_names_for_unmapped)
          all_names <- all_names[all_names != ""]
          if (length(all_names) > 0) {
            formatted_name <- smart_merge_names(all_names)
            tax_tab$common_name[i] <- formatted_name
            assigned_name <- formatted_name
            resolution_method <- "genus_partial"
            resolved <- TRUE
          }
        }
      }

      # Fallback: concatenate or first match
      if (!resolved) {
        if (concatenate_conflicts) {
          unique_names <- unique(matched_conv_names[matched_conv_names != "" & !is.na(matched_conv_names)])
          if (length(unique_names) > 0) {
            formatted_name <- smart_merge_names(unique_names)
            tax_tab$common_name[i] <- formatted_name
            assigned_name <- formatted_name
            resolution_method <- "concatenated"
          } else {
            tax_tab$common_name[i] <- matched_conv_names[1]
            assigned_name <- matched_conv_names[1]
            resolution_method <- "first_match_default"
          }
        } else {
          tax_tab$common_name[i] <- matched_conv_names[1]
          assigned_name <- matched_conv_names[1]
          resolution_method <- "first_match_default"
        }
      }

      all_conflicts[[length(all_conflicts) + 1]] <- data.frame(
        asv_index = i,
        asv_seq = substr(query_seq, 1, 50),
        num_matches = length(matched_indices),
        resolution_method = resolution_method,
        assigned_common_name = assigned_name,
        assigned_taxa = tax_tab$taxa[i],
        all_conventional_names = paste(matched_conv_names, collapse = " | "),
        all_taxa = paste(matched_taxa, collapse = " || "),
        genera = paste(unique(unlist(lapply(matched_taxa, extract_genera))), collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Species-level common name propagation ---
  # For ASVs with no CSV match, inherit the common name from sibling ASVs
  # that share the same species-level (and below) taxonomy, provided all
  # donors for that species agree on the same common name (unambiguous).
  all_cols <- colnames(tax_tab)
  species_col_idx <- which(tolower(all_cols) == "species")
  if (length(species_col_idx) == 1) {
    sp_col <- all_cols[species_col_idx]
    subsp_cols <- all_cols[tolower(all_cols) %in% c("subspecies", "varietas", "forma")]
    group_cols <- c(sp_col, subsp_cols)

    # Normalize NA and empty strings to a common sentinel before building keys
    species_key <- apply(tax_tab[, group_cols, drop = FALSE], 1, function(row) {
      row[is.na(row) | !nzchar(row)] <- "<NA>"
      paste(row, collapse = "|")
    })

    has_species <- !is.na(tax_tab[[sp_col]]) & nzchar(tax_tab[[sp_col]])
    needs_name  <- is.na(tax_tab$common_name) & has_species
    has_name    <- !is.na(tax_tab$common_name) & has_species

    if (any(needs_name) && any(has_name)) {
      propagated <- 0L
      skipped_ambiguous <- 0L
      for (i in which(needs_name)) {
        donors <- which(has_name & species_key == species_key[i])
        if (length(donors) > 0) {
          donor_names <- unique(tax_tab$common_name[donors])
          if (length(donor_names) == 1) {
            tax_tab$common_name[i] <- donor_names
            if (is.na(tax_tab$taxa[i])) {
              donor_taxa <- tax_tab$taxa[donors]
              donor_taxa <- donor_taxa[!is.na(donor_taxa)]
              if (length(donor_taxa) > 0) tax_tab$taxa[i] <- donor_taxa[1]
            }
            propagated <- propagated + 1L
          } else {
            skipped_ambiguous <- skipped_ambiguous + 1L
          }
        }
      }
      if (propagated > 0) {
        message(sprintf("Propagated common names to %d ASVs via species-level matching.", propagated))
      }
      if (skipped_ambiguous > 0) {
        message(sprintf(
          "Skipped %d ASVs with ambiguous species-level common names (multiple different names for same species).",
          skipped_ambiguous
        ))
      }
    }
  }

  # Update the tax_table in phyloseq
  phyloseq::tax_table(physeq) <- as.matrix(tax_tab)

  # Report conflicts if requested
  if (report_conflicts && length(all_conflicts) > 0) {
    conflicts_df <- do.call(rbind, all_conflicts)

    if (report_all_conflicts) {
      message(paste("\nFound", nrow(conflicts_df), "total conflicts (multiple matches):"))
      message("\nBreakdown by resolution method:")
      print(table(conflicts_df$resolution_method))
      message("\n")
      print(conflicts_df)
    } else {
      unresolved <- conflicts_df[conflicts_df$resolution_method %in%
                                   c("unresolved", "first_match_default", "concatenated"), ]
      if (nrow(unresolved) > 0) {
        message(paste("\nFound", nrow(unresolved), "unresolved conflicts:"))
        print(unresolved)
      }
    }

    attr(physeq, "common_name_conflicts") <- conflicts_df
  }

  return(physeq)
}
