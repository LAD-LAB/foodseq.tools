#' @title Find G-to-A / C-to-T Parent-Child ASV Pairs
#'
#' @description Identifies pairs of ASVs consistent with G→A or C→T
#'   nucleotide transitions, which are common DADA2 denoising artefacts
#'   (chimeras or damage-induced variants). A pair is reported when:
#'   \enumerate{
#'     \item The sequences differ at 1–\code{max_subs} positions.
#'     \item Every differing position is a G→A or C→T change (parent→child).
#'     \item The parent read count is at least \code{min_ratio} × the child.
#'   }
#'
#' @details
#' \strong{Modes:}
#' \itemize{
#'   \item \code{"by_length"} (default): processes each group of equal-length
#'     sequences independently. Sequences of different lengths are never
#'     compared.
#'   \item \code{"trim_to_min"}: all sequences are left-trimmed to the minimum
#'     observed length before comparison. Use with caution.
#' }
#'
#' Sequences are retrieved from \code{refseq(ps)} if present; otherwise
#' taxa names are assumed to be the sequences themselves (standard DADA2
#' output).
#'
#' @param ps Phyloseq object.
#' @param min_ratio Numeric; minimum parent/child read count ratio for a pair
#'   to be reported. Default \code{10}.
#' @param max_subs Integer; maximum number of substitutions allowed (1 or 2).
#'   Default \code{2}.
#' @param ignore_N Logical; skip sequences containing non-ACGT characters
#'   (e.g., ambiguous bases). Default \code{TRUE}.
#' @param mode Character; \code{"by_length"} (default) or
#'   \code{"trim_to_min"}. See Details.
#'
#' @return A data frame with one row per identified pair, sorted by descending
#'   ratio, with columns:
#'   \describe{
#'     \item{parent_id}{Parent ASV identifier.}
#'     \item{child_id}{Child ASV identifier.}
#'     \item{parent_seq}{Parent sequence.}
#'     \item{child_seq}{Child sequence.}
#'     \item{n_subs}{Number of substitutions (1 or 2).}
#'     \item{positions}{Comma-separated substitution positions.}
#'     \item{parent_count}{Total reads of the parent ASV.}
#'     \item{child_count}{Total reads of the child ASV.}
#'     \item{ratio}{parent_count / child_count.}
#'   }
#'   Returns an empty data frame if no pairs are found.
#'
#' @export
find_g2a_c2t_pairs <- function(ps,
                                min_ratio = 10,
                                max_subs  = 2,
                                ignore_N  = TRUE,
                                mode      = c("by_length", "trim_to_min")) {

  mode <- match.arg(mode)
  stopifnot(max_subs %in% c(1L, 2L))

  taxa_ids <- phyloseq::taxa_names(ps)

  # Prefer refseq() slot; otherwise assume taxa names are sequences
  rs <- tryCatch(phyloseq::refseq(ps), error = function(e) NULL)
  if (!is.null(rs)) {
    seqs       <- as.character(rs)
    names(seqs) <- names(rs)
    seqs       <- seqs[taxa_ids]
  } else {
    seqs       <- taxa_ids
    names(seqs) <- taxa_ids
  }

  counts       <- phyloseq::taxa_sums(ps)[taxa_ids]

  if (ignore_N) {
    keep <- !grepl("[^ACGT]", seqs, perl = TRUE)
    if (!all(keep)) {
      message("Skipping ", sum(!keep), " sequence(s) with non-ACGT characters.")
      seqs     <- seqs[keep]
      counts   <- counts[keep]
      taxa_ids <- taxa_ids[keep]
    }
  }

  if (!length(seqs)) return(data.frame())

  if (mode == "trim_to_min") {
    Lmin <- min(nchar(seqs))
    if (length(unique(nchar(seqs))) > 1) {
      message("Trimming all sequences to minimum length = ", Lmin, " (left-anchored).")
      seqs <- substr(seqs, 1L, Lmin)
    }
    length_groups <- list(all = rep(TRUE, length(seqs)))
  } else {
    lens          <- nchar(seqs)
    ug            <- sort(unique(lens))
    length_groups <- lapply(ug, function(L) lens == L)
    names(length_groups) <- paste0("L", ug)
  }

  process_group <- function(mask) {
    if (!any(mask)) return(NULL)
    g_seqs   <- seqs[mask]
    g_ids    <- names(g_seqs)
    g_counts <- counts[mask]
    n        <- length(g_seqs)
    if (n < 2) return(NULL)

    # Hash map: sequence -> index within group
    idx_by_seq <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_len(n)) assign(g_seqs[i], i, envir = idx_by_seq)

    set_at <- function(s, pos, ch) {
      paste0(substr(s, 1, pos - 1), ch, substr(s, pos + 1, nchar(s)))
    }

    is_valid_transition <- function(parent, child) {
      p    <- strsplit(parent, "", useBytes = TRUE)[[1]]
      ch   <- strsplit(child,  "", useBytes = TRUE)[[1]]
      diff <- which(p != ch)
      if (!length(diff) || length(diff) > max_subs) return(FALSE)
      all((p[diff] == "G" & ch[diff] == "A") | (p[diff] == "C" & ch[diff] == "T"))
    }

    res_parent <- integer()
    res_child  <- integer()
    res_nsubs  <- integer()
    res_pos    <- character()

    for (i in seq_len(n)) {
      parent_seq <- g_seqs[i]
      parent_ct  <- g_counts[i]
      if (parent_ct == 0) next
      chars  <- strsplit(parent_seq, "", useBytes = TRUE)[[1]]
      gc_pos <- which(chars %in% c("G", "C"))
      if (!length(gc_pos)) next

      # 1-substitution candidates
      for (p1 in gc_pos) {
        mut1  <- if (chars[p1] == "G") "A" else "T"
        cand1 <- set_at(parent_seq, p1, mut1)
        j     <- get0(cand1, envir = idx_by_seq, ifnotfound = NULL)
        if (!is.null(j)) {
          child_ct <- g_counts[j]
          if (child_ct > 0 && parent_ct >= min_ratio * child_ct &&
              is_valid_transition(parent_seq, cand1)) {
            res_parent <- c(res_parent, i); res_child <- c(res_child, j)
            res_nsubs  <- c(res_nsubs,  1L); res_pos  <- c(res_pos,  as.character(p1))
          }
        }
      }

      # 2-substitution candidates
      if (max_subs == 2 && length(gc_pos) >= 2) {
        m <- length(gc_pos)
        for (a in seq_len(m - 1)) {
          p1   <- gc_pos[a]; mut1 <- if (chars[p1] == "G") "A" else "T"
          seq1 <- set_at(parent_seq, p1, mut1)
          for (b in (a + 1):m) {
            p2    <- gc_pos[b]; mut2 <- if (chars[p2] == "G") "A" else "T"
            cand2 <- set_at(seq1, p2, mut2)
            j     <- get0(cand2, envir = idx_by_seq, ifnotfound = NULL)
            if (!is.null(j)) {
              child_ct <- g_counts[j]
              if (child_ct > 0 && parent_ct >= min_ratio * child_ct &&
                  is_valid_transition(parent_seq, cand2)) {
                res_parent <- c(res_parent, i); res_child <- c(res_child, j)
                res_nsubs  <- c(res_nsubs,  2L)
                res_pos    <- c(res_pos,  paste0(p1, ",", p2))
              }
            }
          }
        }
      }
    }

    if (!length(res_parent)) return(NULL)

    data.frame(
      parent_id    = g_ids[res_parent],
      child_id     = g_ids[res_child],
      parent_seq   = g_seqs[res_parent],
      child_seq    = g_seqs[res_child],
      n_subs       = res_nsubs,
      positions    = res_pos,
      parent_count = as.numeric(g_counts[res_parent]),
      child_count  = as.numeric(g_counts[res_child]),
      ratio        = as.numeric(g_counts[res_parent]) /
                     pmax(as.numeric(g_counts[res_child]), 1),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, lapply(length_groups, process_group))
  if (is.null(out) || !nrow(out)) return(data.frame())

  out <- unique(out[order(-out$ratio, out$parent_id, out$child_id), ])
  rownames(out) <- NULL
  out
}
