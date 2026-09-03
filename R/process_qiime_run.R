# ---- unzip all *.qzv in a run dir (idempotent) ----
unzip_qzv_safe <- function(run_dir) {
  qzv <- list.files(run_dir, pattern = "^[0-9].*\\.qzv$", full.names = TRUE, ignore.case = TRUE)
  for (zf in qzv) {
    base <- tools::file_path_sans_ext(basename(zf))
    outdir <- file.path(run_dir, base)
    if (!dir.exists(outdir)) {
      message("Unzipping: ", basename(zf))
      tmpdir <- file.path(run_dir, paste0(".tmp_unpack_", base))
      if (dir.exists(tmpdir)) unlink(tmpdir, recursive = TRUE, force = TRUE)
      utils::unzip(zf, exdir = tmpdir)
      if (dir.exists(outdir)) unlink(outdir, recursive = TRUE, force = TRUE)
      file.rename(tmpdir, outdir)
    }
  }
  # special case
  stats_qzv <- file.path(run_dir, "4_denoised-stats.qzv")
  if (file.exists(stats_qzv)) {
    outdir <- file.path(run_dir, "4_denoised-stats")
    if (!dir.exists(outdir)) {
      tmpdir <- file.path(run_dir, ".tmp_unpack_4_denoised-stats")
      if (dir.exists(tmpdir)) unlink(tmpdir, recursive = TRUE, force = TRUE)
      utils::unzip(stats_qzv, exdir = tmpdir)
      if (dir.exists(outdir)) unlink(outdir, recursive = TRUE, force = TRUE)
      file.rename(tmpdir, outdir)
    }
  }
}

# ---- parse per-sample TSVs into a pipeline tracking table ----
build_track_table <- function(run_dir) {
  tsvs <- list.files(
    run_dir,
    pattern = "(per-sample-fastq-counts\\.tsv|metadata\\.tsv)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  if (length(tsvs) == 0) stop("No relevant TSVs found under: ", run_dir)

  canon_map <- c(
    "1_demultiplexed"   = "raw",
    "2_adapter-trimmed" = "adapter_trim",
    "3_primer-trimmed"  = "primer_trim",
    "filtered"          = "filtered",
    "denoised"          = "denoised",
    "merged"            = "merged",
    "non-chimeric"      = "non_chimeric"
  )
  final_order <- c("sample_ID", "raw", "adapter_trim", "primer_trim",
                   "filtered", "denoised", "merged", "non_chimeric")

  read_step_tsv <- function(f) {
    rel <- sub(paste0("^", normalizePath(run_dir, winslash = "/"), "/?"), "",
               normalizePath(f, winslash = "/"))
    parts <- strsplit(rel, "/")[[1]]
    if (length(parts) < 4) return(NULL)
    step_dir <- parts[length(parts) - 3]

    if (step_dir == "4_denoised-table") return(NULL)

    df <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE))

    if (step_dir == "4_denoised-stats") {
      if (!"sample-id" %in% names(df)) return(NULL)
      out <- df[-1, , drop = FALSE] |>
        dplyr::mutate(sample_ID = .data[["sample-id"]]) |>
        dplyr::select(sample_ID, filtered, denoised, merged, `non-chimeric`) |>
        dplyr::mutate(dplyr::across(-sample_ID, ~ readr::parse_number(as.character(.))))
      return(out)
    }

    names(df) <- gsub(" ", "_", names(df), fixed = TRUE)
    sid <- intersect(c("sample_ID","sample-id","SampleID","sampleid"), names(df))
    if (length(sid) == 0) return(NULL)
    names(df)[names(df) == sid[1]] <- "sample_ID"
    if ("reverse_sequence_count" %in% names(df)) df$reverse_sequence_count <- NULL

    num_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], "sample_ID")
    if (length(num_cols) == 0) return(NULL)
    val_col <- num_cols[1]

    step_norm <- canon_map[[step_dir]]
    if (is.null(step_norm)) step_norm <- step_dir

    out <- df[, c("sample_ID", val_col), drop = FALSE]
    names(out)[2] <- step_norm
    out[[step_norm]] <- readr::parse_number(as.character(out[[step_norm]]))
    out
  }

  pieces <- lapply(tsvs, read_step_tsv)
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) stop("No usable TSVs parsed.")

  track <- NULL
  for (piece in pieces) {
    if (is.null(track)) { track <- piece; next }
    merged <- dplyr::full_join(track, piece, by = "sample_ID")

    step_cols <- setdiff(names(merged), "sample_ID")
    base <- sub("\\.(x|y|\\d+)$", "", step_cols)
    for (b in unique(base)) {
      cols <- step_cols[base == b]
      if (length(cols) > 1) {
        nums <- lapply(cols, function(cn) readr::parse_number(as.character(merged[[cn]])))
        merged[[b]] <- dplyr::coalesce(!!!nums)
        merged[cols] <- NULL
      } else if (cols != b) {
        merged[[b]] <- readr::parse_number(as.character(merged[[cols]]))
        merged[cols] <- NULL
      }
    }
    track <- merged
  }

  step_cols <- setdiff(names(track), "sample_ID")
  track[step_cols] <- lapply(track[step_cols], function(x) readr::parse_number(as.character(x)))

  for (old in names(canon_map)) {
    if (old %in% names(track)) names(track)[names(track) == old] <- canon_map[[old]]
  }

  present <- intersect(final_order, names(track))
  track <- track[, unique(c(present, setdiff(names(track), present))), drop = FALSE]
  track
}

# ---- plot and save QC track ----
plot_and_save_track <- function(track_df, run_dir, title = NULL) {
  stopifnot("sample_ID" %in% names(track_df))
  step_cols <- setdiff(names(track_df), "sample_ID")

  df_long <- track_df |>
    dplyr::rename(sample = sample_ID) |>
    tidyr::pivot_longer(dplyr::all_of(step_cols), names_to = "step", values_to = "count") |>
    dplyr::mutate(
      label = dplyr::if_else(sample == "Undetermined", "Undetermined", "Sample"),
      step  = factor(step, levels = step_cols)
    )

  p <- ggplot2::ggplot(df_long, ggplot2::aes(step, count, group = sample)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::facet_wrap(~label, scales = "free_y") +
    ggplot2::labs(
      x = "Pipeline step", y = "Reads",
      title = if (is.null(title)) paste0(basename(run_dir), " run") else title
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(file.path(run_dir, "QC_track-reads-plot.png"),
                  plot = p, width = 10, height = 6, dpi = 300)
  p
}

#' Process QIIME Run
#'
#' This function optionally unzips files, creates a track table, and outputs a plot of read counts depending on inputs.
#'
#' @param marker The marker gene, either "trnL" or "12S"
#' @param run_dir Path to the QIIME run directory. Defaults to looking up
#'   \code{qiime.dir.<marker>} in the calling environment.
#' @param proj The name of the project
#' @param dl_mode If all files from the QIIME run were downloaded or just the necessary ones, either "full" or "minimal"
#' @return A track table and a plot of read counts, optionally
#' @export
process_qiime_run <- function(marker,
                              run_dir = get(paste0("qiime.dir.", marker), envir = parent.frame()),
                              proj = project,
                              dl_mode = download_mode) {
  plot_title <- paste0(proj, ": ", marker, " run")

  # Files we might use
  track_csv <- list.files(
    run_dir,
    pattern = "(?i)track.*pipeline.*\\.csv$",  # case-insensitive, requires both words
    full.names = TRUE
  )[1]
  qzv_needed <- any(file.exists(
    file.path(run_dir, "1_demultiplexed.qzv"),
    file.path(run_dir, "2_adapter-trimmed.qzv"),
    file.path(run_dir, "3_primer-trimmed.qzv"),
    file.path(run_dir, "4_denoised-stats.qzv")
  ))

  # Helper to try reading an existing CSV
  try_read_track_csv <- function() {
    if (file.exists(track_csv)) {
      tr <- tryCatch(readr::read_csv(track_csv, show_col_types = FALSE),
                     error = function(e) NULL) |>
        dplyr::select(-dplyr::starts_with("..."))
      names(tr)[names(tr) == "sample"] <- "sample_ID"
      if (!is.null(tr) && "sample_ID" %in% names(tr)) return(tr)
    }
    NULL
  }

  # Dispatch by mode
  if (dl_mode == "minimal") {
    tr <- try_read_track_csv()
    if (is.null(tr)) {
      message("[", marker, "] minimal mode: no track-pipeline.csv present -> skipping plot.")
      return(list(track = NULL, plot = NULL, note = "No CSV; minimal mode; skipped."))
    }
    p <- plot_and_save_track(tr, run_dir, title = plot_title)
    return(list(track = tr, plot = p))
  }

  if (dl_mode == "full") {
    if (!qzv_needed) {
      stop("[", marker, "] full mode requires the *.qzv files to build track table, but none were found.")
    }
    unzip_qzv_safe(run_dir)
    tr <- build_track_table(run_dir)
    readr::write_csv(tr, track_csv)
    p <- plot_and_save_track(tr, run_dir, title = plot_title)
    return(list(track = tr, plot = p))
  }

  # dl_mode == "auto"
  # 1) Prefer existing CSV
  tr <- try_read_track_csv()
  if (!is.null(tr)) {
    p <- plot_and_save_track(tr, run_dir, title = plot_title)
    return(list(track = tr, plot = p))
  }
  # 2) Else try to build if .qzv are available
  if (qzv_needed) {
    unzip_qzv_safe(run_dir)
    tr <- build_track_table(run_dir)
    readr::write_csv(tr, track_csv)
    p <- plot_and_save_track(tr, run_dir, title = plot_title)
    return(list(track = tr, plot = p))
  }
  # 3) Else we only have .qza (or nothing) → skip gracefully
  message("[", marker, "] auto mode: no CSV and no *.qzv available -> skipping plot.")
  list(track = NULL, plot = NULL, note = "No CSV/.qzv; skipped.")
}
