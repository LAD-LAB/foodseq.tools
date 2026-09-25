# R/project_pca.R
#
# Projects a new (already ASV-harmonized) phyloseq object into an EXISTING
# PCA's fixed space -- the pca_plot() analog of plan_projection()/
# apply_projection(): pca_plot() fits a fresh PCA (like plan_harmonization()/
# apply_harmonization() fit a fresh merge); project_pca() never re-fits
# anything, it just answers "where does this new data land in a PCA that's
# already been fit."
#
# This only works at all because apply_projection() guarantees $ps_projected
# has ps_ref's EXACT ASV/tax_table space -- the standard "out-of-sample
# projection" technique for PCA is to center new data by the ORIGINAL
# fit's own column means and multiply by its ORIGINAL rotation matrix
# (exactly what stats::predict.prcomp() does), which is only valid when
# the new data's columns are the identical feature set the PCA was fit on.
# Loadings (the biplot arrows) are entirely a property of that original
# fit -- never recomputed here, always passed through from pca_ref.
#
# Deliberately a separate function from pca_plot(), not a mode-switch
# argument on it -- fitting fresh and projecting onto an existing fit are
# different enough operations (very different data flow: one runs
# prcomp() from scratch, the other explicitly never does) that overloading
# one function invites bugs, mirroring why plan_harmonization()/
# apply_harmonization() and plan_projection()/apply_projection() are kept
# separate too. The shared biplot-rendering logic (arrows, quadrant-aware
# text labels) lives in .pca_biplot_layer() (see pca_plot.R) so neither
# function duplicates it.
#
# `together = TRUE`'s $pca.df automatically combines ps_ref's and
# ps_query's own sample metadata (decided explicitly, see Description) --
# usable via `colorVar` exactly as in `together = FALSE`/pca_plot() itself.
# Shared column names merge directly whenever both sides agree on type; a
# genuine type conflict (confirmed on real data: two independently-
# collected batches both had a `seq_date` column, but one held integers
# and the other characters) is never silently coerced or guessed at by
# dplyr::bind_rows() -- each side's conflicting column is instead kept
# separate, renamed `<col>.ref`/`<col>.query`, with a message naming which.
# The built-in `dataset` ("query"/"reference") label is reserved the same
# way `name` already is -- an existing same-named column on either side is
# renamed `dataset.x` first, never silently overwritten.

#' @title Project New Data into an Existing PCA
#'
#' @description Projects a new, already ASV-harmonized phyloseq object
#'   into a PCA that was already fit by \code{\link{pca_plot}} on a
#'   reference dataset -- the \code{pca_plot()} analog of
#'   \code{\link{apply_projection}}. Never re-fits anything: \code{ps_query}
#'   is centered by the reference PCA's own column means and projected onto
#'   its own rotation matrix (via \code{stats::predict.prcomp()}), the
#'   standard technique for placing new samples into an existing PCA's
#'   fixed space. The biplot's loading arrows are the reference's own,
#'   passed through unchanged -- they never differ between the reference's
#'   own biplot and a query projected into it, since both are describing
#'   the same fixed PC-space basis.
#'
#'   This only produces a meaningful result when \code{ps_query} is
#'   already in the reference's exact ASV space -- e.g. \code{$ps_projected}
#'   from \code{\link{apply_projection}}, filtered the same way the
#'   reference was (see Details) -- checked explicitly before anything else
#'   runs.
#'
#' @param ps_ref The reference's CLR-transformed (and filtered) phyloseq
#'   object -- the same one that would be passed as \code{pca_plot()}'s
#'   \code{ps}. Used to label the loading arrows by taxon, and (when
#'   \code{pca_ref} isn't supplied) to fit the PCA itself.
#' @param ps_query New, CLR-transformed and filtered phyloseq object to
#'   project -- must have exactly the reference PCA's fitted ASV set (same
#'   taxa, order doesn't matter). Typically \code{apply_projection()}'s
#'   \code{$ps_projected}, filtered identically to how the reference was
#'   (see Details).
#' @param pca_ref Optional -- an already-fit PCA to project onto, skipping
#'   re-fitting. Accepts \code{\link{pca_plot}}'s full return value, or (for
#'   a much smaller object worth actually persisting long-term) just
#'   \code{list(pca.output = <the prcomp object>)} -- nothing else is ever
#'   read from it. If \code{NULL} (default), the PCA is fit fresh from
#'   \code{ps_ref} internally (equivalent to \code{pca_plot(ps_ref)}, and
#'   just as expensive) -- worth avoiding by supplying \code{pca_ref} when
#'   projecting many query batches onto the same reference repeatedly.
#' @param together If \code{FALSE} (default), the biplot shows only
#'   \code{ps_query}'s projected points, and \code{$pca.df} carries
#'   \code{ps_query}'s own full sample data. If \code{TRUE}, the
#'   reference's own points are shown alongside, and \code{$pca.df}
#'   automatically combines \strong{both} objects' own sample data (see
#'   \code{$pca.df} below for how shared-but-conflicting-type columns are
#'   handled) -- \code{colorVar} et al. below work identically either way,
#'   referencing whichever object's (or both objects' shared) column
#'   directly, exactly as \code{pca_plot()}'s own do. Points default to
#'   coloring by a built-in \code{dataset} ("query"/"reference") label
#'   only when \code{together = TRUE} \emph{and} \code{colorVar} isn't set.
#' @param nTaxa Number of taxa (loading arrows) to display. Default \code{10}.
#' @param colorVar Optional column name to color points by -- identical
#'   meaning to \code{pca_plot()}'s own \code{colorVar}, and usable
#'   regardless of \code{together} (referencing \code{ps_query}'s own
#'   column, or, when \code{together = TRUE}, either object's own or a
#'   column shared by both -- see \code{$pca.df} below). Takes priority
#'   over the built-in \code{dataset} label when \code{together = TRUE}.
#' @param colorName What to display \code{colorVar} as in the legend.
#' @param customColors When \code{colorVar} is set, an optional named
#'   vector of colors for its own levels (as in \code{pca_plot()}),
#'   regardless of \code{together}. When \code{colorVar} isn't set and
#'   \code{together = TRUE}, an optional named length-2 vector
#'   \code{c(query = ..., reference = ...)} overriding the two built-in
#'   dataset colors instead. Ignored (with a warning) when \code{colorVar}
#'   isn't set and \code{together = FALSE} (nothing to color by).
#' @param customGradient When \code{colorVar} is set, an optional data
#'   frame of \code{low}/\code{mid}/\code{high} colors for a continuous
#'   \code{colorVar} gradient -- identical meaning to \code{pca_plot()}'s
#'   own \code{customGradient}, regardless of \code{together}. Ignored
#'   (with a warning) when \code{colorVar} isn't set.
#' @param mid For \code{customGradient} -- midpoint to \code{"median"},
#'   \code{"mean"}, or \code{"middle"} of range. Default \code{"mean"}.
#' @param xPC Principal component for the x-axis. Default \code{1}.
#' @param yPC Principal component for the y-axis. Default \code{2}.
#' @param ellipse Optional centroid ellipse(s) -- one per \code{colorVar}
#'   level when set, one per dataset when \code{together = TRUE} and
#'   \code{colorVar} isn't set, one overall ellipse otherwise.
#' @param bplab Optional tax_table column name (on \code{ps_ref}) to
#'   prioritize for biplot arrow labels -- identical meaning to
#'   \code{pca_plot()}'s own \code{bplab}.
#'
#' @details
#' \strong{Why \code{ps_query} must match exactly}: the projection technique
#' this function uses is only valid when \code{ps_query}'s columns are the
#' identical feature set (same ASVs) the reference PCA was fit on -- a
#' different or partially-overlapping ASV set would silently produce
#' meaningless coordinates.
#'
#' \strong{Getting there}: don't try to reproduce this by filtering
#' \code{ps_query} to match the reference after the fact -- filter the
#' \emph{reference} for taxonomy (removing unassigned ASVs, known
#' controls, host reads, etc.) \emph{before} running
#' \code{\link{plan_harmonization}}/\code{\link{apply_harmonization}} on
#' it, not after. Done that way, \code{\link{apply_projection}}'s
#' \code{$ps_projected} is \emph{already} in the filtered reference's exact
#' final ASV space by construction (it always carries the reference's own
#' final tax_table verbatim) -- pass that straight in as \code{ps_query}
#' with no separate filtering step of its own. Filtering the reference only
#' after harmonizing (or filtering \code{ps_query} independently instead of
#' relying on \code{apply_projection()}'s own output) risks filtering a
#' different ASV universe than the one the reference PCA was actually fit
#' on, even though the filtering criteria themselves stay perfectly
#' consistent -- exactly the mismatch this function's validation above
#' exists to catch.
#'
#' @return A named list: \code{pca.df} (a plain data frame of PC scores plus
#'   \code{name} and \code{dataset}, plus sample metadata -- when
#'   \code{together = FALSE}, \code{ps_query}'s own full sample data,
#'   joined in since there's no second dataset to conflict with; when
#'   \code{together = TRUE}, \strong{both} objects' own sample data,
#'   automatically combined. A shared column name merges directly as long
#'   as both sides agree on type; when they genuinely don't -- two
#'   independently collected batches routinely reuse a column name, e.g.
#'   \code{seq_date}, for incompatibly-typed data \code{dplyr::bind_rows()}
#'   can't silently reconcile -- that one column is kept separate instead,
#'   renamed \code{seq_date.ref}/\code{seq_date.query}, with a
#'   \code{message()} naming which columns this happened to; every other
#'   shared or side-specific column merges normally, \code{NA} on
#'   whichever side didn't have it), \code{pca.biplot}
#'   (the default biplot), and \code{loadings} (the reference's own loading
#'   matrix, passed through unchanged from
#'   \code{pca_ref}). Does \emph{not} return \code{pca.output}/
#'   \code{scree.table}/\code{scree.plot} -- nothing new is fit here, so
#'   those would just duplicate \code{pca_ref}'s own, unchanged; use
#'   \code{pca_ref} directly for those.
#'
#' @examples
#' \dontrun{
#' # filtered_ref was filtered for taxonomy BEFORE ref_harmonized was built
#' # from it (see Details) -- ref_harmonized's final ASV space, and so
#' # $ps_projected below, already reflect that filtering.
#' ref_clr <- microbiome::transform(filtered_ref, "clr")
#'
#' plan <- plan_projection(ref_harmonized, ps)
#' # ... resolve any pairs plan$decisions flags, if needed ...
#' result <- apply_projection(plan, ps)          # ... resolved to "complete" ...
#' query_clr <- microbiome::transform(result$ps_projected, "clr")
#'
#' # Auto-fits the PCA on ref_clr internally -- equivalent to pca_plot(ref_clr).
#' proj <- project_pca(ref_clr, query_clr, colorVar = "study")
#' proj$pca.biplot
#'
#' # Reference and query together, colored by the built-in dataset label
#' # (query/reference) since colorVar isn't set:
#' proj2 <- project_pca(ref_clr, query_clr, together = TRUE)
#' proj2$pca.biplot
#'
#' # Reference and query together, colored instead by a column shared by
#' # both objects' own sample data (auto-combined -- see $pca.df):
#' proj2b <- project_pca(ref_clr, query_clr, together = TRUE, colorVar = "study")
#' proj2b$pca.biplot
#'
#' # Projecting many query batches onto the same reference repeatedly --
#' # fit once, reuse the (much smaller) fit rather than re-fitting each time.
#' fit <- list(pca.output = pca_plot(ref_clr)$pca.output)
#' proj3 <- project_pca(ref_clr, query_clr, pca_ref = fit)
#' }
#'
#' @export
project_pca <- function(ps_ref, ps_query, pca_ref = NULL, together = FALSE,
                        nTaxa = 10, colorVar = NULL, colorName = NULL,
                        customColors = NULL, customGradient = NULL, mid = "mean",
                        xPC = 1, yPC = 2, ellipse = FALSE, bplab = NULL) {

  stopifnot(inherits(ps_ref, "phyloseq"), inherits(ps_query, "phyloseq"))

  # -- Auto-fit when pca_ref isn't supplied -- equivalent to pca_plot(ps_ref)
  # and just as expensive; only worth avoiding by passing pca_ref yourself
  # when projecting many query batches onto the same reference repeatedly.
  # Only $pca.output is ever read from pca_ref -- nothing else (not
  # $pca.df, not the biplot/scree objects) -- so a caller wanting a small,
  # long-term-storable object can pass list(pca.output = <the prcomp
  # object>) instead of pca_plot()'s full return.
  if (is.null(pca_ref)) {
    pca_ref <- list(pca.output = pca_plot(ps_ref)$pca.output)
  }
  stopifnot(is.list(pca_ref), !is.null(pca_ref$pca.output))
  pca_obj <- pca_ref$pca.output
  if (!inherits(pca_obj, "prcomp")) {
    stop("project_pca: `pca_ref$pca.output` must be a prcomp object -- ",
         "pass the list returned by pca_plot() (or one built the same way), ",
         "not something else.")
  }

  if (is.null(colorVar)) {
    if (!is.null(customGradient)) {
      warning("project_pca: `customGradient` is ignored when `colorVar` isn't set -- ",
              if (together)
                paste0("the built-in `dataset` label together = TRUE colors by default is ",
                       "discrete, not continuous (use `customColors` instead to override its ",
                       "two colors).")
              else
                "there's nothing yet to color by.")
    }
    if (!together && !is.null(customColors)) {
      warning("project_pca: `customColors` is ignored when `colorVar` isn't set and ",
              "together = FALSE (there's nothing yet to color by).")
    }
  }

  # Normalize orientation to samples-as-rows / taxa-as-columns, exactly as
  # pca_plot() does -- otherwise a taxa-as-rows ps_query would silently
  # transpose into nonsense once matrix-multiplied against pca_obj$rotation.
  if (phyloseq::taxa_are_rows(ps_query)) {
    otu_samples_rows <- t(methods::as(phyloseq::otu_table(ps_query), "matrix"))
    phyloseq::otu_table(ps_query) <- phyloseq::otu_table(otu_samples_rows, taxa_are_rows = FALSE)
  }

  # -- The precondition this whole function depends on: ps_query's ASVs
  # must be EXACTLY the set the reference PCA was fit on (order doesn't
  # matter -- stats::predict.prcomp() matches and reorders by name itself).
  # Deliberately stricter than predict.prcomp()'s own check, which only
  # requires the fitted ASVs to be PRESENT in newdata and silently ignores
  # anything extra -- a genuinely different (even if superset) ASV set
  # signals ps_query was never actually run through apply_projection()
  # against this reference, or was filtered differently, either of which would
  # otherwise silently produce a technically-computable but meaningless
  # projection. No override, consistent with every other safety-vs-
  # convenience fork already made across this package.
  ref_taxa   <- rownames(pca_obj$rotation)
  query_taxa <- phyloseq::taxa_names(ps_query)
  missing_taxa <- setdiff(ref_taxa, query_taxa)
  extra_taxa   <- setdiff(query_taxa, ref_taxa)
  if (length(missing_taxa) > 0 || length(extra_taxa) > 0) {
    stop("project_pca: ps_query's ASVs do not exactly match the reference ",
         "PCA's fitted feature space -- ",
         if (length(missing_taxa) > 0)
           paste0(length(missing_taxa), " missing (e.g. ",
                  paste(utils::head(missing_taxa, 3), collapse = ", "),
                  if (length(missing_taxa) > 3) ", ..." else "", "). ")
         else "",
         if (length(extra_taxa) > 0)
           paste0(length(extra_taxa), " extra (e.g. ",
                  paste(utils::head(extra_taxa, 3), collapse = ", "),
                  if (length(extra_taxa) > 3) ", ..." else "", "). ")
         else "",
         "ps_query must be apply_projection()'s $ps_projected against this same ",
         "reference, filtered identically to how the reference was before ",
         "its PCA was fit (see ?project_pca Details) -- there is no override.")
  }

  # Reserve `name` (the join key below) and `dataset` (added below,
  # regardless of `together` -- see the dplyr::mutate() a few lines down)
  # the same way pca_plot() already reserves `name` -- an existing
  # same-named column on ps_query is renamed rather than silently
  # overwritten by dplyr::mutate(dataset = "query") later.
  ps_query_samdf_raw <- data.frame(ps_query@sam_data)
  reserved_hit_query <- intersect(c("name", "dataset"), colnames(ps_query_samdf_raw))
  if (length(reserved_hit_query) > 0) {
    phyloseq::sample_data(ps_query) <- ps_query_samdf_raw %>%
      dplyr::rename_with(~ paste0(.x, ".x"), dplyr::all_of(reserved_hit_query))
  }

  samdf_query <- if (!is.null(phyloseq::sample_data(ps_query, errorIfNULL = FALSE))) {
    data.frame(ps_query@sam_data) %>% tibble::rownames_to_column(var = "name")
  } else {
    data.frame(name = phyloseq::sample_names(ps_query), stringsAsFactors = FALSE)
  }

  # -- The actual projection: center by the REFERENCE's own column means
  # and multiply by the REFERENCE's own rotation -- stats::predict.prcomp()
  # matches columns by name (reordering as needed) and does exactly this.
  otu_mat       <- methods::as(phyloseq::otu_table(ps_query), "matrix")
  query_scores  <- stats::predict(pca_obj, newdata = otu_mat)

  eigs         <- pca_obj$sdev^2
  varExplained <- 100 * eigs / sum(eigs)
  names(varExplained) <- paste0('PC', seq_along(varExplained))
  ve.xPC <- as.character(round(varExplained[paste0('PC', xPC)], 3))
  ve.yPC <- as.character(round(varExplained[paste0('PC', yPC)], 3))

  pca.df.query <- data.frame(query_scores) %>%
    tibble::rownames_to_column(var = "name") %>%
    dplyr::left_join(samdf_query, by = "name") %>%
    dplyr::mutate(dataset = "query")

  if (together) {
    # Same `name`/`dataset` reservation as ps_query above.
    ps_ref_samdf_raw  <- data.frame(ps_ref@sam_data)
    reserved_hit_ref  <- intersect(c("name", "dataset"), colnames(ps_ref_samdf_raw))
    if (length(reserved_hit_ref) > 0) {
      phyloseq::sample_data(ps_ref) <- ps_ref_samdf_raw %>%
        dplyr::rename_with(~ paste0(.x, ".x"), dplyr::all_of(reserved_hit_ref))
    }
    samdf_ref <- if (!is.null(phyloseq::sample_data(ps_ref, errorIfNULL = FALSE))) {
      data.frame(ps_ref@sam_data) %>% tibble::rownames_to_column(var = "name")
    } else {
      data.frame(name = phyloseq::sample_names(ps_ref), stringsAsFactors = FALSE)
    }

    pc_cols <- grep("^PC[0-9]+$", colnames(pca_ref$pca.output$x), value = TRUE)
    pca.df.ref <- data.frame(pca_ref$pca.output$x[, pc_cols, drop = FALSE]) %>%
      tibble::rownames_to_column(var = "name") %>%
      dplyr::left_join(samdf_ref, by = "name") %>%
      dplyr::mutate(dataset = "reference")

    # Auto-combine, but never silently coerce a genuine type conflict (see
    # Description/file header) -- a shared column name where ps_ref's and
    # ps_query's own values disagree in class is kept separate instead,
    # each side renamed with its own suffix, rather than left to
    # dplyr::bind_rows() to guess about.
    shared <- setdiff(intersect(colnames(pca.df.ref), colnames(pca.df.query)),
                      c("name", "dataset", pc_cols))
    conflicting <- shared[vapply(shared, function(col)
      !identical(class(pca.df.ref[[col]]), class(pca.df.query[[col]])), logical(1))]
    if (length(conflicting) > 0) {
      message("project_pca: ", length(conflicting), " shared sample-metadata column(s) have ",
              "conflicting types between ps_ref and ps_query -- kept separate as ",
              "'<col>.ref'/'<col>.query' rather than merged: ",
              paste(conflicting, collapse = ", "), ".")
      for (col in conflicting) {
        colnames(pca.df.ref)[colnames(pca.df.ref) == col]     <- paste0(col, ".ref")
        colnames(pca.df.query)[colnames(pca.df.query) == col] <- paste0(col, ".query")
      }
    }

    dupe_names <- intersect(pca.df.query$name, pca.df.ref$name)
    if (length(dupe_names) > 0) {
      message("project_pca: ", length(dupe_names), " sample name(s) appear in ",
              "both ps_query and the reference (e.g. '", dupe_names[1], "') -- ",
              "both are kept (disambiguated by the `dataset` column), but note ",
              "this before joining $pca.df against external metadata by name alone.")
    }
    pca.df <- dplyr::bind_rows(pca.df.ref, pca.df.query)
  } else {
    pca.df <- pca.df.query
  }

  limit <- max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))])) +
    0.05 * max(abs(pca.df[, c(paste0('PC', xPC), paste0('PC', yPC))]))

  if (!is.null(colorVar)) {
    # colorVar, when set, always wins -- regardless of together -- exactly
    # like pca_plot() itself; the built-in `dataset` label is only the
    # together = TRUE default when nothing more specific was asked for.
    pca.plot <- ggplot2::ggplot(pca.df, ggplot2::aes(x = .data[[paste0('PC', xPC)]],
                                                      y = .data[[paste0('PC', yPC)]],
                                                      color = .data[[colorVar]]))
  } else if (together) {
    pca.plot <- ggplot2::ggplot(pca.df, ggplot2::aes(x = .data[[paste0('PC', xPC)]],
                                                      y = .data[[paste0('PC', yPC)]],
                                                      color = .data[["dataset"]]))
  } else {
    pca.plot <- ggplot2::ggplot(pca.df, ggplot2::aes(x = .data[[paste0('PC', xPC)]],
                                                      y = .data[[paste0('PC', yPC)]]))
  }

  pca.plot <- pca.plot +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste0('PC', xPC, ' (', ve.xPC, '%)'),
         y = paste0('PC', yPC, ' (', ve.yPC, '%)')) +
    ggplot2::xlim(-limit, limit) + ggplot2::ylim(-limit, limit) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 1, color = 'black'),
          axis.ticks = ggplot2::element_line(color = 'black'),
          axis.title = ggplot2::element_text(size = 14, face = 'bold', color = 'black'))

  if (!is.null(colorVar)) {
    # customColors means "one color per colorVar level" whenever colorVar
    # is set -- identical to pca_plot(), regardless of together.
    if (!is.null(customColors)) {
      pca.plot <- pca.plot + ggplot2::scale_color_manual(values = customColors)
    }
  } else if (together && !is.null(customColors)) {
    # colorVar unset: customColors instead overrides the two built-in
    # dataset colors (together = TRUE's only meaning for it).
    if (!setequal(names(customColors), c("query", "reference"))) {
      stop("project_pca: `customColors` must be a named vector with exactly ",
           "names 'query' and 'reference' (got: ",
           paste(names(customColors), collapse = ", "), ").")
    }
    pca.plot <- pca.plot + ggplot2::scale_color_manual(values = customColors)
  }

  # Custom color gradient -- meaningful whenever colorVar is set (a
  # continuous column to build a gradient over), regardless of together;
  # already warned as ignored above when colorVar isn't set.
  if (!is.null(colorVar) && !is.null(customGradient)) {
    midpoint <- switch(mid,
                       "middle" = mean(range(pca.df[[colorVar]], na.rm = TRUE)),
                       "median" = stats::median(pca.df[[colorVar]], na.rm = TRUE),
                       "mean" = mean(pca.df[[colorVar]], na.rm = TRUE),
                       stop("'mid' must be one of 'middle', 'median', or 'mean'")
    )

    pca.plot <- pca.plot + ggplot2::scale_color_gradient2(low = paste0(customGradient$low),
                                                 mid = paste0(customGradient$mid),
                                                 high = paste0(customGradient$high),
                                                 midpoint = midpoint)
  }

  if (ellipse) {
    if (!is.null(colorVar)) {
      pca.plot <- pca.plot + ggplot2::stat_ellipse(level = 0.95, ggplot2::aes(group = .data[[colorVar]]), linetype = "dashed")
    } else if (together) {
      pca.plot <- pca.plot + ggplot2::stat_ellipse(level = 0.95, ggplot2::aes(group = .data[["dataset"]]), linetype = "dashed")
    } else {
      pca.plot <- pca.plot + ggplot2::stat_ellipse(level = 0.95, linetype = "dashed")
    }
  }

  biplot_out <- .pca_biplot_layer(pca.plot, pca_obj$rotation, pca_obj$sdev, ps_ref,
                                  xPC, yPC, nTaxa, bplab,
                                  colorName = if (is.null(colorVar) && together) "Dataset" else colorName)

  list(pca.df = pca.df, pca.biplot = biplot_out$pca.biplot, loadings = biplot_out$loadings)
}
