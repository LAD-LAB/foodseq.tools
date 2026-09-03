#' @title Truncate to Folder
#'
#' @description This function allows you to truncate a given path to a
#'   folder of interest.
#'
#' @param path A file path
#' @param folder The folder of interest
#' @return A file path truncated to the folder of interest
#' @export
truncate_to_folder <- function(path, folder) {
  parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
  idx <- match(folder, parts)
  if (is.na(idx)) stop(sprintf("Folder '%s' not found in path: %s", folder, path))
  paste0("/", paste(parts[2:idx], collapse = "/"))
}
