#' @title Join Table Sequences
#'
#' @description This function joins a sequence hash table with a feature
#'   table from QIIME2.
#'
#' @param feature_table A QIIME2-formatted feature table
#' @param sequence_hash A QIIME2-formatted sequence hash table
#' @return A feature table with a sequence hash table joined to it
#' @export
join_table_seqs <- function(feature_table, sequence_hash){
  # feature_table and sequence_hash are the result of reading in QIIME2
  # artifacts with QIIME2R

  # Make dataframe mapping from from hash to ASV
  sequence_hash <-
    data.frame(asv = sequence_hash$data) %>%
    tibble::rownames_to_column(var = 'hash')

  # Substitute hash for ASV in feature table
  feature_table <-
    feature_table$data %>%
    data.frame() %>%
    tibble::rownames_to_column(var = 'hash') %>%
    dplyr::left_join(sequence_hash) %>%
    tibble::column_to_rownames(var = 'asv') %>%
    dplyr::select(-hash)

  # Transform rows and columns and repair plate-well names\
  feature_table <- t(feature_table)

  # Repair names
  row.names(feature_table) <- gsub(pattern = 'X',
                                   replacement = '',
                                   row.names(feature_table))
  row.names(feature_table) <- gsub(pattern = '\\.',
                                   replacement = '-',
                                   row.names(feature_table))

  feature_table
}
