#' @importFrom magrittr %>%
NULL

utils::globalVariables(c(
  ".data","asv","Species","sample_ID","filtered","denoised","non-chimeric",
  "step","count","project","label","Sample","Abundance","species","genus",
  "family","phylum","plate","well","hash",
  "common_name","taxa","conventional_name","genus_conventional_name",
  # NSE column names referenced unquoted inside dplyr/ggplot2 pipes
  ".","Eigenvalue","PC","PC_index","PCx","PCy","Sample_ID",
  "VarianceExplained","X_reads","Y_reads","adj_ang","ang","eigenvalue",
  "name","observed_eig","rand_mean","type","value","measure",
  # Default argument values that intentionally look up FoodSeq pipeline
  # objects (e.g. `qiime.dir.trnL`) from the caller's environment
  "download_mode","qiime.asvtab.12S","qiime.asvtab.trnL",
  "qiime.dir.12S","qiime.dir.trnL","ref.12S","ref.trnL"
))
