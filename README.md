# foodseq.tools

`foodseq.tools` is an R package of helper functions for FoodSeq metabarcoding
workflows: taking raw QIIME 2 output through taxonomy assignment, quality
control, cross-batch ASV harmonization, and ordination. It is developed and
maintained by the David Lab at Duke University.

FoodSeq is a DNA metabarcoding approach for characterizing diet from stool,
food, or environmental samples, typically combining a vertebrate marker
(12SV5) and a plant marker (trnL). `foodseq.tools` provides the R-side
building blocks for that pipeline once QIIME 2 has produced ASVs: assigning
and harmonizing taxonomy, checking quality, reconciling ASVs across
sequencing batches, and running ordinations on the result.

## Who this is for

This package is built for David Lab members and assumes lab-specific
conventions (object naming, file layout, reference files) throughout. Several
functions will not work as intended outside that context.

## Documentation

This README covers installation and a catalog of what's available. For
step-by-step guidance on actually running a FoodSeq analysis with these
functions, see the FoodSeq handbook:
[LAD-LAB/lad-lab.github.io](https://github.com/LAD-LAB/lad-lab.github.io).

## Installation

`foodseq.tools` depends on several Bioconductor packages that aren't
available from CRAN, so install those first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "dada2", "phyloseq", "RcppParallel", "ShortRead"))
```

Then install `foodseq.tools` itself from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("LAD-LAB/foodseq.tools")
```

## What's included

### Sequence processing & track tables

| Function | Description |
|---|---|
| `process_qiime_run()` | Unzips QIIME 2 outputs, builds a track table, and plots read counts through the pipeline. |
| `join_table_seqs()` | Joins a QIIME 2 feature table to a sequence-hash table. |
| `truncate_to_folder()` | Truncates a path down to a folder of interest. |

### Taxonomy assignment

| Function | Description |
|---|---|
| `assignment_12S()` | Builds a taxonomy table for 12SV5 (vertebrate) ASVs. |
| `assignment_trnL()` | Builds a taxonomy table for trnL (plant) ASVs. |
| `update_taxonomy()` | Updates a phyloseq object's taxonomic assignments against a reference. |
| `assign_common_names()` | Resolves human-readable common food names for ASVs from a reference list, with conflict handling and sibling-based propagation. |
| `lowest_level()` | Extracts each ASV's finest resolved taxonomic level. |

### Quality control

| Function | Description |
|---|---|
| `qc_controls()` | QC plots for control samples and possible contamination. |
| `plot_asv_length_hist()` | ASV length-distribution histogram. |
| `percent_assigned_tax()` | Percent of ASVs assigned at each taxonomic rank. |
| `find_g2a_c2t_pairs()` | Flags ASV pairs consistent with G→A/C→T denoising artifacts. |

### ASV harmonization & cross-batch projection

| Function | Description |
|---|---|
| `plan_harmonization()` | Detects cross-batch ASV redundancy and runs interactive review. |
| `apply_harmonization()` | Applies a harmonization plan's decisions and builds the cumulative ASV ledger. |
| `plan_projection()` | Detects correspondences between a new dataset's ASVs and a harmonized reference's fixed ASV space. |
| `apply_projection()` | Applies a projection plan, mapping a query dataset onto a reference's ASV space. |

### Filtering

| Function | Description |
|---|---|
| `filter_phyloseq()` | Subsets a phyloseq object to samples matching a metadata value. |

### Ordination

| Function | Description |
|---|---|
| `clr_transform()` | Centered log-ratio (CLR) transform of a phyloseq count table. |
| `pca_plot()` | Fits a PCA and returns a scree plot, biplot, and loadings. |
| `project_pca()` | Projects new, already-harmonized data into an existing PCA's fixed space. |
| `bstick_pc()` | Broken-stick method for choosing how many PCs to retain. |
| `elbow_pc()` | Elbow method for PC retention. |
| `paran_pc()` | Permutation-based parallel analysis for PC retention. |

### Diversity & taxa relationships

| Function | Description |
|---|---|
| `alpha_diversity()` | Per-sample alpha diversity measures, with an optional grouped boxplot. |
| `plot_taxa_correlation()` | Scatter plot of read counts between two named taxa across samples. |

Each function's full argument list and return value are documented in its
help page (`?function_name`); the handbook linked above covers how they fit
together in a typical analysis.

## Getting help

For questions or issues specific to using this package, reach out within the
David Lab. Bugs and feature requests can be filed as GitHub issues on this
repository.

## License

MIT — see [LICENSE](LICENSE).
