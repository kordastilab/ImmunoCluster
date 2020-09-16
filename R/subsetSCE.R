#' @export
subsetSCE = function(sce, ...) {

  # create dataframes of sce metadata
  md = vapply(sce@metadata, as.character, character(ncol(sce)))
  md = data.frame(k = seq_len(ncol(sce)), md,
                   check.names = FALSE, stringsAsFactors = FALSE)

  # Subset metadata dataframe
  md_filter = dplyr::filter(md, ... )
  idx = md_filter$k
  md_filter = select(md_filter, -"k")

  # convert to factors
  md_filter = mutate_all(md_filter, factor)

  # subset reduced dimensions
  if (length(reducedDims(sce)) > 0) {
    red_dim = reducedDims(sce)
    red_dim = lapply(red_dim, "[", i = idx, j = TRUE)
  } else red_dim = NULL

  # subset assays & return filtered SCE
  asssay = lapply(assays(sce), "[", i = TRUE, j = idx)
  sce_subset = SingleCellExperiment(assays = asssay,
                       reducedDims = red_dim)

  metadata(sce_subset) = md_filter

  return(sce_subset)

  }

