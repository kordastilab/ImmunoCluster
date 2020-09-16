#' @export
subsetSCE = function(sce, ...) {

  # create dataframes of sce metadata
  md = vapply(sce@metadata, as.character, character(ncol(sce)))
  md = data.frame(i = seq_len(ncol(sce)), md,
                   check.names = FALSE, stringsAsFactors = FALSE)

  # Subset metadata dataframe
  mdf = try(dplyr::filter(md, ... ), silent = TRUE)
  if (inherits(mdf, "try-error")) mdf = md
  idx = mdf$i; mdf = select(mdf, -"i")

  # convert to factors
  mdf = mutate_all(mdf, factor)

  # subset reduced dimensions
  if (length(reducedDims(sce)) > 0) {
    red_dim = reducedDims(sce)
    red_dim = lapply(red_dim, "[", i = idx, j = TRUE)
  } else red_dim = NULL

  # subset assays & return filtered SCE
  asssay = lapply(assays(sce), "[", i = TRUE, j = idx)
  sce_subset = SingleCellEsceperiment(assays = asssay,
                       reducedDims = red_dim)

  metadata(sce_subset) = mdf

  return(sce_subset)

  }

