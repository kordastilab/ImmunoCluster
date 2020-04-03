#' @export
subsetSCE = function(x, ...) {
  
  # check validity of input arguments
  stopifnot(is(x, "SingleCellExperiment"))
  
  # construct data.frames of cell & feature metadata
  md = vapply(x@metadata, as.character, character(ncol(x)))
  md = data.frame(i = seq_len(ncol(x)), md, 
                   check.names = FALSE, stringsAsFactors = FALSE)
  
  # filter rows & columns      ...
  mdf <- try(dplyr::filter(md, ... ), silent = TRUE)
  if (inherits(mdf, "try-error")) mdf <- md
  idx <- mdf$i; mdf <- select(mdf, -"i")
  
  # convert to factors
  mdf <- mutate_all(mdf, factor)
  
  # subset reduced dimensions
  if (length(reducedDims(x)) > 0) {
    dr <- reducedDims(x)
    dr <- lapply(dr, "[", i = idx, j = TRUE)
  } else dr <- NULL
  
  # subset assays & returned filtered SCE
  as <- lapply(assays(x), "[", i = TRUE, j = idx)
  x_subset = SingleCellExperiment(assays = as, 
                       reducedDims = dr)
  
  metadata(x_subset) = mdf
  
  return(x_subset)

  }

