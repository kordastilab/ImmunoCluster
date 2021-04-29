#' @export
subsetSCE = function(sce, ...) {

  # create dataframes of sce metadata
  sce_meta = vapply(sce@metadata, as.character, character(ncol(sce)))

  # Initialize the dataframe of metadata
  sce_meta = data.frame(id_row = seq_len(ncol(sce)), sce_meta,
                   check.names = FALSE, stringsAsFactors = FALSE)

  # Subset metadata dataframe
  # Filter based on the input parameters
  sce_meta_filter = dplyr::filter(sce_meta, ... )

  # Extract the indexes to filter
  idx = sce_meta_filter$id_row

  # Select the indexes from the metadata
  sce_meta_filter = dplyr::select(sce_meta_filter, -"id_row")

  # convert the metadaata to factors
  sce_meta_filter = dplyr::mutate_all(sce_meta_filter, factor)

  # extract the assays of the sce and subset the objects
  # based on the filter prarameters
  sce_assay = lapply(assays(sce), "[", i = TRUE, j = idx)

  # subset reduced dimensions of the sce
  # if there are any
  if (length(reducedDims(sce)) > 0) {

    # Extract the reduced dime of the sce
    red_dim = reducedDims(sce)

    # Subset based on the filter parameter
    red_dim = lapply(red_dim, "[", i = idx, j = TRUE)

  } else red_dim = NULL # If no red_dim then return an empty red_dim

  # Remake a new sce object with the assays and the reduced dims
  sce_subset = SingleCellExperiment(assays = sce_assay,
                       reducedDims = red_dim)

  # Add the metadata to the sce
  metadata(sce_subset) = sce_meta_filter

  # Return the sce
  return(sce_subset)

  }

