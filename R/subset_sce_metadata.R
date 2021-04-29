#' @export
subset_sce_metadata = function(sce, ...) {

  # create dataframes of sce metadata
  sce_meta = vapply(sce@metadata, as.character, character(ncol(sce)))

  # Initialize the dataframe of metadata
  sce_meta = data.frame(id_row = seq_len(ncol(sce)),
                        sce_meta,
                        check.names = FALSE,
                        check.rows = FALSE)

  # Subset metadata dataframe
  # Filter based on the input parameters
  sce_meta_filter = dplyr::filter(sce_meta, ... )

  # Extract the indexes to filter
  idx = sce_meta_filter$id_row

  # Drop the ID row variable
  drop_id = which(colnames(sce_meta_filter) == "id_row")
  sce_meta_filter = sce_meta_filter[,-drop_id]

  # extract the assays of the sce and subset the objects
  # based on the filter parameters
  sce_assay = lapply(assays(sce), "[", i = T, j = idx)

  # subset reduced dimensions of the sce
  # if there are any
  if (length(reducedDims(sce)) > 0) {

    # Extract the reduced dime of the sce
    red_dim = reducedDims(sce)

    # Subset based on the filter parameter
    red_dim = lapply(red_dim, "[", i = idx, j = TRUE)

  } else red_dim = NULL # If no red_dim then return an empty red_dim

  # Remake a new sce object with the assays
  sce_subset = SingleCellExperiment(assays = sce_assay)

  # Add the metadata to the sce
  metadata(sce_subset) = sce_meta_filter

  # Add the reduced dims
  reducedDims(sce_subset) = red_dim

  # Return the sce
  return(sce_subset)

  }

