#' @export
performUMAP <- function(
  indata,
  assay = 'scaled',
  downsample = NULL,
  downsample_grouping = 'group',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL)
{
  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')

    if (!is.null(reducedDim)) {
        message('--input data is taken from \'', reducedDim,
          '\' dimensional reduction')
        message('--Dimensions to use: ', paste(dims, collapse = ', '))
        mat <- as.matrix(reducedDims(indata)[[reducedDim]][,dims])

        if (is.null(newDimName)) {
          newDimName <- paste0('UMAP_', reducedDim)
        }
    } else {
      message('--input data is taken from \'', assay, '\' assay slot')
      mat <- t(assay(indata, assay))

      if (is.null(newDimName)) {
        newDimName <- 'UMAP'
      }
    }

    # Downsample
    # downsample # I want to downsample by sample_id
    if (!is.null(downsample)) {
      if (downsample > nrow(mat)) {
        warning('Cannot downsample to ', downsample, ' number of variables as',
                ' there are ', nrow(mat), ' variables currently in the merged ',
                'dataset.')
        message('--Skipping downsampling')
      } else {
        if(is.null(downsample_grouping == TRUE)){
          message('--Cannot downsample without sample grouping')
          message('--Aborting downsampling')
          stop()
        }

        message('--Downsampling to ', downsample, ' variables.')
        message('--Downsampling equally by ', downsample_grouping, ' metadata column')

        ## Data subsampling: create indices by sample
        inds <- split(1:length(indata@metadata[,downsample_grouping]), indata@metadata[,downsample_grouping])
        ## How many cells to downsample per-sample
        down_cells <- pmin(table(indata@metadata[,downsample_grouping]), downsample)
        ## Get subsampled indices
        set.seed(2234)
        down_inds <- lapply(names(inds), function(i){
          s <- sample(inds[[i]], down_cells[i], replace = FALSE)
        })

        # downsample matrix
        down_inds <- unlist(down_inds)
        mat <- mat[down_inds,]

      }
    }

    message('--Performing UMAP...')
    if (is.null(useMarkers)) {
      u <- umap(mat)
    } else if (!is.null(useMarkers) && !is.null(reducedDim)) {
      warning('\'useMarkers\' and \'reducedDim\' are incompatible - ',
        'markers cannot be selected from a reduced dimensional ',
        'component in which they don\'t exist! Dimensions to use ',
        'for UMAP have already been chosen via the \'dims\' parameter')
      u <- umap(mat)
    } else {
      message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers])
    }

    # Create vector of NAs to fill with the inds
    if (!is.null(downsample)) {

      # Create empty vector to fill
      umap_layout = data.frame('UMAP1' = rep(NA, length(colnames(indata))), 'UMAP2' = rep(NA, length(colnames(indata))))

      # Fill the index with the slots
      for(i in 1:length(u$layout[,1])) {

        umap_layout[down_inds[i],] = u$layout[i,]

      }


    } else {
      umap_layout = u$layout

      colnames(umap_layout) <- c('UMAP1', 'UMAP2')
    }

    reducedDim(indata, newDimName) <- umap_layout

    message('--Done')

    return(indata)

  } else {

    message('--input data class is ', class(indata))
    message('Note: all non-SingleCellExperiment objects will be ',
      'coerced to matrix')
    mat <- t(as.matrix(indata))

    message('--Performing UMAP...')
    if (is.null(useMarkers)) {
      u <- umap(mat)
    } else {
      message('Note: only using the following markers for UMAP calculation: ',
        paste(useMarkers, collapse = ', '))
      u <- umap(mat[,useMarkers])
    }

    colnames(u$layout) <- c('UMAP1', 'UMAP2')

    message('--Done')

    return(u$layout)
  }
}
