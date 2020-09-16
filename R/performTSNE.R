#' @rdname performTSNE
#'
#' @title Perform TSNE dimensionality reduction on SingleCellExperiment object
#'
#' @param indata a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay the SCE assay to choose from.
#' @param downsample a numeric of the number of cells to downsample by.
#' @param downsample_grouping a character string indicating which metadata grouping to downsample by.
#' @param reducedDim a character vector represting a dimensionality reduction slot from the SCE to be selected as input data.
#' @param dims numeric of reduced dimensions to select from the reducedDim
#' @param newDimName character string of the name for the resulting TSNE dimensionality reduction in the SCE object.
#' @param useMarkers chacter vector of markers to use as input for the TSNE run.
#' @param perplexity numeric of the perplexity for the TSNE algorithm param.
#' @param theta numeric of the theta for the TSNE algorithm param.
#' @param max_iter numeric of the maximum iterations to perform for the TSNE algorithm param.
#'
#' @author Kevin Blighe, James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{SingleCellExperiment} object.
#'
#' @examples
#' # Download complete ImmunoCluster SCE object from zenodo
#' sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))
#'
#' # Run tSNE and store in sce object
#' sce_gvhd = performTSNE(sce_gvhd)
#'
#' @import Rtsne
#'
#' @export

performTSNE <- function(
  indata,
  assay = 'scaled',
  downsample = NULL,
  downsample_grouping = 'group',
  reducedDim = NULL,
  dims = seq_len(20),
  newDimName = NULL,
  useMarkers = NULL,
  perplexity = 50,
  theta = 0.5,
  max_iter = 1000)
{
  if (is(indata, 'SingleCellExperiment')) {

    message('--input data class is SingleCellExperiment')

    if (!is.null(reducedDim)) {
      message('--input data is taken from \'', reducedDim,
              '\' dimensional reduction')
      message('--Dimensions to use: ', paste(dims, collapse = ', '))
      mat <- as.matrix(reducedDims(indata)[[reducedDim]][,dims])

      if (is.null(newDimName)) {
        newDimName <- paste0('Rtsne_', reducedDim)
      }
    } else {
      message('--input data is taken from \'', assay, '\' assay slot')
      mat <- t(assay(indata, assay))

      if (is.null(newDimName)) {
        newDimName <- 'Rtsne'
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

    message('--Performing Rtsne...')
    if (is.null(useMarkers)) {
      u <- Rtsne(mat, check_duplicates = FALSE, pca = FALSE, perplexity = perplexity, theta = theta, max_iter = max_iter)
    } else if (!is.null(useMarkers) && !is.null(reducedDim)) {
      warning('\'useMarkers\' and \'reducedDim\' are incompatible - ',
              'markers cannot be selected from a reduced dimensional ',
              'component in which they don\'t exist! Dimensions to use ',
              'for Rtsne have already been chosen via the \'dims\' parameter')
      u <- Rtsne(mat, check_duplicates = FALSE, pca = FALSE, perplexity = perplexity, theta = theta, max_iter = max_iter)
    } else {
      message('Note: only using the following markers for Rtsne calculation: ',
              paste(useMarkers, collapse = ', '))
      u <- Rtsne(mat[,useMarkers], check_duplicates = FALSE, pca = FALSE, perplexity = perplexity, theta = theta, max_iter = max_iter)
    }

    # Create vector of NAs to fill with the inds
    if (!is.null(downsample)) {

      # Create empty vector to fill
      Rtsne_layout = data.frame('Rtsne1' = rep(NA, length(colnames(indata))), 'Rtsne2' = rep(NA, length(colnames(indata))))

      # Fill the index with the slots
      for(i in 1:length(u$Y[,1])) {

        Rtsne_layout[down_inds[i],] = u$Y[i,]

      }


    } else {
      Rtsne_layout = u$Y

      colnames(Rtsne_layout) <- c('Rtsne1', 'Rtsne2')
    }

    reducedDim(indata, newDimName) <- Rtsne_layout

    message('--Done')

    return(indata)

  } else {

    message('--input data class is ', class(indata))
    message('Note: all non-SingleCellExperiment objects will be ',
            'coerced to matrix')
    mat <- t(as.matrix(indata))

    message('--Performing Rtsne...')
    if (is.null(useMarkers)) {
      u <- Rtsne(mat, check_duplicates = FALSE, pca = FALSE, perplexity = perplexity, theta = theta, max_iter = max_iter)
    } else {
      message('Note: only using the following markers for Rtsne calculation: ',
              paste(useMarkers, collapse = ', '))
      u <- Rtsne(mat[,useMarkers], check_duplicates = FALSE, pca = FALSE, perplexity = perplexity, theta = theta, max_iter = max_iter)
    }

    colnames(u$Y) <- c('Rtsne1', 'Rtsne2')

    message('--Done')

    return(u$Y)
  }
}
