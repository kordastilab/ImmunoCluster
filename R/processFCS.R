#' @rdname processFCS
#'
#' @title Import .fcs file data into a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object with experiment metadata.
#'
#' @param files a character string representing the filenames of the .fcs files to import.
#' @param assayname a character string to name the assay within the \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param metadata the metadata must be a dataframe where each row represents an individial .fcs file to be incorporated into
#' the SCE object. The number of rows must match the number of files. Each column represents a metadata attribute stored in the SCE.
#' @param filter a boolean of whether to perform noise filtering and normalization on the data.
#' @param bgNoiseThreshold an numeric representing the background noise threshold to filter.
#' @param euclideanNormThreshold an numeric as input to euclidean normalization function.
#' @param transformation a boolean of whether to transform the data.
#' @param TransFun a funtion to perform transformation on the data with, is by default the asinh with cofactor 5.
#' @param asinhFactor a numeric that will be used as the cofactor in asinh transformation.
#' @param downsample a numeric defining the number of cells to downsample to from all cells (if downsample_grouping is NULL).
#' @param downsample_grouping a character string specifying a metadata slot from which to downsample evenly by the numeric specified in downsample.
#' n downsampled cells will be chosen evenly from each downsample_grouping.
#' @param downsampleVar a boolean defining whether to donwsample from all cells by variance. Does not work in combination with downsample_grouping.
#' @param colsDiscard a vector of character strings of column names to discard.
#' @param colsRetain a vector of character strings of column names to retain, inverse of colsDiscard.
#' @param newColnames a vector of character strings to rename the channel parameters from the raw .fcs data. The order of newColnames must match the
#' the exact channels to be renamed in order after filtering. If colsRetain or colsDiscard has been used to filter out certain channels this needs to be taken into
#' account when renaming the subsequent filtered channel list or the new names may be attached to the wrong channels.
#'
#' @author Kevin Blighe, James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#'
#' @examples
#' # Create SCE object with CyTOF data inside
#' require(SingleCellExperiment)
#' sce_gvhd = processFCS(
#'  files = filelist,
#'  metadata = metadata,
#'  filter = FALSE, # Do not perform noise filtering
#'  transformation = TRUE, # Transform data
#'  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
#'  downsampleVar = NULL, # No downsampling by variance
#'  downsample = NULL, # Downsample to n cells per downsample_grouping
#'  downsample_grouping = NULL, # if downsample = n, downsample n cells by the specified metadata grouping.
#'  newColnames = col_names, # Rename columns from panel metadata
#'  colsDiscard = colsDiscard) # Discard columns not selected for downstream analysis
#'
#' @import SingleCellExperiment
#' @import flowCore
#' @export

processFCS <- function(
  files,
  assayname = 'scaled',
  metadata = NULL,
  filter = FALSE,
  bgNoiseThreshold = 1,
  euclideanNormThreshold = 1,
  transformation = TRUE,
  transFun = function (x) asinh(x),
  asinhFactor = 5,
  downsample = 100000,
  downsample_grouping = 'group',
  downsampleVar = NULL,
  colsDiscard = c('Time','Event_length','Center','Offset','Width',
    'Residual','tSNE1','tSNE2','BCKG'),
  colsRetain = NULL,
  newColnames = NULL)
{
  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as filelist
  if (!is.null(metadata)) {
    if(!identical(files, rownames(metadata))) {
      stop("'filelist' is not identical to 'rownames(metadata)'")
    }
  }

  # read in the data to a list
  samples <- list()
  samples <- lapply(files,
    function(x) exprs(read.FCS(x, transformation = FALSE)))
  names(samples) <- files

  # filter markers out
  if (!is.null(colsDiscard)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsDiscard)) > 0) {
        x[,-which(colnames(x) %in% colsDiscard)]} else {return(x)})
  }

  # filter markers in
  if (!is.null(colsRetain)) {
    samples <- lapply(
      samples,
      function(x) if (length(which(colnames(x) %in% colsRetain)) > 0) {
        x[,which(colnames(x) %in% colsRetain)]} else {return(x)})
  }

  # rename markers
  if(!is.null(newColnames)) {
    for(i in seq_len(length(samples))) {
      colnames(samples[[i]]) <- newColnames
    }
  }

  # filter
  if (filter == TRUE) {
    message('--filtering background / noise')

    # Euclidean norm
    samples <- lapply(
      samples,
      function(x)
        x[apply(x, 1, FUN = function(x) sqrt(sum(x^2))) > euclideanNormThreshold,])

    # noise correction
    for(i in seq_len(length(samples))) {
      x <- samples[[i]]
      x[x < bgNoiseThreshold] <- 0
      samples[[i]] <- x
    }
  }

  # transform
  if (transformation == TRUE) {
    message('--transforming data')
    samples <- lapply(
      samples,
      function(x) transFun(x / asinhFactor))
  }

  # load function for downsampling based on variance
  if(!is.null(downsampleVar)) {
    if (downsampleVar > 0) {
      samples <- lapply(
        samples,
        function(x) downsampleByVar(x, varianceFactor = downsampleVar))
    }
  }

  # is there metadata?
  names <- colnames(metadata)
  metanew <- list()
  if (!is.null(metadata)) {
    for (i in seq_len(length(samples))) {
      tmp <- data.frame(row.names = seq_len(nrow(samples[[i]])))
      for (j in seq_len(ncol(metadata))) {
        tmp <- cbind(tmp, rep(metadata[i,j], nrow(samples[[i]])))
      }
      metanew[[i]] <- tmp
    }

    metadata <- do.call(rbind, metanew)
    colnames(metadata) <- names
    rownames(metadata) <- paste0('cell', seq_len(nrow(metadata)))
  }

  # combine all samples
  samples <- do.call(rbind, samples)
  rownames(samples) <- paste0('cell', seq_len(nrow(samples)))

  # downsample # I want to downsample by sample_id
  if (!is.null(downsample)) {
    if (downsample > nrow(samples)) {
      warning('Cannot downsample to ', downsample, ' number of variables as',
        ' there are ', nrow(samples), ' variables currently in the merged ',
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
      inds <- split(1:length(metadata[,downsample_grouping]), metadata[,downsample_grouping])
      ## How many cells to downsample per-sample
      down_cells <- pmin(table(metadata[,downsample_grouping]), downsample)
      ## Get subsampled indices
      set.seed(2234)
      down_inds <- lapply(names(inds), function(i){
        s <- sample(inds[[i]], down_cells[i], replace = FALSE)
      })
      down_inds <- unlist(down_inds)
      samples <- samples[down_inds,]
      metadata <- metadata[down_inds,]

      rownames(metadata) <- paste0('cell', seq_len(nrow(metadata)))
      rownames(samples) <- paste0('cell', seq_len(nrow(samples)))
    }
  }

  # these should be equal
  if (!is.null(metadata)) {
    if (nrow(metadata) != nrow(samples)) {
      stop(paste0('Metadata does not match expression data',
        ' - please check your input.'))
    }
  }

  # return a SingleCellExperiment object
  ret <- list(t(samples))
  names(ret)[1] <- assayname
  ret <- SingleCellExperiment(
    assays = ret)
  metadata(ret) <- metadata
  return(ret)
}
