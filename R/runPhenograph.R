#' @rdname runPhenograph
#'
#' @title run Rphenograph implementation from "JinmiaoChenLab/Rphenograph" clustering on sce expression data
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k a number of desired clusters to cluster cell data into.
#' @param markers a vector of selected markers, representing sce colnames to use for clustering
#' @param set_seed set the random seed for clustering reproducibility
#'
#' @return a \code{SingleCellExperiment} object with a new metadata slot with cluster identities.
#'
#' @import Rphenograph
#'
#' @export
runPhenograph <- function(
  sce,
  k = 30,
  markers = NULL,
  set_seed = 4321
){

  colName = paste0('phenograph_k',k)

  # Need to implement marker selection
  if(is.null(markers) == FALSE){

    data = t(assay(sce))

    # Set the seed
    set.seed(set_seed)

    rphenograph_output = Rphenograph(data[,markers], k = k)

    identity = factor(membership(rphenograph_output[[2]]))

    sce@metadata[,colName] = identity

  return(sce)

  }
  else{

    # Set the seed
    set.seed(set_seed)

    rphenograph_output <- Rphenograph(t(assay(sce)), k = k)

    identity <- factor(membership(rphenograph_output[[2]]))

    sce@metadata[,colName] = identity

    return(sce)

  }

}
