#' @rdname runFlowSOM
#'
#' @title Perform combined FlowSOM clustering followed by consensus clustering as outlined in Nowicka et al. (2017) F1000Research.
#'
#' @param sct a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k a numeric representing the desired maximum number of clusters to group cells into by consensus clustering, the funtion will return all clustering k values from 2-k and these will be stored in the metadata table.
#' @param markers a vector of strings representing the markers to cluster on.
#' @param som_x a numeric indicating the x dimension of the som grid, see FLowSOM documentation.
#' @param som_y a numeric indicating the y dimension of the som grid.
#' @param set_seed a numeric to set the random seed for run reproducibility.
#' @param scale a boolean indicating whether to scale the input data by marker.
#'
#' @author James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{SingleCellExperiment} object.
#'
#' @examples
#' # Download complete ImmunoCluster SCE object from zenodo
#' sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))
#'
#' # Run FlowSOM and store in sce object to clustering k=
#' sce_gvhd = runFlowSOM(sce_gvhd, k = 10, markers = clustering_markers, som_x = 10, som_y = 10)
#'
#' @import FlowSOM
#' @import flowCore
#' @import ConsensusClusterPlus
#'
#' @export
#'
runFlowSOM <- function(
  sct,
  k = 30,
  markers = NULL,
  som_x = 10,
  som_y = 10,
  set_seed = 4321,
  scale = FALSE
){

  if(is.null(markers) == TRUE){

    markers = rownames(sct)

  }
    # Extract data
    data = t(assay(sct))

    # Select markers
    data = data[,markers]

    # create flowset
    ff = writeToFlowSet(data)

    # build igraph object
    fsom = ReadInput(ff, transform = FALSE, scale = scale)

    # Set the seed
    set.seed(set_seed)

    # run flowsom
    som = BuildSOM(fsom, colsToUse = markers, xdim = som_x, ydim = som_y)

    # cell codes
    cell_som_mapping = som$map$mapping[,1]
    sct@metadata$som_codes = cell_som_mapping

    # extract codes for cc_plus
    fsom_codes = som$map$codes

    # Where are the plots going to be put
    cc_plotting_dir = "cc_plus_elbow_plots"

    pdf(NULL)
    # Run consensus cluster plus
    clusters = suppressMessages(ConsensusClusterPlus(t(fsom_codes), maxK = k, reps = 100,
                               pItem = 0.9, pFeature = 1, title = cc_plotting_dir, plot = "png",
                               clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                               distance = "euclidean", seed = 1234))
     dev.off()


    # Create the metadata map of the metaclustering
    for(i in 2:k){

      # Retrieve the clustering ids for each cell for each cluster
      fsom_metaclustering <- clusters[[i]]$consensusClass

      # Retrieve the cell_cluster alignments
      flowsom_cc_clustering <- fsom_metaclustering[cell_som_mapping]

      identity = factor(flowsom_cc_clustering)

      colName = paste0('flowsom_cc_k',i)

      sct@metadata[,colName] = identity

    }

    return(sct)


}
