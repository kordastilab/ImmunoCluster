#' @export
runFlowSOM <- function(
  sct,
  k = 30,
  markers = NULL,
  som_x = 10,
  som_y = 10,
  set_seed = 4321
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
    fsom = ReadInput(ff, transform = FALSE, scale = FALSE)

    # Set the seed
    set.seed(set_seed)

    # run flowsom
    som = BuildSOM(fsom, colsToUse = markers, xdim = som_x, ydim = som_y)

    # cell codes
    cell_som_mapping = som$map$mapping[,1]
    sct@metadata$som_codes = cell_som_mapping

    # extract codes for cc_plus
    codes = som$map$codes

    plot_outdir = "consensusC_plus_plots"

    pdf(NULL)
    # Run consensus clsuter plus
    clusters = suppressMessages(ConsensusClusterPlus(t(codes), maxK = k, reps = 100,
                               pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
                               clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                               distance = "euclidean", seed = 1234))
     dev.off()


    # Create map metaclustering
    for(i in 2:k){

      ## Get cluster ids for each cell for each cluster
      code_clustering <- clusters[[i]]$consensusClass
      flowsom_cc_clustering <- code_clustering[cell_som_mapping]

      identity = factor(flowsom_cc_clustering)

      colName = paste0('flowsom_cc_k',i)

      sct@metadata[,colName] = identity

    }

    dlta = plot_delta_elbow(clusters)

    print(dlta)

    return(sct)


}
