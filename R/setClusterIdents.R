#' @rdname setClusterIdents
#'
#' @title Relabel clustering identities in sce metadata.
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param orig.ident a vector of original cluster names present.
#' @param new.ident a vector of the same length with the new cluster identities to map orig.ident to.
#' @param clustering the name of the clustering metadata slot to map the orig.ident to.
#'
#' @return an sce with a new metadata slot called 'cell_annotation' containing the new cluster mappings.
#'
#' @export
setClusterIdents <- function(
  sce,
  orig.ident,
  new.ident,
  clustering = 'phenograph'
){

  ## Convert to factor with merged clusters in desired order
  # This would mean I have to order all the factors, there are too many
  levels_clusters_merged = unique(new.ident)
  new.ident = factor(new.ident, levels = levels_clusters_merged)

  # Map the clusters to a new annotation
  match_idx = match(sce@metadata[clustering][,1], orig.ident)
  cell_annotation = new.ident[match_idx]

  sce@metadata$cell_annotation = cell_annotation

  return(sce)

}
