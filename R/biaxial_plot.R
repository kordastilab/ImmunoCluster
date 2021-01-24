#' @rdname biaxial_plot
#'
#' @title Generate a biaxial dot plot with contouring of the SCE data 
#'
#' @param indata a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay the SCE assay to choose from.
#' @param marker_x marker to display on the x axis
#' @param marker_y marker to display on the y axis
#' @param clustering metadata slot to colour the plot by
#' @param colkey a list of colours matching the number of conditions found in clustering
#'
#' @author James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#' # Download complete ImmunoCluster SCE object from zenodo
#' sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))
#' 
#' biaxial_plot(sce_gvhd, marker_x = "CD8a", marker_y = "CD4")
#'
#'
#' @export
#'
biaxial_plot <- function(indata,
                         assay = "scaled",
                         marker_x,
                         marker_y,
                         clustering = NULL,
                         colkey = NULL
){
  
  # Extract markers from the sce
  markers = c(marker_x, marker_y)
  
  data = as.data.frame(t(assay(indata))[,markers])
   
  # If clustering isn't null
  if(is.null(clustering) == FALSE){
    
    cell_clust = indata@metadata[,clustering]
    
    data_clust = cbind(data, cell_clust)
    
    plot_obj = ggplot(data, aes_string(marker_x, marker_y, colour = cell_clust)) +
      geom_point(size = 0.01) +
      geom_density2d(binwidth = 0.006) +
      theme_classic()
    
    # sort out custom colour pairing
    if (!is.null(colkey)) {
      plot_obj <- plot_obj + scale_colour_discrete('') +
        scale_color_manual(values = colkey)
    }
    
  }
  else{
  
  # Otherwise do grey/black
  # Plot the biaxial plot
  plot_obj = ggplot(data, aes_string(marker_x, marker_y)) +
              geom_point(colour = "grey", size = 0.01) +
              geom_density2d(color='black', binwidth = 0.006) +
              theme_classic()
  
  }
  
  return(plot_obj)
}
 
  
  
