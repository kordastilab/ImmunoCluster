#' @rdname write_sce_to_fcs
#' @title Export an SCE to.fcs format
#'
#' @description Takes as input a SingleCellExperiment object and writes
#' the data contained into an .fcs file
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param filename the .fcs file name to export, default is export.fcs
#' @param clustering optional metadata clustering (or other) slots to export as parameters in the .fcs file
#' @param reduced_dim optional dimensionality reduction coordinates to export as a parameters in the .fcs file
#'
#' @author James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a .fcs file
#'
#' @examples
#' # Download complete ImmunoCluster SCE object from zenodo
#' sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))
#'
#' # write sct to fcs
#' write_sce_to_fcs(sce_gvhd, filename = "export.fcs")
#'
#' @import flowCore
#' @export
#'

write_sce_to_fcs <- function(sce, filename = "export.fcs", clustering = NULL, reduced_dim = NULL) {

  print('Extracting SCE to write to .fcs file')

  # Extract primary exprs from sce
  exprs = t(assay(sce))

  if(is.null(reduced_dim) == FALSE){

  # Extract and bind Dim_reduc
  dim_reduc = as.data.frame(reducedDim(sce, reduced_dim))
  exprs = cbind(exprs, dim_reduc)

  }

  # Extract metadata and convert it to numeric
  if(is.null(clustering) == FALSE){

    print('Extracting clustering to write to .fcs file as a parameter')
    print('Converting all clustering to numeric')

    # Extract the clustering and bind
    clust_data = sce@metadata[,clustering]

    # convert to numeric and bind
    if(length(clustering) == 1){

      clust_data = as.numeric(clust_data)

    }else{
      for(i in 1:ncol(clust_data)){

        clust_data[,i] = as.numeric(clust_data[,i])

     }
    }

    exprs = cbind(exprs, clust_data)

    if(length(clustering) == 1){
      colnames(exprs)[length(colnames(exprs))] = clustering
    }
  }

  # Combine into one data_frame
  dta = exprs

  # you need to prepare some metadata
  meta = data.frame(name=dimnames(dta)[[2]],
                     desc=paste('export_',dimnames(dta)[[2]])
  )
  meta$range = apply(apply(dta,2,range),2,diff)
  meta$minRange = apply(dta,2,min)
  meta$maxRange = apply(dta,2,max)

  head(meta)

  mat.dta <- as.matrix(dta)

  # All these parameters are required to build the flowframe
  # The flowFrame is the internal representation of a FCS file
  ff = new("flowFrame",
            exprs=mat.dta,
            parameters=AnnotatedDataFrame(meta)
  )

  print('Writing to .fcs file')
  write.FCS(ff, filename = filename)

}

