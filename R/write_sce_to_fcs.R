#' @rdname write_sce_to_fcs
#' @title Export an SCE to.fcs format
#'
#' @description Takes as input a SingleCellExperiment object and writes
#' the data contained into an .fcs file
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param filename the .fcs file name to export
#' @param clustering metadata clustering slots to export as parameters in the .fcs file
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
#' write_sce_to_fcs(sce_gvhd)
#'
#' @import flowCore
#' @export
#'

write_sce_to_fcs <- function(sce, filename = "export.fcs", clustering = NULL) {

  # Extract primary exprs from sce
  exprs = t(assay(sce))

  # Extract metadata and convert it to numeric
  if(is.null(clustering) == FALSE){

    # Extract the clustering
    clust_data = data@metadata[,clustering]

    exprs = cbind(exprs, clust_data)
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

  # all these are required for the following steps to work

  # a flowFrame is the internal representation of a FCS file
  ff = new("flowFrame",
            exprs=mat.dta,
            parameters=AnnotatedDataFrame(meta)
  )

  write.FCS(ff, filename = filename)

  return()
}

