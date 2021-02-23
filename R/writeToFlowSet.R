#' @export
writeToFlowSet <- function(data) {


  dta <- data

  # prepare the metadata flowset export
  meta <- data.frame(name=dimnames(dta)[[2]],
                     desc=paste('export',dimnames(dta)[[2]])
  )
  meta$range <- apply(apply(dta,2,range),2,diff)
  meta$minRange <- apply(dta,2,min)
  meta$maxRange <- apply(dta,2,max)

  head(meta)

  mat.dta <- as.matrix(dta)

  # all these are required for the following steps to work
  # a flowFrame is the internal representation of a FCS file
  ff <- new("flowFrame",
            exprs=mat.dta,
            parameters=AnnotatedDataFrame(meta)
  )

  return(ff)

}
