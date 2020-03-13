#' @export
subsetData <- function(
  indata,
  ident,
  clusters
) 
{
  
  ids = which(indata@metadata[,ident] == clusters)
  indata_sub = indata[,ids]
  indata_sub@metadata = indata_sub@metadata[ids,]
  
  return(indata_sub)
}

