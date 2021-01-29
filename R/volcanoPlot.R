#' @rdname volcanoPlot
#'
#' @title Generate a volcano plot using statistical testing output
#'
#' @param pval_df a dataframe that contains a column name "logfc" and a statistical testing output p values.
#' @param p_val a string identifying the column name from pval_df to plot as p values
#' @param threshold a numeric representing the p_value threshold on which to colour points on the volcano plot
#'
#' @author James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#'
#' @export
#'
volcanoPlot <- function(
  pval_df,
  p_val = "padj",
  threshold = 0.1
){

  # Extract plot data
  plotdf <- data.frame(logfc = pval_df$logfc, pvalue = pval_df[,p_val], cluster = pval_df$cluster)

  # Remove INF rows
  plotdf <- plotdf[!is.infinite(plotdf$logfc),]

  # Create volcano plot
  gg_logFC = ggplot(plotdf, aes(x=logfc, y=-log(pvalue))) +
    geom_point(aes(color=(-log(pvalue)>-log(threshold)), size=1)) +
    theme_bw() +
    scale_color_manual(values = c('grey', 'red2')) +
    geom_label_repel(data = plotdf,
                     aes(x=logfc, y=-log(pvalue),
                         label = ifelse(-log(pvalue) > -log(threshold),
                                        as.character(cluster) ,"")), point.padding = 2, force = 10) +
    xlab("LogFC") + ylab("-log(FDR)") + theme(legend.position = "none")

  return(gg_logFC)

}
