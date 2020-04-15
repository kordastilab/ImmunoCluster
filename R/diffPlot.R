#' @export
diffPlot <- function(
  pval_df,
  p_val = "padj", 
  threshold = 0.1
){
  
  # Create volcano plot
  gg_logFC = ggplot(pval_df, aes(x=logfc, y=-log(pval_df[,p_val]))) +
    geom_point(aes(color=(-log(pval_df[,p_val])>-log(threshold)), size=1.3)) +
    theme_bw() +
    scale_color_manual(values = c('grey', 'red2')) +
    geom_label_repel(data = pval_df,
                     aes(x=logfc, y=-log(pval_df[,p_val]),
                         label = ifelse(-log(pval_df[,p_val]) > -log(threshold),
                                        as.character(cluster) ,"")), point.padding = 2) +
    xlab("LogFC") + ylab("-log(FDR)") + theme(legend.position = "none")
  
  return(gg_logFC)
  
} 
