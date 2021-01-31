#' @import broom
#' @export
#'
diffClust <- function(
  sct,
  group = 'group',
  clustering = 'cell_annotation',
  feature = 'condition',
  p_val = 'padj',
  threshold = 0.1
){

  # Function deprecated
  .Deprecated("stat_test_clust")

  stat_test_clust(sct,
    group = group,
    clustering = clustering,
    feature = feature,
    p_val = p_val,
    threshold = threshold,
    test = "wilcox")
}
