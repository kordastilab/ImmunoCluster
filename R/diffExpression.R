#' @import broom
#' @export
#'
diffExpression <- function(
  indata,
  assay = 'scaled',
  grouping = 'group', # The condensing function
  feature = 'condition', # The contrast
  clusterAssign = 'cell_annotation' # The clustering
){

  # Function deprecated
  .Deprecated("stat_test_expression")

  # Run stat_test_expression
  stat_test_expression(indata,
                       assay = 'scaled',
                       grouping = 'group',
                       feature = 'condition',
                       clusterAssign = 'cell_annotation')

}
