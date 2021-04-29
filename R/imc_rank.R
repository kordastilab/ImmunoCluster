#' @rdname imc_rank
#'
#' @title Creates a ranked heatmap based on SingleCellExperiment object
#'
#' @param indata a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param grouping The metadata grouping to create the median marker exrpression out of
#' @param clusters The specific clusters from this grouping to select
#' @param markers Markers to include for the heatmap
#' @param feature_cols Colors for reature bar annotation
#' @param feature A metadata feature to overlay on the heatmap as an annotation bar
#' @param heat_bar The 'RColorBrewer' heat_bar color ramp designation
#' @param scale_01 A logical determining whether to scale the data
#' @param display_numbers A logical determining whether to display data in heatmap cells
#'
#' @author Kevin Blighe, James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{SingleCellExperiment} object.
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @export

imc_rank <- function(
  indata,
  grouping = NULL,
  clusters = unique(indata@metadata[,grouping]),
  markers = rownames(indata),
  feature = NULL,
  feature_cols = NULL,
  heat_bar = "YlOrRd",
  scale_01 = TRUE,
  display_numbers = TRUE
)
{

  # 01 scale the data
  if(scale_01 == TRUE){

    expr = t(assay(indata))

    rng <- colQuantiles(expr, probs = c(0.01, 0.99))
    expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
    expr01[expr01 < 0] <- 0
    expr01[expr01 > 1] <- 1

    expr = expr01

  } else {

    expr = t(assay(indata))

  }

  # Get the median marker expression per sample
  data.median = data.frame(sample_id = indata@metadata[,grouping], expr)
  colnames(data.median) = c('sample_id', rownames(indata))

  # Subset
  data.median = data.median[which(data.median$sample_id == clusters), c('sample_id', markers)]

  # Find medians
  data.median = data.median %>%
    group_by(sample_id) %>%
    summarize_all(list(median))


  # Create ranking for each marker and cluster
  ranked_data.median = as.data.frame(apply(data.median[,-1], 2, function(x) rank(x)))
  ranked_data.median = cbind("sample_id" = data.median$sample_id, ranked_data.median)

  # reshape to matrix for heatmapping
  data.median.sample = t(ranked_data.median[, -1])
  colnames(data.median.sample) = ranked_data.median$sample_id



  # Using metadata feature slot
  if(is.null(feature) == FALSE){

    # Summarize sampleID to metadata feature
    feature.summary = data.frame(sample_id = data@metadata[,grouping], feat = data@metadata[,feature]) %>%
      distinct()

    # Match feature
    match_idx = match(colnames(data.median.sample), feature.summary$sample_id)
    annotation_col = data.frame(feature = feature.summary$feat[match_idx],
                                row.names = colnames(data.median.sample))

    # If feature Cols is not NULL then do the ggplot hues function
    if(is.null(feature_cols) == FALSE){

      # possible colors
      anno_colors =  feature_cols

    }else{
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }

      anno_colors =  gg_color_hue(length(unique(feature.summary$feat)))

    }

    # Create annotation
    condition = anno_colors[1:length(unique(feature.summary$feat))]
    names(condition) = unique(feature.summary$feat)
    annotation_colors = list(feature = condition)

    # Colors for the heatmap
    color_bar = brewer.pal(n = 9, name = heat_bar)

    ## Somethings not quite lining up here
    # execute heatmap
    heatmap = pheatmap(data.median.sample, color = color_bar,
                       display_numbers = display_numbers,
                       fontsize_number = 9,
                       annotation_col = annotation_col,
                       annotation_colors = annotation_colors,
                       clustering_method = "average")

    return(heatmap)

  } else{


    # Colors for the heatmap
    color_bar = brewer.pal(n = 9, name = heat_bar)

    # execute heatmap
    heatmap = pheatmap(data.median.sample,
                       color = color_bar,
                       display_numbers = display_numbers,
                       clustering_method = "average")

    return(heatmap)

  }




}

