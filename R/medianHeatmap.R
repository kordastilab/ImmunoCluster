#' @export
medianHeatmap <- function(
  data,
  grouping = NULL,
  clusters = unique(data@metadata[,grouping]),
  markers = rownames(data),
  feature = NULL,
  feature_cols = NULL,
  heat_bar = "Greens",
  scale_01 = FALSE
)
{

  # 01 scale the data
  if(scale_01 == TRUE){

    # Extract the expression
    expr <- t(assay(data))

    # Take the upper-lower quantile
    col_quant <- colQuantiles(expr, probs = c(0.01, 0.99))

    # Scale to min max
    scaled_expression_mat <- t((t(expr) - col_quant[, 1]) / (col_quant[, 2] - col_quant[, 1]))

    # set upper lower bounds
    scaled_expression_mat[scaled_expression_mat < 0] <- 0
    scaled_expression_mat[scaled_expression_mat > 1] <- 1

    # Return expression
    expr = scaled_expression_mat

  } else {

  expr = t(assay(data))

  }

  # Get the median marker expression per sample
  data.median = data.frame(sample_id = data@metadata[,grouping], expr)
  colnames(data.median) = c('sample_id', rownames(data))

  # Subset
  data.median = data.median[which(data.median$sample_id == clusters), c('sample_id', markers)]

  # Find medians
  data.median = data.median %>%
    group_by(sample_id) %>%
    summarize_all(list(median))

  # reshape to matrix
  data.median.sample = t(data.median[, -1])
  colnames(data.median.sample) = data.median$sample_id

  # Using metadata feature slot
  if(is.null(feature) == FALSE){

    # Summarize sampleID to metadata feature
    feature.summary = data.frame(sample_id = data@metadata[,grouping], feat = data@metadata[,feature]) %>%
      distinct()

    # Match feature
    mm = match(colnames(data.median.sample), feature.summary$sample_id)
    annotation_col = data.frame(feature = feature.summary$feat[mm],
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
                       display_numbers = F,
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
                       display_numbers = F,
                       clustering_method = "average")

    return(heatmap)

  }




}

