#' @title Create a ggplot object of sce cluster proportion abundances grouped by sample and metadata condition
#'
#' @param indata A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param graph_type Whether to create a boxplot with 'box' or stacked barchart with 'bar'
#' @param group The metadata slot representing sample_id
#' @param clusters The specific clusters to plot within the clustering metatdata slot
#' @param clusterAssign the clustering metadata or other metadata slot by which to group the samples specified in 'group'
#' @param feature A feature representing a metadata slot by which to split the plots
#' @param colkey A colour key by which to manually specify plotting colours
#' @param legendPosition Where to position the legend
#' @param pointSize Dotplot pointsize
#' @param legendLabSize legend label size
#' @param legendIconSize legend key icon size
#' @param legendKeyHeight legend key height
#' @param xlim specify the x limit
#' @param ylim specify the y limit
#' @param xlab specify the x axis label
#' @param xlabAngle specify the x lab angle
#' @param xlabhjust Specify the horizontal justification of the x axis (0 = left jsutified, 1 = right justified)
#' @param xlabvjust Specify the vertical justification of the x axis (0 = left jsutified, 1 = right justified)
#' @param ylab specify the y axis label
#' @param ylabAngle specify the y lab angle
#' @param ylabAngle specify the y lab angle
#' @param ylabhjust Specify the horizontal justification of the y axis (0 = left jsutified, 1 = right justified)
#' @param ylabvjust Specify the vertical justification of the y axis (0 = left jsutified, 1 = right justified)
#' @param axisLabSize specify the axis label size
#' @param stripLabSize specify the facet labels of a plot
#' @param title specify the title text of the plot
#' @param subtitle specify the subtitle text of the plot
#' @param caption specify the caption text of the plot
#' @param titleLabSize specify the title text size
#' @param subtitleLabSize specify the subtitle text size
#' @param captionLabSize specify the caption text size
#'
#' @export
plotAbundance = function(
  indata,
  graph_type = 'box',
  group = 'group',
  clusters = unique(metadata(indata)[['cell_annotation']]),
  clusterAssign = 'cell_annotation',
  feature = 'condition',
  colkey = NULL,
  legendPosition = 'right',
  pointSize = 1,
  legendLabSize = 12,
  legendIconSize = 5.0,
  legendKeyHeight = 2.5,
  xlim = NULL,
  ylim = NULL,
  xlab = 'Population',
  xlabAngle = 90,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = '% of cells',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  stripLabSize = 16,
  title = 'Marker expression per cluster',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  borderWidth = 0.8,
  borderColour = 'black'
){

  if (is(indata, 'SingleCellExperiment')) {
    message('--input data class is SingleCellExperiment')
  } else {
    message('--input data class is ', class(indata))
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background = element_rect(),

      title = element_text(size = legendLabSize),

      plot.title=element_text(angle=0, size=titleLabSize,
                              face='bold', vjust=1),
      plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
                                 face = 'plain', vjust = 1),
      plot.caption=element_text(angle = 0, size = captionLabSize,
                                face = 'plain', vjust = 1),

      axis.text.x=element_text(angle = xlabAngle, size = axisLabSize,
                               hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = axisLabSize,
                               hjust = ylabhjust, vjust = ylabvjust),
      axis.title = element_text(size = axisLabSize),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      legend.title = element_blank(),
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = legendLabSize),
      legend.key.height = unit(legendKeyHeight, 'cm'),

      strip.text.x = element_text(size = stripLabSize,
                                  face = 'bold', margin = margin(b = 5, t = 5)))

  # Create the sample-cluster counts table
  counts_table = table(indata@metadata[,clusterAssign], indata@metadata[,group])

  # Create the proportion table from the counts table
  props_table = t(t(counts_table) / colSums(counts_table)) * 100

  # Select the specified clusters
  props = as.data.frame.matrix(props_table)[clusters,]

  plotobj = melt(data.frame(cluster = rownames(props), props),
               id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")

  plotobj$cluster <- factor(plotobj$cluster)
  plotobj$sample_id = gsub("X", "", plotobj$sample_id)

  if(is.null(feature) == FALSE){ # No contrast

  # Summarize sampleID to metadata feature
  feature.summary = data.frame(Grouping = indata@metadata[,group], feat = indata@metadata[,feature]) %>%
    distinct()

  ## Add condition info
  match_idx <- match(plotobj$sample_id, feature.summary$Grouping)

  # Need to recover the original matrix
  plotobj$condition <- factor(feature.summary$feat[match_idx])

  if(graph_type == 'bar'){

    # Plot group individual bars but do stacked
    # Currently not stacked
    plot = ggplot(plotobj, aes(x = sample_id, y = proportion, fill = cluster)) + th +

      guides(
        fill = guide_legend(),
        shape = guide_legend(),
        alpha = FALSE) +

      geom_bar(stat = "identity") +
      facet_wrap(~ condition, scales = "free_x") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab(xlab) + ylab(ylab)

    # sort out custom colour pairing,
    if (!is.null(colkey)) {
       plot <- plot + scale_colour_discrete('') +
        scale_fill_manual(values = colkey)
    }

    plot

  }

  else{

    plot = ggplot(plotobj, aes(x = cluster, y = proportion, fill = condition)) + th +
      geom_boxplot(outlier.shape = NA) +
      geom_point(aes(color = condition), alpha = 0.5, position = position_jitterdodge(),
                 size = pointSize) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab(xlab) + ylab(ylab)

    # sort out custom colour pairing,
    if (!is.null(colkey)) {
      plot <- plot + scale_colour_manual(values = colkey) +
        scale_fill_manual(values = colkey)
    }

    plot

  }

  }
  else{

    if(graph_type == 'bar'){

      # Plot group individual bars but do stacked
      # Currently not stacked
      plot = ggplot(plotobj, aes(x = sample_id, y = proportion, fill = cluster)) + th +

        guides(
          fill = guide_legend(),
          shape = guide_legend(),
          alpha = FALSE) +

        geom_bar(stat = "identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab(xlab) + ylab(ylab)

      # sort out custom colour pairing,
      if (!is.null(colkey)) {
        plot <- plot + scale_colour_discrete('') +
          scale_fill_manual(values = colkey)
      }



      plot


    }
    else {

      plot = ggplot(plotobj, aes(x = factor(cluster), y = proportion)) + th +
        geom_boxplot(outlier.shape = NA) +
        geom_point(aes(color = cluster),alpha = 0.5,
                   size = pointSize) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab(xlab) + ylab(ylab)


      if (!is.null(colkey)) {
        plot <- plot +
          scale_color_manual(values = colkey)
      }

      plot


    }

  }

}
