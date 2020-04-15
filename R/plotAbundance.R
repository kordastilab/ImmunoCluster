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
  yfixed = FALSE,
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

  ## Create the props table
  counts_table <- table(indata@metadata[,clusterAssign], indata@metadata[,group])

  props_table <- t(t(counts_table) / colSums(counts_table)) * 100
  counts <- as.data.frame.matrix(counts_table)
  props <- as.data.frame.matrix(props_table[clusters,])

  plotobj <- melt(data.frame(cluster = rownames(props), props),
               id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
  plotobj$cluster <- factor(plotobj$cluster)
  plotobj$sample_id = gsub("X", "", plotobj$sample_id)

  if(is.null(feature) == FALSE){ # No contrast

  # Summarize sampleID to metadata feature
  feature.summary = data.frame(Grouping = indata@metadata[,group], feat = indata@metadata[,feature]) %>%
    distinct()

  ## Add condition info
  mm <- match(plotobj$sample_id, feature.summary$Grouping)

  # Need to recover the original matrix
  plotobj$condition <- factor(feature.summary$feat[mm])

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
