#' @export
markerExpressionPerSample <- function(
  indata,
  assay = 'scaled',
  grouping = 'group', # The condensing function
  clusters = NULL, # Which clusters
  feature = 'condition', # The contrast
  clusterAssign = 'cell_annotation', # The clustering
  markers = rownames(indata),
  feature_cols = NULL, # Colors
  ncol = 5,
  nrow = 2,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 5.0,
  legendKeyHeight = 2.5,
  xlim = NULL,
  ylim = NULL,
  yfixed = FALSE,
  xlab = 'Marker',
  xlabAngle = 90,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = 'Expression',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  stripLabSize = 16,
  title = 'Marker expression per sample',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  borderWidth = 0.8,
  borderColour = 'black'){

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

      legend.title = element_blank(),
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = legendLabSize),
      legend.key.height = unit(legendKeyHeight, 'cm'),

      strip.text.x = element_text(size = stripLabSize,
                                  face = 'bold', margin = margin(b = 5, t = 5)))

    idx <- which(rownames(indata) %in% markers)

    plotobj <- data.frame(Cluster = sct@metadata[,clusterAssign],
                          Grouping = sct@metadata[,grouping],
                          as.data.frame(t(as.matrix(assay(indata, assay)[idx,]))))

    plotobj <- plotobj[which(plotobj$Cluster %in% clusters),]

    data.median = plotobj %>% group_by( Grouping, Cluster) %>%
      summarize_all(list(median))

    ## Melt
    expr_median_sample_cluster_melt = melt(data.median,
                                            id.vars = c("Grouping", "Cluster"), value.name = "median_expression",
                                            variable.name = "antigen")

    plotobj = expr_median_sample_cluster_melt


    if(is.null(feature) == FALSE){ # No contrast

      # Summarize sampleID to metadata feature
      feature.summary = data.frame(Grouping = sct@metadata[,grouping], feat = sct@metadata[,feature]) %>%
        distinct()

      mm = match(plotobj$Grouping, feature.summary$Grouping)

      plotobj$Feature <- factor(feature.summary$feat[mm])

      # initialise the plot object
      plot <- ggplot(plotobj, aes(x = antigen, y = median_expression,
                                  color = Feature, fill = Feature)) + th +

        guides(
          fill = guide_legend(),
          shape = guide_legend(),
          alpha = FALSE)

      plot <- plot + geom_boxplot(
        position = position_dodge(),
        alpha = 0.5,
        outlier.color = NA) +
        geom_point(alpha = 0.8, position = position_jitterdodge(),
                   size = 0.7)

      if (yfixed == TRUE) {
        plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol)
      } else {
        plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol,
                                   scales = 'free_y')
      }

      # add elements to the plot for xy labeling and axis limits
      plot <- plot + xlab(xlab) + ylab(ylab)
      if (!is.null(xlim)) {
        plot <- plot + xlim(xlim[1], xlim[2])
      }
      if (!is.null(ylim)) {
        plot <- plot + ylim(ylim[1], ylim[2])
      }

      # add elements to the plot for title, subtitle, caption
      plot <- plot + labs(title = title,
                          subtitle = subtitle, caption = caption)

      # border around plot
      plot <- plot +
        theme(panel.border = element_rect(
          colour = borderColour,
          fill = NA,
          size = borderWidth))

      return(plot)


    }
    else{

      # initialise the plot object
      plot <- ggplot(plotobj, aes(x = antigen, y = median_expression)) + th +

        guides(
          fill = guide_legend(),
          shape = guide_legend(),
          alpha = FALSE)

      plot <- plot + geom_boxplot(
        position = position_dodge(),
        alpha = 0.5,
        outlier.color = NA) +
        geom_point(alpha = 0.8, position = position_jitterdodge(),
                   size = 0.7)

      if (yfixed == TRUE) {
        plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol)
      } else {
        plot <- plot + facet_wrap( ~ Cluster, nrow = nrow, ncol = ncol,
                                   scales = 'free_y')
      }

      # add elements to the plot for xy labeling and axis limits
      plot <- plot + xlab(xlab) + ylab(ylab)
      if (!is.null(xlim)) {
        plot <- plot + xlim(xlim[1], xlim[2])
      }
      if (!is.null(ylim)) {
        plot <- plot + ylim(ylim[1], ylim[2])
      }

      # add elements to the plot for title, subtitle, caption
      plot <- plot + labs(title = title,
                          subtitle = subtitle, caption = caption)

      # border around plot
      plot <- plot +
        theme(panel.border = element_rect(
          colour = borderColour,
          fill = NA,
          size = borderWidth))

      return(plot)

    }

}
