#' @rdname mdsplot
#' @title Create sample-level MDS plot
#'
#' @description Multi-Dimensional Scaling (MDS) plot computed across
#' median marker expression in each sample or metadata grouping.
#'
#' @param data a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param grouping a character string representaing which \code{metadata}
#' slot to group the individual samples by.
#' @param feature a character string represnting a feature to label the heatmap with.
#' @param markers a character string specifying which markers to include.
#' @param colkey a vector of grouping ids the colour strings to indicate the feature labelling.
#' @param legendPosition a string for \code{ggplot2} legend position.
#' @param legendLabSize \code{ggplot2} legend label size
#' @param legendIconSize \code{ggplot2} legend icon size
#' @param xlim \code{ggplot2} x limit specification
#' @param ylim \code{ggplot2} y limit specification
#'
#' @author James Opzoomer \email{james.opzoomer@kcl.ac.uk}
#'
#' @return a \code{ggplot} object.
#'
#' @examples
#' # Download complete ImmunoCluster SCE object from zenodo
#' sce_gvhd = readRDS(url("https://zenodo.org/record/3801882/files/sce_gvhd.rds"))
#'
#' # Generate sample level MDS plot with all markers
#' mdsplot(sce_gvhd, feature = "condition", colkey = c(None = 'royalblue', GvHD = 'red2'))
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom limma plotMDS
#' @export

mdsplot <- function(
  data,
  grouping = 'group',
  feature = NULL,
  markers = NULL,
  colkey = NULL,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 5.0,
  xlim = NULL,
  ylim = NULL,
  celllab = NULL,
  labSize = 3.0,
  labhjust = 1.5,
  labvjust = 0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  xlab = dimColnames[1],
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = dimColnames[2],
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  title = 'Metadata plot',
  subtitle = '',
  caption = 'MDS Plot',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  borderWidth = 0.8,
  borderColour = 'black'
  )
  {


  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background=element_rect(),

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
      axis.title=element_text(size=axisLabSize),

      legend.position=legendPosition,
      legend.key=element_blank(),
      legend.key.size=unit(0.5, 'cm'),
      legend.text=element_text(size=legendLabSize),

      title=element_text(size=legendLabSize),
      legend.title=element_blank())

  if(is.null(markers) == TRUE){

    markers = rownames(data)

  }

  # Get the median marker expression per sample without normalization
  # add marker selection later
  data.median = data.frame(sample_id = data@metadata[,grouping], t(assay(data))[,markers]) %>%
                dplyr::group_by(sample_id) %>%
                dplyr::summarize_all(list(median))

  data.median.sample = t(data.median[, -1])

  colnames(data.median.sample) = data.median$sample_id
  mds = limma::plotMDS(data.median.sample, plot = FALSE)
  mdsdf = data.frame(MDS1 = mds$x, MDS2 = mds$y,
                    sample_id = colnames(data.median.sample))

  # Using metadata feature slot
  if(is.null(feature) == FALSE){

    # Summarize feature
    feature.summary = data.frame(sample_id = data@metadata[,grouping], feature = data@metadata[feature]) %>%
                      dplyr::distinct()

    # Match feature
    match_idx = match(mdsdf$sample_id, feature.summary$sample_id)
    mdsdf[,feature] = as.factor(feature.summary[,feature][match_idx])

    # Plot
    mds = ggplot(mdsdf, aes_string(x = "MDS1", y = "MDS2", color = feature)) +
    geom_point(size = 2, alpha = 0.8) + th +
    geom_label_repel(aes(label = sample_id))

    # sort out custom colour pairing,
    if (!is.null(colkey)) {
      mds <- mds + scale_colour_discrete('') +
        scale_color_manual(values = colkey)
    }


    return(mds)

  } else{

    # Plot
    mds = ggplot(mdsdf, aes(x = MDS1, y = MDS2)) +
        geom_point(size = 2, alpha = 0.8) + th +
        geom_label_repel(aes(label = sample_id))

    return(mds)

  }

  }
