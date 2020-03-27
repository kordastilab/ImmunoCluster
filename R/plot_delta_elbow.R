#' @export
triangle <- function(m) {
  n <- ncol(m)
  nm <- matrix(0, ncol=n, nrow=n)
  fm <- m
  nm[upper.tri(nm)] <- m[upper.tri(m)]
  fm <- t(nm) + nm
  diag(fm) <-  diag(m)
  nm <- fm
  nm[upper.tri(nm)] <- NA
  diag(nm) <- NA
  m[lower.tri(nm)]
}

plot_delta_elbow <- function(mc) {
  # empirical CDF distribution
  maxK <- length(mc)
  v <- lapply(mc[seq_len(maxK)[-1]], function(x) triangle(x$ml))
  h <- lapply(v, function(x) {
    h <- graphics::hist(x, breaks=seq(0, 1, .01), plot=FALSE)
    h$counts <- cumsum(h$counts) / sum(h$counts)
    return(h)
  })

  # calculate area under CDF curve, by histogram method &
  # calculate proportional increase relative to prior k
  areaK <- vapply(h, function(x) cumsum(x$counts * .01)[100], numeric(1))
  deltaK <- c(areaK[1], diff(areaK) / areaK[seq_len(maxK-2)])

  df <- data.frame(k=seq_len(maxK)[-1], y=deltaK)
  y_max <- ceiling(max(df$y)*2)/2

  dlta = ggplot(df, aes_string(x="k", y="y")) +
    theme_classic() + geom_line(color="black") +
    geom_point(size=1, colour = "grey90") +
    scale_x_continuous(breaks=seq(2, maxK, 2), expand=c(0,.5)) +
    scale_y_continuous(limits=c(0, y_max), expand=c(0,.125),
                       breaks=function(x) seq(x[1]+.125, x[2], .5)) +
    ylab("Relative change\nin area under CDF curve") +
    theme(plot.title=element_text(face="bold"),
          axis.text=element_text(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  return(dlta)
}
