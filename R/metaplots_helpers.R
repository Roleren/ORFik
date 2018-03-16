#' Plot the periods of reads around cds start codon
#'
#' @param start_df a data.frame with the periods
#' @return NULL
#' @import ggplot2
plot_periodic_lengths_start <- function(start_df) {
  colnames(start_df) <- c("start_codon", "count", "flength")
  start_df$fill <- factor(rep(c(0, 1, 1), 10))
  start_df$start_codon <- as.character(start_df$start_codon)
  start_df$start_codon <- factor(start_df$start_codon, levels = start_df$start_codon[1:60])
  #avoid check warnings
  start_codon <- NULL
  count <- NULL
  fill <- NULL
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  p <- ggplot(start_df, aes(x=start_codon, y=count, fill=fill)) + geom_bar(stat="identity") +
    scale_fill_manual(values=c("#E72B3C", "#A9A9A9")) + facet_grid( ~ flength) + theme(legend.position="none")
  rect <- data.frame(xmin=30.5, xmax=33.5, ymin=-Inf, ymax=Inf) # rectangle
  p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                     fill="green", alpha=0.25, inherit.aes = FALSE)
  ll <- levels(start_df$start_codon)[c(T, F, F)] # breaks and labels
  p <- p + scale_x_discrete(breaks=ll, labels = ll)
  print(p)
  return(NULL)
}

#' Plot the periods of reads around cds stop codon
#'
#' @param stop_df a data.frame with the periods
#' @return NULL
#' @import ggplot2
plot_periodic_lengths_stop <- function(stop_df) {
  colnames(stop_df) <- c("stop_codon", "count", "flength")
  stop_df$fill <- factor(rep(c(0, 1, 1), 10))
  stop_df$stop_codon <- as.character(stop_df$stop_codon)
  stop_df$stop_codon <- factor(stop_df$stop_codon, levels = stop_df$stop_codon[1:60])
  #avoid check warnings
  stop_codon <- NULL
  count <- NULL
  fill <- NULL
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  p <- ggplot(stop_df, aes(x=stop_codon, y=count, fill=fill)) + geom_bar(stat="identity") +
    scale_fill_manual(values=c("#E72B3C", "#A9A9A9")) + facet_grid( ~ flength) + theme(legend.position="none")
  rect <- data.frame(xmin=30.5, xmax=33.5, ymin=-Inf, ymax=Inf) # rectangle
  p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                     fill="red", alpha=0.25, inherit.aes = FALSE)
  ll <- levels(stop_df$stop_codon)[c(T, F, F)] # breaks and labels
  p <- p + scale_x_discrete(breaks=ll, labels = ll)
  print(p)
  return(NULL)
}

#' Calculate metaplot coverage of reads around input GRanges object.
#'
#' It should create window of window_size around input ranges eg.
#' for granges starting at 100-100: and window_size = 3
#' 97-103
#' @param cdsRanges GRanges object of your CDSs start or stop postions.
#' @param AllRiboReadsRangesResized GRanges object of your reads.
#' You should resize them beforehand to width of 1 to focus on 5' ends of footprints.
#' @return A data.frame with frequencies (Freq) of reads mapped to
#' positions (Position) specified in cdsRanges along with frame (Frame).
#' @export
#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @examples
#' #metaPosition()
#'
metaPosition <- function(cdsRanges, AllRiboReadsRangesResized) {

    window_size <- unique(width(cdsRanges))
    if (length(window_size) != 1) {
        stop("All input GRanges should have the same window i.e. same width()")
    }
    window_size <- (window_size - 1)/2

    # disable warnings
    oldw <- getOption("warn")
    options(warn = -1)

    hits <- findOverlaps(cdsRanges, AllRiboReadsRangesResized)

    # enable warnings
    options(warn = oldw)

    RiboHitsStarts <- start(ranges(AllRiboReadsRangesResized[subjectHits(hits)]))
    CDSHitsStarts <- start(ranges(cdsRanges[queryHits(hits)]))

    position <- RiboHitsStarts - CDSHitsStarts

    # change sign for minus strands
    minusStrands <- as.vector(strand(AllRiboReadsRangesResized[subjectHits(hits)]) == "-")
    position[minusStrands] <- window_size * 2 - position[minusStrands]

    position <- table(as.numeric(position))

    hitMap <- data.frame(rep(0, window_size * 2 + 1), -window_size:window_size)
    names(hitMap) <- c("Freq", "Position")
    hitMap$Freq[as.numeric(names(position)) + 1] <- as.numeric(position)
    hitMap$Frame <- rep(1:3, window_size)[1:(window_size * 2 + 1)]

    return(hitMap)
}
