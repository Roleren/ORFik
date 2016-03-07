#' Creates window around GRanged object.
#'
#' It creates window of window_size around input ranges eg.
#' for GRanges starting at 100-100 and window_size of 3 will give
#' 97-103
#' @param GRanges_obj GRanges object of your CDSs start or stop postions.
#' @param window_size Numeric. Default 30. What size of the window to consider.
#' @return A GRanges object of resized by window_size sites.
#' @export
#' @import GenomicRanges
#' @examples
#' window_resize(GRanges(Rle(c("1"), c(4)),
#'                       IRanges(c(100, 200, 200, 100), width=c(1, 1, 1, 1)),
#'                       Rle(strand(c("+", "+", "-", "-")))),
#'              window_size = 50)

window_resize <- function(GRanges_obj, window_size = 30) {
  GRanges_obj <- promoters(GRanges_obj, upstream = window_size, downstream = window_size  + 1)
  return(GRanges_obj)
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
#' @examples
#' #metaPosition()

meta_position <- function(cdsRanges, AllRiboReadsRangesResized){

  window_size <- unique(width(cdsRanges))
  if(length(window_size) != 1){
    stop("All input GRanges should have the same window i.e. same width()")
  }
  window_size <- (window_size - 1)/2

  #disable warnings
  oldw <- getOption("warn")
  options(warn = -1)

  hits <- findOverlaps(cdsRanges, AllRiboReadsRangesResized)

  #enable warnings
  options(warn = oldw)

  RiboHitsStarts <- start(ranges(AllRiboReadsRangesResized[subjectHits(hits)]))
  CDSHitsStarts <- start(ranges(cdsRanges[queryHits(hits)]))

  position <- RiboHitsStarts - CDSHitsStarts

  #change sign for minus strands
  minusStrands <- as.vector(strand(AllRiboReadsRangesResized[subjectHits(hits)]) == "-")
  position[minusStrands] <- window_size * 2 - position[minusStrands]

  position <- table(as.numeric(position))

  hitMap <- data.frame(rep(0, window_size*2 + 1), -window_size:window_size)
  names(hitMap) <- c("Freq", "Position")
  hitMap$Freq[as.numeric(names(position)) + 1] <- as.numeric(position)
  hitMap$Frame <- rep(1:3, window_size)[1:(window_size*2 + 1)]

  return(hitMap)
}
