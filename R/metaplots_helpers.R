

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
