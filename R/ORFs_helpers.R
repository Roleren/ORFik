#' Defines trailers for ORF.
#'
#' Creates GRanges object as a trailer for ORFranges representing ORF,
#' maintaining restrictions of transcriptRanges. Assumes that ORFranges
#' is on the transcriptRanges, strands and seqlevels are in agreement.
#' When lengthOFtrailer is smaller than space left on the transcript than
#' all available space is returned as trailer.
#'
#' It assumes that ORFranges and transcriptRanges are not sorted when on minus strand. Should be like:
#' [200, 600]
#' [50, 100]
#'
#' @param ORFranges GRanges object of your Open Reading Frame.
#' @param transcriptRanges GRanges object of transtript.
#' @param lengthOftrailer Numeric. Default is 10.
#' @return A GRanges object of trailer.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
#'                      ranges = IRanges(start = c(1, 10, 20),
#'                                       end = c(5, 15, 25)),
#'                      strand = Rle(strand(rep("+", 3))))
#' transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
#'                             ranges = IRanges(start = c(1, 10, 20, 30, 40),
#'                                              end = c(5, 15, 25, 35, 45)),
#'                      strand = Rle(strand(rep("+", 5))))
#' define_trailer(ORFranges, transcriptRanges)
#'
define_trailer <- function(ORFranges, transcriptRanges, lengthOftrailer = 200) {

  strands <- runValue(strand(ORFranges))
  leftSpace <- setdiff(transcriptRanges, ORFranges)
  if (strands == "-") {
    leftSpace <- leftSpace[end(leftSpace) <= start(ORFranges)[1]]
  } else {
    leftSpace <- leftSpace[start(leftSpace) >= end(ORFranges)[1]]
  }

  if (sum(width(leftSpace)) <= lengthOftrailer) {
    return(leftSpace)
  } else {
    widths <- if (strands == "-") rev(cumsum(rev(width(leftSpace))) - lengthOftrailer - 1) else cumsum(width(leftSpace)) - lengthOftrailer - 1
    whichExon <- which(widths >= 0)
    whichExon <- whichExon[if (strands == "-") {length(whichExon)} else 1]

    return(c(resize(leftSpace[whichExon],
                    width(leftSpace[whichExon]) - widths[whichExon] - 1,
                    fix = "start"),
             leftSpace[which(widths < 0)]))
  }
}

#' Creates GRangesList from the results of get_all_ORFs_as_GRangesList and
#'  a list of group indeces
#'
#' @param grl GRangesList. A GRangesList of the original sequences that gave the orfs
#' @param result List. A list of the results of finding uorfs
#' list syntax is: result[1] contain grouping indeces, named index
#' result[2] countains two columns of start and stops,  named orf
#' @return A GRangesList of ORFs.
#' @export
#' @importFrom GenomicFeatures mapFromTranscripts
map_to_GRanges <- function(grl, result) {

  if(class(grl) != "GRangesList") stop("Invalid type of grl, must be GRangesList.")
  if(is.null(names(grl))) stop("grl contains no names")
  if(class(result) != "list") stop("Invalid type of result, must be list.")
  if(length(result) != 2)
    stop("Invalid structure of result, must be list with 2 elements",
         "read info for structure")
  # Check that grl is sorted
  grl <- sortPerGroup(grl, equalSort = F)
  # Create GRanges object from result tx ranges
  gr <- GRanges(seqnames = as.character(names(grl[result$index])),
               ranges = IRanges(start = unlist(result$orf[1]),
                                end = unlist(result$orf[2])),
               strand =as.character(phead(strand(grl[result$index]),1L )))
  names(gr) <- names(grl[result$index])

  # map from transcript, remove duplicates, remove hit columns
  # syntax for mapping:-> Seqnames(gr) == names(grl),
  genomicCoordinates <- mapFromTranscripts(x =  gr, transcripts =  grl)
  genomicCoordinates <-
    genomicCoordinates[names(gr[genomicCoordinates$xHits]) == names(genomicCoordinates)]
  rm(gr)

  genomicCoordinates$xHits <- NULL
  genomicCoordinates$transcriptsHits <- NULL
  names(genomicCoordinates) <- names(grl[result$index])

  newGRL <- split(genomicCoordinates,result$index )
  names(newGRL) <- unique(names(genomicCoordinates))

  # Split by exons and create new exon names
  newGRL <- GrangesSplitByExonSkeleton(newGRL, grl)

  return(newGRL)
}

#' Resizes down ORF to the desired length, removing inside. Preserves exons.
#'
#' @export
#' @param grangesObj A GRanges object of ORF.
#' @param orf_goal_length numeric. Desired length of ORF.
#' @return GRanges object of resized ORF
#'
resize_ORF <- function(grangesObj, orf_goal_length) {

  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("biovizBase needed for this function to work. Please install it.", call. = FALSE)
  }

  length_diff <- (sum(width(grangesObj))/3 - orf_goal_length) * 3
  is_even <- (length_diff %% 2) == 0
  left_diff <- if (is_even) { length_diff/2 } else { (length_diff - 1)/2 }
  right_diff <- if (is_even) { length_diff/2 } else { (length_diff + 1)/2 }

  tiled <- tile(grangesObj, width = 1)
  tiled <- biovizBase::flatGrl(tiled)
  middle <- floor(length(tiled)/2)

  tiled <- tiled[-c((middle - left_diff):(middle + right_diff - 1))]
  tiled <- reduce(tiled)

  return(tiled)
}
