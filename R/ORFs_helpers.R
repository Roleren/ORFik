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


#' Creates list of IRanges with Open Reading Frames.
#'
#' @param ORFdef List of IRanges objects representing found ORFs.
#' @param grangesObj GRanges object to map coordinates back from each IRanges ORF.
#' @param transcriptName String of name that will be added as metadata name column.
#' @return A GRanges objects of ORFs mapped to grangesObj.
#' @export
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' #map_to_GRanges() #rewrite into C++
#'
map_to_GRanges <- function(ORFdef, grangesObj, transcriptName = "") {
  iranges = matrix(data = c(start(ORFdef),end(ORFdef)),ncol = 2)
  txranges = matrix(data = c(start(grangesObj),end(grangesObj)),ncol = 2)
  txstrings = matrix(data = c(seqnames(grangesObj),strand(grangesObj)),ncol = 2)
  map_to_GRangesC(GRanges,IRanges,iranges,txranges,txstrings,transcriptName)
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
