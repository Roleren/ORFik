#' Defines trailers for ORF.
#'
#' Creates GRanges object as a trailer for ORFranges representing ORF,
#' maintaining restrictions of transcriptRanges. Assumes that ORFranges
#' is on the transcriptRanges, strands and seqlevels are in agreement.
#' When lengthOFtrailer is smaller than space left on the transcript than
#' all available space is returned as trailer.
#'
#' It assumes that ORFranges and transcriptRanges are not
#' sorted when on minus strand. Should be like:
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
#' @importFrom S4Vectors runValue
#' @examples
#' ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
#'                      ranges = IRanges(start = c(1, 10, 20),
#'                                       end = c(5, 15, 25)),
#'                      strand = "+")
#' transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
#'                             ranges = IRanges(start = c(1, 10, 20, 30, 40),
#'                                              end = c(5, 15, 25, 35, 45)),
#'                      strand = "+")
#' defineTrailer(ORFranges, transcriptRanges)
#'
defineTrailer <- function(ORFranges, transcriptRanges, lengthOftrailer = 200) {

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
    widths <- if (strands == "-") rev(cumsum(rev(width(leftSpace)))
                          - lengthOftrailer - 1) else cumsum(width(leftSpace))
                                        - lengthOftrailer - 1
    whichExon <- which(widths >= 0)
    whichExon <- whichExon[if (strands == "-") {length(whichExon)} else 1]

    return(c(resize(leftSpace[whichExon],
                    width(leftSpace[whichExon]) - widths[whichExon] - 1,
                    fix = "start"),
             leftSpace[which(widths < 0)]))
  }
}



#' Map orfs to genomic coordinates
#'
#' Creates GRangesList from the results of ORFs_as_List and
#'  the GRangesList used to find the ORFs
#' @param grl A \code{\link[GenomicRanges]{GRangesList}} of the original
#'  sequences that gave the orfs
#' @param result List. A list of the results of finding uorfs
#' list syntax is: result[1] contain grouping indeces, named index
#' result[2] countains two columns of start and stops,  named orf
#' @return A \code{\link[GenomicRanges]{GRangesList}} of ORFs.
#' @export
#' @importFrom GenomicFeatures pmapFromTranscripts
mapToGRanges <- function(grl, result) {

  if (class(grl) != "GRangesList") stop("Invalid type of grl,",
                                        "must be GRangesList.")
  if (is.null(names(grl))) stop("grl contains no names")
  if (class(result) != "list") stop("Invalid type of result, must be list.")
  if (length(result) != 2)
    stop("Invalid structure of result, must be list with 2 elements",
         "read info for structure")
  # Check that grl is sorted
  grl <- sortPerGroup(grl, ignore.strand = T)
  # Create Ranges object from orf scanner result
  ranges = IRanges(start = unlist(result$orf[1], use.names = FALSE),
                   end = unlist(result$orf[2], use.names = FALSE))

  # map transcripts to genomic coordinates, reduce away false hits
  genomicCoordinates <- pmapFromTranscripts(x = ranges,
                                            transcripts = grl[result$index])
  genomicCoordinates <- reduce(genomicCoordinates, drop.empty.ranges = TRUE)

  return(makeORFNames(genomicCoordinates))
}

#' Resizes down ORF to the desired length, removing inside. Preserves exons.
#'
#' @export
#' @param grangesObj A GRanges object of ORF.
#' @param orf_goal_length numeric. Desired length of ORF.
#' @return GRanges object of resized ORF
resizeORF <- function(grangesObj, orf_goal_length) {

  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("biovizBase needed for this function to work.
         Please install it.", call. = FALSE)
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

#' Get transcript names from orf names
#'
#' names must either be a column called names, or the names of the
#' grl object
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} grouped by ORF
#'  or GRanges object
#' @param unique a boolean, if true unique the names,
#'  used if several orfs map to same transcript and you only
#'  want the unique groups
#' @export
#' @return a character vector of transcript names,
#'  without _* naming
OrfToTxNames <- function(grl, unique = F){
  if (!is.gr_or_grl(class(grl))) {
    stop("grl must be GRangesList or GRanges Object")
  }

  if (is.null(names(grl))) {
    if (!is.grl(class(grl))) {
      otherPossibility <- unlist(grl, use.names = F)$names
    } else {
      otherPossibility <- grl$names
    }

    if (is.null(otherPossibility)) {
      stop("grl have no valid orf names to convert")
    } else {
      if (unique) {
        return(gsub("_[0-9]*", "", unique(otherPossibility)))
      }
      return(gsub("_[0-9]*", "", otherPossibility))
    }
  }
  if (unique) {
    return(gsub("_[0-9]*", "", unique(names(grl))))
  }
  return(gsub("_[0-9]*", "", names(grl)))
}

#' get the start sites from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGG, get the position of the A.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param asGR a boolean, return as GRanges object
#' @param keep.names if asGR is False, do you still want
#'  to keep a named vector
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @export
#' @return if asGR is False, a vector, if True a GRanges object
ORFStartSites <- function(grl, asGR = FALSE, keep.names = FALSE,
                          is.sorted = FALSE){
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  posIds <- strandBool(grl)


  startSites <- rep(NA, length(grl))
  startSites[posIds] <- firstStartPerGroup(grl[posIds], FALSE)
  startSites[!posIds] <- firstEndPerGroup(grl[!posIds], FALSE)

  if (asGR) {
    gr <- GRanges(seqnames = seqnamesPerGroup(grl, FALSE),
                  ranges = IRanges(startSites,startSites),
                  strand = strandPerGroup(grl, FALSE))
    names(gr) <- names(grl)
    return(gr)
  }

  if (keep.names) {
    names(startSites) <- names(grl)
    return(startSites)
  }
  return(startSites)
}

#' get the Stop sites from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGC, get the position of the C.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param asGR a boolean, return as GRanges object
#' @param keep.names if asGR is False, do you still want
#'  to keep a named vector
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @export
#' @return if asGR is False, a vector, if True a GRanges object
ORFStopSites <- function(grl, asGR = FALSE, keep.names = FALSE,
                         is.sorted = FALSE){
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  posIds <- strandBool(grl)

  stopSites <- rep(NA, length(grl))
  stopSites[posIds] <- lastExonEndPerGroup(grl[posIds], FALSE)
  stopSites[!posIds] <- lastExonStartPerGroup(grl[!posIds], FALSE)

  if (asGR) {
    gr <- GRanges(seqnames = seqnamesPerGroup(grl, FALSE),
                ranges = IRanges(stopSites,stopSites),
                  strand = strandPerGroup(grl, FALSE))
    names(gr) <- names(grl)
    return(gr)
  }

  if (keep.names) {
    names(stopSites) <- names(grl)
    return(stopSites)
  }
  return(stopSites)
}

#' Get the Start codons(3 bases) from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGC, get the positions ATG.
#' It takes care of exons boundaries, with exons < 3 length.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param is.sorted a boolean, a speedup if you know the ranges are sorted
#' @export
#' @return a GRangesList of start codons, since they might be split on exons
ORFStartCodons <- function(grl, is.sorted = FALSE){
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  firstExons <- firstExonPerGroup(grl)
  widths <- widthPerGroup(firstExons)
  validWidths <- widths >= 3L
  if (!all(validWidths)) { # fix short exons by tiling
    needToFix <- grl[!validWidths]
    tileBy1 <- tile1(needToFix)
    fixedStarts <- reduceKeepAttr(phead(tileBy1, 3L), keep.names = TRUE)
    grl[!validWidths] <- fixedStarts
  }
  # fix the others the easy way
  firstExons <- firstExons[validWidths]
  posIds <- strandBool(firstExons)

  end(firstExons[posIds]) <- start(firstExons[posIds]) + 2L
  start(firstExons[!posIds]) <- end(firstExons[!posIds]) - 2L

  grl[validWidths] <- firstExons

  return(grl)
}


#' Get the Stop codons(3 bases) from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGC, get the positions TGC.
#' It takes care of exons boundaries, with exons < 3 length.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param is.sorted a boolean, a speedup if you know the ranges are sorted
#' @export
#' @return a GRangesList of stop codons, since they might be split on exons
ORFStopCodons <- function(grl, is.sorted = FALSE){
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  lastExons <- lastExonPerGroup(grl)
  widths <- widthPerGroup(lastExons)
  validWidths <- widths >= 3L
  if (!all(validWidths)) { # fix short exons by tiling
    needToFix <- grl[!validWidths]
    tileBy1 <- tile1(needToFix)
    fixedStops <- reduceKeepAttr(ptail(tileBy1, 3L), keep.names = TRUE)
    grl[!validWidths] <- fixedStops
  }
  # fix the others the easy way
  lastExons <- lastExons[validWidths]
  posIds <- strandBool(lastExons)

  start(lastExons[posIds]) <- end(lastExons[posIds]) - 2L
  end(lastExons[!posIds]) <- start(lastExons[!posIds]) + 2L

  grl[validWidths] <- lastExons

  return(grl)
}

#' Get id's for orf
#'
#' These id's can be uniqued by isoform etc,
#' this is not supported by GenomicRanges.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param with.tx a boolean, include transcript names,
#'  if you want unique orfs, so that they dont have multiple
#'  versions on different isoforms, set it to FALSE.
#' @importFrom S4Vectors phead
#' @return a character vector of ids, 1 per orf
orfID <- function(grl, with.tx = FALSE){
  seqnames <- as.character(seqnames(phead(grl,1L)))
  strands <- strandPerGroup(grl,F)

  exonInfo <- paste(start(grl),width(grl))
  exonInfo <- paste(exonInfo, sep = '', collapse = ';')
  names(exonInfo) <- NULL

  uorfID <- paste(seqnames, strands, exonInfo, sep = ",")
  if (with.tx) {
    uorfID <- paste(uorfID, OrfToTxNames(grl))
  }
  return(uorfID)
}

#' Get the unique set of orfs
#'
#' Some orfs might be found several times, from different isoforms
#' If you want to have the unique orfs, not seperated by which
#' isoform it came from, use this function.
#'
#' NB. You will lose the transcript and name information, since
#' they no longer map to a transcript, but are now general.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @return a GRangesList of unique orfs
uniqueORFs <- function(grl){
  ids <- orfID(grl)
  grl <- grl[!duplicated(ids)]
  gr <- unlist(grl, use.names = FALSE)
  names(gr) <- NULL
  gr$names <- NULL
  grl <- relist(gr, grl)
  names(grl) <- seq(1, length(grl))
  return(grl)
}
