#' Creates GRanges object as trailer in ORF.
#'
#' @param ORFranges GRanges object of your Open Reading Frame.
#' @param transcriptRanges GRanges object of transtript.
#' @param lengthOftrailer Numeric. Default is 10.
#' @return A GRanges object of trailer.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' #' #defineTrailer()

define_trailer <- function(ORFranges, transcriptRanges, lengthOftrailer = 10) {
    if (runValue(strand(ORFranges) == "+")) {
        ORFranges <- ORFranges[1]
    } else {
        ORFranges <- ORFranges[length(ORFranges)]
    }
    for (i in lengthOftrailer:1) {
        d <- flank(ORFranges, i, start = F)
        if (overlapsAny(d, transcriptRanges)) {
            return(d)
        }
    }
    return(GRanges())
}


#' Creates list of IRanges with Open Reading Frames.
#'
#' @param fastaSeq DNA sequence to search for Open Reading Frames.
#' @param startCodon Default is 'ATG'.
#' @param stopCodon Default is c('TAA', 'TAG', 'TGA').
#' @param longestORF bolean. Default TRUE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least START + 8*3 [bp] + STOP.
#' @return A List of IRanges objects of ORFs.
#' @export
#' @import IRanges
#' @examples
#' #find_in_frame_ORFs()

find_in_frame_ORFs <- function(fastaSeq, startCodon = c("ATG"), stopCodon = c("TAA", "TAG", "TGA"), longestORF = T, minimumLength = 0) {
    # for specified startCodon it finds all in frame open reading frames on the input sequence
    subchar <- function(string, pos) {
        for (i in pos) {
            string <- sub(paste("^(.{", i - 1, "}).", sep = ""), "\\1N", string, perl = TRUE)
        }
        string
    }
    stop_cod <- paste(stopCodon, collapse = "|")
    start_cod <- paste(startCodon, collapse = "|")

    codpos <- paste("(?:", start_cod, ")(?:[ATGCN]{3}(?<!", stop_cod, ")){", minimumLength, ",}(?:", stop_cod, ")", sep = "")

    frame <- gregexpr(codpos, fastaSeq, perl = T)
    frmat <- lapply(frame, function(x) if (x[1] != -1)
        cbind((as.integer(x)), (attr(x, "match.length"))) else matrix(NA, ncol = 2, nrow = 0))
    gr <- IRanges(start = do.call("rbind", frmat)[, 1], width = do.call("rbind", frmat)[, 2])

    while (frame[[1]][1] != -1 && !longestORF) {
        starts <- c(sapply(frame, as.integer))
        fastaSeq <- subchar(fastaSeq, starts)
        frame <- gregexpr(codpos, fastaSeq, perl = T)
        frmat <- lapply(frame, function(x) if (x[1] != -1)
            cbind((as.integer(x)), (attr(x, "match.length"))) else matrix(NA, ncol = 2, nrow = 0))
        frGR <- IRanges(start = do.call("rbind", frmat)[, 1], width = do.call("rbind", frmat)[, 2])

        gr <- c(gr, frGR)
    }

    gr <- gr[order(start(gr), end(gr))]
    gr
}


#' Creates list of IRanges with Open Reading Frames.
#'
#' @param ORFdef - list of IRanges objects representing found ORFs.
#' @param grangesObj - GRanges object to map coordinates back from each IRanges ORF.
#' @return A List of GRanges objects of ORFs mapped to grangesObj.
#' @export
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' #map_to_GRanges()

map_to_GRanges <- function(ORFdef, grangesObj) {
    returnList <- c()

    if (length(ORFdef) == 0) {
        return(returnList)
    }

    for (i in 1:length(ORFdef)) {
        ORFranges <- IRanges()
        startingPos <- start(ORFdef[i]) - 1
        endingPos <- end(ORFdef[i]) - 1

        j = 1
        while (startingPos > width(grangesObj)[j]) {
            startingPos <- startingPos - width(grangesObj)[j]
            endingPos <- endingPos - width(grangesObj)[j]
            j = j + 1
        }

        while (endingPos > width(grangesObj)[j]) {
            ORFranges <- c(ORFranges, IRanges(start = start(grangesObj)[j] + startingPos, end = end(grangesObj)[j]))
            startingPos <- 0
            endingPos <- endingPos - width(grangesObj)[j]
            j = j + 1
        }
        ORFranges <- c(ORFranges, IRanges(start = start(grangesObj)[j] + startingPos, end = start(grangesObj)[j] + endingPos))

        chrom <- unique(as.character(seqnames(grangesObj)))
        strands <- unique(as.character(strand(grangesObj)))
        if (length(chrom) != 1) {
            stop("Different chromosomes in GRanges object")
        }
        if (length(strands) != 1) {
            stop("Different strands in GRanges object")
        }

        ORFranges <- GRanges(seqnames = Rle(rep(chrom, length(ORFranges))), ranges = ORFranges, strand = Rle(strand(rep(strands,
            length(ORFranges)))))
        GenomeInfoDb::seqlevels(ORFranges) <- GenomeInfoDb::seqlevels(grangesObj)
        returnList <- c(returnList, ORFranges)
    }

    returnList
}

#' Resizes down ORF to the desired length, removing inside. Preserves exons.
#'
#' @export
#' @param grangesObj A GRanges object of ORF.
#' @param orf_goal_length numeric. Desired length of ORF.
#' @return GRanges object of resized ORF
resize_ORF <- function(grangesObj, orf_goal_length){
  length_diff <- (sum(width(grangesObj))/3 - orf_goal_length) * 3
  is_even <- (length_diff %% 2) == 0
  left_diff <- ifelse(is_even, length_diff/2, (length_diff - 1)/2)
  right_diff <- ifelse(is_even, length_diff/2, (length_diff + 1)/2)

  tiled <- tile(grangesObj, width = 1)
  tiled <- flatGrl(tiled)
  middle <- floor(length(tiled)/2)

  tiled <- tiled[-c((middle - left_diff):(middle + right_diff - 1))]
  tiled <- reduce(tiled)

  return(tiled)
}
