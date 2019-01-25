#' Make grouping for exon structures.
#'
#' Either by transcript or by original groupings.
#' Must be ordered, so that same transcripts are ordered together.
#' @param grl a \code{\link{GRangesList}}
#' @param byTranscript if ORfs are by transcript, check duplicates
#' @return an integer vector of indices for exon ranks
#' @importFrom S4Vectors nrun Rle
#'
makeExonRanks <- function(grl, byTranscript = FALSE) {
  validGRL(class(grl))
  g <- groupings(grl)

  if (byTranscript) {
    inds <- rep.int(1L, length(g))
    oldNames <- names(grl)
    oldNames <- data.table::chmatch(oldNames, oldNames)
    for (x in seq.int(2L, length(g))) {
      if (g[x] != g[x - 1L]) {
        if (oldNames[g[x]] == oldNames[g[x] - 1L]) {
          inds[x] <- inds[x - 1L] + 1L
        }
      } else{
        inds[x] <- inds[x - 1L]
      }
    }
    return(inds)
  }
  return(g)
}


#' Make ORF names per orf
#'
#' grl must be grouped by transcript
#' If a list of orfs are grouped by transcripts, but does not have
#' ORF names, then create them and return the new GRangesList
#' @param grl a \code{\link{GRangesList}}
#' @param groupByTx logical (T), should output GRangesList be grouped by
#' transcripts (T) or by ORFs (F)?
#' @return (GRangesList) with ORF names, grouped by transcripts, sorted.
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' makeORFNames(grl)
makeORFNames <- function(grl, groupByTx = TRUE) {
  ranks <- makeExonRanks(grl, byTranscript = TRUE)
  asGR <- unlistGrl(grl)
  mcols(asGR) <- DataFrame(row.names = names(asGR),
                           names = paste0(names(asGR), "_", ranks))
  if (groupByTx) {
    return(groupGRangesBy(asGR))
  } else {
    return(groupGRangesBy(asGR, asGR$names))
  }
}


#' Tile each GRangesList group to 1-base resolution.
#'
#' Will tile a GRangesList into single bp resolution, each group of the list
#' will be splited by positions of 1. Returned values are sorted as the same
#' groups as the original GRangesList, except they are in bp resolutions.
#' This is not supported originally by GenomicRanges.
#' @param grl a \code{\link{GRangesList}} object with names
#' @param sort.on.return logical (T), should the groups be
#'  sorted before return.
#' @param matchNaming logical (T), should groups keep unlisted names
#'  and meta data.(This make the list very big, for > 100K groups)
#' @return a GRangesList grouped by original group, tiled to 1
#' @importFrom S4Vectors DataFrame
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' gr1 <- GRanges("1", ranges = IRanges(start = c(1, 10, 20),
#'                                      end = c(5, 15, 25)),
#'                strand = "+")
#' gr2 <- GRanges("1", ranges = IRanges(start = c(20, 30, 40),
#'                                      end = c(25, 35, 45)),
#'                strand = "+")
#' names(gr1) = rep("tx1_1", 3)
#' names(gr2) = rep("tx1_2", 3)
#' grl <- GRangesList(tx1_1 = gr1, tx1_2 = gr2)
#' tile1(grl)
#'
tile1 <- function(grl, sort.on.return = TRUE, matchNaming = TRUE) {
  if (!is.grl(grl)) stop("grl must be a GRangesList")
  # optimize by removing all unecessary data
  ORFs <- unlist(grl, use.names = FALSE)
  mcols(ORFs) <- NULL
  names(ORFs) <- NULL
  tilex <- tile(ORFs, width =  1L)
  grouping <- groupings(grl)
  names(tilex) <- grouping

  if (matchNaming) {
    if (!is.null(names(grl))) {
      names(tilex) <- names(grl)[as.integer(grouping)]
    }
    tilex <- matchNaming(tilex, grl[1L])
  } else {
    tilex <- unlistGrl(tilex)
    grouping <- names(tilex)
    if (!is.null(names(grl))) {
      grouping <- names(grl)[as.integer(grouping)]
    }
    names(tilex) <- NULL
    tilex <- groupGRangesBy(tilex, grouping)
  }

  # only negative must be sorted
  if (sort.on.return) {
    posIndices <- strandBool(tilex)
    tilex[!posIndices] <- sortPerGroup(tilex[!posIndices])
  }
  return(tilex)
}


#' Map genomic to transcript coordinates by reference
#' @param grl a \code{\link{GRangesList}} of ranges within
#' the reference, grl must have column called names that gives
#' grouping for result
#' @param reference a GrangesList of ranges that include and are bigger or
#' equal to grl ig. cds is grl and gene can be reference
#' @return a GRangesList in transcript coordinates
#' @family ExtendGenomicRanges
#'
asTX <- function(grl, reference) {
  orfNames <- txNames(grl)
  if (sum(orfNames %in% names(reference)) != length(orfNames)) {
    stop("not all references are present, so can not map to transcripts.")
  }
  reference <- reference[orfNames]
  names(reference) <- NULL
  return(pmapToTranscripts(grl, reference))
}

#' Faster pmapFromTranscript
#'
#' This version tries to fix the shortcommings of GenomicFeature's version.
#' Much faster and uses less memory.
#' Implemented as dynamic program optimized c++ code.
#'
#' The length of x must be the same as length of transcripts. Only exception is
#' if x have integer names like (1, 3, 3, 5), so that x[1] maps to 1, x[2] maps
#' to transcript 3 etc.
#' @param x IRangesList/IRanges/GRanges to map to genomic coordinates
#' @param transcripts a GRangesList to map agains
#' @param removeEmpty a logical, remove non hit exons, else they are set
#'  to 0.
#' @return a GRangesList of mapped reads, names from ranges are kept.
#' @export
#' @examples
#' ranges <- IRanges(start = c( 5, 6), end = c(10, 10))
#' seqnames = rep("chr1", 2)
#' strands = rep("-", 2)
#' grl <- split(GRanges(seqnames, IRanges(c(85, 70), c(89, 82)), strands),
#'              c(1, 1))
#' ranges <- split(ranges, c(1,1)) # both should be mapped to transcript 1
#' pmapFromTranscriptF(ranges, grl, TRUE)
#'
pmapFromTranscriptF <- function(x, transcripts, removeEmpty = FALSE) {
  if (is.null(names(x))) {
    if (length(x) != length(transcripts))
      stop("invalid names of ranges, when length(x) != length(transcripts)")
  }
  if (is(x, "GRanges")) x <- ranges(x)
  if (is(x, "IRangesList")) {
    # Create Ranges object from orf scanner result
    x = unlist(x, use.names = TRUE)
    indices <- strtoi(names(x))
    names(x) <- NULL
  } else if (is(x, "IRanges")) {
    if (is.null(names(x))) {
      indices <- seq.int(1,length(x))
    } else {
      indices <- strtoi(names(x))
    }
    names(x) <- NULL
  } else stop("x must either be IRanges, IRangesList or GRanges")

  if (!is.logical(removeEmpty)) stop("removeEmpty must be logical")
  removeEmpty <- as.logical(removeEmpty)
  if (max(indices) > length(transcripts)) stop("invalid names of IRanges")
  if (length(x) != length(indices)) stop("length of ranges != indices")

  tx <- ranges(transcripts)
  names <- names(tx)
  names(tx) <- NULL
  tx <- tx[indices]
  if (!all(width(x) <= as.integer(sum(width(tx))))) {
    stop("Invalid ranges to map, check them. One is bigger than its reference")
  }
  groupings <- groupings(tx)
  tx <- unlist(tx, use.names = FALSE)
  xStrand <- strandBool(transcripts)[indices]
  txStrand <- xStrand[groupings]

  # forward strand
  if (any(txStrand)) {
    pos <- pmapFromTranscriptsCPP(start(x)[xStrand], end(x)[xStrand],
                                  start(tx)[txStrand], end(tx)[txStrand],
                                  groupings[txStrand], '+', removeEmpty)
  } else {
    pos <- list(ranges = list(vector("integer"), vector("integer")),
                     index = vector("integer"))
  }

  # reverse strand
  if (any(!txStrand)) {
    neg <- pmapFromTranscriptsCPP(start(x)[!xStrand], end(x)[!xStrand],
                                  start(tx)[!txStrand], end(tx)[!txStrand],
                                  groupings[!txStrand], '-', removeEmpty)
  } else {
    neg <- list(ranges = list(vector("integer"), vector("integer")),
                     index = vector("integer"))
  }

  if (removeEmpty) {
    txStrand <- order(c(pos$index, neg$index)) <=  length(pos$index)
  }
  temp <- indices
  nIndices <- vector("integer", length(txStrand))
  nIndices[txStrand] <- pos$index
  nIndices[!txStrand] <- neg$index
  indices <- indices[nIndices]

  xStart <- vector("integer", length(txStrand))
  xStart[txStrand] <- unlist(pos$ranges[1], use.names = FALSE)
  xStart[!txStrand] <- unlist(neg$ranges[1], use.names = FALSE)

  xEnd <- vector("integer", length(txStrand))
  xEnd[txStrand] <- unlist(pos$ranges[2], use.names = FALSE)
  xEnd[!txStrand] <- unlist(neg$ranges[2], use.names = FALSE)

  result <- GRanges(seqnames = seqnamesPerGroup(transcripts, FALSE)[indices],
                    ranges = IRanges(xStart, xEnd),
                    strand = strandPerGroup(transcripts, FALSE)[indices])

  names(result) <- names[indices]
  result <- split(result, nIndices)
  names(result) <- names(transcripts)[temp]
  seqlevels(result) <- seqlevels(transcripts)
  return(result)
}

#' Get transcript sequence from a GrangesList and a faFile or BSgenome
#'
#' A small safety wrapper around \code{\link{extractTranscriptSeqs}}
#' @param grl a GRangesList object
#' @param faFile FaFile or BSgenome used to find the transcripts,
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @export
#' @return a DNAStringSet of the transcript sequences
#' @family ExtendGenomicRanges
#'
txSeqsFromFa <- function(grl, faFile, is.sorted = FALSE) {
  if (!(is(faFile, "FaFile") || is(faFile, "BSgenome")))
    stop("only FaFile and BSgenome is valid input class")
  if (!any(seqlevels(grl) %in% seqlevels(faFile)))
    stop("FaFile had no matching seqlevels to ranges object")
  if (!is.sorted) grl <- sortPerGroup(grl)
  return(extractTranscriptSeqs(faFile, transcripts = grl))
}

#' Get window region of tx around point of gr
#'
#' If downstream is 20, it means the window will start 20 downstream of
#' gr start site (-20 in relative transcript coordinates.)
#' If upstream is 20, it means the window will start 20 upstream of
#' gr start site (+20 in relative transcript coordinates.)
#' It will keep exon structure of tx, so if -20 is on next exon, the previous
#' exon is of course deleted.
#'
#' If a region has a part that goes out of bounds,
#' it will set that side to the bound.
#' If region has no hit in bound, a width 0 GRanges object is returned.
#' This is usefull for things like countOverlaps, since 0 hits will then always
#' be returned for the correct object.
#' @param gr a GRanges object (startSites and others, must be single point)
#' @param tx a GRangesList of transcripts or (container region)
#' @param downstream an integer, relative region to get downstream from
#' @param upstream an integer vector, relative region to get upstream from.
#' @return a GRanges, or GRangesList object if any group had > 1 exon.
#' @export
#' @family ExtendGenomicRanges
#' @importFrom data.table chmatch
#' @examples
#' # find 2nd codon of an ORF on a spliced transcript
#' ORF <- GRanges("1", c(3), "+") # start site
#' names(ORF) <- "tx1_1" # ORF 1 on tx1
#' tx <- GRangesList(tx1 = GRanges("1", c(1,3,5,7,9,11,13), "+"))
#' windowPerGroup(ORF, tx, 5, -3) # <- 2nd codon
#'
windowPerGroup <- function(gr, tx, downstream = 0L, upstream = 0L) {
  g <- asTX(gr, tx)

  starts <- pmax(start(g) - upstream, 1L)
  indices <- chmatch(txNames(gr), names(tx))
  if (downstream != 0L) {
    ends <- pmin(pmax(end(g) + downstream, start(g)-1),
                 widthPerGroup(tx[indices], FALSE))
    ranges(g) <- IRanges(starts, ends)
  } else {
    start(g) <- starts
  }
  g <- ranges(g)
  names(g) <- indices
  return(pmapFromTranscriptF(g, tx, TRUE))
}

#' Extend the leaders transcription start sites.
#'
#' Will extend the leaders or transcripts upstream by extension.
#' Remember the extension is general not relative, that means splicing
#' will not be taken into account.
#' Requires the \code{grl} to be sorted beforehand,
#' use \code{\link{sortPerGroup}} to get sorted grl.
#' @param grl a \code{\link{GRangesList}}
#'  of 5' utrs or transcripts.
#' @param extension an integer, how much to extend the leaders.
#' Or a GRangesList where start / stops by strand are the positions
#' to use as new starts.
#' @param cds If you want to extend 5' leaders downstream, to catch
#' upstream ORFs going into cds, include it. It will add first
#' cds exon to grl matched by names.
#' Do not add for transcripts, as they are already included.
#' @return an extended GRangeslist
#' @export
#' @examples
#' library(GenomicFeatures)
#' samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#' txdb <- loadDb(samplefile)
#' fiveUTRs <- fiveUTRsByTranscript(txdb) # <- extract only 5' leaders
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' cds <- cdsBy(txdb,"tx",use.names = TRUE)
#' ## now try(extend upstream 1000, downstream 1st cds exons):
#' extendLeaders(fiveUTRs, extension = 1000, cds)
#'
#' ## when extending transcripts, don't include cds' of course,
#' ## since they are already there
#' extendLeaders(tx, extension = 1000)
#'
extendLeaders <- function(grl, extension = 1000L, cds = NULL) {
  if (is(extension, "numeric") && length(extension) == 1L) {
    posIndices <- strandBool(grl)
    promo <- promoters(unlist(firstExonPerGroup(grl), use.names = FALSE),
                       upstream = extension)
    newStarts <- rep(NA, length(grl))
    newStarts[posIndices] <- as.integer(start(promo[posIndices]))
    newStarts[!posIndices] <- as.integer(end(promo[!posIndices]))
  } else if (is.grl(class(grl))) {
    starts <- startSites(extension)
    changedGRL <-downstreamOfPerGroup(grl[names(extension)], starts)
    return(changedGRL)
  } else {
    stop("extension must either be an integer, or a GRangesList")
  }

  extendedLeaders <- assignFirstExonsStartSite(grl, newStarts)
  if(is.null(cds)) return (extendedLeaders)
  return(addCdsOnLeaderEnds(extendedLeaders, cds))
}

#' Get overlaps and convert to coverage list
#'
#' @param gr a \code{\link{GRanges}}
#'  of 5' utrs or transcripts.
#' @param reads a GAlignment or GRanges object of RiboSeq, RnaSeq etc
#' @param keep.names logical (T), keep names or not.
#' @param type a string (any), argument for countOverlaps.
#' @return a Rle, one list per group with # of hits per position.
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                                 end = c(5, 15, 25)),
#'                strand = "+")
#' names(ORF) <- "tx1"
#' reads <- GRanges("1", IRanges(25, 25), "+")
#' overlapsToCoverage(ORF, reads)
#'
overlapsToCoverage <- function(gr, reads, keep.names = TRUE, type = "any") {
  # could make this more efficient by counting overlaps
  # only on untiled, then tile the ones that hit and count again
  counts <- countOverlaps(gr, reads, type = type)

  names <- names(counts)
  names(counts) <- NULL
  countList <- split(counts, names)
  if (!keep.names) {
    countList <- IRanges::RleList(countList)
    names(countList) <- NULL
    return(countList)
  }
  return(IRanges::RleList(countList))
}


#' Get coverage per group
#'
#' It tiles each GRangesList group, and finds hits per position
#' @param grl a \code{\link{GRangesList}}
#'  of 5' utrs or transcripts.
#' @param reads a GAlignment or GRanges object of RiboSeq, RnaSeq etc
#' @param is.sorted logical (F), is grl sorted.
#' @param keep.names logical (T), keep names or not.
#' @return a Rle, one list per group with # of hits per position.
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                                 end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25), "+")
#' coveragePerTiling(grl, RFP)
#'
coveragePerTiling <- function(grl, reads, is.sorted = FALSE,
                              keep.names = TRUE) {

  if (length(grl) > 10000) { # faster version for big grl
    if (!is.sorted) grl <- sortPerGroup(grl)
    cov <- coverageByTranscript(reads, grl)
    if (!keep.names) names(cov) <- NULL
    return(cov)
  }

  unlTile <- unlistGrl(tile1(grl, matchNaming = FALSE))
  return(overlapsToCoverage(unlTile, reads, keep.names = keep.names))
}


#' Subset GRanges to get desired frame. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’ @export
#' @param x A tiled to size of 1 GRanges object
#' @param frame A numeric indicating which frame to extract
#' @return GRanges object reduced to only first frame
#'
subset_to_frame <- function(x, frame) {
  if (as.vector(strand(x) == "+")[1]) {
    x[seq.int(frame, length(x), 3)]
  } else {
    x[seq.int(length(x) + 1 - frame, 1, -3)]
  }
}


#' Reduce GRanges / GRangesList
#'
#' Extends function \code{\link{reduce}}
#' by trying to keep names and meta columns, if it is a
#' GRangesList. It also does not lose sorting for GRangesList,
#' since original reduce sorts all by ascending.
#' If keep.names == FALSE, it's just the normal GenomicRanges::reduce
#' with sorting negative strands descending for GRangesList.
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @param drop.empty.ranges (FALSE) if a group is empty (width 0), delete it.
#' @param min.gapwidth (1L) how long gap can it be to say they belong together
#' @param with.revmap (FALSE) return info on which mapped to which
#' @param with.inframe.attrib (FALSE) For internal use.
#' @param ignore.strand (FALSE), can different strands be reduced together.
#' @param keep.names (FALSE) keep the names and meta columns of the GRangesList
#' @return A reduced GRangesList
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 2, 3), end = c(1, 2, 3)),
#'                strand = "+")
#' # For GRanges
#' reduceKeepAttr(ORF, keep.names = TRUE)
#' # For GRangesList
#' grl <- GRangesList(tx1_1 = ORF)
#' reduceKeepAttr(grl, keep.names = TRUE)
#'
reduceKeepAttr <- function(grl, keep.names = FALSE,
                           drop.empty.ranges = FALSE, min.gapwidth = 1L,
                           with.revmap = FALSE, with.inframe.attrib = FALSE,
                           ignore.strand = FALSE) {
  if (!is.gr_or_grl(class(grl))) stop("grl must be GRanges or GRangesList")

  reduced <- GenomicRanges::reduce(grl, drop.empty.ranges, min.gapwidth,
                                   with.revmap, with.inframe.attrib,
                                   ignore.strand)
  if (is.grl(class(grl))) {
    reduced <- sortPerGroup(reduced)
  }
  if (keep.names) { # return with names
    if (!is.grl(class(grl))) {
      return(reduced)
    }
    gr <- unlist(reduced, use.names = TRUE)
    if (length(gr) == 0) return(GRangesList())

    return(matchNaming(gr, grl))
  } else { # return original
    return(reduced)
  }
}
