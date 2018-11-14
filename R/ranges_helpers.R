#' Make a meta column with exon ranks
#'
#' Must be ordered, so that same transcripts are ordered together.
#' @param grl a \code{\link{GRangesList}}
#' @param byTranscript if ORfs are by transcript, check duplicates
#' @return an integer vector of indices for exon ranks
#' @importFrom S4Vectors nrun Rle
#'
makeExonRanks <- function(grl, byTranscript = FALSE) {
  validGRL(class(grl))

  if (byTranscript) {
    oldNames <- names(grl)
    names(grl) <- seq_along(grl)
    l <- width(grl) - width(grl)
    t <- unlist(l + seq.int(1,length(grl)), use.names = FALSE)
  } else {
    l <- Rle(names(unlist(grl, use.names = TRUE)))
    t <- unlist(lapply(seq(nrun(l)), function(x) {
      rep(x, runLength(l)[x])
    }))
  }
  if (length(t) == 1) {
    return(1)
  }
  Inds <- rep.int(1, length(t))
  if (!byTranscript) {
    for (x in seq.int(2, length(t))) {
      if (t[x] == t[x - 1]) {
        Inds[x] <- Inds[x - 1] + 1
      }
    }
  } else {
    for (x in seq.int(2, length(t))) {
      if (t[x] != t[x - 1]) {
        if (oldNames[t[x]] == oldNames[t[x] - 1]) {
          Inds[x] <- Inds[x - 1] + 1
        }
      } else{
        Inds[x] <- Inds[x - 1]
      }
    }
  }
  return(Inds)
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


#' Tile a GRangesList by 1
#'
#' Will tile a GRangesList into single bp resolution, each group of the list
#' will be splited by positions of 1. Returned values are sorted.
#' This is not supported originally by GenomicRanges.
#' As a precaution, this function requires the unlisted objects to
#' have names.
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
  ORFs <- unlistGrl(grl)
  if (is.null(names(ORFs))) {
    stop("grl does not have names.")
  }
  ## Try to catch dangerous groupings, that is equally named groups
  if (sum(duplicated(names(ORFs)))) {
    if (!is.null(ORFs$names)) {
      names(ORFs) <- ORFs$names
    } else{
      nametest <- unlist(grl, use.names = TRUE)
      dups <- sum(duplicated(names(nametest)))
      if (dups != sum(duplicated(names(ORFs)))) {
        stop("duplicated ORF names,\n
                need a column called 'names' that are unique,\n
                or change names of groups so they are unique")
      }
      mcols(ORFs) <- DataFrame(row.names = names(ORFs),  mcols(ORFs),
                               names = names(ORFs))
    }
  }
  # special case for only single grouped GRangesList
  # This wastes a lot of memory for big lists!
  if (is.null(ORFs$names)) {
    if(length(ORFs) != length(grl)) {
      stop("wrong naming, could not find unique names")
    }
    mcols(ORFs) <- DataFrame(row.names = names(ORFs),  mcols(ORFs),
                             names = names(grl))
  }

  tilex <- tile(ORFs, width =  1L)
  names(tilex) <- ORFs$names
  if (matchNaming) {
    tilex <- matchNaming(tilex, grl[1])
  } else {
    tilex <- groupGRangesBy(unlistGrl(tilex))
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

#' Faster more secure version of mapFromTranscripts
#'
#' Fixes a bug in function, and should have 10x speedup
#' Also removes hit column for you
#'
#' @param ranges IRanges of ranges within grl
#' @param grl the "transcripts" that contain ranges, GRangesList
#' @param indices integer vector of which index of grl ranges are from:
#' (c(1,1,2)) means first two ranges are from grl[1], third from grl[2])
#' @return A GrangesList of ranges mapped from transcripts
#' @family ExtendGenomicRanges
#'
pmapFromTranscriptF <- function(ranges, grl, indices) {
  if (length(ranges) != length(indices)) stop("length of ranges != indices")
  names <- names(grl)
  names(grl) <- NULL
  genomicCoordinates <- pmapFromTranscripts(x = ranges,
                                            transcripts = grl[indices])
  names(genomicCoordinates) <- names[indices]

  genomicCoordinates <- genomicCoordinates[width(genomicCoordinates) > 0]
  a <- unlistGrl(genomicCoordinates)
  a$hit <- NULL
  return(relist(a, genomicCoordinates))
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
  if(!(is(faFile, "FaFile") || is(faFile, "BSgenome")))
    stop("only FaFile and BSgenome is valid input class")
  if(!is.sorted) grl <- sortPerGroup(grl)
  return(extractTranscriptSeqs(faFile, transcripts = grl))
}

#' Get window region of tx around point of gr
#'
#' If downstreamFrom is 20, it means the window will start -20 downstream of
#' gr start site. It will keep exon structure of tx
#' @param gr a GRanges object (startSites and others, must be single point)
#' @param tx a GRangesList of transcripts or (container region)
#' @param downstream an integer, relative region to get downstream from
#' @param upstream an integer vector, relative region to get upstream from.
#' @return a GRanges/GRangesList object if exon/introns
#' @family ExtendGenomicRanges
#' @importFrom data.table chmatch
#'
windowPerGroup <- function(gr, tx, downstream = 0L, upstream = 0L) {
  g <- asTX(gr, tx)

  start(g) <- pmax(start(g) - downstream, 1L)
  indices <- chmatch(txNames(gr), names(tx))
  if (upstream != 0) {
    end(g) <- pmin(end(g) + upstream, widthPerGroup(tx[indices], FALSE))
  }

  return(pmapFromTranscriptF(g, tx, indices))
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
extendLeaders <- function(grl, extension = 1000, cds = NULL) {
  if (is(extension, "numeric") && length(extension) == 1) {
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
    return(coverageByWindow(reads, grl, is.sorted = is.sorted,
                            keep.names = keep.names))
  }

  unlTile <- unlistGrl(tile1(grl, matchNaming = FALSE))
  if (!is.null(unlTile$names)) { # for orf case
    names(unlTile) <- unlTile$names
  }
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
#' GRangesList. It also does not loose sorting for GRangesList,
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
