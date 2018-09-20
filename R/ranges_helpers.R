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
  }
  l <- Rle(names(unlist(grl, use.names = TRUE)))
  t <- unlist(lapply(seq(nrun(l)), function(x) {
    rep(x, runLength(l)[x])
  }))
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
#' @export
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
asTX <- function(grl, reference) {
  orfNames <- txNames(grl)
  if (sum(orfNames %in% names(reference)) != length(orfNames)) {
    stop("not all references are present, so can not map to transcripts.")
  }
  return(pmapToTranscripts(grl, reference[orfNames]))
}


#' Get transcript sequence from a GrangesList and a faFile or BSgenome
#'
#' A small safety wrapper around GenomicFeatures::extractTranscriptSeqs
#' @param grl a GRangesList object
#' @param faFile FaFile or BSgenome used to find the transcripts,
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @export
#' @return a DNAStringSet of the transcript sequences
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
#'
windowPerGroup <- function(gr, tx, downstream = 0L, upstream = 0L) {
  g <- asTX(gr, tx)

  start(g) <- pmax(start(g) - downstream, 1L)
  if (upstream != 0) {
    end(g) <- pmin(end(g) + upstream, widthPerGroup(tx[names(gr)], FALSE))
  }

  return(pmapFromTranscripts(g, tx[names(gr)]))
}


#' Reassign the start positions of the first exons per group in grl
#' @description make sure your grl is sorted, since start of "-" strand
#' objects should be the
#' max end in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStarts an integer vector of same length as grl, with new start
#' values
#' @return the same GRangesList with new start sites
assignFirstExonsStartSite <- function(grl, newStarts) {
  if (length(grl) != length(newStarts)) stop("length of grl and newStarts ",
                                             "are not equal!")
  posIndices <- strandBool(grl)

  dt <- as.data.table(grl)
  dt[!duplicated(dt$group),]$start[posIndices] <- newStarts[posIndices]
  dt[!duplicated(dt$group),]$end[!posIndices] <- newStarts[!posIndices]

  ngrl <-
    GenomicRanges::makeGRangesListFromDataFrame(dt,
                                                split.field = "group",
                                                names.field = "group_name",
                                                keep.extra.columns = TRUE)
  names(ngrl) <- names(grl)

  return(ngrl)
}


#' Reassign the stop positions of the last exons per group
#' @description make sure your grl is sorted, since stop of "-" strand objects
#' should be the min start in group, use ORFik:::sortPerGroup(grl) to get
#' sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStops an integer vector of same length as grl,
#'  with new start values
#' @return the same GRangesList with new stop sites
#' @importFrom data.table .N .I
#'
assignLastExonsStopSite <- function(grl, newStops) {
  if (length(grl) != length(newStops)) stop("length of grl and newStops ",
                                            "are not equal!")
  posIndices <- strandBool(grl)

  dt <- as.data.table(grl)
  group <- NULL # avoid check warning
  idx = dt[, .I[.N], by = group]
  dt[idx$V1]$end[posIndices] <- newStops[posIndices]
  dt[idx$V1]$start[!posIndices] <- newStops[!posIndices]
  ngrl <-
    GenomicRanges::makeGRangesListFromDataFrame(dt,
                                                split.field = "group",
                                                names.field = "group_name",
                                                keep.extra.columns = TRUE)
  names(ngrl) <- names(grl)

  return(ngrl)
}


#' Get rest of objects downstream (exclusive)
#'
#' Per group get the part downstream of position.
#' downstreamOfPerGroup(tx, stopSites(cds, asGR = TRUE))
#' will return the 3' utrs per transcript as GRangesList,
#' usually used for interesting
#' parts of the transcripts.
#'
#' If you want to include the points given in the region,
#' use downstreamFromPerGroup
#' @param tx a \code{\link{GRangesList}},
#'  usually of Transcripts to be changed
#' @param downstreamOf a vector of integers, for each group in tx, where
#' is the new start point of first valid exon.
#' @return a GRangesList of downstream part
#'
downstreamOfPerGroup <- function(tx, downstreamOf) {
  # Needs speed update!
  posIndices <- strandBool(tx)
  posEnds <- end(tx[posIndices])
  negEnds <- start(tx[!posIndices])
  posDown <- downstreamOf[posIndices]
  negDown <- downstreamOf[!posIndices]
  pos <- posEnds > posDown
  neg <- negEnds < negDown
  posTx <- tx[posIndices][pos]
  negTx <- tx[!posIndices][neg]
  downTx <- tx
  downTx[posIndices] <- posTx
  downTx[!posIndices] <- negTx
  #check if anyone hits boundary, set those to boundary
  if (anyNA(strandPerGroup(downTx, FALSE))) {
    boundaryHits <- which(is.na(strandPerGroup(downTx, FALSE)))
    downTx[boundaryHits] <- firstExonPerGroup(tx[boundaryHits])
    ir <- IRanges(start = downstreamOf[boundaryHits],
                  end = downstreamOf[boundaryHits])
    irl <- split(ir, seq_along(ir))
    names(irl) <- names(tx[boundaryHits])
    ranges(downTx[boundaryHits]) <- irl
  }
  # check boundaries within group exons
  startSites <- startSites(downTx, FALSE, FALSE, TRUE)
  posChecks <- startSites[posIndices] > downstreamOf[posIndices] & any(!pos)
  negChecks <- startSites[!posIndices] < downstreamOf[!posIndices] & any(!neg)
  if (any(posChecks)) {
    downstreamOf[posIndices][posChecks] <- startSites[posIndices][posChecks]
  }
  if (any(negChecks)) {
    downstreamOf[!posIndices][negChecks] <- startSites[!posIndices][negChecks]
  }

  return(assignFirstExonsStartSite(downTx, downstreamOf))
}

#' Get rest of objects downstream (inclusive)
#'
#' Per group get the part downstream of position.
#' downstreamFromPerGroup(tx, startSites(threeUTRs, asGR = TRUE))
#' will return the  3' utrs per transcript as GRangesList,
#' usually used for interesting
#' parts of the transcripts.
#'
#' If you don't want to include the points given in the region,
#' use \code{\link{downstreamOfPerGroup}}
#' @param tx a \code{\link{GRangesList}},
#'  usually of Transcripts to be changed
#' @param downstreamFrom a vector of integers, for each group in tx, where
#' is the new start point of first valid exon.
#' @return a GRangesList of downstream part
#'
downstreamFromPerGroup <- function(tx, downstreamFrom) {
  # Needs speed update!
  posIndices <- strandBool(tx)
  posEnds <- end(tx[posIndices])
  negEnds <- start(tx[!posIndices])
  posDown <- downstreamFrom[posIndices]
  negDown <- downstreamFrom[!posIndices]
  pos <- posEnds >= posDown
  neg <- negEnds <= negDown
  posTx <- tx[posIndices][pos]
  negTx <- tx[!posIndices][neg]
  downTx <- tx
  downTx[posIndices] <- posTx
  downTx[!posIndices] <- negTx
  #check if anyone hits boundary, set those to boundary
  if (anyNA(strandPerGroup(downTx, FALSE))) {
    boundaryHits <- which(is.na(strandPerGroup(downTx, FALSE)))
    downTx[boundaryHits] <- firstExonPerGroup(tx[boundaryHits])
    ir <- IRanges(start = downstreamFrom[boundaryHits],
                  end = downstreamFrom[boundaryHits])
    irl <- split(ir, seq_along(ir))
    names(irl) <- names(tx[boundaryHits])
    ranges(downTx[boundaryHits]) <- irl
  }

  return(assignFirstExonsStartSite(downTx, downstreamFrom))
}


#' Get rest of objects upstream (exclusive)
#'
#' Per group get the part upstream of position
#' upstreamOfPerGroup(tx, startSites(cds, asGR = TRUE))
#' will return the 5' utrs per transcript, usually used for interesting
#' parts of the transcripts.
#'
#' @param tx a \code{\link{GRangesList}},
#'  usually of Transcripts to be changed
#' @param upstreamOf a vector of integers, for each group in tx, where
#'  is the the base after the new stop point of last valid exon.
#' @param allowOutside a logical (T), can upstreamOf extend outside
#'  range of tx, can set boundary as a false hit, so beware.
#' @return a GRangesList of upstream part
#'
upstreamOfPerGroup <- function(tx, upstreamOf, allowOutside = TRUE) {
  posIndices <- strandBool(tx)
  posStarts <- start(tx[posIndices])
  negStarts <- end(tx[!posIndices])
  posGrlStarts <- upstreamOf[posIndices]
  negGrlStarts <- upstreamOf[!posIndices]
  pos <- posStarts < posGrlStarts
  neg <- negStarts > negGrlStarts
  posTx <- tx[posIndices]
  negTx <- tx[!posIndices]

  # need to fix pos/neg with possible cage extensions
  if (allowOutside) {
    outside <- which(sum(pos) == 0)
    pos[outside] = TRUE
    posTx[outside] <- firstExonPerGroup(posTx[outside])
    outside <- which(sum(neg) == 0)
    neg[outside] = TRUE
    negTx[outside] <- firstExonPerGroup(negTx[outside])
  }

  posTx <- posTx[pos]
  negTx <- negTx[neg]
  tx[posIndices] <- posTx
  tx[!posIndices] <- negTx
  nonZero <- widthPerGroup(tx) > 0
  if (all(!nonZero)) { # if no ranges exists
    return(tx)
  }
  upstreamOf <- upstreamOf[nonZero]
  posIndices <- posIndices[nonZero]

  stopSites <- stopSites(tx[nonZero], FALSE, FALSE, TRUE)
  if (any(posIndices)){
    posChecks <- stopSites[posIndices] < upstreamOf[posIndices] &
      any(!pos[nonZero[posIndices]])
  } else {
    posChecks <- FALSE
  }
  if(any(!posIndices)){
    negChecks <- stopSites[!posIndices] > upstreamOf[!posIndices] &
      any(!neg[nonZero[!posIndices]])
  } else {
    negChecks <- FALSE
    }

  if (any(posChecks)) {
    upstreamOf[posIndices][posChecks] <- stopSites[posIndices][posChecks]
  }
  if (any(negChecks)) {
    upstreamOf[!posIndices][negChecks] <- stopSites[!posIndices][negChecks]
  }

  tx[nonZero] <- assignLastExonsStopSite(tx[nonZero], upstreamOf)
  return(tx)
}


#' Extend the leaders transcription start sites.
#'
#' Will extend the leaders or transcripts upstream by extension.
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
#' Extends function [GenomicRanges::reduce()]
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

#' Faster more secure version of mapFromTranscript
#'
#' Fixes a bug in function, and should have 10x speedup
#' Also removes hit column for you
#'
#' @param ranges IRanges of ranges within grl
#' @param grl the "transcripts" that contain ranges, GRangesList
#' @param indices integer vector of which index of grl ranges are from:
#' (c(1,1,2)) means first two ranges are from grl[1], third from grl[2])
#' @return A GrangesList of ranges mapped from transcripts
pmapFromTranscriptF <- function(ranges, grl, indices) {
  names <- names(grl)
  names(grl) <- NULL
  genomicCoordinates <- pmapFromTranscripts(x = ranges,
                                            transcripts = grl[indices])
  names(genomicCoordinates) <- names[indices]

  genomicCoordinates <- genomicCoordinates[width(genomicCoordinates) > 0];
  a <- unlistGrl(genomicCoordinates);
  a$hit <- NULL;
  return(relist(a, genomicCoordinates))
}
