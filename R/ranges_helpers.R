
#' Seqnames cleanup
#'
#' For many datasets, the fa file and the gtf file have different naming
#' This functions tries to fix the naming to the GRanges standard
#' chrX instead of X chr1 instead of 1 etc..
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @return a GRangesList with fixed seqnames
fixSeqnames <- function(grl){
  validGRL(class(grl))

  temp <- unlist(grl)
  seqnamesTransformed <- as.character(seqnames(temp))
  indexes <- which(nchar(seqnamesTransformed) < 6)
  temp <- temp[indexes]
  seqlevels(temp) <- sub(replacement = "chrY", pattern = "Y", seqlevels(temp))
  seqlevels(temp) <- sub(replacement = "chrX", pattern = "X", seqlevels(temp))
  seqlevels(temp) <- as.character(unique(seqnames(unlist(temp))))
  return(groupGRangesBy(temp))
}

#' Make a meta column with exon ranks
#'
#' Must be ordered, so that same transcripts are ordered together.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param byTranscript if ORfs are by transcript, check duplicates
#' @importFrom S4Vectors nrun Rle
#' @return an integer vector of indices for exon ranks
makeExonRanks <- function(grl, byTranscript = F){
  validGRL(class(grl))

  if (byTranscript) {
    oldNames <- names(grl)
    names(grl) <- 1:length(grl)
  }
  l <- Rle(names(unlist(grl, use.names = TRUE)))
  t <- unlist(lapply(1:nrun(l), function(x) {
    rep(x, runLength(l)[x])
  }))

  Inds <- rep(1, length(t))
  if (!byTranscript) {
    for (x in 2:length(t)) {
      if (t[x] == t[x - 1]) {
        Inds[x] <- Inds[x - 1] + 1
      }
    }
  } else {
    for (x in 2:length(t)) {
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @return (GRangesList) with ORF names, grouped by transcripts, sorted.
#'
makeORFNames <- function(grl){
  ranks <- makeExonRanks(grl, byTranscript = TRUE)
  asGR <- unlist(grl, use.names = FALSE)
  if (is.null(names(asGR))) asGR <- unlist(grl, use.names = TRUE)
    asGR$names <- paste0(names(asGR), "_", ranks)
  return(groupGRangesBy(asGR))
}

#' Tile a GRangesList by 1
#'
#' Per group, sepereate the groups and split them on each position.
#' Returned sorted
#' This is not supported originally by GenomicRanges
#' As a precaution, this function requires the unlisted objects to
#' have names.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object with names
#' @examples
#' gr1 <- GRanges("1",
#'            ranges = IRanges(start = c(1, 10, 20),
#'            end = c(5, 15, 25)),
#'            strand = "+")
#'
#' gr2 <- GRanges("1", ranges = IRanges(start = c(20, 30, 40),
#'                                       end = c(25, 35, 45)),
#'                      strand = "+")
#' names(gr1) = rep("tx1_1", 3)
#' names(gr2) = rep("tx1_2", 3)
#' grl <- GRangesList(tx1_1 = gr1, tx1_2 = gr2)
#' tile1(grl)
#' @return a GRangesList grouped by original group, tiled to 1
#' @export
#'
tile1 <- function(grl){
  ORFs <- unlistGrl(grl)
  if (is.null(names(ORFs)) ||
      (length(unique(names(ORFs))) != length(names(grl)))) {
    stop("Unlisted grl has no names.")
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
      ORFs$names <- names(ORFs)
    }
  }
  # special case for only single grouped GRangesList
  if (is.null(ORFs$names)) {
    if(length(ORFs) != length(grl)) {
      stop("wrong naming, could not find unique names")
    }
    ORFs$names <- names(grl)
  }

  tilex <- tile(ORFs, width =  1L)
  names(tilex) <- ORFs$names
  tilex <- matchNaming(tilex, grl[1])
  # only negative must be sorted
  posIndices <- strandBool(tilex)
  tilex[!posIndices] <- sortPerGroup(tilex[!posIndices])
  return(tilex)
}

#' Map genomic to transcript coordinates by reference
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} of ranges within
#'  the reference, grl must have column called names that gives
#'  grouping for result
#' @param reference a GrangesList of ranges
#'  that include and are bigger or equal to grl
#'  ig. cds is grl and gene can be reference
#'  @export
asTX <- function(grl, reference){
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
#' @return a DNAStringSet of the transcript sequences
txSeqsFromFa <- function(grl, faFile, is.sorted = F){
  if(!(class(faFile) == "FaFile" || class(faFile) == "BSgenome"))
    stop("only FaFile and BSgenome is valid input class")
  if(!is.sorted) grl <- sortPerGroup(grl)
  return(extractTranscriptSeqs(faFile, transcripts = grl))
}



#' Reassign the start positions of the first exons per group in grl
#' @description make sure your grl is sorted, since start of "-" strand
#'  objects should be the
#'  max end in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param newStarts an integer vector of same length as grl, with new start values
#' @return the same GRangesList with new start sites
assignFirstExonsStartSite <- function(grl, newStarts){
  if (length(grl) != length(newStarts)) stop("length of grl and newStarts \n
                                            are not equal!")
  posIndices <- strandBool(grl)

  dt <- as.data.table(grl)
  dt[!duplicated(dt$group),]$start[posIndices] <- newStarts[posIndices]
  dt[!duplicated(dt$group),]$end[!posIndices] <- newStarts[!posIndices]

  ngrl <- GenomicRanges::makeGRangesListFromDataFrame(dt,
    split.field = "group", names.field = "group_name", keep.extra.columns = T)
  names(ngrl) <- names(grl)

  return(ngrl)
}
#' Reassign the stop positions of the last exons per group
#' @description make sure your grl is sorted,
#'  since stop of "-" strand objects should be the
#'  min start in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param newStops an integer vector of same length as grl,
#'  with new start values
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @return the same GRangesList with new stop sites
assignLastExonsStopSite <- function(grl, newStops){
  if (length(grl) != length(newStops)) stop("length of grl and newStops \n
                                           are not equal!")
  posIndices <- strandBool(grl)

  dt <- as.data.table(grl)
  group <- NULL # avoid check warning
  idx = dt[, .I[.N], by = group]
  dt[idx$V1]$end[posIndices] <- newStops[posIndices]
  dt[idx$V1]$start[!posIndices] <- newStops[!posIndices]
  ngrl <- GenomicRanges::makeGRangesListFromDataFrame(dt,
    split.field = "group", names.field = "group_name",
      keep.extra.columns = TRUE)
  names(ngrl) <- names(grl)

  return(ngrl)
}
#' Get rest of objects downstream
#'
#' Per group get the part downstream of position
#'  defined in downstreamOf
#'  downstreamOf(tx, ORFik:::stopSites(cds, asGR = F))
#'  will return the 3' utrs per transcript as GRangesList,
#'  usually used for interesting
#'  parts of the transcripts, like upstream open reading frames(uorf).
#'  downstreamOf +/- 1 is start/end site
#'  of transformed tx's, depending on strand
#' @param tx a \code{\link[GenomicRanges]{GRangesList}},
#'  usually of Transcripts to be changed
#' @param downstreamOf a vector of integers, for each group in tx, where
#' is the new start point of first valid exon.
#' @return a GRangesList of downstream part
downstreamOfPerGroup <- function(tx, downstreamOf){
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
    irl <- split(ir, 1:length(ir))
    names(irl) <- names(tx[boundaryHits])
    ranges(downTx[boundaryHits]) <- irl
  }
  return(assignFirstExonsStartSite(downTx, downstreamOf))
}
#' Get rest of objects upstream
#'
#' Per group get the part upstream of position
#'  defined in upstreamOf
#'  upstream(tx, ORFik:::stopSites(cds, asGR = F))
#'  will return the 5' utrs per transcript, usually used for interesting
#'  parts of the transcripts, like upstream open reading frames(uorf).
#'  downstreamOf +/- 1 is start/end site
#'  of transformed tx's, depending on strand
#' @param tx a \code{\link[GenomicRanges]{GRangesList}},
#'  usually of Transcripts to be changed
#' @param upstreamOf a vector of integers, for each group in tx, where
#'  is the new stop point of last valid exon.
#' @return a GRangesList of upstream part
upstreamOfPerGroup <- function(tx, upstreamOf){
  posIndices <- strandBool(tx)
  posStarts <- start(tx[posIndices])
  negStarts <- end(tx[!posIndices])
  posGrlStarts <- upstreamOf[posIndices]
  negGrlStarts <- upstreamOf[!posIndices]
  pos <- posStarts < posGrlStarts
  neg <- negStarts > negGrlStarts
  # need to fix pos/neg with possible cage extensions
  outside <- which(sum(pos) == 0)
  pos[outside] = T
  posTx <- tx[posIndices]
  posTx[outside] <- firstExonPerGroup(posTx[outside])
  outside <- which(sum(neg) == 0)
  neg[outside] = T
  negTx <- tx[!posIndices]
  negTx[outside] <- firstExonPerGroup(negTx[outside])

  posTx <- posTx[pos]
  negTx <- negTx[neg]
  tx[posIndices] <- posTx
  tx[!posIndices] <- negTx

  return(assignLastExonsStopSite(tx, upstreamOf))
}

#' Extend the leaders tss.
#'
#' Either if you have the 5' utr sequences or the
#'  transcript sequences to extend,
#'  any way the 5' part will be extended.
#'  Remember to sort it, if not sorted beforehand.
#'  use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#'  of 5' utrs or transcripts..
#' @param extension an integer, how much to extend the leaders.
#'  Or a GRangesList where start / stops by strand are the positions
#'  to use as new starts.
#' @param cds If you want to extend 5' leaders downstream, to catch
#'  upstream ORFs going into cds, include it. It will add first
#'  cds exon to grl matched by names.
#'  Do not add for transcripts, as they are already included.
#' @examples
#' library(GenomicFeatures)
#' samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#' package = "GenomicFeatures")
#' txdb <- loadDb(samplefile)
#' fiveUTRs <- fiveUTRsByTranscript(txdb) # <- extract only 5' leaders
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' cds <- cdsBy(txdb,"tx",use.names = TRUE)
#' ## now try(extend upstream 1000, downstream 1st cds exons):
#' extendLeaders(fiveUTRs, extension = 1000, cds)
#'
#' ## when extending transcripts, don't include cds' ofcourse,
#' ## since they are already there
#' extendLeaders(tx, extension = 1000)
#'@export
#'@return an extended GRangeslist
extendLeaders <- function(grl, extension = 1000, cds = NULL){
  if (class(extension) == "numeric" && length(extension) == 1) {
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
  return(addFirstCdsOnLeaderEnds(
    extendedLeaders, cds))
}

#' Get coverage per group
#'
#' It tiles each GRangesList group, and finds hits per position
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#'  of 5' utrs or transcripts..
#' @param reads a GAlignment or GRanges object of
#'  ribo-seq, RNA-seq etc
#' @examples
#' ORF <- GRanges(seqnames = "1",
#' ranges = IRanges(start = c(1, 10, 20),
#'                 end = c(5, 15, 25)),
#' strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25),"+")
#' ribosomeStallingScore(grl, RFP)
#' @export
#' @return a Rle, one list per group with # of hits per position.
coveragePerTiling <- function(grl, reads){

  unlTile <- unlistGrl(tile1(grl))

  # could make this more efficient by counting overlaps
  # only on untiled, then tile the ones that hit and count again
  counts <- countOverlaps(unlTile, reads)
  names <- names(counts)
  names(counts) <- NULL
  countList <- split(counts, names)

  return(IRanges::RleList(countList)[names(grl)])
}


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
#' windowResize(GRanges(Rle(c('1'), c(4)),
#'                       IRanges(c(100, 200, 200, 100), width=c(1, 1, 1, 1)),
#'                       Rle(strand(c('+', '+', '-', '-')))),
#'              window_size = 50)
#'
windowResize <- function(GRanges_obj, window_size = 30) {
  GRanges_obj <- promoters(GRanges_obj, upstream = window_size,
                           downstream = window_size + 1)
  return(GRanges_obj)
}


#' Subset GRanges to get desired frame. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’ @export
#' @param x A tiled to size of 1 GRanges object
#' @param frame A numeric indicating which frame to extract
#' @return GRanges object reduced to only first frame
#' @examples
#' #subset_to_frame(x, 1)
#'
subset_to_frame <- function(x, frame){
  if (as.vector(strand(x) == "+")[1]) {
    x[seq(frame, length(x), 3)]
  } else {
    x[seq(length(x) + 1 - frame, 1, -3)]
  }
}

#' Reduce GRanges / GRangesList
#'
#' Extends function \code{\link[GenomicRanges]{reduce}}
#' by trying to keep names and meta columns if it is a
#' GRangesList. It also does not loose sorting for GRangesList,
#' since original reduce sorts all by ascending.
#' If keep.names == F, it's just the normal GenomicRanges::reduce
#' with sorting negative strands descending for GRangesList
#'
#' Only tested for orfik, might not work for other naming conventions.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} or GRanges object
#' @param drop.empty.ranges (FALSE) if a group is empty (width 0), delete it.
#' @param min.gapwidth (1L) how long gap can it be to say they belong together
#' @param with.revmap (FALSE) return info on which mapped to which
#' @param with.inframe.attrib (FALSE) For internal use.
#' @param ignore.strand (FALSE), can different strands be reduced together.
#' @param keep.names (FALSE) keep the names and meta columns of the GRangesList
#' @examples
#' ORF <- GRanges(seqnames = "1",
#' ranges = IRanges(start = c(1, 2, 3),
#'                 end = c(1, 2, 3)),
#' strand = "+")
#' ##For GRanges
#' reduceKeepAttr(ORF, keep.names = TRUE)
#' ## For GRangesList
#' grl <- GRangesList(tx1_1 = ORF)
#' reduceKeepAttr(grl, keep.names = TRUE)
#' ##Only 1 GRanges object in GRangesList returned
#' @export
#' @return A reduced GRangesList
reduceKeepAttr <- function(grl, keep.names = FALSE,
                             drop.empty.ranges = FALSE, min.gapwidth = 1L,
                             with.revmap = FALSE, with.inframe.attrib = FALSE,
                             ignore.strand = FALSE){
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
