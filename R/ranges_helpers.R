#' Make grouping by exons ranks
#'
#' There are two ways to make vector of exon ranking:
#' 1. Iterate per exon in ORF, byTranscript = FALSE
#' 2. Iterate per ORF in transcript, byTranscript = TRUE.
#'
#' Either by transcript or by original groupings.
#' Must be ordered, so that same transcripts are ordered together.
#' @param grl a \code{\link{GRangesList}}
#' @param byTranscript logical (default: FALSE), groups orfs by transcript
#' name or ORF name, if ORfs are by transcript, check duplicates.
#' @return an integer vector of indices for exon ranks
#' @importFrom S4Vectors nrun Rle
#' @keywords internal
makeExonRanks <- function(grl, byTranscript = FALSE) {
  validGRL(class(grl))
  g <- groupings(grl)

  if (byTranscript & length(g) > 1) {
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
#'
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
#' This is not supported originally by GenomicRanges for GRangesList.
#' @param grl a \code{\link{GRangesList}} object with names.
#' @param sort.on.return logical (TRUE), should the groups be
#'  sorted before return (Negative ranges should be in decreasing order).
#'  Makes it a bit slower, but much safer for downstream analysis.
#' @param matchNaming logical (TRUE), should groups keep unlisted names
#'  and meta data.(This make the list very big, for > 100K groups)
#' @param is.sorted logical (TRUE), grl is presorted
#' (negative coordinates are decreasing). Set to FALSE if they are not,
#' else output will most likely be wrong!
#' @return a GRangesList grouped by original group, tiled to 1. Groups with identical names
#' will be merged.
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
tile1 <- function(grl, sort.on.return = TRUE, matchNaming = TRUE,
                  is.sorted = TRUE) {
  if (!is.grl(grl)) stop("grl must be a GRangesList")
  if (!is.sorted) grl <- sortPerGroup(grl)

  # optimize by removing all unecessary data
  ORFs <- unlist(grl, use.names = FALSE)
  mcols(ORFs) <- NULL
  names(ORFs) <- NULL
  tilex <- tile(ORFs, width =  1L)
  grouping <- groupings(grl)
  names(tilex) <- grouping

  # only negative with > 1 exons must be sorted
  if (sort.on.return) {
    negIndices <- !strandBool(tilex) & numExonsPerGroup(tilex) > 1
    tilex[negIndices] <- sortPerGroup(tilex[negIndices], quick.rev = TRUE)
  }

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
  return(tilex)
}

#' Map genomic to transcript coordinates by reference
#'
#' Map range coordinates between features in the genome and
#' transcriptome (reference) space.
#'
#' Similar to GenomicFeatures' pmapToTranscripts, but in this version the
#' grl ranges are compared to reference ranges with same name, not
#' by index. This gives a large speedup, but also requires
#' all objects must be named.
#' @param grl a \code{\link{GRangesList}} of ranges within
#'  the reference, grl must either have names matching, or a meta column called
#' 'names' that gives grouping names. i.e. grl named uORF_1_ENST00001, must then
#'  have a names meta column with ENST00001.
#' @param reference a GRangesList of ranges that include grl as a subset
#' of ranges. Example: cds is grl and mrna can be reference
#' @inheritParams pmapToTranscriptF
#' @return a GRangesList in transcript coordinates
#' @family ExtendGenomicRanges
#' @export
#' @examples
#' seqname <- c("tx1", "tx2", "tx3")
#' seqs <- c("ATGGGTATTTATA", "AAAAA", "ATGGGTAATA")
#' grIn1 <- GRanges(seqnames = "1",
#'                  ranges = IRanges(start = c(21, 10), end = c(23, 19)),
#'                  strand = "-")
#' grIn2 <- GRanges(seqnames = "1",
#'                  ranges = IRanges(start = c(1), end = c(5)),
#'                  strand = "-")
#' grIn3 <- GRanges(seqnames = "1",
#'                  ranges = IRanges(start = c(1010), end = c(1019)),
#'                  strand = "-")
#' grl <- GRangesList(grIn1, grIn2, grIn3)
#' names(grl) <- seqname
#' # Find ORFs
#' test_ranges <- findMapORFs(grl, seqs,
#'                  "ATG|TGG|GGG",
#'                  "TAA|AAT|ATA",
#'                  longestORF = FALSE,
#'                  minimumLength = 0)
#' # Genomic coordinates ORFs
#' test_ranges
#' # Transcript coordinate ORFs
#' asTX(test_ranges, reference = grl)
#' # seqnames will here be index of transcript it came from
#'
asTX <- function(grl, reference,
                 ignore.strand = FALSE,
                 x.is.sorted = TRUE,
                 tx.is.sorted = TRUE) {
  orfNames <- txNames(grl, reference) # Find tx they came from
  if (sum(orfNames %in% names(reference)) != length(orfNames)) {
    orfNames <- names(grl)
    if (is.null(orfNames) || sum(orfNames %in% names(reference)) != length(orfNames)) {
      stop("not all references are present, so can not map to transcripts.",
           " Does your annotation contains  '_'+number naming?",
           "ORFik uses this for a special purpose.")
    }
  }
  reference <- reference[orfNames] # Duplicate needed references
  names(reference) <- NULL
  return(pmapToTranscriptF(grl, reference,
                           ignore.strand = ignore.strand,
                           x.is.sorted = x.is.sorted,
                           tx.is.sorted = tx.is.sorted))
}

#' Faster pmapToTranscript
#'
#' Map range coordinates between features in the transcriptome and
#' genome (reference) space.
#' The length of x must be the same as length of transcripts. Only exception is
#' if x have integer names like (1, 3, 3, 5), so that x[1] maps to 1, x[2] maps
#' to transcript 3 etc.
#'
#' This version tries to fix the shortcommings of GenomicFeature's version.
#' Much faster and uses less memory.
#' Implemented as dynamic program optimized c++ code.
#' @param x GRangesList/GRanges/IRangesList/IRanges to map
#' to transcriptomic coordinates
#' @param transcripts a GRangesList/GRanges/IRangesList/IRanges to
#'  map against (the genomic coordinates). Must be of lower abstraction level
#'  than x. So if x is GRanges, transcripts can not be IRanges etc.
#' @param ignore.strand When ignore.strand is TRUE, strand is ignored in
#'  overlaps operations (i.e., all strands are considered "+") and the
#'  strand in the output is '*'. \cr
#'  When ignore.strand is FALSE (default) strand in the output is taken from the
#'  transcripts argument. When transcripts is a GRangesList, all inner list
#'  elements of a common list element must have the same strand
#'  or an error is thrown. \cr
#'  Mapped position is computed by counting from the transcription start site
#'   (TSS) and is not affected by the value of ignore.strand.
#' @param x.is.sorted if x is a GRangesList object, are "-" strand groups pre-sorted
#' in decreasing order within group, default: TRUE
#' @param tx.is.sorted if transcripts is a GRangesList object,
#' are "-" strand groups pre-sorted in decreasing order within group,
#' default: TRUE
#' @return object of same class as input x, names from ranges are kept.
#' @export
#' @examples
#' library(GenomicFeatures)
#' # Need 2 ranges object, the target region and whole transcript
#' # x is target region
#' x <- GRanges("chr1", IRanges(start = c(26, 29), end = c(27, 29)), "+")
#' names(x) <- rep("tx1_ORF1", length(x))
#' x <- groupGRangesBy(x)
#' # tx is the whole region
#' tx_gr <- GRanges("chr1", IRanges(c(5, 29), c(27, 30)), "+")
#' names(tx_gr) <- rep("tx1", length(tx_gr))
#' tx <- groupGRangesBy(tx_gr)
#' pmapToTranscriptF(x, tx)
#' pmapToTranscripts(x, tx)
#'
#' # Reuse names for matching
#' x <- GRanges("chr1", IRanges(start = c(26, 29, 5), end = c(27, 29, 18)), "+")
#' names(x) <- c(rep("tx1_1", 2), "tx1_2")
#' x <- groupGRangesBy(x)
#' tx1_2 <- GRanges("chr1", IRanges(c(4, 28), c(26, 31)), "+")
#' names(tx1_2) <- rep("tx1", 2)
#' tx <- c(tx, groupGRangesBy(tx1_2))
#'
#' a <- pmapToTranscriptF(x, tx[txNames(x)])
#' b <- pmapToTranscripts(x, tx[txNames(x)])
#' identical(a, b)
#' seqinfo(a)
#' # A note here, a & b only have 1 seqlength, even though the 2 "tx1"
#' # are different in size. This is an artifact of using duplicated names.
#'
#' ## Also look at the asTx for a similar useful function.
pmapToTranscriptF <- function(x, transcripts, ignore.strand = FALSE,
                              x.is.sorted = TRUE, tx.is.sorted = TRUE) {
  if ((length(x) == 0)) return(x)

  if (length(x) != length(transcripts)) { # Recycling
    if (length(x) == 1) {
      x <- rep(x, length(transcripts))
    } else if(length(transcripts) == 1) {
      transcripts <- rep(transcripts, length(x))
    } else stop("recycling is supported when length(x) == 1 or,",
                "length(transcripts) == 1; otherwise the lengths must match")
  }

  # Store original values we need
  xOriginal <- x; oldNames <- names(x); xClass <- class(x)
  xWidths <- width(xOriginal)

  oldTxNames <- names(transcripts)
  txWidths <- if (is.grl(transcripts)) {
    widthPerGroup(transcripts, FALSE)
  } else width(transcripts)

  xStrandOriginal <- if(is.grl(xOriginal)) {
    as.character(unlist(strand(xOriginal), use.names = FALSE))
  } else if (is(xOriginal, "GRanges")) {
    as.character(strand(xOriginal))
  } else NULL

  # subset to ranges and get indices for x
  if (is.grl(x) & !x.is.sorted) x <- sortPerGroup(x, ignore.strand)
  if (is.gr_or_grl(x)) x <- ranges(x)
  if (is(x, "IRangesList")) {
    indices <- groupings(x)
    xWidths <- as.integer(sum(xWidths))
    x <- unlist(x, use.names = FALSE)
    names(x) <- NULL
  } else if (is(x, "IRanges")) {
    indices <- seq.int(1, length(x))
    names(x) <- NULL
  } else stop("x must either be IRanges, IRangesList, GRanges or GRangesList")

  # Sanity tests
  if (!is.logical(ignore.strand)) stop("ignore.strand must be logical")
  ignore.strand <- as.logical(ignore.strand)
  if (max(indices) > length(transcripts)) stop("invalid names of IRanges")
  if (length(x) != length(indices)) stop("length of ranges != indices")
  notEqualSeqnames <- is.gr_or_grl(xOriginal) & is.gr_or_grl(transcripts) &
                         !all(seqlevels(xOriginal) %in% seqlevels(transcripts))
  if (notEqualSeqnames) stop("subscript contains out-of-bounds indices")
  # TODO: add propper test for per row seqnames, not just seqlevels
  if (is.grl(transcripts) & !tx.is.sorted)
    transcripts <- sortPerGroup(transcripts, ignore.strand)
  if (!all(xWidths <= txWidths)) {
    stop("Invalid ranges to map, check them.",
         " One has width bigger than its reference")
  }

  # Unlist tx, if list structure
  tx <- ranges(transcripts)
  names <- names(tx)
  names(tx) <- NULL
  if (is.grl(transcripts) | is(transcripts, "IRangesList")) {
    tx <- unlist(tx, use.names = FALSE)
    groupings <- groupings(transcripts)
    exonN <- lengths(transcripts)
  } else { # not list
    groupings <- seq.int(1, length(transcripts))
    exonN <- seq.int(1, length(transcripts))
  }

  txStrand <- if (ignore.strand) { # If ignore strand, set to '+'
    rep(TRUE, length(transcripts))
  } else strandBool(transcripts)
  xStrand <- if (ignore.strand) { # If ignore strand, set all to '+'
    rep(TRUE, length(indices))
  } else txStrand[indices]

  # Split indices for x into pos / neg-strand
  if (is.grl(xOriginal) | is(xOriginal, "IRangesList")) {
    indicesPos <- groupings(xOriginal[txStrand])
    indicesNeg <- groupings(xOriginal[!txStrand])
  } else {
    indicesPos <- seq_along(xOriginal[txStrand])
    indicesNeg <- seq_along(xOriginal[!txStrand])
  }
  # Make algorithm dynamic, by skipping if you know you can go to next transcript
  exonCumSumPos <- c(0, cumsum(exonN[txStrand]))[indicesPos]
  exonCumSumNeg <- c(0, cumsum(exonN[!txStrand]))[indicesNeg]
  txStrand <- txStrand[groupings]

  # Here is pos and neg direction of the algorithm ->
  # forward strand (c++ code)
  if (any(txStrand)) {
    pos <- pmapToTranscriptsCPP(start(x)[xStrand], end(x)[xStrand],
                                start(tx)[txStrand], end(tx)[txStrand],
                                groupings[txStrand], '+', exonCumSumPos)
  } else {
    pos <- list(ranges = list(vector("integer"), vector("integer")),
                index = vector("integer"))
  }

  # reverse strand (c++ code)
  if (any(!txStrand)) {
    neg <- pmapToTranscriptsCPP(start(x)[!xStrand], end(x)[!xStrand],
                                start(tx)[!txStrand], end(tx)[!txStrand],
                                groupings[!txStrand], '-',
                                exonCumSumNeg)
  } else {
    neg <- list(ranges = list(vector("integer"), vector("integer")),
                index = vector("integer"))
  }
  xStart <- vector("integer", length(indices))
  xStart[xStrand] <- unlist(pos$ranges[1], use.names = FALSE)
  xStart[!xStrand] <- unlist(neg$ranges[1], use.names = FALSE)

  xEnd <- vector("integer", length(indices))
  xEnd[xStrand] <- unlist(pos$ranges[2], use.names = FALSE)
  xEnd[!xStrand] <- unlist(neg$ranges[2], use.names = FALSE)

  result <- IRanges(xStart, xEnd)


  if (is.gr_or_grl(xClass)) { # If it was GR object, recreate
    newSeqnames <- if (!is.null(oldTxNames)) {
      oldTxNames
      # TODO: check that bellow seq.int is valid
    } else seq.int(length(result))
    newStrand <- if (ignore.strand) {
      "*"
    } else if (is.grl(transcripts)) {
      strandPerGroup(transcripts, FALSE)[indices]
    } else as.character(strand(transcripts))[indices]
    result <- GRanges(seqnames = newSeqnames[indices],
                      ranges = result, strand = newStrand)
    # TODO: remove this when not needed
    unmapped <- start(result) == 0
    if (!ignore.strand) # Check strand correction only if not ignore
      unmapped <- unmapped | (strand(result) != xStrandOriginal)
    if (any(unmapped)) {
      strand(result)[unmapped] <- "*"
    }

    seqlevels.used <- seqlevels(result)
    if(is.null(names(transcripts))) {
      seqlevels.used <- as.integer(seqlevels.used)
    }
    seqlengths(result) <- if (is.grl(transcripts)) {
      widthPerGroup(transcripts[seqlevels.used])
    } else as.integer(width(transcripts[seqlevels.used]))
  }

  if (is.grl(xClass) | is(xOriginal, "IRangesList")) {
    result <- split(result, indices)
    names(result) <- oldNames
    result <- reduce(result, drop.empty.ranges = FALSE)
  } else names(result) <- oldNames
  return(result)
}

#' Faster pmapFromTranscript
#'
#' Map range coordinates between features in the transcriptome and
#' genome (reference) space.
#' The length of x must be the same as length of transcripts. Only exception is
#' if x have integer names like (1, 3, 3, 5), so that x[1] maps to 1, x[2] maps
#' to transcript 3 etc.
#'
#' This version tries to fix the short commings of GenomicFeature's version.
#' Much faster and uses less memory.
#' Implemented as dynamic program optimized c++ code.
#' @param x IRangesList/IRanges/GRanges to map to genomic coordinates
#' @param transcripts a GRangesList to map against (the genomic coordinates)
#' @param removeEmpty a logical, remove non hit exons, else they are set
#'  to 0. That is all exons in the reference that the transcript coordinates
#'  do not span.
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
      indices <- seq.int(1, length(x))
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
  seqinfo(result) <- seqinfo(transcripts)
  return(result)
}

#' Get transcript sequence from a GRangesList and a faFile or BSgenome
#'
#' For each GRanges object, find the sequence of it from faFile or BSgenome.
#'
#' A wrapper around \code{\link{extractTranscriptSeqs}} that works for
#' DNAStringSet and ORFik \code{\link{experiment}} input.
#' For debug of errors do:
#' which(!(unique(seqnamesPerGroup(grl, FALSE)) %in% seqlevels(faFile)))
#' This happens usually when the grl contains chromsomes that the fasta
#' file does not have. A normal error is that mitocondrial chromosome is
#' called MT vs chrM even though they have same seqlevelsStyle. The
#' above line will give you which chromosome it is missing.
#' @param grl a \code{\link{GRangesList}} object
#' @inheritParams findFa
#' @param is.sorted a speedup, if you know the grl ranges are sorted
#' @param keep.names a logical, default (TRUE), if FALSE: return as
#' character vector without names.
#' @export
#' @return a \code{\link{DNAStringSet}} of the transcript sequences
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @family ExtendGenomicRanges
#'
txSeqsFromFa <- function(grl, faFile, is.sorted = FALSE,
                         keep.names = TRUE) {
  faFile <- findFa(faFile)
  if (!any(seqlevels(grl) %in% seqlevels(faFile)))
    stop("FaFile had no matching seqlevels to ranges object")
  if (!is.sorted) grl <- sortPerGroup(grl)

  # Check for optimization, if conditions are met
  path <- ifelse(is.character(faFile), faFile, ifelse(is(faFile, "FaFile"), faFile$path, ""))
  if (path != "") {
    if (file.size(path) < 1e7 & length(grl) > 50) {
      # TODO: make failsafe for more special cases
      faFile <- readDNAStringSet(path)
      seqnames_grl <- unique(seqnamesPerGroup(grl, FALSE))
      bad_name_format <- !all(seqnames_grl %in% names(faFile))
      if (bad_name_format) {
        new_names <- try(names(seqinfo(FaFile(path))))
        if (all(seqnames_grl %in% new_names)) names(faFile) <- new_names
      }
    }
  }
  if (is(faFile, "DNAStringSet")) {
    seqs <- faFile[grl]
  } else {
    seqs <- extractTranscriptSeqs(faFile, transcripts = grl)
  }
  if (!keep.names) return(as.character(seqs, use.names = FALSE))
  return(seqs)
}



#' Subset GRanges to get desired frame.
#'
#' Usually used for ORFs to get specific frame (0-2):
#' frame 0, frame 1, frame 2
#'
#' GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#' @param x A tiled to size of 1 GRanges object
#' @param frame A numeric indicating which frame to extract
#' @return GRanges object reduced to only first frame
#' @export
#' @examples
#' subsetToFrame(GRanges("1", IRanges(1:10, width = 1), "+"), 2)
subsetToFrame <- function(x, frame) {
  if (as.vector(strand(x) == "+")[1]) {
    x[seq.int(frame, length(x), 3)]
  } else {
    x[seq.int(length(x) + 1 - frame, 1, -3)]
  }
}


#' Reduce GRanges / GRangesList
#'
#' Reduce away all GRanges elements with 0-width.
#'
#' Extends function \code{\link{reduce}}
#' by trying to keep names and meta columns, if it is a
#' GRangesList. It also does not lose sorting for GRangesList,
#' since original reduce sorts all by ascending position.
#' If keep.names == FALSE, it's just the normal GenomicRanges::reduce
#' with sorting negative strands descending for GRangesList.
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @param drop.empty.ranges (FALSE) if a group is empty (width 0), delete it.
#' @param min.gapwidth (1L) how long gap can it be between two ranges,
#' to merge them.
#' @param with.revmap (FALSE) return info on which mapped to which
#' @param with.inframe.attrib (FALSE) For internal use.
#' @param ignore.strand (FALSE), can different strands be reduced together.
#' @param keep.names (FALSE) keep the names and meta columns of the GRangesList
#' @param min.strand.decreasing (TRUE), if GRangesList, return minus strand group
#' ranges in decreasing order (1-5, 30-50) -> (30-50, 1-5)
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
                           ignore.strand = FALSE, min.strand.decreasing = TRUE) {
  if (!is.gr_or_grl(class(grl))) stop("grl must be GRanges or GRangesList")

  reduced <- GenomicRanges::reduce(grl, drop.empty.ranges, min.gapwidth,
                                   with.revmap, with.inframe.attrib,
                                   ignore.strand)
  if (is.grl(class(grl)) & min.strand.decreasing) {
    reduced <- reverseMinusStrandPerGroup(reduced)
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
