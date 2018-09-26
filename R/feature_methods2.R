#' Get initiation score for a GRangesList of ORFs
#'
#' initiationScore tries to check how much each TIS region resembles, the
#' average of the CDS TIS regions.
#'
#' Since this features uses a distance matrix for scoring, values are
#' distributed like this:
#' As result there is one value per ORF:
#' 0.000: means that ORF had no reads
#' -1.000: means that ORF is identical to average of CDS
#' 1.000: means that orf is maximum different than average of CDS
#' @references doi: 10.1186/s12915-017-0416-0
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @param cds a \code{\link{GRangesList}} object with coding sequences
#' @param tx a GrangesList of transcripts covering grl.
#' @param footprints ribosomal footprints, given as Galignment object or
#'  Granges
#' @param pShifted a logical (TRUE), are riboseq reads p-shifted?
#' @family features
#' @return an integer vector, 1 score per ORF
#' @export
#' @importFrom BiocGenerics Reduce
#' @examples
#' # Good hiting ORF
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(21), end = c(40)),
#'                strand = "+")
#' names(ORF) <- c("tx1")
#' grl <- GRangesList(tx1 = ORF)
#' # 1 width position based
#' RFP <- GRanges("1", IRanges(c(21, 23, 50, 50, 50, 53, 53, 56, 59),
#'  c(21, 23, 50, 50, 50, 53, 53, 56, 59)), "+")
#' score(RFP) <- 28 # original width
#' cds <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(50), end = c(80)),
#'                strand = "+")
#' cds <- GRangesList(tx1 = cds)
#' tx <- GRanges(seqnames = "1",
#'                ranges = IRanges(1,85),
#'                strand = "+")
#' tx <- GRangesList(tx1 = tx)
#'
#' initiationScore(grl, cds, tx, RFP, pShifted = TRUE)
#'
initiationScore <- function(grl, cds, tx, footprints, pShifted = TRUE) {
  # train average cds model
  df <- riboTISCoverageProportion(cds, tx, footprints, average = TRUE,
                                  onlyProportion = FALSE, pShifted = pShifted)
  cdsProp <- split(Rle(df$prop), df$length)
  names(cdsProp) <- NULL

  # get ORF models
  prop <- riboTISCoverageProportion(grl, tx, footprints, average = FALSE,
                                    onlyProportion = TRUE, pShifted = pShifted,
                                    keep.names = TRUE)
  names <- names(prop[[1]])
  names(prop[[1]]) <- NULL
  dif <- lapply(seq.int(1:length(prop)), function(x)
    abs(prop[[x]] - cdsProp[x]))
  dif2 <- lapply(seq.int(1:length(prop)), function(x) sum(dif[[x]]))

  tempAns <- Reduce("+", dif2)/length(dif2) - 1

  ans <- rep.int(0, length(grl))
  ans[names(grl) %in% names] <- tempAns
  return(ans)
}


#' Get ORFscore for a GRangesList of ORFs
#'
#' ORFscore tries to check whether the first frame of the 3 possible frames in
#' an ORF has more reads than second and third frame.
#'
#' Pseudocode:
#' assume rff - is reads fraction in specific frame
#' \preformatted{ORFScore = log(rrf1 + rrf2 + rrf3)}
#' For all ORFs where rrf2 or rrf3 is bigger than rff1,
#' negate the resulting value.
#' \preformatted{ORFScore[rrf1Smaller] <- ORFScore[rrf1Smaller] * -1}
#'
#' As result there is one value per ORF:
#' Positive values say that the first frame have the most reads,
#' negative values say that the first frame does not have the most reads.
#' @references doi: 10.1002/embj.201488411
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @param RFP ribozomal footprints, given as Galignment object,
#'  Granges or GRangesList
#' @param is.sorted logical (F), is grl sorted.
#' @importFrom data.table .SD
#' @importFrom data.table .N
#' @family features
#' @export
#' @return a matrix with 4 columns, the orfscore and score of
#' each of the 3 tiles
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' names(ORF) <- c("tx1", "tx1", "tx1")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25), "+") # 1 width position based
#' score(RFP) <- 28 # original width
#' orfScore(grl, RFP) # negative because more hits on frames 1,2 than 0.
#'
#' # example with positive result, more hits on frame 0 (in frame of ORF)
#' RFP <- GRanges("1", IRanges(c(1, 1, 1, 25), width = 1), "+")
#' score(RFP) <- c(28, 29, 31, 28) # original width
#' orfScore(grl, RFP)
#'
orfScore <- function(grl, RFP, is.sorted = FALSE) {
  if (length(grl) > 50000) { # faster version for big grl
    # only do ORFs that have hits
    validIndices <- hasHits(grl, RFP)
    if (!any(validIndices)) { # no variance in countList, 0 entropy
      frame_zero_RP <- frame_one_RP <- frame_two_RP <-
        ORFScores <- rep(0, length(grl))
      return(data.table(frame_zero_RP, frame_one_RP, frame_two_RP, ORFScores))
    }
    # reduce to unique orfs
    grl <- grl[validIndices]
    reOrdering <- uniqueOrder(grl)
    # find coverage
    cov <- coverageByWindow(RFP, uniqueGroups(grl),
                            is.sorted = is.sorted, keep.names = FALSE)

    countsTile1 <- countsTile2 <- countsTile3 <- rep(0, length(validIndices))
    len <- BiocGenerics::lengths(cov)
    # make the 3 frames
    positionFrame <- lapply(len, function(x){seq.int(1, x, 3)})
    tempTile <- sum(cov[positionFrame])[reOrdering]
    countsTile1[validIndices] <- tempTile # correct order and size
    positionFrame <- lapply(len, function(x){seq.int(2, x, 3)})
    tempTile <- sum(cov[positionFrame])[reOrdering]
    countsTile2[validIndices] <- tempTile
    positionFrame <- lapply(len, function(x){seq.int(3, x, 3)})
    tempTile <- sum(cov[positionFrame])[reOrdering]
    countsTile3[validIndices] <- tempTile
  } else {
    # tile the orfs into a d.t for easy seperation
    dt <- as.data.table(tile1(grl, matchNaming = FALSE))

    group <- NULL
    # seperate the three tiles, by the 3 frames
    tilex1 <- dt[, .SD[seq.int(1, .N, 3)], by = group]
    grl1 <- makeGRangesListFromDataFrame(
      tilex1, split.field = "group")
    tilex2 <- dt[, .SD[seq.int(2, .N, 3)], by = group]
    grl2 <- makeGRangesListFromDataFrame(
      tilex2, split.field = "group")
    tilex3 <- dt[, .SD[seq.int(3, .N, 3)], by = group]
    grl3 <- makeGRangesListFromDataFrame(
      tilex3, split.field = "group")

    countsTile1 <- countOverlaps(grl1, RFP)
    countsTile2 <- countOverlaps(grl2, RFP)
    countsTile3 <- countOverlaps(grl3, RFP)
  }

  RP = countsTile1 + countsTile2 + countsTile3

  Ftotal <- RP/3

  frame1 <- (countsTile1 - Ftotal)^2 / Ftotal
  frame2 <- (countsTile2 - Ftotal)^2 / Ftotal
  frame3 <- (countsTile3 - Ftotal)^2 / Ftotal

  dfORFs <- NULL
  dfORFs$frame_zero_RP <- countsTile1
  dfORFs$frame_one_RP <- countsTile2
  dfORFs$frame_two_RP <- countsTile3

  ORFscore <- log2(frame1 + frame2 + frame3 + 1)
  revORFscore <-  which(frame1 < frame2 | frame1 < frame3)
  ORFscore[revORFscore] <- -1 * ORFscore[revORFscore]
  ORFscore[is.na(ORFscore)] <- 0
  dfORFs$ORFScores <- ORFscore

  return(dfORFs)
}

#' Get distances between ORF Start and TSS of its transcript
#'
#' Matching is done by transcript names.
#' This is applicable practically to the upstream (fiveUTRs) ORFs.
#' If ORF is not within specified search space in fiveUTRs, this function
#' will crash.
#' @references doi: 10.1074/jbc.R116.733899
#' @param ORFs orfs as \code{\link{GRangesList}},
#' names of orfs must be txname_[rank]
#' @param fiveUTRs 5' leaders as \code{\link{GRangesList}}.
#' @return an integer vector, 1 means on TSS, 2 means second base of Tx.
#' @family features
#' @export
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(5, 10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(2, 20), "+"))
#' distToTSS(grl, fiveUTRs)
#'
distToTSS <- function(ORFs, fiveUTRs){
  validGRL(class(ORFs), "ORFs")

  startSites <-  startSites(ORFs, asGR = TRUE, keep.names = TRUE,
                            is.sorted = TRUE)
  return(start(asTX(startSites, fiveUTRs)))
}

#' Get distances between ORF ends and starts of their transcripts cds'.
#'
#' Will calculate distance between each ORF end and begining of the
#' corresponding cds. Matching is done by transcript names.
#' This is applicable practically to the upstream (fiveUTRs) ORFs.
#' The cds start site, will be presumed to be on + 1 of end of fiveUTRs.
#' @references doi: 10.1074/jbc.R116.733899
#' @param ORFs orfs as \code{\link{GRangesList}},
#' names of orfs must be transcript names
#' @param fiveUTRs fiveUTRs as \code{\link{GRangesList}},
#' remember to use CAGE version of 5' if you did CAGE reassignment!
#' @param cds cds' as \code{\link{GRangesList}},
#' only add if you have ORFs going into CDS.
#' @return an integer vector, +1 means one base upstream of cds, -1 means
#' 2nd base in cds, 0 means orf stops at cds start.
#' @family features
#' @export
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1, 10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1, 20), "+"))
#' distToCds(grl, fiveUTRs)
#'
distToCds <- function(ORFs, fiveUTRs, cds = NULL){
  validGRL(class(ORFs), "ORFs")

  cdsStarts <- widthPerGroup(fiveUTRs[
    txNames(ORFs)], FALSE) + 1

  lastExons <- lastExonPerGroup(ORFs)
  if (is.grl(cds)) {
    fiveUTRs <- addCdsOnLeaderEnds(fiveUTRs, cds)
  }
  orfsTx <- asTX(lastExons, fiveUTRs)

  # this is ok, since it is tx not genomic ->
  orfEnds <- lastExonEndPerGroup(orfsTx, FALSE)

  return(cdsStarts - orfEnds)
}


#' Make a score for each ORFs start region by proximity to Kozak
#'
#' The closer the sequence is to the Kozak sequence
#' the higher the score, based on the experimental pwms
#' from article referenced.
#' Minimum score is 0 (worst correlation), max is 1 (the best
#' base per column was chosen).
#' @references doi: https://doi.org/10.1371/journal.pone.0108475
#' @param grl a \code{\link{GRangesList}} grouped by ORF
#' @param faFile a FaFile from the fasta file, see ?FaFile.
#'  Can also be path to fastaFile with fai file in same dir.
#' @param species ("human"), which species to use,
#' currently supports human, zebrafish and mouse (m. musculus).
#' You can also specify a pfm for your own species.
#' Syntax of pfm is an rectangular integer matrix,
#' where all columns must sum to the same value, normally 100.
#' See example for more information.
#' Rows are in order: c("A", "C", "G", "T")
#' @param include.N logical (F), if TRUE, allow N bases to be counted as hits,
#' score will be average of the other bases. If True, N bases will be
#' added to pfm, automaticly, so dont include them if you make your own pfm.
#' @return a numeric vector with values between 0 and 1
#' @return an integer vector, one score per orf
#' @family features
#' @importFrom Biostrings PWM
#' @export
#' @examples
#' # Usually the ORFs are found in orfik, which makes names for you etc.
#' # Here we make an example from scratch
#' seqName <- "Chromosome"
#' ORF1 <- GRanges(seqnames = seqName,
#'                    ranges = IRanges(c(1007, 1096), width = 60),
#'                    strand = c("+", "+"))
#' ORF2 <- GRanges(seqnames = seqName,
#'                     ranges = IRanges(c(400, 100), width = 30),
#'                     strand = c("-", "-"))
#' ORFs <- GRangesList(tx1 = ORF1, tx2 = ORF2)
#' ORFs <- makeORFNames(ORFs) # need ORF names
#' # get faFile for sequences
#' faFile <- FaFile(system.file("extdata", "genome.fasta",
#'   package = "ORFik"))
#' kozakSequenceScore(ORFs, faFile)
#' # For more details see vignettes.
kozakSequenceScore <- function(grl, faFile, species = "human",
                               include.N = FALSE) {
  faFile <- findFa(faFile)
  firstExons <- firstExonPerGroup(grl)
  kozakLocation <- promoters(firstExons, upstream = 9, downstream = 6)

  sequences <- as.character(txSeqsFromFa(kozakLocation,
                                         faFile, is.sorted = TRUE))
  if (!all(nchar(sequences) == 15)) {
    stop("not all ranges had valid kozak sequences length, not 15")
  }

  if(is(species, "matrix")){
    # self defined pfm
    pfm <- species
  } else if (species == "human") {
    # human pfm, see article reference
    pfm <- t(matrix(as.integer(c(20,20,21,21,19,24,46,29,19,22,28,16,
                                 27,33,32,23,32,38,10,38,45,15,39,26,
                                 35,29,28,39,30,26,37,20,28,49,18,37,
                                 18,18,19,17,19,12,7,13,8,14,15,21)),
                    ncol = 4))
  } else if (species == "mouse") {
    # zebrafish pfm, see article reference
    pfm <- t(matrix(as.integer(c(20,19,21,20,18,25,49,28,17,23,28,15,
                                 27,34,31,23,32,38,9,39,47,14,40,26,
                                 34,28,27,39,29,25,36,20,28,49,18,37,
                                 19,19,21,18,21,12,6,13,8,14,14,22)),
                    ncol = 4))
  } else if (species == "zebrafish") {
    # zebrafish pfm, see article reference
    pfm <- t(matrix(as.integer(c(29,26,28,26,22,35,62,39,28,24,27,17,
                                 21,26,24,16,28,32,5,23,35,12,42,21,
                                 25,24,22,33,22,19,28,17,27,47,16,34,
                                 25,24,26,25,28,14,5,21,10,17,15,28)),
                    ncol = 4))
  } else  {
    stop("Either input species as a matrix
         or name of presupported pfm organism")
  }

  bases <- c("A", "C", "G", "T")
  rownames(pfm) <- bases
  pwm <- PWM(pfm)

  if (include.N) {
    bases <- c(bases, "N")
    pwm <- rbind(pwm, colMeans(pwm))
    rownames(pwm) <- bases
  }

  # exclude start codon
  s <- paste0(substr(x = sequences, 1, 9), substr(x = sequences, 13, 15))
  # split strings and relist as letters of 9 rows
  subSplit <- strsplit(s, split = "")
  # this will not when ATG is on start of chr
  mat <- t(matrix(unlist(subSplit, use.names = FALSE), ncol = length(s)))

  scores <- rep(0., length(s))
  for (i in seq(ncol(mat))) {
    for (n in seq_along(bases)) {
      match <- mat[, i] == bases[n]
      scores[match] <- scores[match] + pwm[n, i]
    }
  }
  return(scores)
}


#' Inside/Outside score (IO)
#'
#' Inside/Outside score is defined as
#' \preformatted{(reads over ORF)/(reads outside ORF and within transcript)}
#' A pseudo-count of one was added to both the ORF and outside sums.
#' @references doi: 10.1242/dev.098345
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders, cds', 3' utrs or ORFs
#' @param RFP ribo seq reads as GAlignment, GRanges or GRangesList object
#' @param GtfOrTx if Gtf: a TxDb object of a gtf file that transcripts will be
#' extracted with `exonsBy(Gtf, by = "tx", use.names = TRUE)`, if
#' a GrangesList will use as is
#' @param ds numeric vector (NULL), disengagement score. If you have already
#'  calculated \code{\link{disengagementScore}}, input here to save time.
#' @return a named vector of numeric values of scores
#' @importFrom data.table rbindlist
#' @family features
#' @export
#' @examples
#' # Check inside outside score of a ORF within a transcript
#' ORF <- GRanges("1",
#'                ranges = IRanges(start = c(20, 30, 40),
#'                                   end = c(25, 35, 45)),
#'                strand = "+")
#'
#' grl <- GRangesList(tx1_1 = ORF)
#'
#' tx1 <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
#'                                 end = c(5, 15, 25, 35, 45, 200)),
#'                strand = "+")
#' tx <- GRangesList(tx1 = tx1)
#' RFP <- GRanges(seqnames = "1",
#'                   ranges = IRanges(start = c(1, 4, 30, 60, 80, 90),
#'                                    end = c(30, 33, 63, 90, 110, 120)),
#'                   strand = "+")
#'
#' insideOutsideORF(grl, RFP, tx)
#'
insideOutsideORF <- function(grl, RFP, GtfOrTx, ds = NULL) {

  if (is(GtfOrTx, "TxDb")) {
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = TRUE)
  } else if (is.grl(GtfOrTx)) {
    tx <- GtfOrTx
  } else {
    stop("GtfOrTx is neithter of type TxDb or GRangesList")
  }
  if (length(RFP) > 1e6) {
    RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])
  }

  overlapGrl <- countOverlaps(grl, RFP) + 1
  # find tx with hits
  validIndices <- hasHits(tx, RFP)
  validIndices <- validIndices[data.table::chmatch(txNames(grl), names(tx))]
  if (!any(validIndices)) { # if no hits
    names(overlapGrl) <- NULL
    return(overlapGrl)
  }
  tx <- tx[txNames(grl)][validIndices]
  grl <- grl[validIndices]

  grlStarts <- startSites(grl, asGR = FALSE, is.sorted = TRUE)
  upstreamTx <- upstreamOfPerGroup(tx, grlStarts, allowOutside = FALSE)
  overlapTxOutside <- rep(1, length(validIndices))
  if (!is.null(ds)) { # save time here if ds is defined
    downstreamCounts <- 1/(ds/overlapGrl)
    upstreamCounts <- rep(1, length(validIndices))
    upstreamCounts[validIndices] <- countOverlaps(upstreamTx, RFP)
    overlapTxOutside <- downstreamCounts + upstreamCounts

  } else { # else make ds again
    grlStops <- stopSites(grl, asGR = FALSE, is.sorted = TRUE)
    downstreamTx <- downstreamOfPerGroup(tx, grlStops)

    dtmerge <- data.table::rbindlist(l = list(as.data.table(upstreamTx),
                                              as.data.table(downstreamTx)))
    group <- NULL # for avoiding warning
    txOutside <- makeGRangesListFromDataFrame(
      dtmerge[order(group)], split.field = "group")

    overlapTxOutside[validIndices] <- countOverlaps(txOutside, RFP) + 1
  }

  scores <- overlapGrl / overlapTxOutside
  names(scores) = NULL
  return(scores)
}


#' Find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function
#' [distToCds()]
#'
#' possible outputs:
#' 0: orf is in frame with cds
#' 1: 1 shifted from cds
#' 2: 2 shifted from cds
#'
#' @references doi: 10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @return a logical vector
#' @family features
#' @examples
#' # simple example
#' isInFrame(c(3,6,8,11,15))
#'
#' # GRangesList example
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs, extension = 0)
#' isInFrame <- isInFrame(dist)
#' @export
#'
isInFrame <- function(dists){
  return((dists - 1) %% 3)
}


#' Find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function
#' [distToCds()]
#' @references doi: 10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @family features
#' @return a logical vector
#' @export
#' @examples
#' # simple example
#' isOverlapping(c(-3,-6,8,11,15))
#'
#' # GRangesList example
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs, extension = 0)
#' isOverlapping <- isOverlapping(dist)
isOverlapping <- function(dists) {
  return(dists < 0)
}


#' ORF rank in transcripts
#'
#' @references doi: 10.1074/jbc.R116.733899
#' @description ig. second orf _2 -> 2
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @return a numeric vector of integers
#' @family features
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' grl <- ORFik:::makeORFNames(grl)
#' rankOrder(grl)
rankOrder <- function(grl) {
  gr <- unlist(grl, use.names = FALSE)

  if (is.null(names(grl))) {
    if (is.null(gr$names)) {
      if (is.null(names(gr))) {
        stop("no valid names to find ranks")
      } else {
        orfName <- names(gr)
        if (length(orfName) > length(grl)) {
          orfName <- names(groupGRangesBy(gr, names(gr)))
        }
      }
    } else {
      orfName <- gr$names
      if (length(orfName) > length(grl)) {
        orfName <- names(groupGRangesBy(gr, gr$names))
      }
    }
  } else {
    orfName <- names(grl)
    if (suppressWarnings(anyNA(as.integer(sub(".*_", "", orfName,
                                              perl = TRUE))))) {
      if (!is.null(gr$names)) {
        orfName <- names(groupGRangesBy(gr, gr$names))
      }
    }
  }
  if (length(orfName) > length(grl)) {
    stop("did not find a valid column to find ranks, easiest way to fix is",
         " set grl to: ORFik:::groupGRangesBy(grl, names), ",
         "where names are the orf names with _* in them-")
  }

  if (is.null(orfName)) stop("grl must have column called names")
  ranks <- as.integer(sub(".*_", "", orfName, perl = TRUE))
  if (anyNA(ranks)) {
    stop("no valid names to find ranks, check for orf _* names eg.",
         "tx_1, tx_2.")
  }
  return(ranks)
}
