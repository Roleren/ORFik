#' Get ORF score for a GRangesList of ORFs
#' @description That is, of the 3 possible frame in an ORF
#' Is the first one most important, by how much ?
#' NB! Only support + and - strand, not *
#' See article: 10.1002/embj.201488411
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object with ORFs
#' @param RFP ribozomal footprints, given as Galignment object,
#'  Granges or GRangesList
#' @importFrom data.table .SD
#' @importFrom data.table .N
#' @family features
#' @export
#' @return a matrix with 4 columns, the orfscore and score of
#' each of the 3 tiles
orfScore <- function(grl, RFP){
  # tile the orfs into a d.t for easy seperation
  dt <- as.data.table(tile1(grl))

  group <- NULL
  # seperate the three tiles, by the 3 frames
  tilex1 <- dt[, .SD[seq.int(1, .N, 3)], by = group]
  grl1 <- makeGRangesListFromDataFrame(tilex1,
    split.field = "group", names.field = "group_name")
  tilex2 <- dt[, .SD[seq.int(2, .N, 3)], by = group]
  grl2 <- makeGRangesListFromDataFrame(tilex2,
    split.field = "group", names.field = "group_name")
  tilex3 <- dt[, .SD[seq.int(3, .N, 3)], by = group]
  grl3 <- makeGRangesListFromDataFrame(tilex2,
    split.field = "group", names.field = "group_name")

  countsTile1 <- countOverlaps(grl1, RFP)
  countsTile2 <- countOverlaps(grl2, RFP)
  countsTile3 <- countOverlaps(grl3, RFP)

  RP = countsTile1 + countsTile2 + countsTile3

  Ftotal <- RP/3

  tile1 <- (countsTile1 - Ftotal)^2 / Ftotal
  tile2 <- (countsTile2 - Ftotal)^2 / Ftotal
  tile3 <- (countsTile3 - Ftotal)^2 / Ftotal

  dfORFs <- NULL
  dfORFs$frame_zero_RP <- countsTile1
  dfORFs$frame_one_RP <- countsTile2
  dfORFs$frame_two_RP <- countsTile3

  ORFscore <- log2(tile1 + tile2 + tile3 + 1)
  revORFscore <-  which(tile1 < tile2 | tile1 < tile3)
  ORFscore[revORFscore] <- -1 * ORFscore[revORFscore]
  ORFscore[is.na(ORFscore)] <- 0
  dfORFs$ORFScores <- ORFscore

  return(dfORFs)
}

#' Get distances between orf ends and starts of transcripts cds' belonging to orfs.
#'
#' Matching is done by transcript names.
#' fiveUTRs must be used to make transcript positions possible.
#' The cds start site, will be presumed to be on + 1 of end of fiveUTRs
#' See article:  10.1074/jbc.R116.733899
#' @param ORFs orfs as \code{\link[GenomicRanges]{GRangesList}},
#'  names of orfs must be transcript names
#' @param fiveUTRs fiveUTRs as \code{\link[GenomicRanges]{GRangesList}},
#'  must be original unchanged fiveUTRs
#' @param cds cds' as \code{\link[GenomicRanges]{GRangesList}},
#'  only add if you used cage extension
#' @param extension needs to be set! set to 0 if you did not use cage
#'  if you used cage to change tss' when finding the orfs, standard cage
#'  extension is 1000
#' @family features
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs, extension = 0)
#' ## Now dist contains distances as vector
#' @export
#' @return an integer vector, +1 means one base upstream of cds, -1 means
#'   2nd base in cds, 0 means orf stops at cds start.
distToCds <- function(ORFs, fiveUTRs, cds = NULL, extension = NULL){
  validGRL(class(ORFs), "ORFs")
  if (is.null(extension)) stop("please specify extension, to avoid bugs\n
                              ,if you did not use cage, set it to 0,\n
                              standard cage extension is 1000")

  if (extension > 0) {
    if (class(cds) != "GRangesList") {
      stop("cds must be GRangesList Object,\n
        when extension > 0, cds must be included")
    }
    extendedLeadersWithoutCds <- extendLeaders(fiveUTRs, extension)
    fiveUTRs <- addFirstCdsOnLeaderEnds(
      extendedLeadersWithoutCds, cds)
  }

  lastExons <-  lastExonPerGroup(ORFs)
  orfsTx <- asTX(lastExons, fiveUTRs)
  # this is ok, since it is tx not genomic ->
  orfEnds <- lastExonEndPerGroup(orfsTx, F)
  if (extension > 0) {
    cdsStarts <- widthPerGroup(extendedLeadersWithoutCds[
      txNames(lastExons)], F) + 1
  } else {
      cdsStarts <- widthPerGroup(fiveUTRs[
        txNames(lastExons)], F) + 1
  }
  dists <- cdsStarts - orfEnds

  return(dists)
}

#' Make a score for each ORFs start region
#'
#' The closer the sequence is to the kozak sequence
#' The higher the score, based on the experimental pwms from paper.
#' Minimum score is 0 (worst correlation), max is 1 (the best
#' base per column was chosen)
#'
#' See article: https://doi.org/10.1371/journal.pone.0108475
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} grouped by ORF
#' @param faFile a FaFile from the fasta file, see ?FaFile
#' @param species ("human"), which species to use,
#'  currently only support human and zebrafish. You can also specify a pfm
#'  here for your own, syntax is an rectangular integer matrix,
#'  where all columns must sum to the same value, normally 100.
#' @param include.N logical (F), if TRUE, allow N bases to be counted as hits,
#'  score will be average of the other bases. If True N bases will be
#'  added to pfm, automaticly, so dont include them.
#' @return a numeric vector with values between 0 and 1
#' @family features
#' @importFrom Biostrings PWM
#' @export
#' @return an integer vector, one score per orf
kozakSequenceScore <- function(grl, faFile, species = "human",
                               include.N = FALSE){
  firstExons <- firstExonPerGroup(grl)
  kozakLocation <- promoters(firstExons, upstream = 9, downstream = 6)

  sequences <- as.character(txSeqsFromFa(kozakLocation,
                                                 faFile, is.sorted = TRUE))
  if (!all(nchar(sequences) == 15)) {
    stop("not all ranges had valid kozak sequences length, not 15")
  }

  #template <- sapply(s, function(x) PWMscoreStartingAt(pwm, x))

  if (species == "human") {
    # human pfm, see article reference
    pfm <- t(matrix(as.integer(c(20,20,21,21,19,24,46,29,19,22,28,16,
                                 27,33,32,23,32,38,10,38,45,15,39,26,
                                 35,29,28,39,30,26,37,20,28,49,18,37,
                                 18,18,19,17,19,12,7,13,8,14,15,21)),
                    ncol = 4))
  } else if (species == "zebrafish") {
    # zebrafish pfm, see article reference
    pfm <- t(matrix(as.integer(c(29,26,28,26,22,35,62,39,28,24,27,17,
                                 21,26,24,16,28,32,5,23,35,12,42,21,
                                 25,24,22,33,22,19,28,17,27,47,16,34,
                                 25,24,26,25,28,14,5,21,10,17,15,28)),
                    ncol = 4))
  } else if(class(species) == "matrix"){
    # self defined pfm
    pfm <- species
  } else {
    stop("other species are not supported")
  }

  bases <- c("A", "C", "G", "T")


  rownames(pfm) <- bases
  pwm <- PWM(pfm)

  if(include.N){
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
  for(i in 1:ncol(mat)){
    for(n in 1:length(bases)){
      match <- mat[, i] == bases[n]
      scores[match] <- scores[match] + pwm[n, i]
    }
  }
  return(scores)
}

#' Inside/outside score (IO)
#'
#' is defined as (RPFs over ORF)/(RPFs downstream to tx end).
#' A pseudo-count of one was added to both the ORF and downstream sums.
#'
#' See article: 10.1242/dev.098345
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param GtfOrTx if Gtf: a TxDb object of a gtf file,
#'  if tx: a GrangesList of transcripts, called from:
#'  exonsBy(Gtf, by = "tx", use.names = T)
#' @importFrom data.table rbindlist
#' @family features
#' @export
#' @return a named vector of numeric values of scores
insideOutsideORF <- function(grl, RFP, GtfOrTx){
  overlapGrl <- countOverlaps(grl, RFP) + 1

  if (class(GtfOrTx) == "TxDb") {
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = T)
  } else if (class(GtfOrTx) == "GRangesList") {
    tx <- GtfOrTx
  } else {
    stop("GtfOrTx is neithter of type TxDb or GRangesList")
  }
  tx <- tx[txNames(grl, F)]

  grlStarts <- startSites(grl, asGR = FALSE)
  grlStops <- stopSites(grl, asGR = FALSE)

  downstreamTx <- downstreamOfPerGroup(tx, grlStops)
  upstreamTx <- upstreamOfPerGroup(tx, grlStarts)
  dtd <- as.data.table(downstreamTx)
  dtu <- as.data.table(upstreamTx)
  dtmerge <- data.table::rbindlist(l = list(dtu, dtd))
  group <- NULL # for avoiding warning
  txOutside <- makeGRangesListFromDataFrame(
    dtmerge[order(group)], split.field = "group",
    names.field = "group_name", keep.extra.columns = TRUE)
  names(txOutside) <- names(tx)

  overlapTxOutside <- countOverlaps(txOutside, RFP) + 1
  scores <- overlapGrl / overlapTxOutside
  names(scores) = NULL
  return(scores)
}

#' find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function distToCds
#' possible outputs:
#' 0: orf is in frame with cds
#' 1: 1 shifted from cds
#' 2: 2 shifted from cds
#'
#' See article:  10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @family features
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs, extension = 0)
#' isInFrame <- isInFrame(dist)
#' @export
#' @return a logical vector
isInFrame <- function(dists){

  return((dists - 1) %% 3)
}

#' find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function distToCds
#' See article:  10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @family features
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs, extension = 0)
#' isOverlapping <- isOverlapping(dist)
#' @export
#' @return a logical vector
isOverlapping <- function(dists){
  return(dists < 0)
}

#' Get the orf rank in transcripts
#'
#' See article:  10.1074/jbc.R116.733899
#' @description ig. second orf _2 -> 2
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object with ORFs
#' @family features
#' @export
#' @return a numeric vector of integers
rankOrder <- function(grl){
  gr <- unlist(grl, use.names = FALSE)

  if (is.null(names(grl))) {
    if (is.null(gr$names)) {
      if (is.null(names(gr))) {
        stop("no valid names to find ranks")
      } else {
          orfName <- names(gr)
          if(length(orfName) > length(grl)){
            orfName <- names(groupGRangesBy(grl, names(gr)))
          }
      }
    } else {
        orfName <- gr$names
        if(length(orfName) > length(grl)){
          orfName <- names(groupGRangesBy(grl, gr$names))
        }
    }
  } else {
      orfName <- names(grl)
      if(anyNA(as.integer(gsub(".*_", "", orfName)))){
        if(!is.null(gr$names)){
          orfName <- names(groupGRangesBy(grl, gr$names))
        }
      }
  }
  if(length(orfName) > length(grl)) stop("did not find a valid column\n
    to find ranks, easiest way to fix is set grl to:\n
      ORFik:::groupGRangesBy(grl, names),\n
        where names are the orf names with _* in them-")

  if(is.null(orfName)) stop("grl must have column called names, with orf names")
  ranks <- as.integer(gsub(".*_", "", orfName))
  if (anyNA(ranks)) {
    stop("no valid names to find ranks, check for orf _* names, i.g tx_1, tx_2.")
  }
  return(ranks)
}
