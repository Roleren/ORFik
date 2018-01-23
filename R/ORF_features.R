#' Get ORF score for a GRangesList of ORFs
#' @description That is, of the 3 possible frame in an ORF
#' Is the first one most important, by how much ?
#' NB! Only support + and - strand, not *
#' See article: 10.1002/embj.201488411
#' @param grl a GRangesList object with ORFs
#' @param RFP ribozomal footprints, given as Galignment object,
#'  Granges or GRangesList
#' @return a matrix with 4 columns, the orfscore and score of
#' each of the 3 tiles
ORFScores <- function(grl, RFP){

  sortedTilex <- tile1(grl)

  dt <- as.data.table(sortedTilex)

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
  dfORFs$ORFscore <- ORFscore

  return(dfORFs)
}

#' Get distances between orf ends and starts of transcripts cds' belonging to orfs.
#' @description matching is done by transcript names.
#' fiveUTRs must be used to make transcript positions possible.
#' The cds start site, will be presumed to be on + 1 of end of fiveUTRs
#' @param ORFs orfs as GRangesList, names of orfs must be transcript names
#' @param fiveUTRs fiveUTRs as GRangesList, must be original unchanged fiveUTRs
#' @param cds cds' as GRangesList, only add if you used cage extension
#' @param extension needs to be set! set to 0 if you did not use cage
#'  if you used cage to change tss' when finding the orfs, standard cage
#'  extension is 1000
#'  @return an integer vector, +1 means one base upstream of cds, -1 means
#'   2nd base in cds, 0 means orf stops at cds start.
distOrfToCds <- function(ORFs, fiveUTRs, cds = NULL, extension = NULL){
  if (class(ORFs) != "GRangesList") stop("ORFs must be GRangesList Object")
  if(is.null(extension)) stop("please specify extension, to avoid bugs\n
                              ,if you did not use cage, set it to 0,\n
                              standard cage extension is 1000")

  if(extension > 0){
    if(class(cds) != "GRangesList"){
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
  if(extension > 0){
    cdsStarts <- widthPerGroup(extendedLeadersWithoutCds[
      OrfToTxNames(lastExons)], F) + 1
  } else {
      cdsStarts <- widthPerGroup(fiveUTRs[
        OrfToTxNames(lastExons)], F) + 1
  }
  dists <- cdsStarts - orfEnds

  return(dists)
}

#' Make a score for each ORFs start region
#' @description The closer the sequence is to the kozak sequence
#' The higher the score, based on simplification of PWM
#' score system: 4 upstream, 5 downstream of start
#' CACCATGGC, 1+3+1+2, skip ATG, +2+1 = 10
#' CGCCATGGC, 1+!2+1+2, skip ATG, +2+1 = 9
#' Inspired by experimental bit values for each position
#' @param grl a GRangesList grouped by ORF
#' @param faFile a FaFile from the fasta file, see ?FaFile
#' @param species which species to use, currently only support human
kozakSequenceScore <- function(grl, faFile, species = "human"){
  firstExons <- firstExonPerGroup(grl)
  kozakLocation <- promoters(firstExons, upstream = 4, downstream = 5)

  sequences <- as.character(txSeqsFromFa(grl, faFile))
  names(sequences) <- NULL
  scores <- rep(0, length(sequences))
  if(species == "human"){
    # split strings and relist as letters of 9 rows
    subSplit <- strsplit(sequences, split = "")
    # this will not when ATG is on start of chr
    mat <- matrix(unlist(subSplit), nrow = 9)
    pos1 <- as.character(mat[1,])
    pos2 <- as.character(mat[2,])
    pos3 <- as.character(mat[2,])
    pos4 <- as.character(mat[2,])
    pos8 <- as.character(mat[8,])
    pos9 <- as.character(mat[9,])


    match <- grepl(x = pos1, pattern = "C|G|A")
    scores[match] <- scores[match] + 1

    match <- pos2 == "A"
    scores[match] <- scores[match] + 5
    match <- pos2 == "G"
    scores[match] <- scores[match] + 2
    match <- pos2 == "C"
    scores[match] <- scores[match] + 1

    match <- grepl(x = pos3, pattern = "C|A")
    scores[match] <- scores[match] + 1

    match <- pos4 == "C"
    scores[match] <- scores[match] + 2
    match <- pos4 == "G"
    scores[match] <- scores[match] + 1
    match <- pos4 == "A"
    scores[match] <- scores[match] + 1

    match <- pos8 == "G"
    scores[match] <- scores[match] + 2
    match <- pos8 == "A"
    scores[match] <- scores[match] + 1
    match <- pos8 == "T"
    scores[match] <- scores[match] + 1

    match <- grepl(x = pos9, pattern = "C|A")
    scores[match] <- scores[match] + 1

    return(scores)
  }
  else stop("other species are not supported")
}

#' Inside/outside score (IO)
#'
#' is defined as (RPFs over ORF)/(RPFs downstream to tx end).
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1242/dev.098345
#' @param grl a GRangesList object with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param GtfOrTx if Gtf: a TxDb object of a gtf file,
#'  if tx: a GrangesList of transcripts, called from:
#'  exonsBy(Gtf, by = "tx", use.names = T)
#' @importFrom data.table rbindlist
#' @return a named vector of numeric values of scores
insideOutsideORF <- function(grl, RFP, GtfOrTx){
  overlapGrl <- countOverlaps(grl, RFP) + 1

  if(class(GtfOrTx) == "TxDb"){
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = T)
  } else if(class(GtfOrTx) == "GRangesList") {
    tx <- GtfOrTx
  } else {
    stop("GtfOrTx is neithter of type TxDb or GRangesList")
  }
  tx <- tx[OrfToTxNames(grl, F)]

  grlStarts <- ORFStartSites(grl,asGR = F)
  grlStops <- ORFStopSites(grl,asGR = F)

  downstreamTx <- downstreamOfPerGroup(tx, grlStops)
  upstreamTx <- upstreamOfPerGroup(tx, grlStarts)
  dtd <- as.data.table(downstreamTx)
  dtu <- as.data.table(upstreamTx)
  dtmerge <- data.table::rbindlist(l = list(dtu, dtd))
  txOutside <- makeGRangesListFromDataFrame(
    dtmerge[order(group)], split.field = "group",
    names.field = "group_name", keep.extra.columns = T)
  names(txOutside) <- names(tx)

  overlapTxOutside <- countOverlaps(txOutside, RFP) + 1
  scores <- overlapGrl / overlapTxOutside
  names(scores) = NULL
  return(scores)
}

#' find frame for each orf relative to cds
#' @param dists a vector of distances between ORF and cds
#' Input of this function, is the output of the function distOrfToCds
inFrameWithCDS <- function(dists){
  return(dists %% 3)
}

#' find for each orf if it overlaps cds range
#' @param dists a vector of distances between ORF and cds
#' Input of this function, is the output of the function distOrfToCds
isOverlappingCds <- function(dists){
  return(dists < 0)
}

#' Get the orf rank in transcripts
#' @description ig. second orf _2 -> 2
#' @param grl a GRangesList object with ORFs
OrfRankOrder <- function(grl){
  gr <- unlist(grl, use.names = FALSE)

  if(is.null(names(grl))){
    if(is.null(gr$names)){
      if(is.null(names(gr))){
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
  if(anyNA(ranks)){
    stop("no valid names to find ranks, check for orf _* names, i.g tx_1, tx_2.")
  }
  return(ranks)
}
