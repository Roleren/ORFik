#' Get ORF score for given GRangesList object grouped by orfs
#' That is, of the 3 possible frame in an ORF
#' Is the first one most important, by how much ?
#' NB! Only support + and - strand, not *
#' @param grl a GRangesList object with ORFs
#' @param RFP ribozomal footprints, given as Galignment object
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
#' matching is done by transcript names.
#' fiveUTRs must be used to make transcript positions possible.
#' The cds start site, will be presumed to be on + 1 of end of fiveUTRs
#' @param ORFs orfs as GRangesList, names of orfs must be transcript names
#' @param fiveUTRs fiveUTRs as GRangesList, must be original unchanged fiveUTRs
#' @param cds cds' as GRangesList
#' @param extension if you changed the extension of tss when finding the orfs
distOrfToCds <- function(ORFs, fiveUTRs, cds, extension = 1000){
  if (class(ORFs) != "GRangesList") stop("grl must be GRangesList Object")

  extendedLeadersWithoutCds <- makeGrlAndFilter(
    extendsTSSexons(fiveUTRs, extension = extension), fiveUTRs)
  extendedLeaders <- addFirstCdsOnLeaderEnds(
    extendedLeadersWithoutCds, cds)
  lastExons <-  lastExonPerGroup(ORFs)
  orfsTX <- asTX(lastExons, extendedLeaders)
  orfEnds <- lastExonEndPerGroup(orfsTX, F)
  cdsStarts <- widthPerGroup(extendedLeadersWithoutCds[
    OrfToTxNames(lastExons)], F) + 1
  dists <- cdsStarts - orfEnds

  return(dists)
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

#' Get the orf number in transcripts, ig. second orf _1 -> 2
#' @param grl a GRangesList object with ORFs
OrfRankOrder <- function(grl){
  orfName <- unlist(grl, use.names = F)$names
  if(is.null(orfName)) stop("grl must have column called names, with orf names")
  return(as.integer(gsub(".*_", "", orfName)))
}
