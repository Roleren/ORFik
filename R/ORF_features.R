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

#' get distances between orf ends and starts of transcripts cds' belonging to orfs.
#' matching is done by transcript names, so they must match.
#' @param ORFs orfs as GRangesList
#' @param cds cds as GRangesList
#' @param fiveUTRs fiveUTRs as GRangesList
distOrfToCds = function(ORFs, cds, fiveUTRs){
  stop("not working!")
  cdsFirstExons = firstExonPerGroup(cds)
  namesToUse = as.character(unlist(unique(seqnames(ORFs))))
  cdsToUse = cdsFirstExons[namesToUse]
  c = unlist(cdsToUse, use.names = F)
  cdsToTranscript = mapToTranscripts(c,fiveUTRs)
  c = unlist(cdsToUse)
  cdsToUse = cdsToTranscript[names(c[cdsToTranscript$xHits]) == seqnames(cdsToTranscript)]

  ORFsPos = ORFs[as.character(strand(ORFs)) == "+"]
  ORFsMin = ORFs[as.character(strand(ORFs)) == "-"]

  distPos = start(cdsToUse[as.character(strand(cdsToUse)) == "+"]) - as.integer(end(endsPos))
  distMin = as.integer(start(endsMin)) - end(cdsToUse[as.character(strand(cdsToUse)) == "-"])
  dists = rep(NA,length(names(unlist(ends, use.names = F))))
  dists[as.character(strand(unlist(ends))) == "+"] = distPos
  dists[as.character(strand(unlist(ends))) == "-"] = distMin
  return(dists)
}
