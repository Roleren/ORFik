#' Reassign the start positions of the first exons per group in grl
#' @description make sure your grl is sorted, since start of "-" strand
#' objects should be the
#' max end in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStarts an integer vector of same length as grl, with new start
#' values
#' @return the same GRangesList with new start sites
#' @family GRangesPositions
#'
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
#' @family GRangesPositions
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
#' @family GRangesPositions
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
#' @family GRangesPositions
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
#' @family GRangesPositions
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
