#' Reassign the start positions of the first exons per group in grl
#'
#' Per group in GRangesList, assign the most upstream site.
#'
#' make sure your grl is sorted, since start of "-" strand
#' objects should be the
#' max end in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStarts an integer vector of same length as grl, with new start
#' values (absolute coordinates, not relative)
#' @return the same GRangesList with new start sites
#' @family GRanges
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
#'
#' Per group in GRangesList, assign the most upstream site.
#'
#' make sure your grl is sorted, since stop of "-" strand objects
#' should be the min start in group, use ORFik:::sortPerGroup(grl) to get
#' sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStops an integer vector of same length as grl,
#'  with new start values (absolute coordinates, not relative)
#' @return the same GRangesList with new stop sites
#' @importFrom data.table .N .I
#' @family GRanges
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
#' is the new start point of first valid exon. Can also be a GRangesList,
#' then stopsites will be used.
#' @return a GRangesList of downstream part
#' @family GRanges
#'
downstreamOfPerGroup <- function(tx, downstreamOf) {
  if (is.grl(downstreamOf)) {
    downstreamOf <- stopSites(downstreamOf, asGR = TRUE, is.sorted = TRUE,
                              keep.names = TRUE)
  } else if (is.numeric(downstreamOf) & length(downstreamOf) == length(tx)) {
    downstreamOf <- IRanges(downstreamOf, width  = 1)
    if (is.null(names(tx))) names(tx) <- seq_along(tx)
    names(downstreamOf) <- names(tx)
  } else stop("downstreamOf must be GRangesList, or numeric of equal",
              "size to tx")
  return(windowPerGroup(downstreamOf, tx, upstream = -1,
                      downstream = max(widthPerGroup(tx, FALSE))))
}
# downstreamOfPerGroup <- function(tx, downstreamOf) {
#   # Needs speed update!
#   posIndices <- strandBool(tx)
#   posEnds <- end(tx[posIndices])
#   negEnds <- start(tx[!posIndices])
#   posDown <- downstreamOf[posIndices]
#   negDown <- downstreamOf[!posIndices]
#   pos <- posEnds > posDown
#   neg <- negEnds < negDown
#   posTx <- tx[posIndices][pos]
#   negTx <- tx[!posIndices][neg]
#   downTx <- tx
#   downTx[posIndices] <- posTx
#   downTx[!posIndices] <- negTx
#   #check if anyone hits boundary, set those to boundary
#   if (anyNA(strandPerGroup(downTx, FALSE))) {
#     boundaryHits <- which(is.na(strandPerGroup(downTx, FALSE)))
#     downTx[boundaryHits] <- firstExonPerGroup(tx[boundaryHits])
#     ir <- IRanges(start = downstreamOf[boundaryHits],
#                   end = downstreamOf[boundaryHits])
#     irl <- split(ir, seq_along(ir))
#     names(irl) <- names(tx[boundaryHits])
#     ranges(downTx[boundaryHits]) <- irl
#   }
#   # check boundaries within group exons
#   startSites <- startSites(downTx, FALSE, FALSE, TRUE)
#   posChecks <- startSites[posIndices] > downstreamOf[posIndices] & any(!pos)
#   negChecks <- startSites[!posIndices] < downstreamOf[!posIndices] & any(!neg)
#   if (any(posChecks)) {
#     downstreamOf[posIndices][posChecks] <- startSites[posIndices][posChecks]
#   }
#   if (any(negChecks)) {
#     downstreamOf[!posIndices][negChecks] <- startSites[!posIndices][negChecks]
#   }
#
#   return(assignFirstExonsStartSite(downTx, downstreamOf))
# }

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
#' @family GRanges
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
#' @family GRanges
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

  # Usually from pos/neg with possible cage extensions
  if (allowOutside) {
    outside <- which(sum(pos) == 0)
    pos[outside] = TRUE
    posTx[outside] <- firstExonPerGroup(posTx[outside])
    outside <- which(sum(neg) == 0)
    neg[outside] = TRUE
    negTx[outside] <- firstExonPerGroup(negTx[outside])
  }
  upTx <- tx
  upTx[posIndices] <- posTx[pos]
  upTx[!posIndices] <- negTx[neg]

  nonZero <- widthPerGroup(upTx) > 0
  if (all(!nonZero)) { # if no ranges exists
    return(upTx)
  }
  upstreamOf <- upstreamOf[nonZero]
  oldPosIndices <- posIndices
  posIndices <- posIndices[nonZero]

  stopSites <- stopSites(upTx[nonZero], FALSE, FALSE, TRUE)
  # check boundaries within group exons
  if (any(posIndices)){
    posChecks <- stopSites[posIndices] < upstreamOf[posIndices] &
      any(!pos[nonZero[oldPosIndices]])
  } else {
    posChecks <- FALSE
  }
  if(any(!posIndices)){
    negChecks <- stopSites[!posIndices] > upstreamOf[!posIndices] &
      any(!neg[nonZero[!oldPosIndices]])
  } else {
    negChecks <- FALSE
  }

  if (any(posChecks)) {
    upstreamOf[posIndices][posChecks] <- stopSites[posIndices][posChecks]
  }
  if (any(negChecks)) {
    upstreamOf[!posIndices][negChecks] <- stopSites[!posIndices][negChecks]
  }

  upTx[nonZero] <- assignLastExonsStopSite(upTx[nonZero], upstreamOf)
  return(upTx)
}

#' Get rest of objects upstream (inclusive)
#'
#' Per group get the part upstream of position.
#' upstreamFromPerGroup(tx, stopSites(fiveUTRs, asGR = TRUE))
#' will return the  5' utrs per transcript as GRangesList,
#' usually used for interesting
#' parts of the transcripts.
#'
#' If you don't want to include the points given in the region,
#' use \code{\link{upstreamOfPerGroup}}
#' @param tx a \code{\link{GRangesList}},
#'  usually of Transcripts to be changed
#' @param upstreamFrom a vector of integers, for each group in tx, where
#' is the new start point of first valid exon.
#' @return a GRangesList of upstream part
#' @family GRanges
#'
upstreamFromPerGroup <- function(tx, upstreamFrom) {
  posIndices <- strandBool(tx)
  posStarts <- start(tx[posIndices])
  negStarts <- end(tx[!posIndices])
  posGrlStarts <- upstreamFrom[posIndices]
  negGrlStarts <- upstreamFrom[!posIndices]
  pos <- posStarts <= posGrlStarts
  neg <- negStarts >= negGrlStarts
  upTx <- tx
  upTx[posIndices] <- upTx[posIndices][pos]
  upTx[!posIndices] <- upTx[!posIndices][neg]
  # check if any hits boundary, set those to boundary
  if (anyNA(strandPerGroup(upTx, FALSE))) {
    boundaryHits <- which(is.na(strandPerGroup(upTx, FALSE)))
    upTx[boundaryHits] <- firstExonPerGroup(tx[boundaryHits])
    ir <- IRanges(start = upstreamFrom[boundaryHits],
                  end = upstreamFrom[boundaryHits])
    irl <- split(ir, seq_along(ir))
    names(irl) <- names(tx[boundaryHits])
    ranges(upTx[boundaryHits]) <- irl
  }

  return(assignLastExonsStopSite(upTx, upstreamFrom))
}
