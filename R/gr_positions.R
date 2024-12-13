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
#' @inheritParams extendLeaders
#' @return the same GRangesList with new start sites
#' @family GRanges
#' @keywords internal
assignFirstExonsStartSite <- function(grl, newStarts, is.circular =
                                        all(isCircular(grl) %in% TRUE)) {
  if (length(grl) != length(newStarts)) stop("length of grl and newStarts ",
                                             "are not equal!")
  posIndices <- strandBool(grl)
  if (!is.circular) {
    if (any(newStarts < 1)) {
      message("Transcript found that would be extended below coordinate position 0, setting to 1.")
      newStarts <- pmax(newStarts, 1)
    }
  }

  group <- NULL
  if (!is.null(grl@unlistData$group)) {
    group <- grl@unlistData$group
    grl@unlistData$group <- NULL
  }
  dt <- as.data.table(grl)
  dt[!duplicated(group) & posIndices[group], start := newStarts[posIndices]]
  if (is.circular) { # For negative strand, make failsafe on seqlength
    dt[!duplicated(group) & !posIndices[group], end := newStarts[!posIndices]]
  } else { # Not circular, check if seqlengths exist
    if (all(!is.na(seqlengths(grl)))) { # All seqlengths seq
      seqlengths.per <- seqlengths(grl)[seqnamesPerGroup(grl[!posIndices], FALSE)]
      dt[!duplicated(group) & !posIndices[group], end := pmin(newStarts[!posIndices],
                                                       seqlengths.per)]
    } else {
      dt[!duplicated(group) & !posIndices[group], end := newStarts[!posIndices]]
    }
  }

  ngrl <-
    GenomicRanges::makeGRangesListFromDataFrame(dt,
                                                split.field = "group",
                                                names.field = "group_name",
                                                keep.extra.columns = TRUE,
                                                seqinfo = seqinfo(grl))
  names(ngrl) <- names(grl)
  if (!is.null(group)) ngrl@unlistData$group <- group
  return(ngrl)
}


#' Reassign the stop positions of the last exons per group
#'
#' Per group in GRangesList, assign the most downstream site.
#'
#' make sure your grl is sorted, since stop of "-" strand objects
#' should be the min start in group, use ORFik:::sortPerGroup(grl) to get
#' sorted grl.
#' @param grl a \code{\link{GRangesList}} object
#' @param newStops an integer vector of same length as grl,
#'  with new start values (absolute coordinates, not relative)
#' @inheritParams extendLeaders
#' @return the same GRangesList with new stop sites
#' @importFrom data.table .N .I
#' @family GRanges
#' @keywords internal
assignLastExonsStopSite <- function(grl, newStops, is.circular =
                                      all(isCircular(grl) %in% TRUE)) {
  if (length(grl) != length(newStops)) stop("length of grl and newStops ",
                                            "are not equal!")
  posIndices <- strandBool(grl)
  if (!is.circular) {
    if (any(newStops < 1)) {
      message("Transcript found that would be extended below coordinate position 0, setting to 1.")
      newStops <- pmax(newStops, 1)
    }
  }
  group <- NULL
  if (!is.null(grl@unlistData$group)) {
    group <- grl@unlistData$group
    grl@unlistData$group <- NULL
  }
  dt <- as.data.table(grl)
  idx = dt[, .I[.N], by = group]


  dt[idx$V1]$start[!posIndices] <- newStops[!posIndices]
  if (is.circular) { # For positive strand, make failsafe on seqlength
    dt[idx$V1]$end[posIndices] <- newStops[posIndices]
  } else { # Not circular, check if seqlengths exist
    if (all(!is.na(seqlengths(grl)))) { # All seqlengths seq
      seqlengths.per <- seqlengths(grl)[seqnamesPerGroup(grl[posIndices], FALSE)]
      dt[idx$V1]$end[posIndices] <- pmin(newStops[posIndices], seqlengths.per)
    } else {
      dt[idx$V1]$end[posIndices] <- newStops[posIndices]
    }
  }
  ngrl <-
    GenomicRanges::makeGRangesListFromDataFrame(dt,
                                                split.field = "group",
                                                names.field = "group_name",
                                                keep.extra.columns = TRUE,
                                                seqinfo = seqinfo(grl))
  names(ngrl) <- names(grl)
  if (!is.null(group)) ngrl@unlistData$group <- group

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
#' @keywords internal
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
#' @inheritParams extendLeaders
#' @return a GRangesList of downstream part
#' @family GRanges
#' @keywords internal
downstreamFromPerGroup <- function(tx, downstreamFrom, is.circular =
                                     all(isCircular(tx) %in% TRUE)) {
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

  return(assignFirstExonsStartSite(downTx, downstreamFrom, is.circular))
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
#' @inheritParams extendLeaders
#' @return a GRangesList of upstream part
#' @family GRanges
#' @keywords internal
upstreamOfPerGroup <- function(tx, upstreamOf, allowOutside = TRUE,
                               is.circular = all(isCircular(tx) %in% TRUE)) {
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

  upTx[nonZero] <- assignLastExonsStopSite(upTx[nonZero], upstreamOf,
                                           is.circular = is.circular)
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
#' @keywords internal
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
