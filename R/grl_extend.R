#' Extend the leaders transcription start sites.
#'
#' Will extend the leaders or transcripts upstream (5' end) by extension.
#' The extension is general not relative, that means splicing
#' will not be taken into account.
#' Requires the \code{grl} to be sorted beforehand,
#' use \code{\link{sortPerGroup}} to get sorted grl.
#' @param grl usually a \code{\link{GRangesList}} of 5' utrs or transcripts.
#' Can be used for any extension of groups.
#' @param extension an integer, how much to max extend upstream (5' end).
#' Either single value that will apply for all, or same as length of grl
#' which will give 1 update value per grl object.
#' Or a GRangesList where start / stops by strand are the positions
#' to use as new starts.
#' Will not cross the chromosome boundary for non circular chromosomes.
#' @param cds a \code{\link{GRangesList}} of coding sequences,
#' If you want to extend 5' leaders downstream, to catch
#' upstream ORFs going into cds, include it. It will add first
#' cds exon to grl matched by names.
#' Do not add for transcripts, as they are already included.
#' @param is.circular logical, default FALSE if not any is: all(isCircular(grl) %in% TRUE).
#' Where grl is the ranges checked. If TRUE, allow ranges to extend
#' below position 1 on chromosome. Since circular genomes can have negative coordinates.
#' @return an extended GRangeslist
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' library(GenomicFeatures)
#' samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#' txdb <- loadDb(samplefile)
#' fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE) # <- extract only 5' leaders
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' cds <- cdsBy(txdb,"tx",use.names = TRUE)
#' ## extend leaders upstream 1000
#' extendLeaders(fiveUTRs, extension = 1000)
#' ## now try(extend upstream 1000, add all cds exons):
#' extendLeaders(fiveUTRs, extension = 1000, cds)
#'
#' ## when extending transcripts, don't include cds' of course,
#' ## since they are already there
#' extendLeaders(tx, extension = 1000)
#' ## Circular genome (allow negative coordinates)
#' circular_fives <- fiveUTRs
#' isCircular(circular_fives) <- rep(TRUE, length(isCircular(circular_fives)))
#' extendLeaders(circular_fives, extension = 32672841L)
#'
extendLeaders <- function(grl, extension = 1000L, cds = NULL,
                          is.circular = all(isCircular(grl) %in% TRUE)) {
  if (is(extension, "numeric") && length(extension) %in% c(1L, length(grl))) {
    posIndices <- strandBool(grl)
    promo <- promoters(unlist(firstExonPerGroup(grl), use.names = FALSE),
                       upstream = extension)
    newStarts <- rep(NA, length(grl))
    newStarts[posIndices] <- as.integer(start(promo[posIndices]))
    newStarts[!posIndices] <- as.integer(end(promo[!posIndices]))
  } else if (is.grl(class(grl))) {
    starts <- startSites(extension)
    changedGRL <- downstreamFromPerGroup(grl[names(extension)], starts, is.circular)
    return(changedGRL)
  } else {
    stop("extension must either be an integer, or a GRangesList")
  }

  extendedLeaders <- assignFirstExonsStartSite(grl, newStarts, is.circular)
  if(is.null(cds)) return (extendedLeaders)
  return(addCdsOnLeaderEnds(extendedLeaders, cds))
}

#' Extend the Trailers transcription stop sites
#'
#' Will extend the trailers or transcripts downstream (3' end) by extension.
#' The extension is general not relative, that means splicing
#' will not be taken into account.
#' Requires the \code{grl} to be sorted beforehand,
#' use \code{\link{sortPerGroup}} to get sorted grl.
#' @param grl usually a \code{\link{GRangesList}} of 3' utrs or transcripts.
#' Can be used for any extension of groups.
#' @param extension an integer, how much to max extend downstream (3' end).
#' Either single value that will apply for all, or same as length of grl
#' which will give 1 update value per grl object.
#' Or a GRangesList where start / stops sites by strand are the positions
#' to use as new starts.
#' Will not cross the chromosome boundary for non circular chromosomes.
#' @inheritParams extendLeaders
#' @return an extended GRangeslist
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' library(GenomicFeatures)
#' samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#' txdb <- loadDb(samplefile)
#' threeUTRs <- threeUTRsByTranscript(txdb) # <- extract only 5' leaders
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' ## now try(extend downstream 1000):
#' extendTrailers(threeUTRs, extension = 1000)
#' ## Or on transcripts
#' extendTrailers(tx, extension = 1000)
#' ## Circular genome (allow negative coordinates)
#' circular_three <- threeUTRs
#' isCircular(circular_three) <- rep(TRUE, length(isCircular(circular_three)))
#' extendTrailers(circular_three, extension = 126200008L)[41] # <- negative stop coordinate
#'
extendTrailers <- function(grl, extension = 1000L, is.circular =
                             all(isCircular(grl) %in% TRUE)) {
  if (is(extension, "numeric") && length(extension) %in% c(1L, length(grl))) {
    posIndices <- strandBool(grl)
    promo <- flank(unlist(lastExonPerGroup(grl), use.names = FALSE),
                   width = extension, start = FALSE)
    newEnds <- rep(NA, length(grl))
    newEnds[posIndices] <- as.integer(end(promo[posIndices]))
    newEnds[!posIndices] <- as.integer(start(promo[!posIndices]))
  } else if (is.grl(class(extension))) {
    starts <- startSites(extension)
    changedGRL <-upstreamOfPerGroup(grl[names(extension)], starts,
                                    allowOutside = TRUE, is.circular)
    return(changedGRL)
  } else {
    stop("extension must either be an integer, or a GRangesList")
  }
  return(assignLastExonsStopSite(grl, newEnds, is.circular = is.circular))
}

#' Distance to following range group
#'
#' Follow means downstream GRangesList element, not including itself.
#' @param grl a GRangesList
#' @param grl2 a GRangesList, default 'grl'. The list that defines
#' restrictions on extension. Can also be another set, which is
#' used as 'roadblocks' for extension.
#' @param ignore.strand logical, default FALSE
#' @return numeric vector of distance
distanceToFollowing <- function(grl, grl2 = grl, ignore.strand = FALSE) {
  stops <- unlistGrl(stopRegion(grl, downstream = 0, upstream = 0))
  starts <- unlistToExtremities(grl2)
  if (!ignore.strand){
    followers <- precede(stops, starts)
  } else {
    minus_strand <- which(stops@strand == "-")
    plus_strand <- which(stops@strand == "+")
    followers_plus <- precede(stops[plus_strand], starts, ignore.strand = TRUE)
    followers_minus <- follow(stops[minus_strand], starts, ignore.strand = TRUE)
    followers <- c()
    if (length(minus_strand > 0)) followers[minus_strand] <- followers_minus
    if (length(plus_strand > 0)) followers[plus_strand] <- followers_plus
  }
  nas <- is.na(followers)
  dists <- distance(stops[!is.na(followers)], starts[followers[!is.na(followers)]], ignore.strand = ignore.strand)
  out <- seq_along(grl)
  out[!nas] <- dists
  out[nas] <- NA
  out
}

#' Distance to preceding range group
#'
#' Preceding means upstream GRangesList element, not including itself.
#' @param grl a GRangesList
#' @param grl2 a GRangesList, default 'grl'. The list that defines
#' restrictions on extension. Can also be another set, which is
#' used as 'roadblocks' for extension.
#' @param ignore.strand logical, default FALSE
#' @return numeric vector of distance
distanceToPreceding <- function(grl, grl2 = grl, ignore.strand = FALSE) {
  stops <- unlistToExtremities(grl2)
  starts <- unlistGrl(startRegion(grl, downstream = 0,upstream = 0))
  if (!ignore.strand){
    preceders <- follow(starts, stops)
  } else {
    minus_strand <- which(starts@strand == "-")
    plus_strand <- which(starts@strand == "+")
    preceders_plus <- follow(starts[plus_strand], stops, ignore.strand = TRUE)
    preceders_minus <- precede(starts[minus_strand], stops, ignore.strand=TRUE)
    preceders <- c()
    if (length(minus_strand > 0)) preceders[minus_strand] <- preceders_minus
    if (length(plus_strand > 0)) preceders[plus_strand] <- preceders_plus
  }
  nas <- is.na(preceders)
  dists <- distance(starts[!is.na(preceders)], stops[preceders[!is.na(preceders)]], ignore.strand = ignore.strand)
  out <- seq_along(grl)
  out[!nas] <- dists
  out[nas] <- NA
  out

}

#' Extend Trailers Until
#'
#' Extend trailers until a restriction group / position. This makes you extend until
#' you hit another gene boundary etc.
#' @inheritParams distanceToFollowing
#' @inheritParams extendTrailers
#' @param until numeric, default 200. The nearest you can go to the boundary.
#' #' Defined as boundary hit + 1, so if hit on 22, and until is 2, will set to 22+2+1 = 25.
#' Usually if Leaders/trailers are not defined,
#' this makes a good pseudo leader boundary around your other genes.
#' @param min_ext numeric, default 25. What is the minimum extension, even though it crosses a boundary.
#' Will not cross the chromosome boundary for non circular chromosomes.
#' @param ... Arguments sent to distanceToFollowing
#' @return a GRangesList of extended grl input
#' @export
#' @examples
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(c(10, 15), c(13, 20)), "+"), tx2 = GRanges("1", IRanges(30, 50), "+"))
#' extendTrailersUntil(grl, min_ext = 5)
#' extendTrailersUntil(grl, min_ext = 5, until = 1)
#' extendTrailersUntil(grl, min_ext = 5, until = 1)
#' extendTrailersUntil(grl, min_ext = 5, until = 1, extension = 4)
extendTrailersUntil <- function(grl, grl2=grl, extension = 500, until = 200, min_ext = 25,
                                is.circular = all(isCircular(grl) %in% TRUE), ...) {
  dists <- distanceToFollowing(grl, grl2, ...)
  dists[is.na(dists)] <- extension + until + 1
  diff <- pmax(dists - until, min_ext)
  extension <- pmin(extension, diff)

  return(extendTrailers(grl, extension, is.circular))
}

#' Extend Leaders Until
#'
#' Extend leaders until a restriction group / position. This makes you extend until
#' you hit another gene boundary etc.
#' @inheritParams distanceToFollowing
#' @inheritParams extendLeaders
#' @param until numeric, default 200. The nearest you can go to the neighbour
#' boundaries of grl2 (the "other" genes).
#' Defined as boundary hit + 1, so if hit other gene with distance 22, and
#' 'until' argument is 2, will set final extension to 22-2-1 = 19.
#' Usually if Leaders/trailers are not defined,
#' this makes a good pseudo leader boundary around your other genes.
#' @param min_ext numeric, default 25. What is the minimum extension, even though it crosses a boundary.
#' Will not cross the chromosome boundary for non circular chromosomes.
#' @param ... Arguments sent to distanceToPreceding
#' @return a GRangesList of extended grl input
#' @export
#' @examples
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(c(10, 15), c(13, 20)), "+"), tx2 = GRanges("1", IRanges(30, 50), "+"))
#' extendLeadersUntil(grl, min_ext = 5)
#' extendLeadersUntil(grl, min_ext = 5, until = 1)
#' extendLeadersUntil(grl, min_ext = 5, until = 1)
#' extendLeadersUntil(grl, min_ext = 5, until = 1, extension = 4)
extendLeadersUntil <- function(grl, grl2=grl, extension = 500, until = 200, min_ext = 25,
                               is.circular = all(isCircular(grl) %in% TRUE), ...) {
  dists <- distanceToPreceding(grl, grl2, ...)
  dists[is.na(dists)] <- extension + until + 1
  diff <- pmax(dists - until, min_ext)
  extension <- pmin(extension, diff)

  return(extendLeaders(grl, extension, is.circular = is.circular))
}
