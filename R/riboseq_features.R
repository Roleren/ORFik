#' Create normalizations of overlapping read counts.
#'
#' FPKM is short for "Fragments Per Kilobase of transcript per Million
#' fragments in library". When calculating RiboSeq data FPKM over ORFs,
#' use ORFs as `grl`.
#' When calculating RNASeq data FPKM, use full transcripts as
#' `grl`. It is equal to RPKM given that you do not have paired end reads.
#'
#' Note also that you must consider if you will use the whole read
#' library or just the reads overlapping `grl`.
#' To only overlap do:
#' reads <- reads[countOverlaps(reads, grl) > 0]
#' @references doi: 10.1038/nbt.1621
#' @param grl a \code{\link{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs as a special case (uORFs, potential new cds' etc). If
#'  regions are not spliced you can send a \code{\link{GRanges}} object.
#' @param reads a GAlignment, GRanges or GRangesList object,
#'  usually of RiboSeq, RnaSeq, CageSeq, etc.
#' @param pseudoCount an integer, by default is 0, set it to 1 if you want to
#' avoid NA and inf values.
#' @return a numeric vector with the fpkm values
#' @export
#' @family features
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25),"+")
#' fpkm(grl, RFP)
#'
fpkm <- function(grl, reads, pseudoCount = 0) {
  if (is.gr_or_grl(grl)) {
    if(is.grl(grl)) {
      grl_len <- widthPerGroup(grl, FALSE)
    } else grl_len <- width(grl)
  } else stop("grl must be GRangesList or GRanges")

  overlaps <- countOverlaps(grl, reads)
  librarySize <- length(reads)
  return(fpkm_calc(overlaps, grl_len, librarySize) + pseudoCount)
}

#' Calucalte entropy value of overlapping input reads per GRanges.
#'
#' Calculates entropy of the `reads` coverage over each `grl` group.
#' The entropy value per group is a real number in the interval (0:1),
#' where 0 indicates no variance in reads over group.
#' For example c(0,0,0,0) has 0 entropy, since no reads overlap.
#' @param grl a \code{\link{GRangesList}} that the reads will
#' be overlapped with
#' @param reads a GAlignment, GRanges or GRangesList object,
#'  usually of RiboSeq, RnaSeq, CageSeq, etc.
#' @return A numeric vector containing one entropy value per element in
#' `grl`
#' @family features
#' @export
#' @examples
#' # a toy example with ribo-seq p-shifted reads
#' ORF <- GRanges("1", ranges = IRanges(start = c(1, 12, 22),
#'                                      end = c(10, 20, 32)),
#'                strand = "+",
#'                names = rep("tx1_1", 3))
#' names(ORF) <- rep("tx1", 3)
#' grl <- GRangesList(tx1_1 = ORF)
#' reads <- GRanges("1", IRanges(c(25, 35), c(25, 35)), "+")
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' entropy(grl, reads)
#' # or on cds
#' cdsORF <- GRanges("1", IRanges(35, 44), "+", names = "tx1")
#' names(cdsORF) <- "tx1"
#' cds <-  GRangesList(tx1 = cdsORF)
#' entropy(cds, reads)
#'
entropy <- function(grl, reads) {
  # Optimize: Get count list of only groups with hits
  validIndices <- hasHits(grl, reads)
  if (!any(validIndices)) { # no variance in countList, 0 entropy
    return(rep(0, length(validIndices)))
  }
  grl <- grl[validIndices]
  reOrdering <- uniqueOrder(grl)

  # entropy function, interval 0:1 real number
  # Xi is the ratio of hits per postion per group
  Xi <- codonSumsPerGroup(grl, reads)

  validXi <- Xi$codonSums > 0 # avoid log2(0)
  Xi[, `:=` (Hx = rep(0, nrow(Xi)))]
  Xi[validXi, Hx := codonSums * log2(codonSums)] # Hx: The codon sum part
  Xi <- Xi[, .(Hx = sum(Hx)), by = genes]

  codons <- numCodons(grl)
  MHx <- 1/codons
  Xi[, MHx := MHx * log2(MHx) * codons] # MHx: The length part
  Xi[, entropy := Hx / MHx] # entropy is read sums over lengths

  entropy <- rep(0.0, length(validIndices))
  # non 0 entropy values set to HX / MHX
  Xi[is.na(entropy), entropy := 0.]
  tempEntro <- Xi$entropy[reOrdering] # order back from unique
  entropy[validIndices] <- tempEntro # order back from hits
  return(entropy)
}

#' Fragment Length Organization Similarity Score
#'
#' This feature is usually calcualted only for RiboSeq reads. For reads of
#' width between `start` and `end`,
#' sum the fraction of RiboSeq reads (per widths)
#' that overlap ORFs and normalize by CDS.
#'
#' Pseudo explanation of the function:
#' \preformatted{
#' SUM[start to stop]((grl[start:end][name]/grl) / (cds[start:end][name]/cds))
#' }
#' Please read more in the article.
#' @references doi: 10.1016/j.celrep.2014.07.045
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @param RFP ribosomal footprints, given as Galignment or GRanges object,
#' must be already shifted and resized to the p-site
#' @param cds a \code{\link{GRangesList}} of coding sequences,
#' cds has to have names as grl so that they can be matched
#' @param start usually 26, the start of the floss interval
#' @param end usually 34, the end of the floss interval
#' @return a vector of FLOSS of length same as grl
#' @family features
#' @export
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 12, 22),
#'                end = c(10, 20, 32)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' # RFP is 1 width position based GRanges
#' RFP <- GRanges("1", IRanges(c(1, 25, 35, 38), width = 1), "+")
#' score(RFP) <- c(28, 28, 28, 29) # original width in score col
#' cds <-  GRangesList(tx1 = GRanges("1", IRanges(35, 44), "+"))
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' floss(grl, RFP, cds)
#' # or change ribosome start/stop, more strict
#' floss(grl, RFP, cds, 28, 28)
#'
floss <- function(grl, RFP, cds, start = 26, end = 34){

  if (start > end) stop("start is bigger than end")
  if (is.grl(class(RFP))) {
    stop("RFP must be either GAlignment or GRanges type")
  }
  # for orfs
  overlaps <- findOverlaps(grl, RFP)
  rfpWidth <- readWidths(RFP[to(overlaps)])
  rfpPassFilter <- (rfpWidth >= start) & (rfpWidth <= end)
  rfpValidMatch <- rfpWidth[rfpPassFilter]
  ORFGrouping <- from(overlaps)[rfpPassFilter]
  if (sum(as.numeric(ORFGrouping)) == 0) {
    return(as.numeric(rep(0, length(grl))))
  }
  whichNoHit <- NULL # which ribo-seq did not hit grl
  if (length(unique(ORFGrouping)) != length(grl)) {
    whichNoHit <- S4Vectors::setdiff.Vector(
      seq_along(grl), unique(ORFGrouping))
  }
  orfFractions <- split(rfpValidMatch, ORFGrouping)
  listing<- IRanges::RleList(orfFractions)
  tableFracs <- table(listing)
  colnames(tableFracs) <- NULL
  orfFractions <- lapply(seq_along(listing), function(x) {
    tableFracs[x,] / sum(tableFracs[x,])
  })

  # for cds
  overlapsCds <- findOverlaps(cds, RFP)
  rfpWidth <- readWidths(RFP[to(overlapsCds)])
  rfpPassFilterCDS <- ((rfpWidth >= start) & (rfpWidth <= end))
  rfpValidMatchCDS <- rfpWidth[rfpPassFilterCDS]
  cdsFractions <- split(rfpValidMatchCDS, rfpValidMatchCDS)
  totalLength <- length(rfpValidMatchCDS)
  cdsFractions <- vapply(cdsFractions, FUN.VALUE = c(1.0), FUN = function(x) {
    length(x) / totalLength
  })
  cdsFractions <- as.double(cdsFractions)
  # floss score ->
  score <- vapply(seq_along(orfFractions), FUN.VALUE = c(1.0),
                  FUN = function(x) {
                    sum(abs(orfFractions[[x]] - cdsFractions)) * 0.5
  })
  if (!is.null(whichNoHit)) {
    tempScores <- as.numeric(rep(NA, length(grl)))
    tempScores[unique(ORFGrouping)] <- score
    tempScores[whichNoHit] <- 0.
    score <- tempScores
  }
  if (length(score) != length(grl) || anyNA(score)) {
    stop("could not find floss-score for all objects, most",
         "likely objects are wrongly annotated.")
  }
  return(score)
}

#' Translational efficiency
#'
#' Uses RnaSeq and RiboSeq to get translational efficiency of every element in
#' `grl`. Translational efficiency is defined as:
#' \preformatted{
#' (density of RPF within ORF) / (RNA expression of ORFs transcript)
#' }
#' @references doi: 10.1126/science.1168978
#' @param grl a \code{\link{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs as a special case (uORFs, potential new cds' etc).
#' @param RNA RnaSeq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param RFP RiboSeq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param tx a GRangesList of the transcripts. If you used cage data, then
#' the tss for the the leaders have changed, therefor the tx lengths have
#' changed. To account for that call:
#' `
#' translationalEff(grl, RNA, RFP, tx = extendLeaders(tx, cageFiveUTRs))
#' ` where cageFiveUTRs are the reannotated by CageSeq data leaders.
#' @param with.fpkm logical F, if true return the fpkm values together with
#' translational efficiency
#' @param pseudoCount an integer, 0, set it to 1 if you want to avoid NA and
#' inf values. It also helps against bias from low depth libraries.
#' @return a numeric vector of fpkm ratios, if with.fpkm is TRUE, return a
#' data.table with te and fpkm values
#' @export
#' @importFrom data.table data.table
#' @family features
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25), "+")
#' RNA <- GRanges("1", IRanges(1, 50), "+")
#' tx <-  GRangesList(tx1 = GRanges("1", IRanges(1, 50), "+"))
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' te <- translationalEff(grl, RNA, RFP, tx, with.fpkm = TRUE, pseudoCount = 1)
#' te$fpkmRFP
#' te$te
#'
translationalEff <- function(grl, RNA, RFP, tx, with.fpkm = FALSE,
                             pseudoCount = 0) {
  tx <- tx[txNames(grl)]
  #normalize by tx lengths
  fpkmRNA <- fpkm(tx, RNA, pseudoCount)
  #normalize by grl lengths
  fpkmRFP <- fpkm(grl, RFP, pseudoCount)
  if (with.fpkm) {
    return(data.table(fpkmRFP = fpkmRFP, fpkmRNA = fpkmRNA,
                      te = fpkmRFP / fpkmRNA))
  }
  return(fpkmRFP / fpkmRNA)
}

#' Disengagement score (DS)
#'
#' Disengagement score is defined as
#' \preformatted{(RPFs over ORF)/(RPFs downstream to transcript end)}
#' A pseudo-count of one is added to both the ORF and downstream sums.
#' @references doi: 10.1242/dev.098344
#' @param grl a \code{\link{GRangesList}} object
#' with usually either leaders, cds', 3' utrs or ORFs.
#' @param RFP RiboSeq reads as GAlignment, GRanges
#' or GRangesList object
#' @param GtfOrTx If it is \code{\link{TxDb}} object
#'  transcripts will be extracted using
#'  \code{exonsBy(Gtf, by = "tx", use.names = TRUE)}.
#'  Else it must be \code{\link{GRangesList}}
#' @param RFP.sorted logical (F), an optimizer, have you ran this line:
#' \code{RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])}
#' Normally not touched, for internal optimization purposes.
#' @return a named vector of numeric values of scores
#' @export
#' @family features
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' tx <- GRangesList(tx1 = GRanges("1", IRanges(1, 50), "+"))
#' RFP <- GRanges("1", IRanges(c(1,10,20,30,40), width = 3), "+")
#' disengagementScore(grl, RFP, tx)
#'
disengagementScore <- function(grl, RFP, GtfOrTx, RFP.sorted = FALSE) {
  tx <- loadRegion(GtfOrTx)

  # exclude non hits and set them to 0
  validIndices <- hasHits(tx, RFP)
  validIndices <- validIndices[data.table::chmatch(txNames(grl), names(tx))]
  if (!any(validIndices)) { # if no hits
    score <- countOverlaps(grl, RFP) + 1
    names(score) <- NULL
    return(score)
  }
  overlapDownstream <- rep(1, length(grl))

  grlStops <- stopSites(grl[validIndices], asGR = FALSE, is.sorted = TRUE)
  downstreamTx <- downstreamOfPerGroup(tx[txNames(grl)][validIndices],
                                       grlStops)
  # check for big lists
  if (length(downstreamTx) > 5e5) {
    if(!RFP.sorted){
      RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])
    }
    ordering <- uniqueOrder(downstreamTx)
    downstreamTx <- uniqueGroups(downstreamTx)
    overlapDownstream[validIndices] <- countOverlaps(downstreamTx,
                                                     RFP)[ordering] + 1
  } else {
    overlapDownstream[validIndices] <- countOverlaps(downstreamTx, RFP) + 1
  }

  overlapGrl <- countOverlaps(grl, RFP) + 1

  score <- overlapGrl / overlapDownstream
  names(score) <- NULL
  return(score)
}

#' Inside/Outside score (IO)
#'
#' Inside/Outside score is defined as
#' \preformatted{(reads over ORF)/(reads outside ORF and within transcript)}
#' A pseudo-count of one is added to both the ORF and outside sums.
#' @references doi: 10.1242/dev.098345
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders, cds', 3' utrs or ORFs
#' @param RFP ribo seq reads as GAlignment, GRanges or GRangesList object
#' @param GtfOrTx if Gtf: a TxDb object of a gtf file that transcripts will be
#' extracted with `exonsBy(Gtf, by = "tx", use.names = TRUE)`, if
#' a GrangesList will use as is
#' @param ds numeric vector (NULL), disengagement score. If you have already
#'  calculated \code{\link{disengagementScore}}, input here to save time.
#' @param RFP.sorted logical (F), have you ran this line:
#' \code{RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])}
#' Normally not touched, for internal optimization purposes.
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
insideOutsideORF <- function(grl, RFP, GtfOrTx, ds = NULL,
                             RFP.sorted = FALSE) {
  tx <- loadRegion(GtfOrTx)

  if (!RFP.sorted) RFP <- optimizeReads(tx, RFP)

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
    downstreamCounts <- 1 / (ds / overlapGrl)
    upstreamCounts <- rep(1, length(validIndices))
    upstreamCounts[validIndices] <- countOverlaps(upstreamTx, RFP)
    overlapTxOutside <- downstreamCounts + upstreamCounts

  } else { # else make ds again
    grlStops <- stopSites(grl, asGR = FALSE, is.sorted = TRUE)
    downstreamTx <- downstreamOfPerGroup(tx, grlStops)

    overlapTxOutside[validIndices] <- countOverlaps(upstreamTx, RFP) +
      countOverlaps(downstreamTx, RFP) + 1
  }

  scores <- overlapGrl / overlapTxOutside
  names(scores) = NULL
  return(scores)
}

#' Ribosome Release Score (RRS)
#'
#' Ribosome Release Score is defined as
#' \preformatted{(RPFs over ORF)/(RPFs over 3' utrs)} and
#' additionaly normalized by lengths.
#' If RNA is added as argument, it will normalize by RNA counts
#' to justify location of 3' utrs.
#' It can be understood as a ribosome stalling feature.
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' @references doi: 10.1016/j.cell.2013.06.009
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders,
#'  cds', 3' utrs or ORFs.
#' @param RFP RiboSeq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param GtfOrThreeUtrs if Gtf: a TxDb object of a gtf file transcripts is
#'  called from: `threeUTRsByTranscript(Gtf, use.names = TRUE)`,
#'  if object is GRangesList, it is presumed to be the 3' utrs
#' @param RNA RnaSeq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param overlapGrl an integer vector of overlaps
#' (default: countOverlaps(grl, RFP) + 1)
#' @return a named vector of numeric values of scores, NA means that
#' no 3' utr was found for that transcript.
#' @export
#' @family features
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' threeUTRs <- GRangesList(tx1 = GRanges("1", IRanges(40, 50), "+"))
#' RFP <- GRanges("1", IRanges(25, 25), "+")
#' RNA <- GRanges("1", IRanges(1, 50), "+")
#' ribosomeReleaseScore(grl, RFP, threeUTRs, RNA)
#'
ribosomeReleaseScore <- function(grl, RFP, GtfOrThreeUtrs, RNA = NULL,
                                 overlapGrl = countOverlaps(grl, RFP) + 1){
  threeUTRs <- loadRegion(GtfOrThreeUtrs, part = "trailer")
  # check that naming is correct, else change it.
  orfNames <- txNames(grl, FALSE)
  validNamesThree <- names(threeUTRs) %in% orfNames
  validNamesGRL <- orfNames %in% names(threeUTRs)
  rrs <- rep(NA,length(grl))
  if (sum(validNamesGRL) != length(grl)) {
    threeUTRs <- threeUTRs[validNamesThree]
    grl <- grl[validNamesGRL]
  }

  threeUTRs <- threeUTRs[orfNames[validNamesGRL]]
  overlapThreeUtrs <- countOverlaps(threeUTRs, RFP) + 1

  rrs[validNamesGRL] <- (overlapGrl / widthPerGroup(grl)) /
    (overlapThreeUtrs / widthPerGroup(threeUTRs))

  if (!is.null(RNA)) { # normalize by rna ratio
    rnaRatio <- (countOverlaps(grl, RNA) + 1) /
      (countOverlaps(threeUTRs, RNA) + 1)
    rrs[validNamesGRL] <- rrs[validNamesGRL] / rnaRatio
  }
  names(rrs) <- NULL
  return(rrs)
}

#' Ribosome Stalling Score (RSS)
#'
#' Is defined as \preformatted{(RPFs over ORF stop sites)/(RPFs over ORFs)}
#' and normalized by lengths
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' @references doi: 10.1016/j.cels.2017.08.004
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders,
#'  cds', 3' utrs or ORFs.
#' @param RFP RiboSeq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param overlapGrl an integer vector of overlaps
#' (default: countOverlaps(grl, RFP))
#' @return a named vector of numeric values of RSS scores
#' @export
#' @family features
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25), "+")
#' ribosomeStallingScore(grl, RFP)
#'
ribosomeStallingScore <- function(grl, RFP,
                                  overlapGrl = countOverlaps(grl, RFP)) {
  grl_len <- widthPerGroup(grl, FALSE)
  stopCodons <- stopCodons(grl, is.sorted = TRUE)
  overlapStop <- countOverlaps(stopCodons, RFP)

  rss <- ((overlapStop + 1) / 3) / ((overlapGrl + 1) / grl_len)
  names(rss) <- NULL
  return(rss)
}

#' Start region coverage
#'
#' Get the number of reads in the start region of each ORF. If you want the
#' start codon coverage only, set upstream = 0. Standard is 2 upstream
#' and 2 downstream, a width 5 window centered at start site. since
#' p-shifting is not 100% accurate, this window is usually the reads from the
#' start site.
#'
#' If tx is null, then upstream will be force to 0 and downstream to
#' a maximum of grl width. Since there is no reference for splicing.
#' @param RFP ribo seq reads as GAlignment, GRanges or GRangesList object
#' @inheritParams startRegion
#' @family features
#' @return a numeric vector of counts
startRegionCoverage <- function(grl, RFP, tx = NULL, is.sorted = TRUE,
                                upstream = 2L, downstream = 2L) {
  region <- startRegion(grl, tx, is.sorted, upstream, downstream)
  return(countOverlaps(region, RFP))
}

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
#' @param reads ribosomal footprints, given as Galignment object or
#'  Granges
#' @param pShifted a logical (TRUE), are riboseq reads p-shifted?
#' @family features
#' @return an integer vector, 1 score per ORF, with names of grl
#' @export
#' @importFrom BiocGenerics Reduce
#' @examples
#' # Good hiting ORF
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(21, 40),
#'                strand = "+")
#' names(ORF) <- c("tx1")
#' grl <- GRangesList(tx1 = ORF)
#' # 1 width p-shifted reads
#' reads <- GRanges("1", IRanges(c(21, 23, 50, 50, 50, 53, 53, 56, 59),
#'                             width = 1), "+")
#' score(reads) <- 28 # original width
#' cds <- GRanges(seqnames = "1",
#'                ranges = IRanges(50, 80),
#'                strand = "+")
#' cds <- GRangesList(tx1 = cds)
#' tx <- GRanges(seqnames = "1",
#'                ranges = IRanges(1, 85),
#'                strand = "+")
#' tx <- GRangesList(tx1 = tx)
#'
#' initiationScore(grl, cds, tx, reads, pShifted = TRUE)
#'
initiationScore <- function(grl, cds, tx, reads, pShifted = TRUE) {
  if (length(grl) == 0) stop("grl must have length > 0")
  # meta coverage of cds
  cdsMeta <- windowPerReadLength(cds, tx, reads, pShifted = pShifted)

  # coverage per ORF
  prop <- windowPerReadLength(grl, tx, reads, pShifted = pShifted,
                              scoring = "fracPos")

  # find a better scoring pattern
  prop[, `:=` (dif = abs(score - cdsMeta$score))]
  len <- length(unique(prop$fraction))
  ans <- prop[, .(difPer = sum(dif)), by = list(fraction, genes)]
  ans <- ans[, .(score = sum(difPer)/len - 1), by = list(genes)]$score

  ans[is.na(ans) | is.nan(ans)] <- 0
  names(ans) <- names(grl)
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
#' NOTE: If reads are not of width 1, then a read from 1-4 on range of 1-4,
#' will get scores frame1 = 2, frame2 = 1, frame3 = 1. What could be logical
#' is that only the 5' end is important, so that only frame1 = 1,
#' to get this, you first resize reads to 5'end only.
#'
#' NOTE: p shifting is not exact, so some functional ORFs will get a
#' bad ORF score.
#' @references doi: 10.1002/embj.201488411
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @param RFP ribosomal footprints, given as Galignment object,
#'  Granges or GRangesList
#' @param is.sorted logical (F), is grl sorted.
#' @importFrom data.table .SD
#' @importFrom data.table .N
#' @family features
#' @export
#' @return a data.table with 4 columns, the orfscore (ORFScores) and score of
#' each of the 3 tiles (frame_zero_RP, frame_one_RP, frame_two_RP)
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
  if (any(widthPerGroup(grl, FALSE) < 3)) stop("width < 3 ORFs not allowed")

  counts <- coveragePerTiling(grl, RFP, is.sorted, as.data.table = TRUE,
                              withFrames = TRUE)
  total <- coverageScorings(counts, scoring = "frameSum")
  countsTile1 <- total[frame == 0,]$score
  countsTile2 <- total[frame == 1,]$score
  countsTile3 <- total[frame == 2,]$score

  RP = countsTile1 + countsTile2 + countsTile3
  Ftotal <- RP/3

  frame1 <- (countsTile1 - Ftotal)^2 / Ftotal
  frame2 <- (countsTile2 - Ftotal)^2 / Ftotal
  frame3 <- (countsTile3 - Ftotal)^2 / Ftotal

  dfORFs <- data.table(frame_zero_RP = countsTile1)
  dfORFs[, frame_one_RP := countsTile2]
  dfORFs[, frame_two_RP := countsTile3]

  ORFscore <- log2(frame1 + frame2 + frame3 + 1)
  revORFscore <-  which(frame1 < frame2 | frame1 < frame3)
  ORFscore[revORFscore] <- -1 * ORFscore[revORFscore]
  ORFscore[is.na(ORFscore)] <- 0
  dfORFs$ORFScores <- ORFscore
  dfORFs[] # for print
  return(dfORFs)
}
