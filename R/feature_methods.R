#' Create normalizations of overlapping read counts.
#'
#' FPKM is short for "Fragments Per Kilobase of transcript per Million
#' fragments". When calculating RiboSeq data FPKM over ORFs use ORFs as
#' `grl`. When calculating RNASeq data FPKM use full transcripts as
#' `grl`.
#' @references doi: 10.1038/nbt.1621
#' @param grl a \code{\link{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs as a special case (uORFs, potential new cds' etc).
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
  grl_len <- widthPerGroup(grl, FALSE)
  overlaps <- countOverlaps(grl, reads)
  librarySize <- length(reads)
  return(fpkm_calc(overlaps, grl_len, librarySize) + pseudoCount)
}


#' Subset GRanges to get coverage.
#'
#' GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#' @param cov A coverage object from coverage()
#' @param y GRanges object for which coverage should be extracted
#' @return numeric vector of coverage of input GRanges object
#' @family features
#'
subsetCoverage <- function(cov, y) {
  cov1 <- cov[[as.vector(seqnames(y)[1])]]
  return(as.vector(cov1[ranges(y)]))
}

#' Calucalte entropy value of overlapping input reads.
#'
#' Calculates entropy of the `reads` coverage over each `grl` group.
#' The entropy value per group is a real number in the interval (0:1),
#' where 0 indicates no variance in reads over group.
#' For example c(0,0,0,0) has 0 entropy, since no reads overlap.
#' @param grl a \code{\link{GRangesList}} that the reads will
#' be overlapped with
#' @param reads a GAlignment object or GRanges or GRangesList, usualy data from
#' RiboSeq or RnaSeq
#' @return A numeric vector containing one entropy value per element in
#' `grl`
#' @family features
#' @export
#' @examples
#' ORF <- GRanges("1", ranges = IRanges(start = c(1, 12, 22),
#'                                      end = c(10, 20, 32)),
#'                strand = "+",
#'                names = rep("tx1_1", 3))
#' names(ORF) <- rep("tx1", 3)
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(c(25, 35), c(25, 35)), "+")
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' entropy(grl, RFP)
#' # or on cds
#' cdsORF <- GRanges("1", IRanges(35, 44), "+", names = "tx1")
#' names(cdsORF) <- "tx1"
#' cds <-  GRangesList(tx1 = cdsORF)
#' entropy(cds, RFP)
#'
entropy <- function(grl, reads) {
  # Get count list of groups with hits
  validIndices <- hasHits(grl, reads)
  if (!any(validIndices)) { # no variance in countList, 0 entropy
    return(rep(0, length(validIndices)))
  }
  # get coverage per group
  grl <- grl[validIndices]
  reOrdering <- uniqueOrder(grl)
  countList <- coveragePerTiling(uniqueGroups(grl), reads, is.sorted = TRUE,
                                 keep.names = FALSE)

  # generate the entropy variables
  # Number of hits per group
  N <- sum(countList)
  # length per group
  L <- BiocGenerics::lengths(countList)
  reg_len <- rep.int(3, length(L)) # tuples per orf, start on 3

  # which ORFs should not have tuplet length of 3
  # is hits < length
  bool <- N <= L
  floors <- floor(L / N)
  reg_len[bool] <- floors[bool]
  reg_counts <- as.integer(L / reg_len) # number of tuplets per orf
  runLengths <- reg_counts

  tripletSums <- codonSumsPerGroup(countList, reg_len, runLengths)

  # expand N for easy vectorization
  N <- rep.int(N, reg_counts)

  # entropy function, interval 0:1 real number
  # Xi is the ratio of hits per postion per group
  Xi <-  tripletSums / N
  validXi <- !is.nan(Xi * log2(Xi)) # avoid log2(0)
  Hx <- rep(0, length(validXi))
  Hx[validXi] <- Xi[validXi] * log2(Xi[validXi])

  MHx <- rep(0, length(validXi))

  reg_counts <- rep.int(reg_counts, runLengths)
  MHx <- 1/reg_counts * log2(1 / reg_counts)

  # sum the mhx to groups

  grouping <- rep.int(seq_along(runLengths), runLengths)

  Hx <- sum(NumericList(split(Hx, grouping)))
  MHx <- sum(NumericList(split(MHx, grouping)))

  entropy <- rep(0.0, length(validIndices))
  # non 0 entropy values set to HX / MHX
  tempEntro <- Hx / MHx
  tempEntro[is.na(tempEntro)] <- 0.
  tempEntro <- tempEntro[reOrdering] # order back from unique
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
    tableFracs[x,]/sum(tableFracs[x,])
  })

  # for cds
  overlapsCds <- findOverlaps(cds, RFP)
  rfpWidth <- readWidths(RFP[to(overlapsCds)])
  rfpPassFilterCDS <- ((rfpWidth >= start) & (rfpWidth <= end))
  rfpValidMatchCDS <- rfpWidth[rfpPassFilterCDS]
  cdsFractions <- split(rfpValidMatchCDS, rfpValidMatchCDS)
  totalLength <- length(rfpValidMatchCDS)
  cdsFractions <- vapply(cdsFractions, FUN.VALUE = c(1.0), FUN = function(x) {
    length(x)/totalLength
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


#' Fraction Length
#' @description Fraction Length is defined as
#' \preformatted{(lengths of grl)/(length of tx_len)}
#' so that each group in
#' the grl is divided by the corresponding transcript.
#' @references doi: 10.1242/dev.098343
#' @param grl a \code{\link{GRangesList}} object
#' with usually either leaders,
#' cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param tx_len the transcript lengths of the transcripts,
#' a named (tx names) vector of integers.
#' If you have the transcripts as GRangesList,
#' call `ORFik:::widthPerGroup(tx, TRUE)`.
#'
#' If you used CageSeq to reannotate leaders, then the tss for the the leaders
#' have changed, therefore the tx lengths have changed. To account for that
#' call: `tx_len <- widthPerGroup(extendLeaders(tx, cageFiveUTRs))`
#' and calculate graction length using `fractionLength(grl, tx_len)`.
#' @return a numeric vector of ratios
#' @family features
#' @export
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' tx <-  GRangesList(tx1 = GRanges("1", IRanges(1, 50), "+"))
#' fractionLength(grl, ORFik:::widthPerGroup(tx, keep.names = TRUE))
#'
fractionLength <- function(grl, tx_len){
  grl_len <- widthPerGroup(grl, FALSE)
  tx_len <- tx_len[txNames(grl)]
  names(tx_len) <- NULL
  return(grl_len / tx_len)
}


#' Disengagement score (DS)
#'
#' Disengagement score is defined as
#' \preformatted{(RPFs over ORF)/(RPFs downstream to tx end)}
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
disengagementScore <- function(grl, RFP, GtfOrTx){

  if (is(GtfOrTx,"TxDb")) {
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = TRUE)
  } else if (is.grl(GtfOrTx)) {
    tx <- GtfOrTx
  } else {
    stop("GtfOrTx is neithter of type TxDb or GRangesList")
  }
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
    RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])
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
ribosomeReleaseScore <- function(grl, RFP, GtfOrThreeUtrs, RNA = NULL){
  if (is(GtfOrThreeUtrs,"TxDb")) {
    threeUTRs <- threeUTRsByTranscript(GtfOrThreeUtrs, by = "tx",
                                       use.names = TRUE)
  } else if (is.grl(GtfOrThreeUtrs)) {
    threeUTRs <- GtfOrThreeUtrs
  } else {
    stop("GtfOrThreeUtrs is neithter of type TxDb or GRangesList")
  }
  # check that naming is correct, else change it.
  orfNames <- txNames(grl, FALSE)
  validNamesThree <- names(threeUTRs) %in% orfNames
  validNamesGRL <- orfNames %in% names(threeUTRs)
  rrs <- rep(NA,length(grl))
  if (sum(validNamesGRL) != length(grl)) {
    threeUTRs <- threeUTRs[validNamesThree]
    grl <- grl[validNamesGRL]
  }
  overlapGrl <- countOverlaps(grl, RFP) + 1
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
ribosomeStallingScore <- function(grl, RFP){
  grl_len <- widthPerGroup(grl, FALSE)
  overlapGrl <- countOverlaps(grl, RFP)
  stopCodons <- stopCodons(grl, is.sorted = TRUE)
  overlapStop <- countOverlaps(stopCodons, RFP)

  rss <- ((overlapStop + 1) / 3) / ((overlapGrl + 1) / grl_len)
  names(rss) <- NULL
  return(rss)
}

#' Get all possible features in ORFik
#'
#' If you want to get all the features easily, you can use this function.
#' Each feature have a link to an article describing its creation and idea
#' behind it. Look at the functions in the feature family to see all of them.
#'
#' If you used CageSeq to reannotate your leaders your txDB object, must
#' contain the reassigned leaders. In the future release reasignment will
#' create txdb objects for you, but currently this is not supported,
#' therefore be carefull.
#'
#' @param grl a \code{\link{GRangesList}} object
#'  with usually ORFs, but can also be either leaders, cds', 3' utrs, etc.
#' @param RFP RiboSeq reads as GAlignment, GRanges or GRangesList object
#' @param RNA RnaSeq reads as GAlignment, GRanges or GRangesList object
#' @param Gtf a TxDb object of a gtf file,
#' @param faFile a FaFile or BSgenome from the fasta file, see ?FaFile
#' @param riboStart usually 26, the start of the floss interval, see ?floss
#' @param riboStop usually 34, the end of the floss interval
#' @param orfFeatures a logical, is the grl a list of orfs?
#' @param includeNonVarying a logical, if TRUE, include all features not
#' dependent on RiboSeq data and RNASeq data, that is: Kozak,
#' fractionLengths, distORFCDS, isInFrame, isOverlapping and rankInTx
#' @param grl.is.sorted logical (F), a speed up if you know argument grl
#'  is sorted, set this to TRUE.
#' @return a data.table with scores, each column is one score type, name of
#' columns are the names of the scores, i.g [floss()]
#' or [fpkm()]
#' @importFrom data.table data.table
#' @export
#' @family features
#' @examples
#' # Usually the ORFs are found in orfik, which makes names for you etc.
#' # Here we make an example from scratch
#' gtf <- system.file("extdata", "annotations.gtf",
#' package = "ORFik") ## location of the gtf file
#' suppressWarnings(txdb <-
#'                   GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf"))
#' # use cds' as ORFs for this example
#' ORFs <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
#' ORFs <- makeORFNames(ORFs) # need ORF names
#' # make Ribo-seq data,
#' RFP <- unlistGrl(firstExonPerGroup(ORFs))
#' suppressWarnings(computeFeatures(ORFs, RFP, Gtf = txdb))
#' # For more details see vignettes.
#'
computeFeatures <- function(grl, RFP, RNA = NULL,  Gtf = NULL, faFile = NULL,
                            riboStart = 26, riboStop = 34, orfFeatures = TRUE,
                            includeNonVarying = TRUE, grl.is.sorted = FALSE) {
  #### Check input and load data ####
  validGRL(class(grl), "grl")
  checkRFP(class(RFP))
  checkRNA(class(RNA))
  if(!is(Gtf,"TxDb")) stop("gtf must be TxDb object")

  # get transcript parts
  fiveUTRs <- fiveUTRsByTranscript(Gtf, use.names = TRUE)
  cds <- cdsBy(Gtf, by = "tx", use.names = TRUE)
  threeUTRs <- threeUTRsByTranscript(Gtf, use.names = TRUE)
  tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)

  if (!grl.is.sorted) {
    grl <- sortPerGroup(grl)
  }

  #### Get all features, append 1 at a time, to save memory ####
  scores <- data.table(floss = floss(grl, RFP, cds, riboStart, riboStop))
  scores[, entropyRFP := entropy(grl, RFP)]
  scores[, disengagementScores := disengagementScore(grl, RFP, tx)]
  scores[, RRS := ribosomeReleaseScore(grl, RFP, threeUTRs, RNA)]
  scores[, RSS := ribosomeStallingScore(grl, RFP)]

  if (includeNonVarying) {
    scores[, fractionLengths := fractionLength(grl, widthPerGroup(tx, TRUE))]
  }

  if (!is.null(RNA)) { # if rna seq is included
    TE <- translationalEff(grl, RNA, RFP, tx, with.fpkm = TRUE)
    scores[, te := TE$te]
    scores[, fpkmRFP := TE$fpkmRFP]
    scores[, fpkmRNA := TE$fpkmRNA]
  } else {
    scores[, fpkmRFP := fpkm(grl, RFP)]
  }
  if (orfFeatures) { # if features are found for orfs
    scores[, ORFScores := orfScore(grl, RFP, grl.is.sorted)$ORFScores]
    scores[, ioScore := insideOutsideORF(grl, RFP, tx,
                                         scores$disengagementScores)]

    if (includeNonVarying) {

      if (is(faFile, "FaFile") || is(faFile, "BSgenome")) {
        scores$kozak <- kozakSequenceScore(grl, faFile)
      } else {
        message("faFile not included, skipping kozak sequence score")
      }
      # switch five with tx, is it possible to use ?
      distORFCDS <- distToCds(grl, fiveUTRs, cds)
      scores[, distORFCDS := distORFCDS]
      scores[, inFrameCDS := isInFrame(distORFCDS)]
      scores[, isOverlappingCds := isOverlapping(distORFCDS)]
      scores[, rankInTx := rankOrder(grl)]
    }
  } else {
    message("orfFeatures set to False, dropping all orf features.")
  }
  return(scores)
}

#' Get all possible features in ORFik
#'
#' If you have a txdb with correct lists, use:
#' [computeFeatures()]
#'
#' A specialized version if you don't have a correct txdb, for example with
#' CAGE reassigned leaders while txdb is not updated.
#' It is 2x faster for tested data.
#' If you do have a txdb with e.g. cage reassignments, use computeFeatures
#' instead.
#' The point of this function is to give you the ability to input
#' transcript etc directly into the function, and not load them from txdb.
#' Each feature have a link to an article describing feature,
#' try ?floss
#' @inheritParams computeFeatures
#' @param tx a GrangesList of transcripts,
#'  normally called from: exonsBy(Gtf, by = "tx", use.names = T)
#'  only add this if you are not including Gtf file
#'  You do not need to reassign these to the cage peaks, it will do it for you.
#' @param fiveUTRs fiveUTRs as GRangesList, if you used cage-data to
#'  extend 5' utrs, remember to input CAGE assigned version and not original!
#' @param cds a GRangesList of coding sequences
#' @param threeUTRs  a GrangesList of transcript 3' utrs,
#'  normally called from: threeUTRsByTranscript(Gtf, use.names = T)
#' @importFrom data.table data.table
#' @family features
#' @return a data.table with scores, each column is one score type, name of
#'  columns are the names of the scores, i.g [floss()]
#'  or [fpkm()]
#' @examples
#'  # a small example without cage-seq data:
#'  # we will find ORFs in the 5' utrs
#'  # and then calculate features on them
#'  \dontrun{
#'  if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'   library(GenomicFeatures)
#'   # Get the gtf txdb file
#'   txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'   package = "GenomicFeatures")
#'   txdb <- loadDb(txdbFile)
#'
#'   # Extract sequences of fiveUTRs.
#'   fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)[1:10]
#'   faFile <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#'   # need to suppress warning because of bug in GenomicFeatures, will
#'   # be fixed soon.
#'   tx_seqs <- suppressWarnings(extractTranscriptSeqs(faFile, fiveUTRs))
#'
#'   # Find all ORFs on those transcripts and get their genomic coordinates
#'   fiveUTR_ORFs <- findMapORFs(fiveUTRs, tx_seqs)
#'   unlistedORFs <- unlistGrl(fiveUTR_ORFs)
#'   # group GRanges by ORFs instead of Transcripts
#'   fiveUTR_ORFs <- groupGRangesBy(unlistedORFs, unlistedORFs$names)
#'
#'   # make some toy ribo seq and rna seq data
#'   starts <- unlistGrl(ORFik:::firstExonPerGroup(fiveUTR_ORFs))
#'   RFP <- promoters(starts, upstream = 0, downstream = 1)
#'   score(RFP) <- rep(29, length(RFP)) # the original read widths
#'
#'   # set RNA seq to duplicate transcripts
#'   RNA <- unlistGrl(exonsBy(txdb, by = "tx", use.names = TRUE))
#'
#'   computeFeaturesCage(grl = fiveUTR_ORFs, orfFeatures =  TRUE, RFP = RFP,
#'    RNA = RNA, Gtf = txdb, faFile = faFile)
#'
#' }
#' # See vignettes for more examples
#' }
#'
computeFeaturesCage <- function(grl, RFP, RNA = NULL,  Gtf = NULL, tx = NULL,
                                fiveUTRs = NULL, cds = NULL, threeUTRs = NULL,
                                faFile = NULL, riboStart = 26, riboStop = 34,
                                orfFeatures = TRUE, includeNonVarying = TRUE,
                                grl.is.sorted = FALSE) {
  #### Check input and load data ####
  validGRL(class(grl))
  checkRFP(class(RFP))
  checkRNA(class(RNA))
  if (is.null(Gtf)) {
    validGRL(c(class(fiveUTRs), class(cds), class(threeUTRs)),
             c("fiveUTRs", "cds", "threeUTRs"))
    if (!is.grl(tx)) { stop("if Gtf is not given,\n
                            tx must be specified as a GRangesList")
    }
    } else {
      if(!is(Gtf,"TxDb")) stop("gtf must be TxDb object")

      notIncluded <- validGRL(c(class(fiveUTRs), class(cds),
                                class(threeUTRs), class(tx)),
                              c("fiveUTRs","cds", "threeUTRs", "tx"), TRUE)
      if (notIncluded[1]) {
        fiveUTRs <- fiveUTRsByTranscript(Gtf, use.names = TRUE)
      }
      if (notIncluded[2]) {
        cds <- cdsBy(Gtf, by = "tx", use.names = TRUE)
      }
      if (notIncluded[3]) {
        threeUTRs <- threeUTRsByTranscript(Gtf, use.names = TRUE)
      }
      if (notIncluded[4]) {
        tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
      }
  }
  if (!grl.is.sorted) {
    grl <- sortPerGroup(grl)
  }

  #### Get all features, append 1 at a time, to save memory ####
  scores <- data.table(floss = floss(grl, RFP, cds, riboStart, riboStop))
  scores[, entropyRFP := entropy(grl, RFP)]
  scores[, disengagementScores := disengagementScore(grl, RFP, tx)]
  scores[, RRS := ribosomeReleaseScore(grl, RFP, threeUTRs, RNA)]
  scores[, RSS := ribosomeStallingScore(grl, RFP)]

  if (includeNonVarying) {
    scores[, fractionLengths := fractionLength(grl, widthPerGroup(tx, TRUE))]
  }

  if (!is.null(RNA)) { # if rna seq is included
    TE <- translationalEff(grl, RNA, RFP, tx, with.fpkm = TRUE)
    scores[, te := TE$te]
    scores[, fpkmRFP := TE$fpkmRFP]
    scores[, fpkmRNA := TE$fpkmRNA]
  } else {
    scores[, fpkmRFP := fpkm(grl, RFP)]
  }
  if (orfFeatures) { # if features are found for orfs
    scores[, ORFScores := orfScore(grl, RFP, grl.is.sorted)$ORFScores]
    scores[, ioScore := insideOutsideORF(grl, RFP, tx,
                                         scores$disengagementScores)]

    if (includeNonVarying) {

      if (is(faFile, "FaFile") || is(faFile, "BSgenome")) {
        scores$kozak <- kozakSequenceScore(grl, faFile)
      } else {
        message("faFile not included, skipping kozak sequence score")
      }
      # switch five with tx, is it possible to use ?
      distORFCDS <- distToCds(grl, fiveUTRs, cds)
      scores[, distORFCDS := distORFCDS]
      scores[, inFrameCDS := isInFrame(distORFCDS)]
      scores[, isOverlappingCds := isOverlapping(distORFCDS)]
      scores[, rankInTx := rankOrder(grl)]
    }
  } else {
    message("orfFeatures set to False, dropping all orf features.")
  }
  return(scores)
}
