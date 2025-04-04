#' Create normalizations of overlapping read counts.
#'
#' FPKM is short for "Fragments Per Kilobase of transcript per Million
#' fragments in library". When calculating RiboSeq data FPKM over ORFs,
#' use ORFs as `grl`.
#' When calculating RNASeq data FPKM, use full transcripts as
#' `grl`. It is equal to RPKM given that you do not have paired end reads.
#'
#' Note also that you must consider if you will use the whole read
#' library or just the reads overlapping `grl` for library size.
#' A normal question here is, does it make sense to include rRNA in library
#' size ?
#' If you only want overlapping grl, do:
#' librarySize = "overlapping"
#' @references doi: 10.1038/nbt.1621
#' @param grl a \code{\link{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs as a special case (uORFs, potential new cds' etc). If
#'  regions are not spliced you can send a \code{\link{GRanges}} object.
#' @param reads a \code{\link{GAlignments}}, \code{\link{GRanges}} or
#' \code{\link{GRangesList}} object, usually of RiboSeq, RnaSeq, CageSeq, etc.
#' @param pseudoCount a numeric, default 0, set it to 1 if you want to
#' avoid NA and inf values.
#' @param librarySize either numeric value or character vector.
#' Default ("full"), number of alignments in library (reads).
#' If you just have a subset, you can give the value by
#' librarySize = length(wholeLib) or sum(wholeLib$score),
#' if you want lib size to be only number of reads overlapping grl, do:
#' librarySize = "overlapping"
#' sum(countOverlaps(reads, grl) > 0),
#' if reads[1] has 3 hits in grl, and reads[2] has 2 hits,
#' librarySize will be 2, not 5.
#' You can also get the inverse overlap,
#' if you want lib size to be total number of overlaps, do:
#' librarySize = "DESeq"
#' This is standard fpkm way of DESeq2::fpkm(robust = FALSE)
#' sum(countOverlaps(grl, reads))
#' if grl[1] has 3 reads and grl[2] has 2 reads, librarySize is 5, not 2.
#' @inheritParams getWeights
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
#' # With weights (10 reads at position 25)
#' RFP <- GRanges("1", IRanges(25, 25),"+", score = 10)
#' fpkm(grl, RFP, weight = "score")
#'
fpkm <- function(grl, reads, pseudoCount = 0, librarySize = "full",
                 weight = 1L) {
  if (is.gr_or_grl(grl)) {
    if(is.grl(grl)) {
      grl_len <- widthPerGroup(grl, FALSE)
    } else grl_len <- width(grl)
  } else stop("grl must be GRangesList or GRanges")

  overlaps <- countOverlapsW(grl, reads, weight)

  if (librarySize == "full") {
    librarySize <-
      if (is(reads, "covRle")) {
        sum(sum(reads))
      } else sum(getWeights(reads, weight))
  } else if (is.numeric(librarySize)) {
    # Use value directly
  } else if (librarySize == "overlapping") {
    librarySize <- sum(overlaps > 0)
    # librarySize <- sum((countOverlaps(reads, grl) > 0) *
    #                      getWeights(reads, weight))
  } else if (librarySize == "DESeq") {
    librarySize <- sum(overlaps)
  } else stop("librarySize must be numeric or full, overlapping, DESeq!")
  return(fpkm_calc(overlaps, grl_len, librarySize) + pseudoCount)
}

#' Percentage of maximum entropy
#'
#' Calculates percentage of maximum entropy of the `reads`
#' coverage over each ORF in `grl` group.
#' The entropy value per group is a real number in the interval (0:1),
#' where 0 indicates no variance in reads over all codons of group
#' For example c(0,0,0,0) has 0 entropy, since no reads overlap.\cr
#' Interval: [0]: No reads or all reads in 1 place \cr
#' Interval: [0.01-0.99]: >= 2 positions covered \cr
#' Interval: [1]: all positions covered perfectly in frame\cr
#'
#' @inheritParams fpkm
#' @param is.sorted logical (FALSE), is grl sorted. That is + strand groups in
#' increasing ranges (1,2,3), and - strand groups in decreasing ranges (3,2,1)
#' @param overlapGrl an integer, (default: NULL),
#' if defined must be countOverlaps(grl, RFP),
#' added for speed if you already have it.
#' @return A numeric vector containing one entropy value per element in
#' `grl`
#' @family features
#' @export
#' @examples
#' # a toy example with ribo-seq p-shifted reads
#' ORF <- GRangesList(tx1 = GRanges("1", IRanges(1, width = 9), "+"))
#' entropy(ORF, GRanges()) # 0
#' entropy(ORF, GRanges("1", IRanges(c(1)), "+")) # 0
#' entropy(ORF, GRanges("1", IRanges(c(1,4,6,7)), "+")) # 0.94
#' entropy(ORF, GRanges("1", IRanges(c(1,4,7)), "+", score = c(1,2,1)),
#'         weight = "score") # 0.94
#' entropy(ORF, GRanges("1", IRanges(c(1,4,7)), "+")) # Perfect = 1
entropy <- function(grl, reads, weight = 1L, is.sorted = FALSE,
                    overlapGrl = NULL) {
  # Optimize: Get count list of only groups with hits
  validIndices <- hasHits(grl, reads, overlaps = overlapGrl)
  if (!any(validIndices)) { # no variance in countList, 0 entropy
    return(rep(0, length(validIndices)))
  }
  grl <- grl[validIndices]
  #reOrdering <- uniqueOrder(grl)

  # entropy function, interval 0:1 real number
  # Xi is the ratio of hits per postion per group
  Xi <- codonSumsPerGroup(grl, reads, weight, is.sorted)

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
  entropy[validIndices] <- Xi$entropy # order back from hits
  return(entropy)
}

#' Fragment Length Organization Similarity Score
#'
#' This feature is usually calcualted only for RiboSeq reads. For reads of
#' width between `start` and `end`,
#' sum the fraction of RiboSeq reads (per read widths)
#' that overlap ORFs and normalize by CDS read width fractions.
#' So if all read length are width 34 in ORFs and CDS, value is 1.
#' If width is 33 in ORFs and 34 in CDS, value is 0.
#' If width is 33 in ORFs and 50/50 (33 and 34) in CDS, values will be
#' 0.5 (for 33).
#'
#' Pseudo explanation of the function:
#' \preformatted{
#' SUM[start to stop]((grl[start:end][name]/grl) / (cds[start:end][name]/cds))
#' }
#' Where 'name' is transcript names.\cr
#' Please read more in the article.
#' @references doi: 10.1016/j.celrep.2014.07.045
#' @inheritParams fpkm
#' @param RFP ribosomal footprints, given as \code{\link{GAlignments}}
#' or \code{\link{GRanges}} object,
#' must be already shifted and resized to the p-site. Requires a $size column
#' with original read lengths.
#' @param cds a \code{\link{GRangesList}} of coding sequences,
#' cds has to have names as grl so that they can be matched
#' @param start usually 26, the start of the floss interval (inclusive)
#' @param end usually 34, the end of the floss interval (inclusive)
#' @return a vector of FLOSS of length same as grl, 0 means no RFP reads
#' in range, 1 is perfect match.
#' @family features
#' @importFrom BiocGenerics weights
#' @export
#' @examples
#' ORF1 <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 12, 22),
#'                end = c(10, 20, 32)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF1)
#' # RFP is 1 width position based GRanges
#' RFP <- GRanges("1", IRanges(c(1, 25, 35, 38), width = 1), "+")
#' RFP$size <- c(28, 28, 28, 29) # original width in size col
#' cds <-  GRangesList(tx1 = GRanges("1", IRanges(35, 44), "+"))
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' floss(grl, RFP, cds)
#' # or change ribosome start/stop, more strict
#' floss(grl, RFP, cds, 28, 28)
#'
#' # With repeated alignments in score column
#' ORF2 <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(12, 22, 36),
#'                end = c(20, 32, 38)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF1, tx1_2 = ORF2)
#' score(RFP) <- c(5, 10, 5, 10)
#' floss(grl, RFP, cds, weight = "score")
#'
floss <- function(grl, RFP, cds, start = 26, end = 34, weight = 1L){

  if (start > end) stop("start is bigger than end")
  if (end < start) stop("end is smaller than start")
  if (is.grl(class(RFP))) {
    stop("RFP must be either GAlignment or GRanges type")
  }
  if (is(RFP, "covRle")) {
    return(rep(0, length(grl)))
  }
  # Get all weights
  dt <- data.table(weights = getWeights(RFP, weight),
                   widths  = readWidths(RFP))

  # for orfs
  overlaps <- findOverlaps(grl, RFP)
  dt.rfp <- dt[to(overlaps), ] # keep only matches to grl
  dt.rfp[, ORFGrouping := from(overlaps)]
  dt.rfp <-  dt.rfp[(widths >= start) & (widths <= end),] # matches to widths
  if (nrow(dt.rfp) == 0) { # if no valid hits
    return(as.numeric(rep(0, length(grl))))
  }
  # Get total reads per read length per valid ORF
  orfFractions <- dt.rfp[, .(weights = sum(weights)),
                         by = .(ORFGrouping, widths)]
  orfFractions <- orfFractions[, .(widths, fraction = weights / sum(weights)), by = .(ORFGrouping)]

  # for cds
  overlaps <- findOverlaps(cds, RFP)
  dt.cds <- dt[to(overlaps), ] # keep only matches to grl
  dt.cds[, CDSGrouping := from(overlaps)]
  dt.cds <-  dt.cds[(widths >= start) & (widths <= end),] # matches to widths
  if (nrow(dt.cds) == 0) { # if no valid hits
    warning("No valid reads over CDS, check annotation")
    return(as.numeric(rep(0, length(grl))))
  }
  # Get total reads per read length per valid ORF
  cdsFractions <- dt.cds[, .(weights = sum(weights)), by = .(widths)]
  cdsFractions <- cdsFractions[, .(widths, fraction = weights / sum(weights))]
  merged <- data.table::merge.data.table(orfFractions, cdsFractions, by = "widths")
  score <- merged[, sum(abs(fraction.x - fraction.y))*0.5, by = ORFGrouping]

  zero <- rep(0, length(grl)) # Put 0 for all without hits
  zero[score$ORFGrouping] <- score$V1
  return(zero)
}

#' Translational efficiency
#'
#' Uses RnaSeq and RiboSeq to get translational efficiency of every element in
#' `grl`. Translational efficiency is defined as:
#' \preformatted{
#' (density of RPF within ORF) / (RNA expression of ORFs transcript)
#' }
#' @references doi: 10.1126/science.1168978
#' @param RNA RnaSeq reads as \code{\link{GAlignments}},
#' \code{\link{GRanges}} or GRangesList object
#' @param RFP RiboSeq reads as \code{\link{GAlignments}},
#' \code{\link{GRanges}} or GRangesList object
#' @param tx a GRangesList of the transcripts. If you used cage data, then
#' the tss for the the leaders have changed, therefor the tx lengths have
#' changed. To account for that call:
#' `
#' translationalEff(grl, RNA, RFP, tx = extendLeaders(tx, cageFiveUTRs))
#' ` where cageFiveUTRs are the reannotated by CageSeq data leaders.
#' @param with.fpkm logical, default: FALSE, if true return the fpkm
#' values together with translational efficiency as a data.table
#' @inheritParams fpkm
#' @param weight.RFP a vector (default: 1L). Can also be character name of
#' column in RFP. As in translationalEff(weight = "score") for:
#' GRanges("chr1", 1, "+", score = 5), would mean score column tells
#' that this alignment region was found 5 times.
#' @param weight.RNA Same as weightRFP but for RNA weights.
#' (default: 1L)
#' @return a numeric vector of fpkm ratios, if with.fpkm is TRUE, return a
#' data.table with te and fpkm values (total 3 columns then)
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
                             pseudoCount = 0, librarySize = "full",
                             weight.RFP = 1L, weight.RNA = 1L) {
  tx <- tx[txNames(grl, tx)]
  #normalize by tx lengths
  fpkmRNA <- fpkm(tx, RNA, pseudoCount, librarySize, weight.RNA)
  #normalize by grl lengths
  fpkmRFP <- fpkm(grl, RFP, pseudoCount, librarySize, weight.RFP)
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
#' @param RFP RiboSeq reads as GAlignments, GRanges
#' or GRangesList object
#' @param GtfOrTx If it is \code{\link{TxDb}} object
#'  transcripts will be extracted using
#'  \code{exonsBy(Gtf, by = "tx", use.names = TRUE)}.
#'  Else it must be \code{\link{GRangesList}}
#' @param RFP.sorted logical (FALSE), an optimizer, have you ran this line:
#' \code{RFP <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])}
#' Normally not touched, for internal optimization purposes.
#' @param overlapGrl an integer, (default: NULL),
#' if defined must be countOverlaps(grl, RFP),
#' added for speed if you already have it
#' @inheritParams fpkm
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
disengagementScore <- function(grl, RFP, GtfOrTx, RFP.sorted = FALSE,
                               weight = 1L, overlapGrl = NULL) {
  tx <- loadRegion(GtfOrTx)
  # Optimize
  if(!RFP.sorted) {
    RFP <- optimizeReads(tx, RFP)
  }
  if (is.null(overlapGrl))
    overlapGrl <- countOverlapsW(grl, RFP, weight)
  overlapGrl <- overlapGrl + 1 # Pseudo count

  # exclude non hits and set them to 0
  validIndices <- hasHits(tx, RFP)
  txNames_grl <- txNames(grl, tx)
  validIndices <- validIndices[data.table::chmatch(txNames_grl, names(tx))]

  if (!any(validIndices)) { # if no hits
    overlapGrl
    names(score) <- NULL
    return(score)
  }
  overlapDownstream <- rep(1, length(grl))
  downstreamTx <- downstreamOfPerGroup(tx[txNames_grl][validIndices],
                                       grl[validIndices])

  # # check for big lists
  # if (length(downstreamTx) > 5e5) {
  #   ordering <- uniqueOrder(downstreamTx)
  #   downstreamTx <- uniqueGroups(downstreamTx)
  #   overlapDownstream[validIndices] <- countOverlapsW(downstreamTx,
  #                                                     RFP, weight)[ordering] + 1
  # } else {
  #   overlapDownstream[validIndices] <- countOverlapsW(downstreamTx, RFP,
  #                                                     weight) + 1
  # }
  overlapDownstream[validIndices] <- countOverlapsW(downstreamTx, RFP,
                                                    weight) + 1
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
#' @inheritParams disengagementScore
#' @param ds numeric vector (NULL), disengagement score. If you have already
#'  calculated \code{\link{disengagementScore}}, input here to save time.
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
                             RFP.sorted = FALSE,
                             weight = 1L, overlapGrl = NULL) {
  tx <- loadRegion(GtfOrTx)
  # Optimize
  if(!RFP.sorted) RFP <- optimizeReads(tx, RFP)
  if (is.null(overlapGrl))
    overlapGrl <- countOverlapsW(grl, RFP, weight)
  overlapGrl <- overlapGrl + 1 # Pseudo count

  # find tx with hits
  validIndices <- hasHits(tx, RFP)
  txNames_grl <- txNames(grl, tx)
  validIndices <- validIndices[data.table::chmatch(txNames_grl, names(tx))]
  if (!any(validIndices)) { # if no hits
    names(overlapGrl) <- NULL
    return(overlapGrl)
  }
  tx <- tx[txNames_grl][validIndices]
  grl <- grl[validIndices]

  grlStarts <- startSites(grl, asGR = FALSE, is.sorted = TRUE,
                          keep.names = FALSE)
  upstreamTx <- upstreamOfPerGroup(tx, grlStarts, allowOutside = FALSE)
  overlapTxOutside <- rep(1, length(validIndices))
  if (!is.null(ds)) { # save time here if ds is defined
    downstreamCounts <- 1 / (ds / overlapGrl)
    upstreamCounts <- rep(1, length(validIndices))
    upstreamCounts[validIndices] <- countOverlapsW(upstreamTx, RFP,
                                                   weight)
    overlapTxOutside <- downstreamCounts + upstreamCounts

  } else { # else make ds again
    downstreamTx <- downstreamOfPerGroup(tx, grl)

    overlapTxOutside[validIndices] <- countOverlaps(upstreamTx, RFP, weight) +
      countOverlapsW(downstreamTx, RFP, weight) + 1
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
#' @inheritParams insideOutsideORF
#' @inheritParams translationalEff
#' @param GtfOrThreeUtrs if Gtf: a TxDb object of a gtf file transcripts is
#'  called from: `threeUTRsByTranscript(Gtf, use.names = TRUE)`,
#'  if object is GRangesList, it is presumed to be the 3' utrs
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
                                 weight.RFP = 1L, weight.RNA = 1L,
                                 overlapGrl = NULL) {
  threeUTRs <- loadRegion(GtfOrThreeUtrs, part = "trailer")
  # Optimize
  if (is.null(overlapGrl))
    overlapGrl <- countOverlapsW(grl, RFP, weight.RFP)
  overlapGrl <- overlapGrl + 1 # Pseudo count

  # check that naming is correct, else change it.
  orfNames <- txNames(grl, threeUTRs)
  validNamesThree <- names(threeUTRs) %in% orfNames
  validNamesGRL <- orfNames %in% names(threeUTRs)
  rrs <- rep(NA,length(grl))
  if (sum(validNamesGRL) != length(grl)) {
    threeUTRs <- threeUTRs[validNamesThree]
    grl <- grl[validNamesGRL]
  }
  overlapGrl <- overlapGrl[validNamesGRL]
  threeUTRs <- threeUTRs[orfNames[validNamesGRL]]
  overlapThreeUtrs <- countOverlapsW(threeUTRs, RFP, weight.RFP) + 1

  rrs[validNamesGRL] <- (overlapGrl / widthPerGroup(grl)) /
    (overlapThreeUtrs / widthPerGroup(threeUTRs))

  if (!is.null(RNA)) { # normalize by rna ratio
    rnaRatio <- (countOverlapsW(grl, RNA, weight.RNA) + 1) /
      (countOverlapsW(threeUTRs, RNA, weight.RNA) + 1)
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
#' @inheritParams disengagementScore
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
ribosomeStallingScore <- function(grl, RFP, weight = 1L, overlapGrl = NULL) {
  # Optimize
  if (is.null(overlapGrl))
    overlapGrl <- countOverlapsW(grl, RFP, weight)
  overlapGrl <- overlapGrl + 1 # Pseudo count

  grl_len <- widthPerGroup(grl, FALSE)
  stopCodons <- stopCodons(grl, is.sorted = TRUE)
  overlapStop <- countOverlapsW(stopCodons, RFP, weight) + 1

  rss <- (overlapStop / 3) / (overlapGrl / grl_len)
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
#' @param RFP ribo seq reads as GAlignments, GRanges or GRangesList object
#' @inheritParams startRegion
#' @inheritParams getWeights
#' @family features
#' @return a numeric vector of counts
#' @export
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(21, 40),
#'                strand = "+")
#' names(ORF) <- c("tx1")
#' grl <- GRangesList(tx1 = ORF)
#' tx <- extendLeaders(grl, 20)
#' # 1 width p-shifted reads
#' reads <- GRanges("1", IRanges(c(21, 23, 50, 50, 50, 53, 53, 56, 59),
#'                             width = 1), "+")
#' score(reads) <- 28 # original width
#' startRegionCoverage(grl, reads, tx)
startRegionCoverage <- function(grl, RFP, tx = NULL, is.sorted = TRUE,
                                upstream = 2L, downstream = 2L, weight = 1L) {
  region <- startRegion(grl, tx, is.sorted, upstream, downstream)
  return(countOverlapsW(region, RFP, weight))
}

#' Get initiation score for a GRangesList of ORFs
#'
#' initiationScore tries to check how much each TIS region resembles, the
#' average of the CDS TIS regions.
#'
#' Since this features uses a distance matrix for scoring, values are
#' distributed like this:\cr
#' As result there is one value per ORF:\cr
#' 0.000: means that ORF had no reads\cr
#' -1.000: means that ORF is identical to average of CDS\cr
#' 1.000: means that orf is maximum different than average of CDS\cr
#'
#' If a score column is defined, it will use it as weights,
#' see \code{\link{getWeights}}
#' @references doi: 10.1186/s12915-017-0416-0
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @param reads ribo seq reads as \code{\link{GAlignments}},
#' GRanges or GRangesList object
#' @param cds a \code{\link{GRangesList}} object with coding sequences
#' @param tx a GRangesList of transcripts covering grl.
#' @param pShifted a logical (TRUE), are riboseq reads p-shifted?
#' @inheritParams getWeights
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
initiationScore <- function(grl, cds, tx, reads, pShifted = TRUE,
                            weight = "score") {
  if (length(grl) == 0) stop("grl must have length > 0")
  # meta coverage of cds
  cdsMeta <- windowPerReadLength(cds, tx, reads, pShifted = pShifted,
                                 weight = weight)

  # coverage per ORF
  prop <- windowPerReadLength(grl, tx, reads, pShifted = pShifted,
                              scoring = "fracPos", weight = weight)

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
#' an ORF has more reads than second and third frame. IMPORTANT: Only use
#' p-shifted libraries, see (\code{\link{detectRibosomeShifts}}).
#' Else this score makes no sense.
#'
#' Pseudocode:
#' assume rff - is reads fraction in specific frame
#' \preformatted{ORFScore = log(rff1 + rff2 + rff3)}
#' If rff2 or rff3 is bigger than rff1,
#' negate the resulting value.
#' \preformatted{ORFScore[rff1Smaller] <- ORFScore[rff1Smaller] * -1}
#'
#' As result there is one value per ORF:
#' - Positive values say that the first frame have the most reads,
#' - zero values means it is uniform:
#' (ORFscore between -2.5 and 2.5 can be considered close to uniform),
#' - negative values say that the first frame does not have the most reads.
#' NOTE non-pshifted reads:
#' If reads are not of width 1, then a read from 1-4 on range of 1-4,
#' will get scores frame1 = 2, frame2 = 1, frame3 = 1. What could be logical
#' is that only the 5' end is important, so that only frame1 = 1,
#' to get this, you first resize reads to 5'end only.
#'
#' General NOTES:
#' 1. p shifting is not exact, so some functional ORFs will get a
#' bad ORF score. \cr
#' 2. If a score column is defined, it will use it as weights, set
#' to weight = 1L if you don't have weight, and score column is
#' something else.
#' 3. If needed a test for significance and critical values,
#' use chi-squared. There are 3 degrees of freedom (3 frames),
#' so critical 0.05 (3-1 degrees of freedm = 2), value is: log2(6) = 2.58
#' see \code{\link{getWeights}}
#' @references doi: 10.1002/embj.201488411
#' @inheritParams coveragePerTiling
#' @inheritParams floss
#' @inheritParams entropy
#' @param coverage a data.table from coveragePerTiling of length same as 'grl' argument.
#' Save time if you have already computed it.
#' @param stop3 logical, default TRUE. Stop if any input is of width < 3.
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
orfScore <- function(grl, RFP, is.sorted = FALSE, weight = "score",
                     overlapGrl = NULL, coverage = NULL, stop3 = TRUE) {
  widths <- widthPerGroup(grl, FALSE)
  valid <- widths > 2
  if (stop3 & any(!valid))
    stop("ORFs with width < 3  not allowed")

  coverage_was_given <- !is.null(coverage)
  counts <- if (coverage_was_given) {
    stopifnot(is(coverage, "data.table"))
    hasHits <- valid
    stopifnot(length(unique((coverage$genes))) == sum(widths > 0))
    coverage[genes %in% which(hasHits),]
  } else {
    hasHits <- valid & hasHits(grl, RFP, overlaps = overlapGrl)
    grl <- grl[hasHits]
    coveragePerTiling(grl, RFP, is.sorted, as.data.table = TRUE,
                      withFrames = TRUE, weight = weight)
  }

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
  revORFscore <-  which(countsTile1 < countsTile2 | countsTile1 < countsTile3)
  ORFscore[revORFscore] <- -1 * ORFscore[revORFscore]
  ORFscore[is.na(ORFscore)] <- 0
  dfORFs$ORFScores <- ORFscore

  # insert back empty ones
  temp <- data.table(frame_zero_RP = rep.int(0, length(hasHits)),
                     frame_one_RP = rep.int(0, length(hasHits)),
                     frame_two_RP = rep.int(0, length(hasHits)),
                     ORFScores = rep(0, length(hasHits)))
  temp[hasHits, ] <- dfORFs
  temp[] # for print
  return(temp)
}
