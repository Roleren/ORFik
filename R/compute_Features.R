
#' Get all possible features in ORFik
#'
#' If you want to get all the features easily, you can use this function.
#' Each feature have a link to an article describing its creation and idea
#' behind it. Look at the functions in the feature family to see all of them.
#'
#' If you used CageSeq to reannotate your leaders, your txDB object must
#' contain the reassigned leaders. Use [reassignTxDbByCage()] to get the txdb.
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
#' # Here we make an example from scratch
#' # Usually the ORFs are found in orfik, which makes names for you etc.
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
        scores$kozak <- kozakSequenceScore(grl, tx, faFile)
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
  scores[] # for print
  return(scores)
}

#' Get all possible features in ORFik
#'
#' If you have a txdb with correctly reassigned transcripts, use:
#' [computeFeatures()]
#'
#' A specialized version if you don't have a correct txdb, for example with
#' CAGE reassigned leaders while txdb is not updated.
#' It is 2x faster for tested data.
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
  if (length(RFP) > 1e6) { # speedup on big riboseq libraries
    rfp <- sort(RFP[countOverlaps(RFP, tx, type = "within") > 0])
  } else {
    rfp <- RFP
  }

  #### Get all features, append 1 at a time, to save memory ####
  scores <- data.table(floss = floss(grl, RFP, cds, riboStart, riboStop))
  scores[, entropyRFP := entropy(grl, rfp)]
  scores[, disengagementScores := disengagementScore(grl, rfp, tx, TRUE)]
  scores[, RRS := ribosomeReleaseScore(grl, rfp, threeUTRs, RNA)]
  scores[, RSS := ribosomeStallingScore(grl, rfp)]

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
    scores[, ORFScores := orfScore(grl, rfp, grl.is.sorted)$ORFScores]
    scores[, ioScore := insideOutsideORF(grl, rfp, tx,
                                         scores$disengagementScores, TRUE)]

    if (includeNonVarying) {

      if (is(faFile, "FaFile") || is(faFile, "BSgenome")) {
        scores$kozak <- kozakSequenceScore(grl, tx, faFile)
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
  scores[] # for print
  return(scores)
}
