
#' Get all possible features in ORFik
#'
#' If you want to get all the NGS and/or sequence features easily,
#' you can use this function.
#' Each feature have a link to an article describing its creation and idea
#' behind it. Look at the functions in the feature family to see all of them.
#' Example, if you want to know what the "te" column is, check out:
#' ?translationalEff.
#' \cr\cr
#' If you used CageSeq to reannotate your leaders, your txDB object must
#' contain the reassigned leaders. Use [reassignTxDbByCage()] to get the txdb.
#'
#' As a note the library is reduced to only reads overlapping 'tx', so the
#' library size in fpkm calculation is done on this subset. This will help
#' remove rRNA and other contaminants.\cr
#' Also if you have only unique reads with a weight column, explaining the
#' number of duplicated reads, set weights to make calculations correct.
#' See \code{\link{getWeights}}
#' @param grl a \code{\link{GRangesList}} object
#'  with usually ORFs, but can also be either leaders, cds', 3' utrs, etc.
#'  This is the regions you want to score.
#' @param RFP RiboSeq reads as \code{\link{GAlignments}} , \code{\link{GRanges}}
#'  or \code{\link{GRangesList}} object
#' @param RNA RnaSeq reads as \code{\link{GAlignments}} , \code{\link{GRanges}}
#'  or \code{\link{GRangesList}} object
#' @param Gtf a TxDb object of a gtf file or path to gtf, gff .sqlite etc.
#' @param faFile a path to fasta indexed genome, an open \code{\link{FaFile}},
#' a BSgenome, or path to ORFik \code{\link{experiment}} with valid genome.
#' @param riboStart usually 26, the start of the floss interval, see ?floss
#' @param riboStop usually 34, the end of the floss interval
#' @param sequenceFeatures a logical, default TRUE, include all sequence
#' features, that is: Kozak, fractionLengths, distORFCDS, isInFrame,
#' isOverlapping and rankInTx
#' @param grl.is.sorted logical (F), a speed up if you know argument grl
#'  is sorted, set this to TRUE.
#' @inheritParams translationalEff
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
computeFeatures <- function(grl, RFP, RNA = NULL,  Gtf, faFile = NULL,
                            riboStart = 26, riboStop = 34,
                            sequenceFeatures = TRUE, grl.is.sorted = FALSE,
                            weight.RFP = 1L, weight.RNA = 1L) {
  #### Check input and load data ####
  validGRL(class(grl), "grl")
  checkRFP(class(RFP))
  checkRNA(class(RNA))
  Gtf <- loadTxdb(Gtf)

  # get transcript parts
  fiveUTRs <- loadRegion(Gtf, "leaders")
  cds <- loadRegion(Gtf, "cds")
  threeUTRs <- loadRegion(Gtf, "trailers")
  tx <- loadRegion(Gtf)

  return(allFeaturesHelper(grl, RFP, RNA, tx, fiveUTRs, cds, threeUTRs, faFile,
                           riboStart, riboStop, sequenceFeatures,
                           grl.is.sorted, weight.RFP, weight.RNA))
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
#'  If you are using CAGE, you do not need to reassign these to the cage
#'  peaks, it will do it for you.
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
#'  \donttest{
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
#'   ORFik:::computeFeaturesCage(grl = fiveUTR_ORFs, RFP = RFP,
#'    RNA = RNA, Gtf = txdb, faFile = faFile)
#'
#' }
#' # See vignettes for more examples
#' }
#'
computeFeaturesCage <- function(grl, RFP, RNA = NULL, Gtf = NULL, tx = NULL,
                                fiveUTRs = NULL, cds = NULL, threeUTRs = NULL,
                                faFile = NULL, riboStart = 26, riboStop = 34,
                                sequenceFeatures = TRUE, grl.is.sorted = FALSE,
                                weight.RFP = 1L, weight.RNA = 1L) {
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
      Gtf <- loadTxdb(Gtf)

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
        tx <- loadRegion(Gtf)
      }
    }
  return(allFeaturesHelper(grl, RFP, RNA, tx, fiveUTRs, cds, threeUTRs, faFile,
                           riboStart, riboStop, sequenceFeatures,
                           grl.is.sorted, weight.RFP, weight.RNA))
}

#' Calculate the features in computeFeatures function
#'
#' Not used directly, calculates all features internally for computeFeatures.
#' @inheritParams computeFeaturesCage
#' @param st (NULL), if defined must be: st = startRegion(grl, tx, T, -3, 9)
#' @return a data.table with features
allFeaturesHelper <- function(grl, RFP, RNA, tx, fiveUTRs, cds , threeUTRs,
                              faFile, riboStart, riboStop,
                              sequenceFeatures, grl.is.sorted,
                              weight.RFP = 1L, weight.RNA = 1L,
                              st = NULL) {
  # Clean and optimize
  if (!grl.is.sorted){
    grl <- sortPerGroup(grl)
    grl.is.sorted <- TRUE
  }
  #type <- ifelse(all(unique(width(RFP)) == 1), "within", "any")
  tx <- tx[names(tx) %in% txNames(grl, tx, unique = TRUE)]
  RFP <- optimizeReads(tx, RFP)
  tx_old <- tx
  tx <- tx[txNames(grl)] # Subset tx to only those in grl.
  weight.RFP <- getWeights(RFP, weight.RFP)


  #### Get all features, append 1 at a time, to save memory ####
  scores <- data.table(countRFP = countOverlapsW(grl, RFP, weight.RFP))
  if (!is.null(RNA)) { # if rna seq is included
    TE <- translationalEff(grl, RNA, RFP, tx, with.fpkm = TRUE,
                           weight.RFP = weight.RFP, weight.RNA = weight.RNA)
    scores[, te := TE$te]
    scores[, fpkmRFP := TE$fpkmRFP]
    scores[, fpkmRNA := TE$fpkmRNA]
  } else {
    scores[, fpkmRFP := fpkm_calc(countRFP, widthPerGroup(grl),
                                  sum(weight.RFP))]
  }
  scores[, floss := floss(grl, RFP, cds, riboStart, riboStop, weight.RFP)]
  scores[, entropyRFP := entropy(grl, RFP, weight.RFP, grl.is.sorted,
                                 overlapGrl = countRFP)]
  scores[, disengagementScores := disengagementScore(grl, RFP, tx_old, TRUE,
                                                     weight.RFP, countRFP)]
  scores[, RRS := ribosomeReleaseScore(grl, RFP, threeUTRs, RNA,
                                       weight.RFP, weight.RNA, countRFP)]
  scores[, RSS := ribosomeStallingScore(grl, RFP, weight.RFP, countRFP)]
  scores[, ORFScores := orfScore(grl, RFP, grl.is.sorted,
                                 weight.RFP, overlapGrl = countRFP)$ORFScores]
  scores[, ioScore := insideOutsideORF(grl, RFP, tx_old,
                                       scores$disengagementScores, TRUE,
                                       weight.RFP, countRFP)]
  scores[, startCodonCoverage := startRegionCoverage(grl, RFP, tx,
                                                     weight = weight.RFP)]

  if (is.null(st)) st <- startRegion(grl, tx, TRUE, -3, 9)
  st <- countOverlapsW(st, RFP, weight.RFP)  /
           (pmax(widthPerGroup(st), 1) / 5) # normalize to same size as startCodon
  scores[, startRegionCoverage := st]
  scores[, startRegionRelative := (startCodonCoverage + 1) /
           (startRegionCoverage + 1)] # Relative score

  if (sequenceFeatures) { # sequence features
    if (!is.null(faFile)) {
      faFile <- findFa(faFile)
      scores[, kozak := kozakSequenceScore(grl, tx, faFile)]
      scores[, gc := gcContent(grl, faFile)]

      # Start and stop codons
      starts <- startCodons(grl, is.sorted = TRUE)
      stops <- stopCodons(grl, is.sorted = TRUE)
      scores[, StartCodons := txSeqsFromFa(starts, faFile, TRUE, FALSE)]
      scores[, StopCodons := txSeqsFromFa(stops, faFile, TRUE, FALSE)]
    } else {
      message("Notification: faFile not included,",
      "skipping features dependent fasta genome")
    }
    # switch five with tx, is it possible to use ?
    scores[, fractionLengths := fractionLength(grl, widthPerGroup(tx, TRUE))]
    scores[, distORFCDS := distToCds(grl, fiveUTRs, cds)]
    scores[, inFrameCDS := isInFrame(distORFCDS)]
    scores[, isOverlappingCds := isOverlapping(distORFCDS)]
    scores[, rankInTx := rankOrder(grl)]
  } else {
    message("Notification: sequenceFeatures set to False,",
    "dropping all sequenceFeatures features.")
  }

  scores[] # for print
  return(scores)
}
