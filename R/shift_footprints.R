#' Shift footprints by selected offsets
#'
#' Function shifts footprints (GRanges) using specified offsets for every of
#' the specified lengths. Reads that do not conform to the specified lengths
#' are filtered out and rejected. Reads are resized to single base in 5' end
#' fashion, treated as p site.
#' This function takes account for junctions in cigars of the reads. Length of
#' the footprint is saved in size' parameter of GRanges output. Footprints are
#' also sorted according to their genomic position, ready to be saved as a
#' bed file.
#'
#' The two columns in shift are:
#' - fraction Numeric vector of lengths of footprints you select
#' for shifting.
#' - offsets_start Numeric vector of shifts for coresponding
#' selected_lengths. eg. c(10, -10) with selected_lengths of c(31, 32) means
#' length of 31 will be shifted left by 10. Footprints of length 32 will be
#' shifted right by 10.
#'
#' NOTE: It will remove softclips from valid width, the CIGAR 3S30M is qwidth
#' 33, but will remove 3S so final read width is 30 in ORFik.
#' @param footprints (GAlignments) object of RiboSeq reads
#' @param shifts a data.frame / data.table with minimum 2 columns,
#' selected_lengths and selected_shifts.
#' Output from \code{\link{detectRibosomeShifts}}
#' @return A GRanges object of shifted footprints, sorted and resized to 1bp of
#' p-site, with metacolumn "size" indicating footprint size before shifting and
#' resizing, sorted in increasing order.
#' @family pshifting
#' @export
#' @examples
#' \dontrun{
#' # input path to gtf, or load it as TxDb.
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' # load reads
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' footprints <- GenomicAlignments::readGAlignments(
#'   riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
#'     isDuplicate = FALSE, isSecondaryAlignment = FALSE)))
#' # detect the shifts automagically
#' shifts <- detectRibosomeShifts(footprints, gtf_file)
#' # shift the RiboSeq footprints
#' shiftedReads <- shiftFootprints(footprints, shifts)
#' }
shiftFootprints <- function(footprints, shifts) {
  if (!is(shifts, "data.frame")) stop("shifts must be data.frame/data.table")
  selected_lengths <- shifts$fraction
  selected_shifts <- -1 * shifts$offsets_start
  allFootrpintsShifted <- GRanges()

  if (length(selected_lengths) != length(selected_shifts)) {
      stop("Incorrect input. Not equal number of elements in",
           " selected_lengths and selected_shifts!")
  }
  if (sum(selected_lengths > abs(selected_shifts)) !=
      length(selected_shifts)) {
      stop("Incorrect input. selected_shifts cant be bigger",
           " than selected_lengths!")
  }

  for (i in seq_along(selected_lengths)) {
      message("Shifting footprints of length ", selected_lengths[i])
      riboReadsW <- footprints[readWidths(footprints) == selected_lengths[i]]
      if (length(riboReadsW) == 0) {
          next
      }
      is_cigar <- width(riboReadsW) != readWidths(riboReadsW)
      cigar_strings <- cigar(riboReadsW[is_cigar])
      sizes <- readWidths(riboReadsW)

      riboReadsW <- granges(riboReadsW, use.mcols = TRUE)
      riboReadsW$size <- sizes  #move footprint length to size
      riboReadsW <- resize(riboReadsW, 1L)  #resize to 5' only

      cigars <- riboReadsW[is_cigar]
      notCigars <- riboReadsW[!is_cigar]

      # Not Cigars - shift 5' ends, + shift right, - shift left
      if (length(notCigars) != 0) {
          is_plus <- as.vector(strand(notCigars) == "+")
          shiftedNotCigarsPlus <- shift(notCigars[is_plus],
                                        selected_shifts[i])
          shiftedNotCigarsMinus <- shift(notCigars[!is_plus],
                                         -1 * selected_shifts[i])
          allFootrpintsShifted <- c(allFootrpintsShifted,
                                    shiftedNotCigarsPlus,
                                    shiftedNotCigarsMinus)
      }
      # Cigars
      if (length(cigars) != 0) {
          is_plus <- as.vector(strand(cigars) == "+")
          shift_by <- rep(selected_shifts[i], length(cigars))
          shift_by <- mapply(parseCigar, cigar_strings, shift_by, is_plus)
          shift_by[!is_plus] <- -1 * shift_by[!is_plus]
          shifted_cigars <- shift(cigars, shift_by)
          allFootrpintsShifted <- c(allFootrpintsShifted, shifted_cigars)
      }
  }

  message("Sorting shifted footprints...")
  allFootrpintsShifted <- sortSeqlevels(allFootrpintsShifted)
  allFootrpintsShifted <- sort(allFootrpintsShifted)
  return(allFootrpintsShifted)
}

#' Detect ribosome shifts
#'
#' Utilizes periodicity measurement (fourier transform) and change point
#' analysis to detect ribosomal footprint shifts for each of the ribosomal
#' read lengths. Returns subset of read lengths and their shifts for which
#' top covered transcripts follow periodicity measure. Each shift value
#' assumes 5' anchoring of the reads, so that output offsets values will
#' shift 5' anchored footprints to be on the p-site of the ribosome.
#'
#' Check out vignette for the examples of plotting RiboSeq metaplots over start
#' and stop codons, so that you can verify visually whether this function
#' detects correct shifts.
#'
#' NOTE: It will remove softclips from valid width, the CIGAR 3S30M is qwidth
#' 33, but will remove 3S so final read width is 30 in ORFik.
#' @param footprints (GAlignments) object of RiboSeq reads - footprints, can
#' also be path to the file.
#' @inheritParams loadTxdb
#' @param start (logical) Whether to include predictions based on the start
#' codons. Default TRUE.
#' @param stop (logical) Whether to include predictions based on the stop
#' codons. Default FASLE. Only use if there exists 3' UTRs for the annotation.
#' @param top_tx (integer) Specify which % of the top covered by RiboSeq reads
#' transcripts to use for estimation of the shifts. By default we take top 10%
#' top covered transcripts as they represent less noisy dataset. This is only
#' applicable when there are more than 1000 transcripts.
#' @inheritParams filterTranscripts
#' @param firstN (integer) Represents how many bases of the transcripts
#' downstream of start codons to use for initial estimation of the
#' periodicity.
#' @param tx a GRangesList, if you do not have 5' UTRs in annotation, send
#' your own version. Example: extendLeaders(tx, 30)
#' Where 30 bases will be new "leaders". Since each original transcript was
#' either only CDS or non-coding (filtered out).
#' @return a data.table with lengths of footprints and their predicted
#' coresponding offsets
#' @family pshifting
#' @export
#' @examples
#' \dontrun{
#' # Transcriptome annotation ->
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' # The ribo seq file, usually .bam file ->
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' footprints <- GenomicAlignments::readGAlignments(
#'   riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
#'     isDuplicate = FALSE, isSecondaryAlignment = FALSE)))
#'
#' detectRibosomeShifts(footprints, gtf_file, stop = TRUE)
#'
#' # Without 5' Annotation
#' library(GenomicFeatures)
#'
#' txdb <- loadTxdb(gtf_file)
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' tx <- extendLeaders(tx, 30)
#' # Now run function, without 5' and 3' UTRs
#' detectRibosomeShifts(footprints, txdb, start = TRUE, minFiveUTR = NULL,
#'                      minCDS = 150L, minThreeUTR = NULL, firstN = 150L,
#'                      tx = tx)
#'
#' }
#'
detectRibosomeShifts <- function(footprints, txdb, start = TRUE, stop = FALSE,
  top_tx = 10L, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = 30L,
  firstN = 150L, tx = NULL) {
  txdb <- loadTxdb(txdb)
  # Filters for cds and footprints
  txNames <- filterTranscripts(txdb, minFiveUTR = minFiveUTR, minCDS = minCDS,
                               minThreeUTR = minThreeUTR)
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)[txNames]
  footprints <- fimport(footprints, cds)

  # reduce data-set to only matching seqlevels
  seqMatch <- validSeqlevels(cds, footprints)
  cds <- keepSeqlevels(cds, seqMatch, pruning.mode = "coarse")
  footprints <- keepSeqlevels(footprints, seqMatch, pruning.mode = "coarse")
  if (length(cds) == 0 | length(footprints) == 0) {
    stop("txdb and footprints did not have any matched seqnames")
  }
  txNames <- txNames[txNames %in% names(cds)]
  if (is.null(tx)) tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
  tx <- tx[txNames]

  # find periodic read lengths
  footprints <- convertToOneBasedRanges(footprints, addSizeColumn = TRUE)
  periodicity <- windowPerReadLength(grl = cds, tx = tx, reads = footprints,
                      pShifted = FALSE, upstream = 0, downstream = 149,
                      zeroPosition = 0, scoring = "periodic")
  validLengths <- periodicity[score == TRUE,]$fraction

  # find shifts
  offset <- data.table()
  if (start) {
    rw <- windowPerReadLength(grl = cds, tx = tx, reads = footprints,
                              pShifted = FALSE, upstream = 30, downstream = 29,
                              acceptedLengths = validLengths)
    offset <- rw[, .(offsets_start = changePointAnalysis(score)),
                 by = fraction]
  }
  if (stop & !is.null(minThreeUTR)) {
    threeUTRs <- threeUTRsByTranscript(txdb, use.names = TRUE)[txNames]
    rw <- windowPerReadLength(grl = threeUTRs, tx = tx, reads = footprints,
                              pShifted = FALSE, upstream = 30, downstream = 29,
                              acceptedLengths = validLengths)
    if (nrow(offset) == 0) {
      offset <- rw[, .(offsets_stop = changePointAnalysis(score, "stop")),
                   by = fraction]
    } else offset$offsets_stop <- rw[, .(changePointAnalysis(score, "stop")),
                                     by = fraction]$V1
  }

  return(offset)
}
