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
#' @param footprints (GAlignments) object of RiboSeq reads
#' @param selected_lengths Numeric vector of lengths of footprints you select
#' for shifting.
#' @param selected_shifts Numeric vector of shifts for coresponding
#' selected_lengths. eg. c(10, -10) with selected_lengths of c(31, 32) means
#' length of 31 will be shifted left by 10. Footprints of length 32 will be
#' shifted right by 10.
#' @return A GRanges object of shifted footprints, sorted and resized to 1bp of
#' p-site, with metacolumn "size" indicating footprint size before shifting and
#' resizing.
#' @export
#' @examples
#' \dontrun{
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' footprints <- GenomicAlignments::readGAlignments(
#'   riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
#'     isDuplicate = FALSE, isSecondaryAlignment = FALSE)))
#' # detect the shifts automagically
#' shifts <- detectRibosomeShifts(footprints, txdb)
#' # shift the RiboSeq footprints
#' shiftedReads <- shiftFootprints(footprints, shifts$fragment_length,
#'                                 shifts$offsets_start)
#' }
shiftFootprints <- function(footprints, selected_lengths, selected_shifts) {

    selected_shifts <- -1 * selected_shifts
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
        riboReadsW <- footprints[qwidth(footprints) == selected_lengths[i]]
        if (length(riboReadsW) == 0) {
            next
        }
        is_cigar <- width(riboReadsW) != qwidth(riboReadsW)
        cigar_strings <- cigar(riboReadsW[is_cigar])
        sizes <- qwidth(riboReadsW)

        riboReadsW <- granges(riboReadsW, use.mcols = TRUE)
        riboReadsW$size <- sizes  #move footprint length to size
        riboReadsW <- resize(riboReadsW, 1)  #resize to 5' only

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
#' @param footprints (GAlignments) object of RiboSeq reads - footprints
#' @param txdb a txdb object from a gtf file
#' @param start (logical) Whether to include predictions based on the start
#' codons. Default TRUE.
#' @param stop (logical) Whether to include predictions based on the stop
#' codons. Default FASLE.
#' @param top_tx (integer) Specify which % of the top covered by RiboSeq reads
#' transcripts to use for estimation of the shifts. By default we take top 10%
#' top covered transcripts as they represent less noisy dataset. This is only
#' applicable when there are more than 1000 transcripts.
#' @inheritParams txNamesWithLeaders
#' @param firstN (integer) Represents how many bases of the transcripts
#' downstream of start codons to use for initial estimation of the
#' periodicity.
#' @return a data.frame with lengths of footprints and their predicted
#' coresponding offsets
#' @importFrom BiocGenerics Reduce
#' @export
#' @examples
#' \dontrun{
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' footprints <- GenomicAlignments::readGAlignments(
#'   riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
#'     isDuplicate = FALSE, isSecondaryAlignment = FALSE)))
#'
#' detectRibosomeShifts(footprints, txdb, stop = TRUE)
#' }
#'
detectRibosomeShifts <- function(
  footprints, txdb, start = TRUE, stop = FALSE,
  top_tx = 10L, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = 30L,
  firstN = 150L) {
  window_size = 30L


  txNames <- txNamesWithLeaders(txdb, minFiveUTR = minFiveUTR,
                                minCDS = minCDS, minThreeUTR = minThreeUTR)
  if (length(txNames) == 0) stop("No transcript has leaders and trailers of",
                                 " specified minFiveUTR, minCDS, minThreeUTR")
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)[txNames]

  # reduce data set to only matching seqlevels
  cds <- keepSeqlevels(cds, unique(seqnames(footprints))[unique(seqnames(
                       footprints)) %in% unique(seqnamesPerGroup(cds, FALSE))],
                        pruning.mode = "coarse")

  txNames <- txNames[txNames %in% names(cds)]
  footprints <- keepSeqlevels(footprints, unique(seqnamesPerGroup(cds, FALSE)),
                              pruning.mode = "coarse")
  if (length(cds) == 0 | length(footprints) == 0) {
    stop("txdb and footprints did not have any matched seqnames")
  }
  ## start stop windows
  ss <- getStartStopWindows(txdb, txNames, start = start, stop = stop,
                            window_size = window_size, cds)
  cds <- reduceKeepAttr(downstreamN(cds, firstN = firstN))
  rWidth <- readWidths(footprints)
  all_lengths <- sort(unique(rWidth))
  selected_lengths <- c()
  offsets_start <- c()
  offsets_stop <- c()
  footprints <- resize(granges(footprints), 1L)
  for (l in all_lengths) {
    ends_uniq <- footprints[rWidth == l]

    # get top_tx of covered tx
    counts <- countOverlaps(cds, ends_uniq)
    counts <- counts[counts > 1]
    if (length(counts) > 1000) {
      counts <- sort(counts, decreasing = TRUE)
      counts <- counts[seq_len(floor(top_tx * length(counts) / 100))]
    }
    if (length(counts) == 0) next
    # This is the slow line, we need to speed this up! ->
    cvgCDS <- coverageByWindow(ends_uniq, cds[names(counts)], is.sorted = TRUE,
                               keep.names = FALSE)
    cvgCDS <- Reduce(`+`, cvgCDS)
    if (isPeriodic(as.vector(cvgCDS))) {
      selected_lengths <- c(selected_lengths, l)
      if (start) {
        start_meta <- metaWindow(ends_uniq, ss$starts)
        offset <- changePointAnalysis(start_meta$avg_counts, feature = "start")
        offsets_start <- c(offsets_start, offset)
      }
      if (stop) {
        stop_meta <- metaWindow(ends_uniq, ss$stops)
        offset <- changePointAnalysis(stop_meta$avg_counts, feature = "stop")
        offsets_stop <- c(offsets_stop, offset)
      }
    }
  }

  shifts <- data.frame(fragment_length = selected_lengths)
  if (start) shifts$offsets_start <- offsets_start
  if (stop) shifts$offsets_stop <- offsets_stop
  return(shifts)
}
