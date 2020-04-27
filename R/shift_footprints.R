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
#' @param footprints \code{\link{GAlignments}} object of RiboSeq reads
#' @param shifts a data.frame / data.table with minimum 2 columns,
#' fraction (selected_lengths) and selected_shifts (relative position).
#' Output from \code{\link{detectRibosomeShifts}}
#' @return A \code{\link{GRanges}} object of shifted footprints, sorted and
#' resized to 1bp of p-site,
#' with metacolumn "size" indicating footprint size before shifting and
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
  if (nrow(shifts) == 0) stop("No shifts in data.frame")

  selected_lengths <- shifts$fraction
  selected_shifts <- shifts$offsets_start

  lengthsAll <- readWidths(footprints, along.reference = TRUE)
  validLengths <- lengthsAll %in% selected_lengths
  lengthsAll <- lengthsAll[validLengths]
  footprints <- footprints[validLengths]
  cigarsAll <- cigar(footprints)

  shiftsAll <- -selected_shifts[chmatch(as.character(lengthsAll),
                          as.character(selected_lengths))]

  # Shift cigar and map back to corrected GRanges.
  is_pos <- strandBool(footprints)
  starts <- rep(NA, length(cigarsAll))
  ends <- starts
  starts[is_pos] <- shiftsAll[is_pos] + 1
  ends[!is_pos] <- lengthsAll[!is_pos] - shiftsAll[!is_pos]

  shiftedCigar <- cigarNarrow(cigarsAll, starts, ends)

  pos <- start(footprints) + unlist(attributes(shiftedCigar),
                                    use.names = FALSE)

  shifted <- GAlignments(seqnames(footprints), pos = pos,
                         cigar = shiftedCigar, strand = strand(footprints))
  shifted <- convertToOneBasedRanges(shifted)
  shifted$size <- lengthsAll

  message("Sorting shifted footprints...")
  shifted <- sortSeqlevels(shifted)
  shifted <- sort(shifted)
  return(shifted)
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
#' 33, but will remove 3S so final read width is 30 in ORFik. This is standard
#' for ribo-seq.
#' @param footprints (GAlignments) object of RiboSeq reads - footprints, can
#' also be path to the file.
#' @inheritParams loadTxdb
#' @param start (logical) Whether to include predictions based on the start
#' codons. Default TRUE.
#' @param stop (logical) Whether to include predictions based on the stop
#' codons. Default FASLE. Only use if there exists 3' UTRs for the annotation.
#' If peridicity around stop codon is stronger than at the start codon, use
#' stop instead of start region for p-shifting.
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
#' @param min_reads default (1000), how many reads must a read-length have to
#' be considered for periodicity.
#' @param accepted.lengths accepted readlengths, default 1:1000, usually ribo-seq
#' is between 26:34.
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
  firstN = 150L, tx = NULL, min_reads = 1000, accepted.lengths = 1:1000) {
  txdb <- loadTxdb(txdb)
  # Filters for cds and footprints
  txNames <- filterTranscripts(txdb, minFiveUTR = minFiveUTR, minCDS = minCDS,
                               minThreeUTR = minThreeUTR)
  cds <- cdsBy(txdb, by = "tx", use.names = TRUE)[txNames]
  footprints <- fimport(footprints, cds)

  # reduce data-set to only matching seqlevels
  seqMatch <- validSeqlevels(cds, footprints)
  cds <- keepSeqlevels(cds, seqMatch, pruning.mode = "coarse")
  footprints <- keepSeqlevels(footprints, seqMatch, pruning.mode = "coarse")
  if (length(cds) == 0 | length(footprints) == 0) {
    stop("txdb and footprints did not have any matched seqnames")
  }
  txNames <- txNames[txNames %in% names(cds)]
  if (is.null(tx)) tx <- loadRegion(txdb)
  tx <- tx[txNames]

  # find periodic read lengths
  footprints <- convertToOneBasedRanges(footprints, addSizeColumn = TRUE,
                                        addScoreColumn = TRUE,
                                        along.reference = TRUE)
  # Filter if < 1000 counts read size or not in accepted.lengths
  lengths <- data.table(score = footprints$score, size = footprints$size)
  tab <- lengths[, .(counts = sum(score)), by = size]
  tab <- tab[(counts >= min_reads) & (size %in% accepted.lengths),]
  if (nrow(tab) == 0) stop("No valid read lengths found with",
                            "accepted.lengths and counts > min_reads")

  periodicity <- windowPerReadLength(grl = cds, tx = tx, reads = footprints,
                      pShifted = FALSE, upstream = 0, downstream = 149,
                      zeroPosition = 0, scoring = "periodic",
                      acceptedLengths = tab$size)
  validLengths <- periodicity[score == TRUE,]$fraction

  # find shifts
  offset <- data.table()
  if (start) {
    rw <- windowPerReadLength(grl = cds, tx = tx, reads = footprints,
                              pShifted = FALSE, upstream = 30, downstream = 29,
                              acceptedLengths = validLengths)
    offset <- rw[, .(offsets_start = changePointAnalysis(score)),
                 by = fraction]
    # Figure how this is possible ->
    offset <- offset[offsets_start < 0, ]
  }
  if (stop & !is.null(minThreeUTR)) {
    threeUTRs <- loadRegion(txdb, "trailers")[txNames]
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

#' Shift footprints of each file in experiment
#'
#' Saves files to a specified location as .bed, it will include a score column
#' containing read width.
#'
#' For more details, see: \code{\link{detectRibosomeShifts}}
#' @param df an ORFik \code{\link{experiment}}
#' @param out.dir output directory for files,
#' default: dirname(df$filepath[1]), making a /pshifted
#' folder at that location
#' @inheritParams detectRibosomeShifts
#' @param output_format default (bed), use export.bed or ORFik optimized
#' (bedo) using \code{\link{export.bedo}} ?
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return NULL (Objects are saved to out.dir/pshited/"name_pshifted.bed"
#' or .bedo)
#' @importFrom rtracklayer export.bed
#' @family pshifting
shiftFootprintsByExperiment <- function(df,
                                        out.dir = pasteDir(dirname(
                                          df$filepath[1]), "/pshifted/"),
                                        start = TRUE, stop = FALSE,
                                        top_tx = 10L, minFiveUTR = 30L,
                                        minCDS = 150L, minThreeUTR = 30L,
                                        firstN = 150L, min_reads = 1000,
                                        accepted.lengths = 1:1000,
                                        output_format = "bed",
                                        BPPARAM = bpparam()) {
  path <- out.dir
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(path)) stop(paste("out.dir", out.dir, "does not exist!"))
  message(paste("Saving", output_format, "files to:", out.dir))
  message(paste0("Shifting reads in experiment:", df@experiment))

  txdb <- loadTxdb(df)
  rfpFiles <- filepath(df, "default")
  bplapply(rfpFiles, FUN = function(file, path, txdb, start, stop,
                                    top_tx, minFiveUTR, minCDS, minThreeUTR,
                                    firstN, min_reads, accepted.lengths,
                                    output_format) {
    message(file)
    rfp <- fimport(file)
    shifts <- detectRibosomeShifts(rfp, txdb, start = start, stop = stop,
                                   top_tx = top_tx, minFiveUTR = minFiveUTR,
                                   minCDS = minCDS, minThreeUTR = minThreeUTR,
                                   firstN = firstN, min_reads = min_reads,
                                   accepted.lengths = accepted.lengths)
    shifted <- shiftFootprints(rfp, shifts)
    name <- paste0(path, remove.file_ext(df$filepath[i], basename = TRUE))
    if (output_format == "bed") {
      shifted$score <- shifted$size
      export.bed(shifted, paste0(name, "_pshifted.bed"))
    } else if (output_format == "bedo") {
      shifted <- convertToOneBasedRanges(shifted, addScoreColumn = TRUE,
                                         addSizeColumn = TRUE)
      export.bedo(shifted, paste0(name, "_pshifted.bedo"))
    } else stop("output_format must be bed or bedo")
    return(invisible(NULL))
  }, path = path, txdb = txdb, start = start, stop = stop,
      top_tx = top_tx, minFiveUTR = minFiveUTR,
      minCDS = minCDS, minThreeUTR = minThreeUTR,
      firstN = firstN, min_reads = min_reads,
      accepted.lengths = accepted.lengths, output_format = output_format,
      BPPARAM = BPPARAM)
  return(invisible(NULL))
}
