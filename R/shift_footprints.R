
#' Shifts footprints
#'
#' Function shifts footprints from bam file and loads them to memory as GRanges.
#' Resizes them to single base in 5' end fashion,
#' treated as p site. Takes account for junctions in cigars. Length of footprint is saved in
#' 'size' parameter of GRanges output. Footprints are also sorted according to their genomic
#' position, ready for saving e.g. as bed file.
#' @param footprints a bam or bed file path, or a GAlignment object
#'  of ribo-seq reads
#' @param selected_lengths Numeric vector of lengths of footprints you select for shifting.
#' @param selected_shifts Numeric vector of shifts for coresponding selected_lengths.
#' eg. c(10, -10) with selected_lengths of c(31, 32) means length of 31 will be shifted left by 10.
#' Footprints of length 32 will be shifted right by 10.
#' @return A GRanges object of shifted footprints.
#' @export
#' @import GenomicAlignments
#' @import GenomeInfoDb
#' @examples
#' gtf <- system.file("extdata", "example.gtf",
#'                    package = "ORFik") ## location of the gtf file
#' suppressWarnings(txdb <-
#'                    GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf"))
#' bam <- system.file("extdata", "example.bam",
#'                    package = "ORFik") ## location of the bam file
#' footprints <-readFootprints(bam)
#' shifts <- detectRibosomeShifts(footprints, txdb, start=TRUE,
#'                                stop=FALSE, offset_plots=FALSE, top_tx=10)
#'
#' shiftedReads <- shiftFootprints(footprints, shifts$fragment_length,
#'                                 shifts$offsets_start)
shiftFootprints <- function(footprints, selected_lengths, selected_shifts) {

    selected_shifts <- -1 * selected_shifts
    allFootrpintsShifted <- GRanges()

    riboReads <- readFootprints(footprints)

    # check for input correctness
    if (length(selected_lengths) != length(selected_shifts)) {
        stop("Incorrect input. Not equal number of elements in
             selected_lengths and selected_shifts!")
    }
    if (sum(selected_lengths > abs(selected_shifts)) !=
        length(selected_shifts)) {
        stop("Incorrect input. selected_shifts cant be bigger
             than selected_lengths!")
    }

    for (i in 1:length(selected_lengths)) {
        message("Shifting footprints of length ", selected_lengths[i])
        riboReadsW <- riboReads[qwidth(riboReads) == selected_lengths[i]]
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
            shift_by <- mapply(parse_cigar, cigar_strings, shift_by, is_plus)
            shift_by[!is_plus] <- -1 * shift_by[!is_plus]
            shifted_cigars <- shift(cigars, shift_by)
            allFootrpintsShifted <- c(allFootrpintsShifted, shifted_cigars)
        }
    }

    message("Sorting shifted footprints...")
    allFootrpintsShifted <- sortSeqlevels(allFootrpintsShifted)
    allFootrpintsShifted <- sort(allFootrpintsShifted)
    message("Process completed successfuly.")
    return(allFootrpintsShifted)
}

#' Shift ribosome footprints
#'
#' Given a TxDb object and a bam file path, get the predicted shifts
#' of the footprints to the p-site.
#' @param footprints a bam or bed file path, or a GAlignment object
#'  of ribo-seq reads
#' @param txdb a txdb object from a gtf file
#' @param start logical (T), include starts
#' @param stop logical (F), include stops
#' @param offset_plots logical (F), plot shifting for validation
#' @param top_tx numeric(10), which top percentage of transcripts to use.
#' @return a data.frame with shift lengths and offsets
#' @export
#' @examples
#' gtf <- system.file("extdata", "example.gtf",
#'         package = "ORFik") ## location of the gtf file
#' suppressWarnings(txdb <-
#'  GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf"))
#' bam <- system.file("extdata", "example.bam",
#'         package = "ORFik") ## location of the bam file
#' detectRibosomeShifts(bam, txdb, start=TRUE, stop=FALSE,
#'                      offset_plots=FALSE, top_tx=10)
#'
detectRibosomeShifts <- function(footprints, txdb, start = TRUE, stop = FALSE,
                                 offset_plots = FALSE, top_tx = 10) {
  alignments <- readFootprints(footprints)
  txNames <- get_longest_tx(txdb)
  cds150 <- get_cds150(txdb, txNames)
  ## start stop windows
  ss <- get_start_stop(txdb, txNames, start=start, stop=stop)

  ## footprint lengths
  if(class(alignments) == "GAlignments"){
    all_lengths <- sort(unique(qwidth(alignments)))
  } else {
    all_lengths <- sort(unique(width(alignments)))
  }

  selected_lengths <- c()
  offsets_start <- c()
  offsets_stop <- c()

  if (start) {
    start_df <- data.frame(start_codon=c(), count=c(), flength=c())
  }
  if (stop) {
    stop_df <- data.frame(stop_codon=c(), count=c(), flength=c())
  }

  unlCds150 <- unlist(cds150, use.names = FALSE)
  ## for each ribo-seq width, i.g. 28, 29, 30, 31 etc
  for(riboLength in all_lengths) {
    ### 5'ends
    ends_uniq <- get_5ends(alignments, riboLength)
    periodic <- periodicity(unlCds150, ends_uniq, top_tx)
    if (periodic) {
      selected_lengths <- c(selected_lengths, riboLength)
      if (start) {
        start_meta <- start_ribo(ss[[1]], ends_uniq, top_tx)
        df <- data.frame(codon=c(-30:30)[-31], count=start_meta,
                         flength=rep(riboLength, 60))
        start_df <- rbind(start_df, df)
        ### change point analysis
        offset <- change_point_analysis(df, feature="start")
        offsets_start <- c(offsets_start, offset)
      }
      if (stop) {
        stop_meta <- stop_ribo(ss[[2]], ends_uniq, top_tx)
        df <- data.frame(codon=c(-27:33)[-28], count=stop_meta,
                         flength=rep(riboLength, 60))
        stop_df <- rbind(stop_df, df)
        ### change point analysis
        offset <- change_point_analysis(df, feature="stop")
        offsets_stop <- c(offsets_stop, offset)
      }
    }
  }

  ### print plots if TRUE
  if (offset_plots) {
    if (start) {
      plot_periodic_lengths_start(start_df)
    }
    if (stop) {
      plot_periodic_lengths_stop(stop_df)
    }
  }

  ## return lengths and offsets
  shifts <- data.frame(fragment_length = selected_lengths)
  if (start) {
    shifts$offsets_start <-  offsets_start
  }
  if (stop) {
    shifts$offsets_stop <-  offsets_stop
  }
  return(shifts)
}
