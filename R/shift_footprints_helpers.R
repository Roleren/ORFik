#' Find if there is a periodicity of 3 in the vector
#'
#' It uses Fourier transform + periodogram for finding periodic
#' vectors on the transcript
#' normalized counts over all CDS regions from position 0 (TIS) to 149
#' (or other max position if increased by the user.\cr
#' Checks if there is a periodicity and if the periodicity is 3,
#' more precisely between 2.9 and 3.1.\cr
#'
#' Input data:\cr
#' Transcript normalized means per CDS TIS region, count reads per position,
#' divide that number per position by the total of that transcript, then sum
#' up these numbers per position for all transcripts.\cr
#' Detection method:\cr
#' The maximum dominant Fourier frequencies is found by finding which
#' period has the highest spectrum density (using a 10% cosine taper).
#' @param x (numeric) Vector of values to detect periodicity of 3 like in
#' RiboSeq data.
#' @param info specify read length if wanted for verbose output.
#' @param verbose logical, default FALSE.
#' Report details of analysis/periodogram. Good if you are not sure
#' if the analysis was correct.
#' @param strict.fft logical, TRUE. Use a FFT without noise filter.
#' This means keep only reads lengths that are "periodic for the human eye".
#' If you want more coverage, set to FALSE, to also get read lengths
#' that are "messy", but the noise filter detects the periodicity of 3.
#' This should only be done when you do not need high quality
#' periodic reads! Example would be differential translation analysis by
#' counts over each ORF.
#' @return a logical, if it is periodic.
#' @keywords internal
#' @importFrom stats fft spec.pgram
#'
isPeriodic <- function(x, info = NULL, verbose = FALSE, strict.fft = TRUE) {
  if (sum(x) == 0) return(FALSE)
  amplitudes <- abs(fft(x))
  amp <- amplitudes[2 : (length(amplitudes) / 2 + 1)]
  specter <- spec.pgram(x = x, plot = FALSE)
  periods <- 1 / specter$freq
  if (verbose) {

    dt <- data.table(periods,
                     amp, spec = specter$spec)[order(spec, decreasing = TRUE)][1:10,]
    if (!is.null(info)) message("Read length: ", info)
    message("Top 10 periods from spectrogram, look for period > 2.9 & < 3.1:")
    print(dt)
  }
  if (strict.fft) {
    return((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1))
  } else {
    return((periods[which.max(specter$spec)] > 2.9) &
             (periods[which.max(specter$spec)] < 3.1))
  }
}

#' Get the offset for specific RiboSeq read width
#'
#' Creates sliding windows of transcript normalized counts per position
#' and check which window has most in upstream window vs downstream window.
#' Pick the position with highest absolute value maximum of the window difference.
#' Checks windows with split sites between positions -17 to -7, where 0 is TIS.
#' Normally you expect the shift around -12 for Ribo-seq, in TCP-seq / RCP-seq
#' it is usually a bit higher, usually because of cross-linking variations.
#'
#' For visual explanation, see the supl. data of ORFik paper:
#' Transcript normalized means per CDS TIS region, count reads per position,
#' divide that number per position by the total of that transcript, then sum
#' up these numbers per position for all transcripts.
#' @param x a vector with count per position to analyse, assumes the zero
#' position (TIS) is in the middle + 1 (position 0). Default it is size 60,
#' from -30 to 29 in p-shifting
#' @param feature (character) either "start" or "stop"
#' @param max.pos integer, default 40L, subset x to go from index 1 to max.pos,
#'  if tail is not relevant.
#' @param interval integer vector , default seq.int(14L, 24L).
#'  The possible shift locations, default
#'  Seperation points for upstream and downstream windows.
#'  That is (+/- 5 from -12) position.
#' @param center.pos integer, default 12. Centering position for likely p-site.
#' A first qualified guess to save time. 12 means 12 bases before TIS.
#' @param info specify read length if wanted for verbose output.
#' @param verbose logical, default FALSE. Report details of change point analysis.
#' @return a single numeric offset, -12 would mean p-site is 12 bases upstream
#' @keywords internal
#' @family pshifting
#'
changePointAnalysis <- function(x, feature = "start", max.pos = 40L,
                                interval = seq.int(14L, 24L),
                                center.pos = 12, info = NULL, verbose = FALSE) {
  if (length(x) == 0) stop("Length of x is 0")
  if (max.pos > length(x)) stop("Can not subset max.pos > length of x")
  if (!all(interval %in% seq.int(x)))
    stop("interval vector must be subset of x indices!")
  meta <- x[seq.int(max.pos)] # First 40 or other specified
  pos <- -(length(x)/2):(length(x)/2 - 1) # The positions
  frames <- rep(c(0,1,2), 20)[seq.int(max.pos)] # Frame vector
  # Find best frame
  max.frame = which.max(c(sum(meta[frames == 0]), sum(meta[frames == 1]),
                        sum(meta[frames == 2]))) - 1
  # subset to best frame
  interval <- interval[c(0,1,2) %% 3 == max.frame] - (pos[interval[1]] %% 3)
  if (feature == "start") {
    means <- c()
    downs <- c()
    ups <- c()
    # upstream window vs downstream window, Check counts in area: pos -18 to -6
    for (j in interval) {
      down <- meta[seq.int(j, max.pos, by = 3)]
      downstream <- mean(down) # down window
      up <- meta[seq.int(1, j - 1, by = 3)]
      upstream <- mean(up) # up window
      m <- downstream - upstream
      downs <- c(downs, downstream); ups <- c(ups, upstream)
      means <- c(means, m)
    }
    # Find center position of frame:
    center.pos <- c(center.pos, center.pos -1, center.pos+2)[max.frame + 1]
    # New scaler, punishes rare extreme regions far away from center position.
    scaled_means <- means / (abs(pos[interval] + center.pos) + 1)^0.25

    # Debug information:
    if (verbose) {
      dt <- data.table(ups, downs, means, scaled_means, pos = pos[interval])
      message("Info / Read length: ", info)
      message("Possible change points of the max frame with sliding windows:")
      print(dt)
    }
    #
    offset <- pos[interval[which.max(abs(scaled_means))]]
  } else if (feature == "stop") {
    shift <- which.max(meta)
    offset <- pos[shift] + 6
  }
  return(offset)
}

#' Convert percentage to ratio of 1
#'
#' 50 -> 0.5 etc, if length cds > minimum.cds
#' @param top_tx numeric
#' @param cds GRangesList object
#' @param minimum.cds numeric, default 1000
#' @keywords internal
#' @return numeric, as ratio of 1
percentage_to_ratio <- function(top_tx, cds, minimum.cds = 1000) {
  if (!is.numeric(top_tx) |top_tx < 0 | top_tx > 100)
    stop("top_tx must be numeric between 0 and 100 %")
  top_tx <- if (length(cds) > minimum.cds) { # as ratio of 1
    (100 - top_tx) / 100
  } else 0
  return(top_tx)
}

# Shift cigar and map back to corrected GRanges.
shift_narrow <- function(footprints, shiftsAll) {
  is_pos <- strandBool(footprints)
  shiftsAll[!is_pos] <- -shiftsAll[!is_pos]
  qnarrow(narrow(footprints), shiftsAll, shiftsAll)
}

#' Pre shifting plot analysis
#'
#' For internal use only!
#' @param rw a data.table of position, score and fraction (read length)
#' of either TIS or TES (translation end site, around 3' UTR)
#' @param heatmap a logical or character string, default FALSE.
#' If TRUE, will plot heatmap of
#' raw reads before p-shifting to console, to see if shifts given make sense.
#' You can also set a filepath to save the file there.
#' @param region a character string, default "start of CDS"
#' @return invisible(NULL)
#' @keywords internal
footprints.analysis <- function(rw, heatmap, region = "start of CDS") {
  region <- paste("Position relative to", region)
  if (heatmap != FALSE) {
    gg <- coverageHeatMap(rw, scoring = "transcriptNormalized",
                          xlab = region)
    plot(gg)
    if (is.character(heatmap)) {
      ggsave(heatmap, gg)
    }
  }
  return(invisible(NULL))
}

#' Get periodogram data per read length
#'
#' A data.table of periods and amplitudes, great to detect
#' ribosomal read lengths. Uses 5' end of reads to detect periodicity.
#' Works both before and after p-shifting. Plot results with ribo_fft_plot.
#' @inheritParams detectRibosomeShifts
#' @param footprints Ribosome footprints in either \code{\link{GAlignments}} or
#' \code{\link{GRanges}}
#' @param cds a \code{\link{GRangesList}} of coding sequences. Length must match
#' length of argument mrna, and all must have length > arugment firstN.
#' @param read_lengths integer vector, default: 26:34,
#'  which read length to check for. Will exclude all read_lengths that
#'  does not exist for footprints.
#' @return a data.table with read_length, amplitude and periods
#' @export
#' @examples
#' ## Note, this sample data is not intended to be strongly periodic.
#' ## Real data should have a cleaner peak for x = 3 (periodicity)
#' # Load sample data
#' df <- ORFik.template.experiment()
#' # Load annotation
#' loadRegions(df, "cds", names.keep = filterTranscripts(df))
#' # Select a riboseq library
#' df <- df[df$libtype == "RFP", ]
#' footprints <- fimport(filepath(df[1,], "default"))
#' fft_dt <-ribo_fft(footprints, cds)
#' ribo_fft_plot(fft_dt)
ribo_fft <- function(footprints, cds, read_lengths = 26:34, firstN = 150) {
  stopifnot(all(widthPerGroup(cds, FALSE) >= firstN))
  # 5' ends only, to detect periodicity
  footprints <- convertToOneBasedRanges(footprints, addSizeColumn = TRUE,
                                        addScoreColumn = TRUE,
                                        along.reference = TRUE)

  # Get a fixed size coverage window
  cov <- windowPerReadLength(cds, cds, footprints,
                             pShifted = FALSE, upstream = 0,
                             downstream = firstN - 1,
                             zeroPosition = 0, scoring = "transcriptNormalized",
                             acceptedLengths = read_lengths,
                             drop.zero.dt = TRUE, append.zeroes = TRUE)
  # Draw spectrogram of FFT
  fft_dt<- data.table()
  read_lengths <- read_lengths[read_lengths %in% readWidths(footprints)]
  for (i in read_lengths) {
    spec <- spec.pgram(x = cov[fraction == i,]$score, plot = FALSE)
    fft_dt <- rbindlist(list(fft_dt, data.table(read_length = i, amplitude = spec$spec, periods = 1 / spec$freq)))
  }
  fft_dt[]
  return(fft_dt)
}

pshifts_export <- function(shifted, name, df, output_format) {
  if ("bedo" %in% output_format) {
    export.bedo(convertToOneBasedRanges(shifted, addScoreColumn = TRUE,
                                        addSizeColumn = TRUE),
                paste0(name, "_pshifted.bedo"))
  }
  if ("ofst" %in% output_format) {
    export.ofst(convertToOneBasedRanges(shifted, addScoreColumn = TRUE,
                                        addSizeColumn = TRUE),
                paste0(name, "_pshifted.ofst"))
  }
  if ("bed" %in% output_format) {
    export.bed(convertToOneBasedRanges(shifted, addScoreColumn = TRUE,
                                       addSizeColumn = FALSE),
               paste0(name, "_pshifted.bed"))
  }
  if ("wig" %in% output_format) {
    export.wiggle(shifted, paste0(name, "_pshifted.wig"))
  }
  if ("bigWig" %in% output_format) {
    if (anyNA(seqlengths(shifted))) {
      seqinfo(shifted) <- seqinfo(df)[seqlevels(shifted),]
    }
    if (anyNA(seqlengths(shifted))) {
      seqinfo(shifted) <- seqinfo(loadTxdb(df))[seqlevels(shifted),]
    }
    if (!anyNA(seqlengths(shifted))) {
      export.bigWig(shifted, paste0(name, "_pshifted.bigWig"))
    } else warning("Could not export bigWig, reads does not have defined seqlengths!")
  }
}

#' Get periodogram plot per read length
#' @param fft_dt a data.table with read_length, amplitude and periods
#' @param period_window x axis limits, default c(0,6)
#' @return a ggplot, geom_line plot facet by read length.
#' @export
#' @examples
#' ## Note, this sample data is not intended to be strongly periodic.
#' ## Real data should have a cleaner peak for x = 3 (periodicity)
#' # Load sample data
#' df <- ORFik.template.experiment()
#' # Load annotation
#' cds <- loadRegion(df, "cds", names.keep = filterTranscripts(df))
#' # Select a riboseq library
#' df <- df[df$libtype == "RFP", ]
#' footprints <- fimport(filepath(df[1,], "default"))
#' fft_dt <-ribo_fft(footprints, cds)
#' ribo_fft_plot(fft_dt)
ribo_fft_plot <- function(fft_dt, period_window = c(0, 6)) {
  stopifnot(is(fft_dt, "data.table"))
  ggplot(fft_dt, aes(x = periods, y = amplitude)) +
    geom_line() +
    facet_wrap( ~ read_length, scales = "free") +
    scale_x_continuous(limits = period_window, breaks = c(1,3,5)) +
    theme_classic()
}

#' Load the shifts from experiment
#'
#' When you p-shift using the function shiftFootprintsByExperiment,
#' you will get a list of shifts per library. To automatically load them, you
#' can use this function. Defaults to loading pshifts, if you made a-sites or
#' e-sites, change the path argument to ashifted/eshifted folder instead.
#' @inheritParams shiftFootprintsByExperiment
#' @param path path, default file.path(libFolder(df), "pshifted", "shifting_table.rds").
#' Path to .rds file containing the shifts as a list,
#' one list element per shifted bam file.
#' @return a list of the shifts, one list element per shifted bam file.
#' @family pshifting
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' # subset on Ribo-seq
#' df <- df[df$libtype == "RFP",]
#' #shiftFootprintsByExperiment(df)
#' #shifts_load(df)
shifts_load <- function(df, path = file.path(libFolder(df), "pshifted",
                                            "shifting_table.rds")) {
  return(readRDS(file = path))
}

#' @inherit shifts_load
shifts.load <- shifts_load

#' Save shifts for Ribo-seq
#'
#' Should be stored in pshifted folder relative to default files
#' @param shifts a list of data.table/data.frame objects.
#' Must be named with the full path to ofst/bam files that defines the shifts.
#' @param folder directory to save file,
#' Usually: file.path(libFolder(df), "pshifted"), where df is the ORFik
#' experiment / or your path of default file types.
#' It will be named file.path(folder, "shifting_table.rds").
#' For ORFik to work optimally,
#' the folder should be the /pshifted/ folder relative to default files.
#' @return invisible(NULL), file saved to disc as "shifting_table.rds".
#' @family pshifting
#' @export
#' @examples
#' df <- ORFik.template.experiment.zf()
#' shifts <- shifts_load(df)
#' original_shifts <- file.path(libFolder(df), "pshifted", "shifting_table.rds")
#' # Move to temp
#' new_shifts_path <- file.path(tempdir(), "shifting_table.rds")
#' new_shifts <- c(shifts, shifts)
#' names(new_shifts)[2] <- file.path(tempdir(), "RiboSeqTemp.ofst")
#' saveRDS(new_shifts, new_shifts_path)
#' new_shifts[[1]][1,2] <- -10
#' # Now update the new shifts, here we input only first
#' shifts_save(new_shifts[1], tempdir())
#' readRDS(new_shifts_path) # You still get 2 outputs
#'
shifts_save <- function(shifts, folder) {
  stopifnot(is(shifts, "list"))
  stopifnot(is(shifts[[1]], "data.frame"))
  if (length(shifts) == 0) {
    warning("Tried to save shift table of length 0, returning without saving!")
  }
  stopifnot(length(unique(names(shifts))) == length(shifts) & !anyNA(names(shifts)))
  names(shifts) <- pasteDir(names(shifts))
  folder <- pasteDir(folder)
  shift_table_path <- file.path(folder, "shifting_table.rds")
  if (file.exists(shift_table_path)) {
    old_shifts <- readRDS(shift_table_path)
    identical_libs <- identical(names(shifts), names(old_shifts))
    is_subset <- all(names(shifts) %in% names(old_shifts))
    if (!identical_libs & is_subset) {
      old_shifts[names(shifts)] <- shifts
      shifts <- old_shifts
    }
  }
  saveRDS(shifts, file = shift_table_path)
  return(invisible(NULL))
}




