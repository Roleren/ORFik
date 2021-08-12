#' Shift footprints by selected offsets
#'
#' Function shifts footprints (GRanges) using specified offsets for every of
#' the specified lengths. Reads that do not conform to the specified lengths
#' are filtered out and rejected. Reads are resized to single base in 5' end
#' fashion, treated as p site.
#' This function takes account for junctions in cigars of the reads. Length of
#' the footprint is saved in size' parameter of GRanges output. Footprints are
#' also sorted according to their genomic position, ready to be saved as a
#' ofst, bed or wig file.
#'
#' The two columns in the shift data.frame/data.table argument are:\cr
#' - fraction Numeric vector of lengths of footprints you select
#' for shifting.\cr
#' - offsets_start Numeric vector of shifts for corresponding
#' selected_lengths. eg. c(-10, -10) with selected_lengths of c(31, 32) means
#' length of 31 will be shifted left by 10. Footprints of length 32 will be
#' shifted right by 10.
#'
#' NOTE: It will remove softclips from valid width, the CIGAR 3S30M is qwidth
#' 33, but will remove 3S so final read width is 30 in ORFik.
#' @param footprints \code{\link{GAlignments}} object of RiboSeq reads
#' @param shifts a data.frame / data.table with minimum 2 columns,
#' fraction (selected read lengths) and offsets_start (relative position in nt).
#' Output from \code{\link{detectRibosomeShifts}}.\cr
#' Run \code{ORFik::shifts.load(df)[[1]]} for an example of input.
#' @param sort logical, default TRUE. If False will keep original order of
#' reads, and not sort output reads in increasing genomic location per
#' chromosome and strand.
#' @return A \code{\link{GRanges}} object of shifted footprints, sorted and
#' resized to 1bp of p-site,
#' with metacolumn "size" indicating footprint size before shifting and
#' resizing, sorted in increasing order.
#' @family pshifting
#' @references https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4912-6
#' @export
#' @examples
#' ## Basic run
#' # Transcriptome annotation ->
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' # Ribo seq data ->
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' \dontrun{
#' footprints <- readBam(riboSeq_file)
#'
#' # detect the shifts automagically
#' shifts <- detectRibosomeShifts(footprints, gtf_file)
#' # shift the RiboSeq footprints
#' shiftedReads <- shiftFootprints(footprints, shifts)
#' }
shiftFootprints <- function(footprints, shifts, sort = TRUE) {
  if (!is(shifts, "data.frame")) stop("shifts must be data.frame/data.table")
  if (nrow(shifts) == 0) stop("No shifts found in data.frame")
  if (is.null(shifts$fraction) | is.null(shifts$offsets_start))
    stop("Either fraction (read lengths) or offsets_start (shifts by nt) column in shifts is not set!")

  selected_lengths <- shifts$fraction
  selected_shifts <- shifts$offsets_start

  lengthsAll <- readWidths(footprints, along.reference = TRUE)
  validLengths <- lengthsAll %in% selected_lengths
  lengthsAll <- lengthsAll[validLengths]
  footprints <- footprints[validLengths]
  shiftsAll <- -selected_shifts[match(lengthsAll,
                          selected_lengths)] + 1

  # Shift cigar and map back to corrected GRanges.
  is_pos <- strandBool(footprints)

  shiftsAll[!is_pos] <- -shiftsAll[!is_pos]
  shifted <- qnarrow(footprints, shiftsAll, shiftsAll)
  shifted <- resize(GRanges(shifted), width = 1)
  shifted$size <- lengthsAll

  shifted <- sortSeqlevels(shifted)
  if (sort) {
    message("Sorting shifted footprints...")
    shifted <- sort(shifted)
  }
  return(shifted)
}

#' Detect ribosome shifts
#'
#' Utilizes periodicity measurement (Fourier transform), and change point
#' analysis to detect ribosomal footprint shifts for each of the ribosomal
#' read lengths. Returns subset of read lengths and their shifts for which
#' top covered transcripts follow periodicity measure. Each shift value
#' assumes 5' anchoring of the reads, so that output offsets values will
#' shift 5' anchored footprints to be on the p-site of the ribosome. The
#' E-site will be shift + 3 and A site will be shift - 3. So update to these,
#' if you rather want those.
#'
#' Check out vignette for the examples of plotting RiboSeq metaplots over start
#' and stop codons, so that you can verify visually whether this function
#' detects correct shifts.
#'
#' For how the Fourier transform works, see: \code{\link{isPeriodic}}\cr
#' For how the changepoint analysis works, see: \code{\link{changePointAnalysis}}\cr
#'
#' NOTE: It will remove softclips from valid width, the CIGAR 3S30M is qwidth
#' 33, but will remove 3S so final read width is 30 in ORFik. This is standard
#' for ribo-seq.
#'
#' @param footprints \code{\link{GAlignments}} object of RiboSeq reads -
#' footprints, can also be path to the .bam /.ofst file. If GAlignment object
#' has a meta column called "score", this will be used as replicate numbering
#' for that read. So be careful if you have custom files with score columns,
#' with another meaning.
#' @inheritParams loadTxdb
#' @param start (logical) Whether to include predictions based on the start
#' codons. Default TRUE.
#' @param stop (logical) Whether to include predictions based on the stop
#' codons. Default FASLE. Only use if there exists 3' UTRs for the annotation.
#' If peridicity around stop codon is stronger than at the start codon, use
#' stop instead of start region for p-shifting.
#' @param top_tx (integer), default 10. Specify which \% of the top covered by RiboSeq
#' reads transcripts to use for estimation of the shifts. By default we take top 10%
#' top covered transcripts as they represent less noisy dataset. This is only
#' applicable when there are more than 1000 transcripts.
#' @inheritParams filterTranscripts
#' @param txNames a character vector of subset of CDS to use. Default:
#' txNames = filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR)\cr
#' Example:
#' c("ENST1000005"), will use only that transcript (You should use at least 100!).
#' Remember that top_tx argument, will by default specify to use top 10 \%
#' of those CDSs. Set that to 100, to use all these specified transcripts.
#' @param firstN (integer) Represents how many bases of the transcripts
#' downstream of start codons to use for initial estimation of the
#' periodicity.
#' @param tx a GRangesList, if you do not have 5' UTRs in annotation, send
#' your own version. Example: extendLeaders(tx, 30)
#' Where 30 bases will be new "leaders". Since each original transcript was
#' either only CDS or non-coding (filtered out).
#' @param min_reads default (1000), how many reads must a read-length have to
#' be considered for periodicity.
#' @param accepted.lengths accepted readlengths, default 26:34, usually ribo-seq
#' is strongest between 27:32.
#' @inheritParams footprints.analysis
#' @param must.be.periodic logical TRUE, if FALSE will not filter on
#' periodic read lengths. (The Fourier transform filter will be skipped).
#' This is useful if you are not going to do periodicity analysis, that is:
#' for you more coverage depth (more read lengths)
#' is more important than only keeping the high quality periodic read lengths.
#' @inheritParams isPeriodic
#' @return a data.table with lengths of footprints and their predicted
#' coresponding offsets
#' @family pshifting
#' @references https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4912-6
#' @importFrom IRanges quantile
#' @export
#' @examples
#' ## Basic run
#' # Transcriptome annotation ->
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' # Ribo seq data ->
#' riboSeq_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' \dontrun{
#' footprints <- readBam(riboSeq_file)
#' ## Using CDS start site as reference point:
#' detectRibosomeShifts(footprints, gtf_file)
#' ## Using CDS start site and stop site as 2 reference points:
#' #detectRibosomeShifts(footprints, gtf_file, stop = TRUE)
#' ## Debug and detailed information for accepted reads lengths and p-site:
#' detectRibosomeShifts(footprints, gtf_file, heatmap = TRUE, verbose = TRUE)
#' ## Debug why read length 31 was not accepted or wrong p-site:
#' #detectRibosomeShifts(footprints, gtf_file, must.be.periodic = FALSE,
#' #              accepted.lengths = 31, heatmap = TRUE, verbose = TRUE)
#'
#' ## Subset bam file
#' param = ScanBamParam(flag = scanBamFlag(
#'                        isDuplicate = FALSE,
#'                        isSecondaryAlignment = FALSE))
#' footprints <- readBam(riboSeq_file, param = param)
#' detectRibosomeShifts(footprints, gtf_file, stop = TRUE)
#'
#' ## Without 5' Annotation
#' library(GenomicFeatures)
#'
#' txdb <- loadTxdb(gtf_file)
#' tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#' tx <- extendLeaders(tx, 30)
#' ## Now run function, without 5' and 3' UTRs
#' detectRibosomeShifts(footprints, txdb, start = TRUE, minFiveUTR = NULL,
#'                      minCDS = 150L, minThreeUTR = NULL, firstN = 150L,
#'                      tx = tx)
#'
#' }
#'
detectRibosomeShifts <- function(footprints, txdb, start = TRUE, stop = FALSE,
  top_tx = 10L, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = 30L,
  txNames = filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR),
  firstN = 150L, tx = NULL, min_reads = 1000, accepted.lengths = 26:34,
  heatmap = FALSE, must.be.periodic = TRUE, strict.fft = TRUE, verbose = FALSE) {

  txdb <- loadTxdb(txdb)
  cds <- loadRegion(txdb, part = "cds", names.keep = txNames)
  footprints <- fimport(footprints, cds)

  # reduce data-set to only matching seqlevels
  seqMatch <- validSeqlevels(cds, footprints)
  cds <- keepSeqlevels(cds, seqMatch, pruning.mode = "coarse")
  footprints <- keepSeqlevels(footprints, seqMatch, pruning.mode = "coarse")
  if (length(cds) == 0 | length(footprints) == 0) {
    stop("txdb and footprints did not have any matched seqnames")
  }
  txNames <- txNames[txNames %in% names(cds)]
  if (is.null(tx)) tx <- loadRegion(txdb, part = "tx")
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
                            " accepted.lengths and counts > min_reads")
  cds <- cds[countOverlapsW(cds, footprints, "score") > 0]
  if (verbose) message("Number of CDSs used for p-site detection: ", length(cds))
  top_tx <- percentage_to_ratio(top_tx, cds)
  if (must.be.periodic) {
    periodicity <- windowPerReadLength(cds, tx, footprints,
                                       pShifted = FALSE, upstream = 0,
                                       downstream = firstN - 1,
                                       zeroPosition = 0, scoring = "transcriptNormalized",
                                       acceptedLengths = tab$size,
                                       drop.zero.dt = TRUE, append.zeroes = TRUE)
    periodicity <- periodicity[, .(score = isPeriodic(score, unique(fraction), verbose, strict.fft)),
                               by = fraction]
    validLengths <- periodicity[score == TRUE,]$fraction
  } else validLengths <- accepted.lengths
  if (length(validLengths) == 0)
    stop(paste("Library contained no periodic or accepted read-lengths,",
         "check your library. Are you using the correct genome?",
         "Is this Ribo-seq?"))
  # find shifts
  if (start) {
    rw <- windowPerReadLength(cds, tx, footprints, pShifted,
                              upstream = 30, downstream = 29,
                              acceptedLengths = validLengths,
                              scoring = NULL)
    rw[, sum.count := sum(count), by = genes]
    rw <- rw[sum.count >= quantile(sum.count, top_tx), ]
    rw <- coverageScorings(rw, scoring = "sum")
    footprints.analysis(rw, heatmap)
    if (verbose) {
      message("-------------------")
      message("Change point analysis (start): ups (upstream window score), downs (downstream window score)")
    }
    offset <- rw[, .(offsets_start = changePointAnalysis(score, info = unique(fraction),
                                                         verbose = verbose)),
                 by = fraction]
  }
  if (stop & !is.null(minThreeUTR)) {
    threeUTRs <- loadRegion(txdb, "trailers", names.keep = txNames)
    rw <- windowPerReadLength(threeUTRs, tx, footprints, FALSE,
                              upstream = 30, downstream = 29,
                              acceptedLengths = validLengths,
                              scoring = "sum")
    footprints.analysis(rw, heatmap, region = "start of 3' UTR")
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
#' A function that combines the steps of periodic read length detection,
#' p-site shift detection and p-shifting into 1 function.
#' For more details, see: \code{\link{detectRibosomeShifts}}\cr
#' Saves files to a specified location as .ofst and .wig,
#' The .ofst file will include a score column containing read width. \cr
#' The .wig files, will be saved in pairs of +/- strand, and score column
#' will be replicates of reads starting at that position,
#' score = 5 means 5 reads.\cr
#' Remember that different species might have different default Ribosome
#' read lengths, for human, mouse etc, normally around 27:30.
#' @inheritParams detectRibosomeShifts
#' @param df an ORFik \code{\link{experiment}}
#' @param out.dir output directory for files,
#' default: dirname(df$filepath[1]), making a /pshifted
#' folder at that location
#' @param output_format default c("ofst", "wig"), use export.ofst or
#' wiggle format (wig) using \code{\link{export.wiggle}} ? Default is both.
#' The wig format version can be used in IGV, the score column is counts of that
#' read with that read length, the cigar reference width is lost,
#' ofst is much faster to save and load in R, and retain cigar reference width,
#' but can not be used in IGV. \cr You can also do bedoc format, bed format
#' keeping cigar: \code{\link{export.bedoc}}. bedoc is usually not used for
#' p-shifting.
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @param log logical, default (TRUE), output a log file with parameters used and
#' a .rds file with all shifts per library
#' (can be loaded with \code{\link{shifts.load}})
#' @param shift.list default NULL, or a list containing named data.frames / data.tables
#' with minimum 2 columns, fraction (selected read lengths) and
#' offsets_start (relative position in nt). 1 named data.frame / data.table per library.
#' Output from \code{\link{detectRibosomeShifts}}.\cr
#' Run \code{ORFik::shifts.load(df)} for an example of input. The names of the list must
#' be the file.paths of the Ribo-seq libraries. Use this to edit the shifts, if
#' you suspect some of them are wrong in an experiment.
#' @return NULL (Objects are saved to out.dir/pshited/"name_pshifted.ofst",
#' wig, bedo or .bedo)
#' @importFrom rtracklayer export.bed
#' @importFrom utils packageVersion
#' @family pshifting
#' @references https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4912-6
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[3,] #lets only p-shift RFP sample at index 3
#' ## If you want to check it in IGV do:
#' shiftFootprintsByExperiment(df)
#' # Then use the .wig files that are created, which are readable in IGV.
#' # If you only need in R, do: (then you get no .wig files)
#' #shiftFootprintsByExperiment(df, output_format = "ofst")
#' ## With debug info:
#' #shiftFootprintsByExperiment(df, verbose = TRUE)
#' ## Re-shift, if you think some are wrong
#' ## Here we update library 1, third read length to shift 12
#' shift.list <- shifts.load(df)
#' shift.list[[1]]$offsets_start[3] <- -12
#' #shiftFootprintsByExperiment(df, shift.list = shift.list)
shiftFootprintsByExperiment <- function(df,
                                        out.dir = pasteDir(dirname(
                                          df$filepath[1]), "/pshifted/"),
                                        start = TRUE, stop = FALSE,
                                        top_tx = 10L, minFiveUTR = 30L,
                                        minCDS = 150L, minThreeUTR = 30L,
                                        firstN = 150L, min_reads = 1000,
                                        accepted.lengths = 26:34,
                                        output_format = c("ofst", "wig"),
                                        BPPARAM = bpparam(), tx = NULL,
                                        shift.list = NULL,
                                        log = TRUE, heatmap = FALSE,
                                        must.be.periodic = TRUE, strict.fft = TRUE,
                                        verbose = FALSE) {
  path <- out.dir
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(path)) stop(paste("out.dir", out.dir, "does not exist!"))
  if (!any(c("bed", "bedo", "wig", "ofst") %in% output_format))
    stop("output_format allowed bed, bedo, wig or ofst")
  rfpFiles <- filepath(df, "ofst") # If ofst file not present, uses bam file
  if (!is.null(shift.list)) {
    if (!all(names(shift.list) %in% rfpFiles))
      stop("shift.list does not contain all files to be shifted!")
  }
  for (out.form in output_format)
    message(paste("Saving", out.form, "files to:", out.dir))
  message(paste("Shifting reads in experiment:", df@experiment))

  txdb <- loadTxdb(df)
  tx <- loadRegion(txdb, part = "mrna")
  txNames <- filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR)


  shifts <- bplapply(rfpFiles,
           FUN = function(file, path, df, start, stop,
                          top_tx, minFiveUTR, minCDS, minThreeUTR,
                          firstN, min_reads, accepted.lengths,
                          output_format, heatmap, tx, shift.list,
                          must.be.periodic, txNames, strict.fft,
                          verbose = verbose
                          ) {
    message(file)
    rfp <- fimport(file)
    if (!is.null(shift.list)) { # Pre defined shifts
      shifts <- shift.list[file][[1]]
    } else {
      shifts <- detectRibosomeShifts(rfp, txdb = loadTxdb(df), start = start,
                                     stop = stop, top_tx = top_tx,
                                     minFiveUTR = minFiveUTR,
                                     minCDS = minCDS, minThreeUTR = minThreeUTR,
                                     txNames = txNames, firstN = firstN,
                                     min_reads = min_reads,
                                     accepted.lengths = accepted.lengths,
                                     heatmap = heatmap, tx = tx,
                                     must.be.periodic = must.be.periodic,
                                     strict.fft = strict.fft, verbose = verbose)
    }

    shifted <- shiftFootprints(rfp, shifts)
    name <- paste0(path, remove.file_ext(file, basename = TRUE))

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

    return(shifts)
  }, path = path, df = df, start = start, stop = stop,
     top_tx = top_tx, minFiveUTR = minFiveUTR,
     minCDS = minCDS, minThreeUTR = minThreeUTR,
     firstN = firstN, min_reads = min_reads,
     accepted.lengths = accepted.lengths, output_format = output_format,
     heatmap = heatmap, must.be.periodic = must.be.periodic,
     strict.fft = strict.fft, verbose = verbose, tx = tx,
     shift.list = shift.list, txNames = txNames, BPPARAM = BPPARAM)

  if (log) {
    fileConn<-file(paste0(path, "/pshifting_arguments.txt"), "w")
    cat(paste("From ORFik version:", packageVersion("ORFik"), "\n"),
        file = fileConn)
    cat("All arguments not specificed below are default:\n", file = fileConn)
    cat(paste(as.character(sys.call()), "\n"), file = fileConn)
    close(fileConn)
    # Save shifts
    names(shifts) <- rfpFiles
    saveRDS(shifts, file = file.path(path, "shifting_table.rds"))
  }
  if (verbose) {
    message("Shifting done, detected shifts per file:")
    print(shifts)
  }

  return(invisible(NULL))
}

#' Plot shifted heatmaps per library
#'
#' Around CDS TISs, plot coverage.
#' A good validation for you p-shifting, to see shifts are corresponding
#' and close to the CDS TIS.
#' @inheritParams shiftFootprintsByExperiment
#' @inheritParams windowPerReadLength
#' @param output name to save file, full path. (Default NULL) No saving.
#' Sett to "auto" to save to QC_STATS folder of experiment named:
#' "pshifts_barplots.png" or "pshifts_heatmaps.png" depending on type argument.
#' Folder must exist!
#' @param scoring which scoring scheme to use for heatmap, default
#' "transcriptNormalized". Some alternatives: "sum", "zscore".
#' @param type character, default "bar". Plot as faceted bars,
#' gives more detailed information of read lengths,
#' but harder to see patterns over multiple read lengths.
#' Alternative: "heatmaps", better overview of patterns over
#' multiple read lengths.
#' @param title Title for top of plot, default "Ribo-seq".
#' A more informative name could be "Ribo-seq zebrafish Chew et al. 2013"
#' @param addFracPlot logical, default TRUE, add positional sum plot on top
#' per heatmap.
#' @param plot.ext default ".pdf". Alternative ".png". Only added if output is
#' "auto".
#' @importFrom gridExtra grid.arrange
#' @return a ggplot2 grob object
#' @family pshifting
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[3,] #lets only p-shift RFP sample at index 3
#' #shiftFootprintsByExperiment(df, output_format = "bedo)
#' #grob <- shiftPlots(df, title = "Ribo-seq Human ORFik et al. 2020")
#' #plot(grob) #Only plot in RStudio for small amount of files!
shiftPlots <- function(df, output = NULL, title = "Ribo-seq",
                       scoring = "transcriptNormalized",
                       pShifted = TRUE,
                       upstream = if (pShifted) 5 else 20,
                       downstream = if (pShifted) 20 else 5,
                       type = "bar",
                       addFracPlot = TRUE,
                       plot.ext = ".pdf",
                       BPPARAM = bpparam()) {
  stopifnot(plot.ext %in% c(".png", ".pdf"))
  if (!(type %in% c("bar", "heatmap")))
    stop("The 'type' argument must be bar or heatmap")
  txdb <- loadTxdb(df)
  txNames <- filterTranscripts(txdb, upstream, downstream + 1, 0)
  cds <-  loadRegion(txdb, part = "cds", names.keep = txNames)
  mrna <- loadRegion(txdb, part = "mrna", names.keep = txNames)
  style <- seqlevelsStyle(cds)
  plots <- bplapply(seq(nrow(df)),
                    function(x, cds, mrna, style, paths, df, upstream,
                             downstream, type) {
    miniTitle <- gsub("_", " ", bamVarName(df, skip.experiment = TRUE)[x])
    hitMap <- windowPerReadLength(cds, mrna,  fimport(paths[x], style),
                                  upstream = upstream, downstream = downstream)
    if (type == "heatmap") {
      coverageHeatMap(hitMap, scoring = scoring, addFracPlot = addFracPlot,
                      title = miniTitle)
    } else {
      hitMap[, frame := position %% 3]
      pSitePlot(hitMap, scoring = scoring,
                facet = TRUE, frameSum = TRUE, title = miniTitle)
    }
  }, cds = cds, mrna = mrna, style = style, BPPARAM = BPPARAM,
     paths = filepath(df, "pshifted"), df = df, upstream = upstream,
     downstream = downstream, type = type)
  res <- do.call("arrangeGrob", c(plots, ncol=1, top = title))
  if (!is.null(output)) {
    if (type == "heatmap") {
      if (output == "auto") {
        dir.to.save <- file.path(dirname(df$filepath[1]), "QC_STATS")
        output <- file.path(dir.to.save, paste0("pshifts_heatmaps", plot.ext))
      }
      ggsave(output, res,
             width = 225, height = (length(res) -1)*65,
             units = "mm", dpi = 300, limitsize = FALSE)
    } else {
      if (output == "auto") {
        dir.to.save <- file.path(dirname(df$filepath[1]), "QC_STATS")
        output <- file.path(dir.to.save, paste0("pshifts_barplots", plot.ext))
      }
      dpi <- ifelse(nrow(df) < 22, 300, 200)
      ggsave(output, res,
             width = 225, height = (length(res) -1)*95,
             units = "mm", dpi = dpi, limitsize = FALSE)
    }
    message("Saved pshift plots to location: ",
            output)
  }

  return(res)
}
