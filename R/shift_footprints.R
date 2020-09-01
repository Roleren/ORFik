#' Shift footprints by selected offsets
#'
#' Function shifts footprints (GRanges) using specified offsets for every of
#' the specified lengths. Reads that do not conform to the specified lengths
#' are filtered out and rejected. Reads are resized to single base in 5' end
#' fashion, treated as p site.
#' This function takes account for junctions in cigars of the reads. Length of
#' the footprint is saved in size' parameter of GRanges output. Footprints are
#' also sorted according to their genomic position, ready to be saved as a
#' bed or wig file.
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
#' fraction (selected_lengths) and selected_shifts (relative position).
#' Output from \code{\link{detectRibosomeShifts}}
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
#' # Basic run
#' #shiftFootprints(footprints, shifts)
#' # Full example
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
shiftFootprints <- function(footprints, shifts, sort = TRUE) {
  if (!is(shifts, "data.frame")) stop("shifts must be data.frame/data.table")
  if (nrow(shifts) == 0) stop("No shifts found in data.frame")
  if (is.null(shifts$fraction) | is.null(shifts$offsets_start))
    stop("Either fraction or offsets_start column in shifts is not set!")

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
  shifted <- GRanges(shifted)
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
#' shift 5' anchored footprints to be on the p-site of the ribosome.
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
#' for that read.
#' @inheritParams loadTxdb
#' @param start (logical) Whether to include predictions based on the start
#' codons. Default TRUE.
#' @param stop (logical) Whether to include predictions based on the stop
#' codons. Default FASLE. Only use if there exists 3' UTRs for the annotation.
#' If peridicity around stop codon is stronger than at the start codon, use
#' stop instead of start region for p-shifting.
#' @param top_tx (integer), default 10. Specify which % of the top covered by RiboSeq
#' reads transcripts to use for estimation of the shifts. By default we take top 10%
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
#' @param accepted.lengths accepted readlengths, default 26:34, usually ribo-seq
#' is strongest between 27:32.
#' @inheritParams footprints.analysis
#' @param must.be.periodic logical TRUE, if FALSE will not filter on
#' periodic read lengths. (The Fourier transform filter will be skipped).
#' @return a data.table with lengths of footprints and their predicted
#' coresponding offsets
#' @family pshifting
#' @references https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4912-6
#' @importFrom IRanges quantile
#' @export
#' @examples
#' # Basic run
#' #detectRibosomeShifts(footprints, txdb)
#' # Full example
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
  firstN = 150L, tx = NULL, min_reads = 1000, accepted.lengths = 26:34,
  heatmap = FALSE, must.be.periodic = TRUE) {

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
  scoreExists <- "score" %in% colnames(mcols(a))
  footprints <- convertToOneBasedRanges(footprints, addSizeColumn = TRUE,
                                        addScoreColumn = !scoreExists,
                                        along.reference = TRUE)
  # Filter if < 1000 counts read size or not in accepted.lengths
  lengths <- data.table(score = footprints$score, size = footprints$size)
  tab <- lengths[, .(counts = sum(score)), by = size]
  tab <- tab[(counts >= min_reads) & (size %in% accepted.lengths),]
  if (nrow(tab) == 0) stop("No valid read lengths found with",
                            "accepted.lengths and counts > min_reads")
  cds <- cds[countOverlapsW(cds, footprints, "score") > 0]
  top_tx <- percentage_to_ratio(top_tx, cds)
  if (must.be.periodic) {
    periodicity <- windowPerReadLength(cds, tx, footprints,
                                       pShifted = FALSE, upstream = 0,
                                       downstream = firstN - 1,
                                       zeroPosition = 0, scoring = "periodic",
                                       acceptedLengths = tab$size)
    validLengths <- periodicity[score == TRUE,]$fraction
  } else validLengths <- accepted.lengths

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
    offset <- rw[, .(offsets_start = changePointAnalysis(score)),
                 by = fraction]
    # Figure if we want to keep this.
    offset <- offset[offsets_start < 0, ]
  }
  if (stop & !is.null(minThreeUTR)) {
    threeUTRs <- loadRegion(txdb, "trailers")[txNames]
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
#' Saves files to a specified location as .ofst and .wig,
#' .ofst will include a score column containing read width. \cr
#' For more details, see: \code{\link{detectRibosomeShifts}}
#'
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
#' @param log logical, default (TRUE), output a log file with parameters used.
#' @return NULL (Objects are saved to out.dir/pshited/"name_pshifted.ofst",
#' wig, bedo or .bedo)
#' @importFrom rtracklayer export.bed
#' @family pshifting
#' @references https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4912-6
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[3,] #lets only p-shift RFP sample at index 3
#' # If you want to check it in IGV do:
#' shiftFootprintsByExperiment(df)
#' # Then use the bed files that are created, which are readable in IGV.
#' # If you only need in R, do: (then you get no .wig files)
#' #shiftFootprintsByExperiment(df, output_format = "ofst")
shiftFootprintsByExperiment <- function(df,
                                        out.dir = pasteDir(dirname(
                                          df$filepath[1]), "/pshifted/"),
                                        start = TRUE, stop = FALSE,
                                        top_tx = 10L, minFiveUTR = 30L,
                                        minCDS = 150L, minThreeUTR = 30L,
                                        firstN = 150L, min_reads = 1000,
                                        accepted.lengths = 26:34,
                                        output_format = c("ofst", "wig"),
                                        BPPARAM = bpparam(),
                                        log = TRUE) {
  path <- out.dir
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(path)) stop(paste("out.dir", out.dir, "does not exist!"))
  if (!any(c("bed", "bedo", "wig", "ofst") %in% output_format))
    stop("output_format allowed bed, bedo, wig or ofst")
  for (out.form in output_format)
    message(paste("Saving", out.form, "files to:", out.dir))
  message(paste("Shifting reads in experiment:", df@experiment))

  txdb <- loadTxdb(df)
  rfpFiles <- filepath(df, "ofst") # If ofst file not present, use bam file
  bplapply(rfpFiles, FUN = function(file, path, df, start, stop,
                                    top_tx, minFiveUTR, minCDS, minThreeUTR,
                                    firstN, min_reads, accepted.lengths,
                                    output_format) {
    message(file)
    rfp <- fimport(file)
    shifts <- detectRibosomeShifts(rfp, txdb = loadTxdb(df), start = start,
                                   stop = stop, top_tx = top_tx,
                                   minFiveUTR = minFiveUTR,
                                   minCDS = minCDS, minThreeUTR = minThreeUTR,
                                   firstN = firstN, min_reads = min_reads,
                                   accepted.lengths = accepted.lengths)
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
      export.bed(convertToOneBasedRanges(shifted, addScoreColumn = TRUE,
                                         addSizeColumn = FALSE),
                 paste0(name, "_pshifted.wig"))
    }

    return(invisible(NULL))
  }, path = path, df = df, start = start, stop = stop,
      top_tx = top_tx, minFiveUTR = minFiveUTR,
      minCDS = minCDS, minThreeUTR = minThreeUTR,
      firstN = firstN, min_reads = min_reads,
      accepted.lengths = accepted.lengths, output_format = output_format,
      BPPARAM = BPPARAM)

  if (log) {
    fileConn<-file(paste0(path, "/pshifting_arguments.txt"))
    writeLines("All arguments not specificed below are default:", fileConn)
    writeLines(as.character(sys.call()), fileConn)
    close(fileConn)
  }
  return(invisible(NULL))
}

#' Plot shifted heatmaps per library
#'
#' A good validation for you p-shifting
#' @inheritParams shiftFootprintsByExperiment
#' @param output name to save file, full path. (Default NULL) No saving.
#' @param scoring which scoring scheme to use for heatmap, default
#' "transcriptNormalized". Some alternatives: "sum", "zscore".
#' @param title Title for top of plot, default "Ribo-seq".
#' A more informative name could be "Ribo-seq zebrafish Chew et al. 2013"
#' @param addFracPlot logical, default TRUE, add positional sum plot on top
#' per heatmap.
#' @importFrom gridExtra grid.arrange
#' @return a ggplot2 grob object
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[3,] #lets only p-shift RFP sample at index 3
#' #shiftFootprintsByExperiment(df, output_format = "bedo)
#' #shiftPlots(df, title = "Ribo-seq Human ORFik et al. 2020")
shiftPlots <- function(df, output = NULL, title = "Ribo-seq",
                       scoring = "transcriptNormalized",
                       addFracPlot = TRUE,
                       BPPARAM = bpparam()) {
  txNames <- filterTranscripts(df, 20, 21, 1)
  txdb <- loadTxdb(df)
  cds <-  loadRegion(txdb, part = "cds", names.keep = txNames)
  mrna <- loadRegion(txdb, part = "mrna", names.keep = txNames)
  style <- seqlevelsStyle(cds)
  plots <- bplapply(seq(nrow(df)),
                    function(x, cds, mrna, style, paths, df) {
    miniTitle <- gsub("_", " ", bamVarName(df, skip.experiment = TRUE)[x])
    hitMap <- windowPerReadLength(cds, mrna,  fimport(paths[x], style),
                                  pShifted = TRUE)
    coverageHeatMap(hitMap, scoring = scoring, addFracPlot = addFracPlot,
                    title = miniTitle)
  }, cds = cds, mrna = mrna, style = style, BPPARAM = BPPARAM,
  paths = filepath(df, "pshifted"), df = df)
  res <- do.call("grid.arrange", c(plots, ncol=1, top = title))
  if (!is.null(output))
    ggsave(output, res,
           width = 225, height = (length(res) -1)*65,
           units = "mm", dpi = 300)
  return(res)
}
