#' A post Alignment quality control of reads
#'
#' The ORFik QC uses the aligned files (usually bam files),
#' fastp and STAR log files
#' combined with annotation to create relevant statistics.\cr\cr
#' This report consists of several steps:\cr
#' 1. Convert bam file / Input files to ".ofst" format, if not already done.
#' This format is around 400x faster to use in R than the bam format.
#' Files are also outputted to R environment specified by \code{envExp(df)}\cr
#' 2. From this report you will get a summary csv table, with distribution of
#' aligned reads and overlap counts over transcript regions like:
#' leader, cds, trailer, lincRNAs, tRNAs, rRNAs, snoRNAs etc. It will be called
#' STATS.csv. And can be imported with \code{\link{QCstats}} function.\cr
#' 3. It will also make correlation plots and meta coverage plots,
#' so you get a good understanding of how good the quality of your NGS
#' data production + aligner step were.\cr
#' 4. Count tables are produced, similar to HTseq count tables.
#' Over mrna, leader, cds and trailer separately. This tables
#' are stored as \code{\link{SummarizedExperiment}}, for easy loading into
#' DEseq, conversion to normalized fpkm values,
#' or collapsing replicates in an experiment.
#' And can be imported with \code{\link{countTable}} function.\cr\cr
#' Everything will be outputed in the directory of your NGS data,
#' inside the folder ./QC_STATS/, relative to data location in 'df'.
#' You can specify new out location with out.dir if you want.\cr
#' To make a ORFik experiment, see ?ORFik::experiment \cr
#' To see some normal mrna coverage profiles of different RNA-seq protocols:
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4310221/figure/F6/}
#' @inheritParams outputLibs
#' @param out.dir character, output directory, default:
#' \code{resFolder(df)}.
#' Will make a folder within this called "QC_STATS" with all results in this directory.
#' Warning: If you assign not default path, you will have a hazzle to load files later.
#' Much easier to load count tables, statistics, ++ later with default. Update
#' resFolder of df instead if needed.
#' @param plot.ext character, default: ".pdf". Alternatives: ".png" or ".jpg".
#' Note that in pdf format the complex correlation plots become very slow to load!
#' @param create.ofst logical, default TRUE. Create ".ofst" files from the input
#' libraries, ofst is much faster to load in R, for later use. Stored
#' in ./ofst/ folder relative to experiment main folder.
#' @param complex.correlation.plots logical, default TRUE. Add in addition
#' to simple correlation plot two computationally heavy dots + correlation plots.
#' Useful for deeper analysis, but takes longer time to run, especially on low-quality
#' gpu computers. Set to FALSE to skip these.
#' @return invisible(NULL) (objects are stored to disc)
#' @family QC report
#' @importFrom utils write.csv
#' @export
#' @examples
#' # Load an experiment
#' df <- ORFik.template.experiment()
#' # Run QC
#' #QCreport(df, tempdir())
#' # QC on subset
#' #QCreport(df[9,], tempdir())
QCreport <- function(df, out.dir = resFolder(df),
                     plot.ext = ".pdf", create.ofst = TRUE,
                     complex.correlation.plots = TRUE,
                     BPPARAM = bpparam()) {
  # Check input
  validateExperiments(df)
  stopifnot(plot.ext %in% c(".pdf", ".png"))
  stopifnot(create.ofst %in% c(TRUE, FALSE))
  stopifnot(is.character(out.dir))

  message("Started ORFik QC report for experiment: ", df@experiment)
  stats_folder <- pasteDir(out.dir, "/QC_STATS/")
  if (!dir.create(stats_folder, recursive = TRUE)) {
    if (!dir.exists(stats_folder)) stop("Could not create output directory!")
  }
  if (create.ofst) {
    convertLibs(df, reassign.when.saving = TRUE)
  }
  message("--------------------------")
  message("- Creating read length tables:")
  dt_read_lengths <- readLengthTable(df, output.dir = stats_folder)
  # Get count tables
  QC_count_tables(df, out.dir, force = FALSE,  BPPARAM = BPPARAM)
  # Alignment statistcs
  finals <- alignmentFeatureStatistics(df, force = FALSE, BPPARAM = BPPARAM)
  # Do trimming detection
  finals <- trim_detection(df, finals)
  # Save file
  write.csv(finals, file = pasteDir(stats_folder, "STATS.csv"))
  # Get plots
  QCplots(df, "mrna", stats_folder, plot.ext = plot.ext,
          complex.correlation.plots = complex.correlation.plots,
          force = FALSE, BPPARAM = BPPARAM)

  message("--------------------------")
  message(paste("Everything done, saved QC to:", stats_folder))
  return(invisible(NULL))
}

# Keep for legacy purpose for now

#' @inherit QCreport
#' @export
ORFikQC <- QCreport

#' Correlation and coverage plots for ORFikQC
#'
#' Correlation plots default to mRNA covering reads.
#' Meta plots defaults to leader, cds, trailer.\cr
#' Output will be stored in same folder as the
#' libraries in df.\cr
#' Correlation plots will be fpkm correlation and
#' log2(fpkm + 1) correlation between samples.
#'
#' Is part of \code{\link{QCreport}}
#' @inheritParams QCreport
#' @inheritParams outputLibs
#' @param region a character (default: mrna), make raw count matrices of
#' whole mrnas or one of (leaders, cds, trailers)
#' @param stats_folder directory to save, default:
#' \code{QCfolder(df)}
#' @return invisible(NULL) (objects stored to disc)
#' @family QC report
#' @importFrom AnnotationDbi metadata
#' @keywords internal
QCplots <- function(df, region = "mrna",
                    stats_folder = QCfolder(df),
                    plot.ext = ".pdf",
                    complex.correlation.plots = TRUE,
                    force = TRUE,
                    BPPARAM) {
  message("--------------------------")
  message("Making QC plots:")
  message("- Annotation to NGS libraries plot:")
  QCstats.plot(stats_folder, stats_folder, plot.ext = plot.ext)

  correlation.plots(df, stats_folder, region, plot.ext = plot.ext,
                    complex.correlation.plots = complex.correlation.plots)
  message("- PCA outlier plot:")
  pcaExperiment(df, stats_folder, plot.ext = plot.ext)
  # window coverage over mRNA regions
  message("- Meta coverage plots")
  txdb <- loadTxdb(df)
  txNames <- filterTranscripts(txdb, 100, 100, 100, longestPerGene = TRUE,
                               stopOnEmpty = FALSE)
  if (length(txNames) == 0) { # No valid tx to plot
    warning("No 5' UTRs or 3' of significant length defined, UTR metacoverage plots",
    " can not be made, check your annotation file. In case no UTRs exist in your annotation,
    you can add pseudo UTRs, to also see coverage profiles over those areas.")
    # Check if CDS exists
    txNames <- filterTranscripts(txdb, 0, 100, 0, longestPerGene = TRUE,
                                 stopOnEmpty = FALSE)
    if (length(txNames) == 0) {
      warnings("No CDS of length 100 detected, skipping meta coverage completely!")
      return(invisible(NULL))
    }
    message("  - Metacoverage of CDS region only")
    transcriptWindow(GRangesList(), loadRegion(txdb, "cds", txNames),
                     GRangesList(), df = df, outdir = stats_folder,
                     scores = c("sum", "transcriptNormalized"),
                     is.sorted = TRUE, windowSize = 100, plot.ext = plot.ext,
                     verbose = FALSE, force = force,
                     BPPARAM = BPPARAM)
    return(invisible(NULL))
  }


  loadRegions(txdb, parts = c("leaders", "cds", "trailers"),
              names.keep = txNames)
  # Plot seperated by leader, cds & trailer
  message("  - seperated into 5' UTR, CDS and 3' UTR regions")
  transcriptWindow(leaders, get("cds", mode = "S4"),
                   trailers, df = df, outdir = stats_folder,
                   scores = c("sum", "transcriptNormalized"),
                   is.sorted = TRUE, plot.ext = plot.ext,
                   verbose = FALSE, force = force,
                   BPPARAM = BPPARAM)
  return(invisible(NULL))
}
