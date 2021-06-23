#' A post Alignment quality control of reads
#'
#' The ORFik QC uses the aligned files (usually bam files),
#' fastp and STAR log files
#' combined with annotation to create relevant statistics.\cr\cr
#' This report consists of several steps:\cr
#' 1. Convert bam file / Input files to ".ofst" format, if not already done.
#' This format is around 400x faster to use in R than the bam format.\cr
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
#' @param out.dir optional output directory, default:
#' \code{dirname(df$filepath[1])}.
#' Will make a folder called "QC_STATS" with all results in this directory.
#' @param plot.ext character, default: ".pdf". Alternatives: ".png" or ".jpg".
#' @return invisible(NULL) (objects are stored to disc)
#' @family QC report
#' @importFrom utils write.csv
#' @export
#' @examples
#' # Load an experiment
#' df <- ORFik.template.experiment()
#' # Run QC
#' # QCreport(df)
QCreport <- function(df, out.dir = dirname(df$filepath[1]),
                     plot.ext = ".pdf", BPPARAM = bpparam()) {
  # When experiment is ready, everything down from here is automatic
  message("Started ORFik QC report:")
  validateExperiments(df)
  stats_folder <- pasteDir(out.dir, "/QC_STATS/")
  if (!dir.create(stats_folder, recursive = TRUE)) {
    if (!dir.exists(stats_folder)) stop("Could not create output directory!")
  }
  message("- Converting input files to .ofst")
  convertLibs(df, reassign.when.saving = TRUE)
  message("- Creating read length tables:")
  dt_read_lengths <- readLengthTable(df, output.dir = stats_folder)
  # Get count tables
  finals <- QC_count_tables(df, out.dir, BPPARAM)
  # Do trimming detection
  finals <- trim_detection(df, finals, out.dir)
  # Save file
  write.csv(finals, file = pasteDir(stats_folder, "STATS.csv"))
  # Get plots
  QCplots(df, "mrna", stats_folder, plot.ext = plot.ext, BPPARAM = BPPARAM)

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
#' @param region a character (default: mrna), make raw count matrices of
#' whole mrnas or one of (leaders, cds, trailers)
#' @param stats_folder directory to save, default:
#' paste0(dirname(df$filepath[1]), "/QC_STATS/")
#' @return invisible(NULL) (objects stored to disc)
#' @family QC report
#' @importFrom GGally ggpairs
#' @importFrom AnnotationDbi metadata
QCplots <- function(df, region = "mrna",
                    stats_folder = paste0(dirname(df$filepath[1]),
                                          "/QC_STATS/"),
                    plot.ext = ".pdf", BPPARAM) {
  message("--------------------------")
  message("Making QC plots:")
  message("- Annotation to NGS libraries plot:")
  QCstats.plot(df, stats_folder, plot.ext = plot.ext)

  correlation.plots(df, stats_folder, region, plot.ext = plot.ext)
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
                     scores = c("sum", "zscore", "transcriptNormalized"),
                     is.sorted = TRUE, BPPARAM = BPPARAM, windowSize = 100)
    return(invisible(NULL))
  }


  loadRegions(txdb, parts = c("leaders", "cds", "trailers"),
              names.keep = txNames)
  # Plot seperated by leader, cds & trailer
  message("  - seperated into 5' UTR, CDS and 3' UTR regions")
  transcriptWindow(leaders, get("cds", mode = "S4"),
                   trailers, df = df, outdir = stats_folder,
                   scores = c("sum", "zscore", "transcriptNormalized"),
                   is.sorted = TRUE, BPPARAM = BPPARAM)
  # Plot all transcripts as 1 region
  # TODO: Make this safe enough to include for 32GB computers
  # message("  - whole transcripts")
  # transcriptWindow1(df = df, outdir = stats_folder,
  #                   scores = c("sum", "zscore", "transcriptNormalized"),
  #                   BPPARAM = BPPARAM)
  return(invisible(NULL))
}
