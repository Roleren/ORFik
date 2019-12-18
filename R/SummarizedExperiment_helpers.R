#' Make a count matrix from a libraries
#'
#' Make a summerizedExperiment / matrix object from bam files
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' NOTE: If the file called saveName exists, it will the load file,
#' not remake it!
#' @param df an ORFik \code{\link{experiment}}
#' @param saveName a character (default NULL),
#' if set save experiment to path given. Always saved as .rds.,
#' it is optional to add .rds, it will be added for you if not present.
#' @param longestPerGene a logical (default TRUE), if FALSE all transcript
#' isoforms per gene.
#' @param geneOrTxNames a character (default "gene"), should row names
#' keep trancripts names ("tx") or change to gene names ("gene")
#' @param region a character vector (default: "mrna"), make raw count matrices
#'  of whole mrnas or one of (leaders, cds, trailers). Can also be a
#' \code{\link{GRangesList}}, then it uses this region directly.
#' @param type default: "count" (raw counts matrix), alternative is "fpkm",
#' "log2fpkm" or "log10fpkm"
#' @import SummarizedExperiment
#' @export
#' @return a \code{\link{SummarizedExperiment}} object or data.table if
#' "type" is not "count, with rownames as transcript / gene names.
#' @examples
#' # 1. Pick directory
#' dir <- system.file("extdata", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               viewTemplate = FALSE)
#' template$X5[6] <- "heart" # <- fix non unique row
#' # read experiment
#' df <- read.experiment(template)
#'
#' # makeSummarizedExperimentFromBam(df)
#' # Only cds (coding sequences):
#' # makeSummarizedExperimentFromBam(df, region = "cds")
#' # FPKM instead of raw counts on whole mrna regions
#' # makeSummarizedExperimentFromBam(df, type = "fpkm")
makeSummarizedExperimentFromBam <- function(df, saveName = NULL,
                                            longestPerGene = TRUE,
                                            geneOrTxNames = "gene",
                                            region = "mrna", type = "count") {


  if(!is.null(saveName)) {
    if (file_ext(saveName) != "rds") saveName <- paste0(saveName,".rds")
    if (file.exists(saveName)) return(readRDS(saveName))
  }
  validateExperiments(df)

  if (is(region, "character")) {
    txdb <- loadTxdb(df@txdb)
    tx <- loadRegion(txdb, region)
  } else tx <- region

  if (geneOrTxNames == "gene") {
    if (!is(region, "character")) txdb <- loadTxdb(df@txdb)
    names(tx) <- txNamesToGeneNames(names(tx), txdb)
  }

  varNames <- bamVarName(df)
  outputLibs(df, tx)

  rawCounts <- data.table(matrix(0, ncol = length(varNames),
                                 nrow = length(tx)))
  for (i in 1:length(varNames)) { # For each sample
    print(varNames[i])
    co <- countOverlaps(tx, get(varNames[i]))
    rawCounts[, (paste0("V",i)) := co]
  }
  mat <- as.matrix(rawCounts);colnames(mat) <- NULL

  colData <- DataFrame(SAMPLE = bamVarName(df, TRUE),
                       row.names=varNames)
  if (!is.null(df$rep)) colData$replicate <- df$rep

  res <- SummarizedExperiment(assays=list(counts=mat), rowRanges=tx,
                              colData=colData)
  if (type %in% c("fpkm", "log2fpkm", "log10fpkm")) {
    res <- as.data.table(scoreSummarizedExperiment(res, score = type))
    rownames(res) <- names(tx)
  }
  if(!is.null(saveName)) {
    if (file_ext(saveName) != "rds") saveName <- paste0(saveName,".rds")
    saveRDS(res, file = saveName)
  }
  return(res)
}

#' Helper function for makeSummarizedExperimentFromBam
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' @param final ranged summarized experiment object
#' @param score default: "transcriptNormalized"
#' (row normalized raw counts matrix),
#' alternative is "fpkm", "log2fpkm" or "log10fpkm"
#' @param collapse a logical/character (default FALSE), if TRUE all samples
#' within the group SAMPLE will be collapsed to one. If "all", all
#' groups will be merged into 1 column called merged_all. Collapse is defined
#' as rowSum(elements_per_group) / ncol(elements_per_group)
#' @import SummarizedExperiment DESeq2
#' @export
#' @return a DEseq summerizedExperiment object (transcriptNormalized)
#'  or matrix (if fpkm input)
scoreSummarizedExperiment <- function(final, score = "transcriptNormalized",
                                      collapse = FALSE) {
  if (is.factor(final$SAMPLE)) {
    lvls <- levels(final$SAMPLE) %in% unique(colData(final)$SAMPLE)
    final$SAMPLE <- factor(final$SAMPLE, levels = levels(final$SAMPLE)[lvls])
  }
  if (collapse %in% c(TRUE, "all")) {
    if (collapse == TRUE) {
      collapsedAll <- collapseReplicates(final, final$SAMPLE)
    } else {
      collapsedAll <- collapseReplicates(final, rep("merged_all",
                                                    ncol(final)))
    }
    # Number of samples per group as matrix
    nlibs <- t(matrix(as.double(table(colData(final)$SAMPLE)),
                      ncol = nrow(assay(collapsedAll)) ,
                      nrow = length(unique(colData(final)$SAMPLE))))
    assay(collapsedAll) <- ceiling(assay(collapsedAll) / nlibs)
  } else collapsedAll <- final


  dds <- DESeqDataSet(collapsedAll, design = ~ SAMPLE)
  fpkmCollapsed <- DESeq2::fpkm(dds)
  if (score == "transcriptNormalized") {
    normalization <- matrix(rep(rowSums2(fpkmCollapsed),
                                ncol(fpkmCollapsed)),
                            ncol = ncol(fpkmCollapsed))
    fpkmTranscriptNormalized <- fpkmCollapsed / normalization
    assay(dds) <- fpkmTranscriptNormalized
    return(dds)
  }
  if (score == "fpkm") {
    return(fpkmCollapsed)
  } else if (score == "log2fpkm") {
    return(log2(fpkmCollapsed))
  } else if (score == "log10fpkm") {
    return(log10(fpkmCollapsed))
  }
  return(dds)
}

#' Extract count table directly from experiment
#'
#' Extracts by getting /QC_STATS directory, and searching for region
#' Requires \code{\link{ORFikQC}} to have been run on experiment!
#' @param df an ORFik \code{\link{experiment}} or path to folder with
#' countTable, use path if not same folder as experiment libraries.
#' @param region a character vector (default: "mrna"), make raw count matrices
#'  of whole mrnas or one of (leaders, cds, trailers). Can also be a
#' @param type default: "count" (raw counts matrix), alternative is "fpkm",
#' "log2fpkm" or "log10fpkm"
#' @param collapse a logical/character (default FALSE), if TRUE all samples
#' within the group SAMPLE will be collapsed to one. If "all", all
#' groups will be merged into 1 column called merged_all. Collapse is defined
#' as rowSum(elements_per_group) / ncol(elements_per_group)
#' @return a DEseq summerizedExperiment object (count)
#'  or matrix (if fpkm input)
#' @export
#' @examples
#' # 1. Pick directory
#' dir <- system.file("extdata", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               viewTemplate = FALSE)
#' template$X5[6] <- "heart" # <- fix non unique row
#' # read experiment
#' df <- read.experiment(template)
#' # Make QC report to get counts ++
#' # ORFikQC(df)
#'
#' # Get count Table of mrnas
#' # countTable(df, "mrna")
#' # Get count Table of cds
#' # countTable(df, "cds")
#' # Get count Table of mrnas as fpkm values
#' # countTable(df, "mrna", type = "count")
#' # Get count Table of mrnas with collapsed replicates
#' # countTable(df, "mrna", collapse = TRUE)
countTable <- function(df, region = "mrna", type = "count",
                       collapse = FALSE) {
  if (is(df, "experiment")) {
    dir = dirname(df$filepath[1])
    df <- paste0(dir, "/QC_STATS")
  }
  if (is(df, "character")) {
    if (dir.exists(df)) {
      df <- list.files(path = df, pattern = paste0(region, ".rds"),
                       full.names = T)
    }
    if (length(df) == 1) {
      res <- readRDS(df)
      ress <- as.data.table(scoreSummarizedExperiment(res, type, collapse))
      rownames(ress) <- names(ranges(res))

      return(ress)
    }
  }
  message(paste("Invalid count table:", df))
  stop("df must be filepath to dir, table or ORFik experiment!")
}
