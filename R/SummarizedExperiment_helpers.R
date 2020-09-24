#' Make a count matrix from a library or experiment
#'
#' Make a summerizedExperiment / matrix object from bam files
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' NOTE: If the file called saveName exists, it will then load file,
#' not remake it!
#' @param df an ORFik \code{\link{experiment}}
#' @param saveName a character (default NULL),
#' if set save experiment to path given. Always saved as .rds.,
#' it is optional to add .rds, it will be added for you if not present.
#' Also used to load existing file with that name.
#' @param longestPerGene a logical (default TRUE), if FALSE all transcript
#' isoforms per gene.
#' @param geneOrTxNames a character vector (default "tx"), should row names
#' keep trancript names ("tx") or change to gene names ("gene")
#' @param region a character vector (default: "mrna"), make raw count matrices
#' of whole mrnas or one of (leaders, cds, trailers).
#' Can also be a \code{\link{GRangesList}}, then it uses this region directly.
#' @param type default: "count" (raw counts matrix), alternative is "fpkm",
#' "log2fpkm" or "log10fpkm"
#' @param weight numeric or character, a column to score overlaps by. Default "score",
#' will check for a metacolumn called "score" in libraries. If not found,
#' will not use weights.
#' @import SummarizedExperiment
#' @export
#' @return a \code{\link{SummarizedExperiment}} object or data.table if
#' "type" is not "count, with rownames as transcript / gene names.
#' @examples
#' ##Make experiment
#' df <- ORFik.template.experiment()
#' # makeSummarizedExperimentFromBam(df)
#' # Only cds (coding sequences):
#' # makeSummarizedExperimentFromBam(df, region = "cds")
#' # FPKM instead of raw counts on whole mrna regions
#' # makeSummarizedExperimentFromBam(df, type = "fpkm")
makeSummarizedExperimentFromBam <- function(df, saveName = NULL,
                                            longestPerGene = TRUE,
                                            geneOrTxNames = "tx",
                                            region = "mrna", type = "count",
                                            weight = "score") {
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
  for (i in seq(length(varNames))) { # For each sample
    print(varNames[i])
    if (is.character(weight) & length(weight) == 1) {
      if (!(weight %in% colnames(mcols(get(varNames[i])))))
        weight <- NULL
    }

    co <- countOverlapsW(tx, get(varNames[i]), weight = weight)
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
    if (!dir.exists(dirname(saveName))) {
      dir.create(dirname(saveName), showWarnings = FALSE, recursive = TRUE)
    }
    if (file_ext(saveName) != "rds") saveName <- paste0(saveName,".rds")
    saveRDS(res, file = saveName)
  }
  return(res)
}

#' Helper function for makeSummarizedExperimentFromBam
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' For FPKM values, DESeq2::fpkm(robust = FALSE) is used
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
      nlibs <- t(matrix(as.double(table(colData(final)$SAMPLE)),
                        ncol = nrow(assay(collapsedAll)) ,
                        nrow = length(unique(colData(final)$SAMPLE))))
    } else {
      collapsedAll <- collapseReplicates(final, rep("merged_all",
                                                    ncol(final)))
      nlibs <- ncol(final)
    }
    # Number of samples per group as matrix

    assay(collapsedAll) <- ceiling(assay(collapsedAll) / nlibs)
  } else collapsedAll <- final
  if (collapse == "all") {
    dds <- DESeqDataSet(collapsedAll, design = ~ 1)
  } else {
    dds <- DESeqDataSet(collapsedAll, design = ~ SAMPLE)
  }

  if (score %in% c("transcriptNormalized", "fpkm", "log2fpkm", "log10fpkm")) {
    fpkmCollapsed <- DESeq2::fpkm(dds, robust = FALSE)
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
  }
  return(dds)
}

#' Extract count table directly from experiment
#'
#' Used to quickly load read count tables to R. \cr
#' If df is experiment:
#' Extracts by getting /QC_STATS directory, and searching for region
#' Requires \code{\link{ORFikQC}} to have been run on experiment!
#'
#' If df is path to folder:
#' Loads the the file in that directory with the regex region.rds,
#' where region is what is defined by argument.
#' @param df an ORFik \code{\link{experiment}} or path to folder with
#' countTable, use path if not same folder as experiment libraries.
#' @param region a character vector (default: "mrna"), make raw count matrices
#'  of whole mrnas or one of (leaders, cds, trailers).
#' @param type default: "count" (raw counts matrix),
#' "summarized" (SummarizedExperiment object),
#' "deseq" (Deseq2 experiment, design will be all valid non-unique
#' columns except replicates, change by using DESeq2::design,
#' normalization alternatives are: "fpkm", "log2fpkm" or "log10fpkm"
#' @param collapse a logical/character (default FALSE), if TRUE all samples
#' within the group SAMPLE will be collapsed to one. If "all", all
#' groups will be merged into 1 column called merged_all. Collapse is defined
#' as rowSum(elements_per_group) / ncol(elements_per_group)
#' @return a data.table of columns as counts per library, column name
#' is name of library. Rownames must be unique for now. Might change.
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom stats as.formula
#' @export
#' @examples
#' # Make experiment
#' ORFik.template.experiment()
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
#' # Get count Table of mrnas as summarizedExperiment
#' # countTable(df, "mrna", type = "summarized")
#' # Get count Table of mrnas as DESeq2 object,
#' # for differential expression analysis
#' # countTable(df, "mrna", type = "deseq")
countTable <- function(df, region = "mrna", type = "count",
                       collapse = FALSE) {
  # TODO fix bug if deseq!
  if (is(df, "experiment")) {
    dir = dirname(df$filepath[1])
    df <- paste0(dir, "/QC_STATS")
  }
  if (is(df, "character")) {
    if (dir.exists(df)) {
      df <- list.files(path = df, pattern = paste0(region, ".rds"),
                       full.names = TRUE)
    }
    if (length(df) == 1) {
      res <- readRDS(df)
      if (type == "summarized") return(res)
      if (type == "deseq") {
        # remove replicate from formula
        formula <- colnames(colData(res))
        if ("replicate" %in% formula)
          formula <- formula[-grep("replicate", formula)]
        formula <- as.formula(paste(c("~", paste(formula,
                                    collapse = " + ")), collapse = " "))
        return(DESeqDataSet(res, design = formula))
      }
      ress <- scoreSummarizedExperiment(res, type, collapse)
      if (is(ress, "matrix")) {
        ress <- as.data.table(ress)
      } else { # is deseq
        ress <- as.data.table(assay(ress))
      }
      rownames(ress) <- names(ranges(res))
      return(ress)
    }
  }
  message(paste("Invalid count table:", df))
  stop("df must be filepath to directory with countTable, the path
       to the countTable or ORFik experiment with a QC_STATS folder!")
}

#' Make a list of count matrices from experiment
#'
#' @inheritParams makeSummarizedExperimentFromBam
#' @inheritParams QCreport
#' @param regions a character vector, default:
#'  c("mrna", "leaders", "cds", "trailers"), make raw count matrices
#' of whole regions specified.
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a list of data.table, 1 data.table per region. The regions
#' will be the names the list elements.
countTable_regions <- function(df, out.dir = dirname(df$filepath[1]),
                               longestPerGene = TRUE,
                               geneOrTxNames = "tx",
                               regions = c("mrna", "leaders", "cds",
                                           "trailers"),
                               type = "count",
                               weight = "score",
                               BPPARAM = bpparam()) {

  stats_folder <- pasteDir(out.dir, "/QC_STATS/")
  countDir <- paste0(stats_folder, "countTable_")
  libs <- bplapply(
    regions,
    function(region, countDir, df, geneOrTxNames, longestPerGene) {
     message("Making count tables for region:")
     message(region)
     path <- paste0(countDir, region)
     makeSummarizedExperimentFromBam(df, region = region,
                                     geneOrTxNames = "tx",
                                     longestPerGene = FALSE,
                                     saveName = path)
    },
    countDir = countDir, df = df,
    geneOrTxNames = geneOrTxNames,
    longestPerGene = longestPerGene, BPPARAM = BPPARAM
  )
  names(libs) <- regions
  return(libs)
}
