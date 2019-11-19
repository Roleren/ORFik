#' Make a summerizedExperiment object from bam files
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' @param df an ORFik experiment, to make it, see: ?experiment
#' @param saveName a character (default NULL),
#' if set save experiment to path given. Always saved as .rds.,
#' it is optional to add .rds, it will be added for you if not present.
#' @param longestPerGene a logical (default TRUE), if FALSE all transcript
#' isoforms per gene.
#' @param geneOrTxNames a character (default gene), if tx use with tx names
#' @param region (default: mrna), make raw count matrices of
#' whole mrnas or one of (leaders, cds, trailers)
#' @param type default: "count" (raw counts matrix), alternative is "fpkm",
#' "log2fpkm" or "log10fpkm"
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @return a summerizedExperiment object
makeSummarizedExperimentFromBam <- function(df, saveName = NULL,
                                            longestPerGene = TRUE,
                                            geneOrTxNames = "gene",
                                            region = "mrna", type = "count") {
  if(!is.null(saveName) && file.exists(saveName)) {
    return(readRDS(saveName))
  }
  libTypes <- libraryTypes(df)
  validateExperiments(df)

  if (is(region, "character")) {
    txdb <- loadTxdb(df@txdb)
    tx <- loadRegion(txdb, region)
  } else tx <- region

  if (geneOrTxNames == "gene"){
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
  if (type %in% c("fpkm", "log2fpkm", "log10fpkm"))
    res <- as.data.table(scoreSummarizedExperiment(res, score = type))
  if(!is.null(saveName)) {
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
#' @param collapse a logical (default FALSE), if TRUE all samples
#' within group
#' @importFrom DESeq2 DESeqDataSet
#' @export
#' @return a DEseq summerizedExperiment object
scoreSummarizedExperiment <- function(final, score = "transcriptNormalized",
                                      collapse = FALSE) {
  if (is.factor(final$SAMPLE))
    lvls <- levels(final$SAMPLE) %in% unique(colData(final)$SAMPLE)
    final$SAMPLE <- factor(final$SAMPLE, levels = levels(final$SAMPLE)[lvls])

  if (collapse) {
    collapsedAll <- collapseReplicates(final, final$SAMPLE)
    assay(collapsedAll) <- ceiling(assay(collapsedAll) /
                              t(matrix(as.double(table(colData(final)$SAMPLE)),
                                       ncol = nrow(assay(collapsedAll)) ,
                                       nrow = length(unique(
                                         colData(final)$SAMPLE)))))
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
  } else if (score == "fpkm") {
    return(fpkmCollapsed)
  }
  else if (score == "fpkm") {
    return(fpkmCollapsed)
  }
  else if (score == "log2fpkm") {
    return(log2(fpkmCollapsed))
  }
  else if (score == "log10fpkm") {
    return(log10(fpkmCollapsed))
  }
  return(dds)
}
