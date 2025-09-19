#' Make a count matrix from a library or experiment
#'
#' Make a summerizedExperiment / matrix object from bam files or
#' other library formats sepcified by lib.type argument. Works
#' like HTSeq, to give you count tables per library.
#'
#' If txdb or gtf path is added, it is a rangedSummerizedExperiment
#' NOTE: If the file called saveName exists, it will then load file,
#' not remake it!\cr
#' There are different ways of counting hits on transcripts, ORFik does
#' it as pure coverage (if a single read aligns to a region with 2 genes, both
#' gets a count of 1 from that read).
#' This is the safest way to avoid false negatives
#' (genes with no assigned hits that actually have true hits).
#' @inheritParams QCreport
#' @inheritParams outputLibs
#' @param saveName a character (default NULL),
#' if set save experiment to path given. Always saved as .rds.,
#' it is optional to add .rds, it will be added for you if not present.
#' Also used to load existing file with that name.
#' @param longestPerGene a logical (default FALSE), if FALSE all transcript
#' isoforms per gene. Ignored if "region" is not a character of either:
#' "mRNA","tx", "cds", "leaders" or "trailers".
#' @param geneOrTxNames a character vector (default "tx"), should row names
#' keep trancript names ("tx") or change to gene names ("gene")
#' @param region a character vector (default: "mrna"), make raw count matrices
#' of whole mrnas or one of (leaders, cds, trailers).
#' Can also be a \code{\link{GRangesList}}, then it uses this region directly.
#' Can then be uORFs or a subset of CDS etc.
#' @param type default: "count" (raw counts matrix), alternative is "fpkm",
#' "log2fpkm" or "log10fpkm"
#' @param lib.type a character(default: "default"), load files in experiment
#' or some precomputed variant, either "ofst", "pshifted" or "cov"
#' These are made with ORFik:::convertLibs() or shiftFootprintsByExperiment().
#' Can also be custom user made folders inside the experiments bam folder.
#' Format "cov" (i.e. covRle format) is by far the fastest to use if existing.
#' @param weight numeric or character, a column to score overlaps by. Default "score",
#' will check for a metacolumn called "score" in libraries. If not found,
#' will not use weights.
#' @param forceRemake logical, default FALSE. If TRUE, will not look for existing file count table files.
#' @param libraries The call to output libraries, the input is not used! Default:
#' outputLibs(df, chrStyle = seqinfo(df), type = lib.type, force = force,
#'  library.names = library.names, BPPARAM = BPPARAM)
#' @param format character, default "qs", alternative: "rds". Which format to save summarizedExperiment.
#' @param BPPARAM how many cores/threads to use? default: BiocParallel::SerialParam()
#' @import SummarizedExperiment
#' @export
#' @return a \code{\link{SummarizedExperiment}} object or data.table if
#' "type" is not "count, with rownames as transcript / gene names.
#' @examples
#' ##Make experiment
#' df <- ORFik.template.experiment()
#' # makeSummarizedExperimentFromBam(df)
#' ## Only cds (coding sequences):
#' # makeSummarizedExperimentFromBam(df, region = "cds")
#' ## FPKM instead of raw counts on whole mrna regions
#' # makeSummarizedExperimentFromBam(df, type = "fpkm")
#' ## Make count tables of pshifted libraries over uORFs
#' uorfs <- GRangesList(uorf1 = GRanges("chr23", 17599129:17599156, "-"))
#' #saveName <- file.path(dirname(df$filepath[1]), "uORFs", "countTable_uORFs")
#' #makeSummarizedExperimentFromBam(df, saveName, region = uorfs)
#' ## To load the uORFs later
#' # countTable(df, region = "uORFs", count.folder = "uORFs")
makeSummarizedExperimentFromBam <- function(df, saveName = NULL,
                                            longestPerGene = FALSE,
                                            geneOrTxNames = "tx",
                                            region = "mrna", type = "count",
                                            lib.type = "ofst",
                                            weight = "score", forceRemake = FALSE,
                                            force = TRUE, library.names = bamVarName(df),
                                            libraries = outputLibs(df, chrStyle = seqinfo(df),
                                                                   paths = filepath(df, lib.type, suffix_stem = c("", "_pshifted")),
                                                                   type = lib.type, force = force,
                                                                   library.names = library.names,
                                                                   BPPARAM = BPPARAM),
                                            format = "qs",
                                            BPPARAM = BiocParallel::SerialParam()) {

  if(!is.null(saveName)) {
    if (file_ext(saveName) != format) saveName <- paste0(saveName,".", format)
    if (file.exists(saveName) & !forceRemake) {
      message("Loading existing count table, set forceRemake=TRUE if you want to remake")
      return(read_RDSQS(saveName))
    }
  }
  validateExperiments(df, library.names)
  tx <- loadRegionCustom(df, region, geneOrTxNames, longestPerGene)


  rawCounts <- data.table(matrix(0, ncol = length(library.names),
                                 nrow = length(tx)))
  colnames(rawCounts) <- library.names
  message("    - Counting overlaps")
  envir <- envExp(df)
  force(libraries)
  for (lib in library.names) { # For each sample
    message(lib)
    if (is.character(weight) & length(weight) == 1) {
      has_no_weights <- !(weight %in% colnames(mcols(get_lib_from_env(lib, envir))))
      if (has_no_weights) weight <- NULL
    }
    rawCounts[, (paste0(lib)) := countOverlapsW(tx, get_lib_from_env(lib, envir), weight = weight)]
  }
  res <- SummarizedExperimentExp(df, rawCounts, tx, library.names)

  if (type %in% c("fpkm", "log2fpkm", "log10fpkm")) {
    res <- as.data.table(scoreSummarizedExperiment(res, score = type))
    rownames(res) <- names(tx)
  }
  if(!is.null(saveName)) {
    if (!dir.exists(dirname(saveName))) {
      dir.create(dirname(saveName), showWarnings = FALSE, recursive = TRUE)
    }
    if (file_ext(saveName) != format) saveName <- paste0(saveName, ".", format)
    save_RDSQS(res, file = saveName)
  }
  return(res)
}

SummarizedExperimentExp <- function(df, rawCounts, rowRanges, library.names = bamVarName(df)) {
  mat <- as.matrix(rawCounts);colnames(mat) <- NULL
  return(SummarizedExperiment(assays=list(counts=mat), rowRanges=rowRanges,
                              colData=colDataFromExp(df, library.names)))
}

colDataFromExp <- function(df, library.names = bamVarName(df)) {
  colData <- DataFrame(SAMPLE = as.factor(bamVarName(df, TRUE)),
                       row.names=library.names)
  # Add sample columns
  if (!is.null(df$rep)) colData$replicate <- as.factor(df$rep)
  if (!is.null(df$stage)) colData$stage <- as.factor(df$stage)
  if (!is.null(df$libtype)) colData$libtype <- as.factor(df$libtype)
  if (!is.null(df$condition)) colData$condition <- as.factor(df$condition)
  if (!is.null(df$fraction)) colData$fraction <- as.factor(as.character(df$fraction))
  return(colData)
}

get_lib_from_env <- function(lib_variable_name, envir = .GlobalEnv) {
  return(get(lib_variable_name, envir))
}

loadRegionCustom <- function(df, region, geneOrTxNames = "tx", longestPerGene= TRUE) {
  stopifnot(length(geneOrTxNames) == 1)
  stopifnot(geneOrTxNames %in% c("tx", "gene"))
  if (is(region, "character")) {
    txdb <- loadTxdb(df)
    if (longestPerGene) {
      longestTxNames <- filterTranscripts(txdb, 0, 0, 0, longestPerGene = TRUE)
      tx <- loadRegion(txdb, region, names.keep = longestTxNames)
    } else tx <- loadRegion(txdb, region)
  } else tx <- region

  if (geneOrTxNames == "gene") {
    if (!is(region, "character")) txdb <- loadTxdb(df)
    names(tx) <- txNamesToGeneNames(names(tx), txdb)
  }
  return(tx)
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
    } else { # all
      collapsedAll <- collapseReplicates(final, rep("merged_all",
                                                    ncol(final)))
      nlibs <- ncol(final)
    }
    # Number of samples per group as matrix

    assay(collapsedAll) <- ceiling(assay(collapsedAll) / nlibs)
  } else collapsedAll <- final

  only.one.group <- (length(unique(collapsedAll$SAMPLE)) == 1) |
    (ncol(collapsedAll) == 1)
  if ((collapse == "all") | only.one.group) {
    dds <- DESeqDataSet(collapsedAll, design = ~ 1)
  } else {
    dds <- DESeqDataSet(collapsedAll, design = ~ SAMPLE)
  }

  if (score %in% c("transcriptNormalized", "fpkm", "log2fpkm", "log10fpkm")) {
    fpkmCollapsed <- DESeq2::fpkm(dds, robust = FALSE)
    fpkmCollapsed[is.nan(fpkmCollapsed)] <- 0
    if (score == "transcriptNormalized") {
      normalization <- matrix(rep(rowSums2(fpkmCollapsed),
                                  ncol(fpkmCollapsed)),
                              ncol = ncol(fpkmCollapsed))
      fpkmTranscriptNormalized <- fpkmCollapsed / normalization
      assay(dds) <- fpkmTranscriptNormalized
      return(dds)
    } else if (score == "fpkm") {
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
#' Used to quickly load pre-created read count tables to R. \cr
#' If df is experiment:
#' Extracts by getting /QC_STATS directory, and searching for region
#' Requires \code{\link{ORFikQC}} to have been run on experiment,
#' to get default count tables!
#'
#' If df is path to folder:
#' Loads the the file in that directory with the regex region.rds,
#' where region is what is defined by argument, if multiple exist,
#' see if any start with "countTable_", if so, subset. If loaded as SummarizedExperiment
#' or deseq, the colData will be made from ORFik.experiment information.
#' @param df an ORFik \code{\link{experiment}} or path to folder with
#' countTables, use path if not same folder as experiment libraries. Will subset to
#' the count tables specified if df is experiment. If experiment has 4 rows and you subset it
#' to only 2, then only those 2 count tables will be outputted.
#' @param region a character vector (default: "mrna"), make raw count matrices
#'  of whole mrnas or one of (leaders, cds, trailers).
#' @param type character, default: "count" (raw counts matrix).
#' Which object type and normalization do you want ?
#' "summarized" (SummarizedExperiment object),
#' "deseq" (Deseq2 experiment, design will be all valid non-unique
#' columns except replicates, change by using DESeq2::design,
#' normalization alternatives are: "fpkm", "log2fpkm" or "log10fpkm".
#' @param collapse a logical/character (default FALSE), if TRUE all samples
#' within the group SAMPLE will be collapsed to one. If "all", all
#' groups will be merged into 1 column called merged_all. Collapse is defined
#' as rowSum(elements_per_group) / ncol(elements_per_group)
#' @param count.folder character, default "auto" (Use count tables from
#' original bam files stored in "QC_STATS", these are like HTseq count tables).
#' To load your custome count tables from pshifted reads, set to "pshifted"
#' (remember to create the pshifted tables first!). If you
#' have custom ranges, like reads over uORFs stored in a folder called
#' "/uORFs" relative to the bam files, set to "uORFs". Always create these
#' custom count tables with \code{\link{makeSummarizedExperimentFromBam}}.
#' Always make the location of the folder directly
#' inside the bam file directory!
#' @param full_path Full path to countTable, default: countTablePath(df, region, count.folder)
#' @return a data.table/SummarizedExperiment/DESeq object
#' of columns as counts / normalized counts per library, column name
#' is name of library. Rownames must be unique for now. Might change.
#' @importFrom DESeq2 DESeqDataSet
#' @importFrom stats as.formula
#' @export
#' @family countTable
#' @examples
#' # Make experiment
#' df <- ORFik.template.experiment()
#' # Make QC report to get counts ++ (not needed for this template)
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
                       collapse = FALSE,
                       count.folder = "default",
                       full_path = countTablePath(df, region, count.folder)) {
  df.temp <- attr(full_path, "experiment")
  if (length(full_path) == 1) {
    res <- read_RDSQS(full_path)

    # Subset to samples wanted
    if (!is.null(df.temp)) {
      if ((ncol(res) != nrow(df.temp))) {
        res <- subset_count_table(res, df.temp)
      }
    }

    is_ribo <- any(c("RFP", "RPF", "LSU","80S") %in% colData(res)$libtype, na.rm = TRUE)
    if(count.folder != "pshifted" & is_ribo)
      message("Loading default 80S counts, update count.folder to pshifted if wanted?")
    if (type == "count") return(as.data.table(assay(res)))

    res <- metadata_count_table(res, df.temp, type)
    # Give important sanity check info:

    # Decide output format
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
  } else if (length(full_path) > 1) {
    message(paste("More than 1 count table: ", df, collapse = ", "))
    stop("Folder contains multiple count tables for the same region, ORFik does not
         know which to pick. Delete or move the one that is not supposed to be there!")
  }
}

countTablePath <- function(df, region = "mrna", count.folder = "default") {
  # TODO fix bug if deseq!
  full_path <- experiment <- NULL
  if (is(df, "experiment")) {
    if (nrow(df) == 0) stop("df experiment has 0 rows (samples)!")
    experiment <- df
    df <-
      if (count.folder == "default") {
        QCfolder(df)
      } else paste0(libFolder(df), "/", count.folder)
  }
  if (is(df, "character")) {
    if (dir.exists(df)) {
      full_path <- list.files(path = df, pattern = paste0(region, "\\.qs$"),
                       full.names = TRUE)
      if (length(full_path) == 0) {
        full_path <- list.files(path = df, pattern = paste0(region, "\\.rds$"),
                         full.names = TRUE)
      }
      if (length(full_path) > 1) {
        hits <- grep("^countTable_", basename(full_path))
        if (length(hits) == 1) {
          full_path <- full_path[hits]
        }
      }
    }
  }

  if (length(full_path) == 0) {
    message(paste("Invalid count table directory:", df))
    stop("Table not found!",
         " Must be either: filepath to directory with defined countTable of region, the full path
       to the countTable, run ORFikQC to get default countTables!")
  }

  attr(full_path, "experiment") <- experiment
  return(full_path)
}

#' Make a list of count matrices from experiment
#'
#' By default will make count tables over mRNA, leaders, cds and trailers for
#' all libraries in experiment. Saved as "qs" or "rds" format files.
#'
#' @inheritParams makeSummarizedExperimentFromBam
#' @inheritParams QCreport
#' @param regions a character vector, default:
#'  c("mrna", "leaders", "cds", "trailers"), make raw count matrices
#' of whole regions specified. Can also be a custom GRangesList of
#' for example uORFs or a subset of cds etc.
#' @param rel.dir relative output directory for out.dir, default:
#' "QC_STATS". For pshifted, write "pshifted".
#' @param path_prefix the prefix names of tables, default:
#' if (!is.null(out.dir) {pasteDir(file.path(out.dir, rel.dir, "countTable_"))} else NULL,
#' i.e. directory + countTable_ or NULL if out.dir is NULL.
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a list of data.table, 1 data.table per region. The regions
#' will be the names the list elements.
#' @family countTable
#' @export
#' @examples
#' ##Make experiment
#' df <- ORFik.template.experiment()
#' ## Create count tables for all default regions
#' countTable_regions(df, NULL)
#' ## Pshifted reads (first create pshiftead libs)
#' # countTable_regions(df, lib.type = "pshifted", rel.dir = "pshifted")
countTable_regions <- function(df, out.dir = libFolder(df),
                               longestPerGene = FALSE,
                               geneOrTxNames = "tx",
                               regions = c("mrna", "leaders", "cds",
                                           "trailers"),
                               type = "count", lib.type = "ofst",
                               weight = "score",
                               rel.dir = "QC_STATS", forceRemake = FALSE,
                               library.names = bamVarName(df),
                               format = "qs",
                               path_prefix = if (!is.null(out.dir)) {pasteDir(file.path(out.dir, rel.dir, "countTable_"))} else {NULL},
                               libraries = outputLibs(df, chrStyle = seqinfo(df),
                                                      paths = filepath(df, lib.type, suffix_stem = c("", "_pshifted")),
                                                      type = lib.type, force = FALSE,
                                                      library.names = library.names,
                                                      BPPARAM = BiocParallel::SerialParam()),
                               BPPARAM = bpparam()) {

  libs <- bplapply(
    regions,
    function(region, path, df, geneOrTxNames, longestPerGene, forceRemake,
             library.names, libraries) {
     message("- Creating read count tables for region:")
     message("  - ", region)
     if (!is.null(path)) path <- paste0(path, region)
     makeSummarizedExperimentFromBam(df, region = region,
                                     geneOrTxNames = geneOrTxNames,
                                     longestPerGene = longestPerGene,
                                     saveName = path, lib.type = lib.type,
                                     library.names = library.names,
                                     forceRemake = forceRemake, force = FALSE,
                                     libraries = libraries,
                                     format = format)
    },
    path = path_prefix, df = df,
    geneOrTxNames = geneOrTxNames, library.names = library.names,
    longestPerGene = longestPerGene, libraries = libraries,
    forceRemake = forceRemake, BPPARAM = BPPARAM
  )
  names(libs) <- regions
  return(libs)
}

subset_count_table <- function(res, df.temp) {
  subset <- df.temp$index
  if (length(subset) > ncol(res)) {
    # Fall back to old method, useful if you made some edits and forgot to
    # update count table
    subset <- if (sum(colnames(res) %in% bamVarName(df.temp, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, skip.experiment = FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, skip.experiment = FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, TRUE, FALSE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, TRUE, FALSE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, FALSE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, FALSE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, FALSE, FALSE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, FALSE, FALSE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, TRUE, FALSE, TRUE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, TRUE, FALSE, TRUE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, TRUE, TRUE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, TRUE, TRUE, FALSE)
    } else if (sum(colnames(res) %in%
                   bamVarName(df.temp, FALSE, FALSE, TRUE, FALSE)) == nrow(df.temp)) {
      colnames(res) %in% bamVarName(df.temp, FALSE, FALSE, TRUE, FALSE)
    } else stop("No valid names for count tables found from experiment")
  }
  return(res[, subset])
}

metadata_count_table <- function(res, df.temp, type) {
  # Add all sample columns if not existing and it is possible
  if (is.null(colData(res)$stage)) {
    if (length(unique(df.temp$stage)) > 1) {
      colData(res)$stage <- as.factor(df.temp$stage)
      if (anyNA(colData(res)$stage))
        colData(res)$stage[is.na(colData(res)$stage)] <- ""
    }
    if ((length(unique(df.temp$libtype)) > 0) & (type != "deseq")) {
      colData(res)$libtype <- as.factor(df.temp$libtype)
      if (anyNA(colData(res)$libtype))
        colData(res)$libtype[is.na(colData(res)$libtype)] <- ""
    }
    if (length(unique(df.temp$condition)) > 1) {
      colData(res)$condition <- as.factor(df.temp$condition)
      if (anyNA(colData(res)$condition))
        colData(res)$condition[is.na(colData(res)$condition)] <- ""
    }
    if (length(unique(df.temp$fraction)) > 1) {
      colData(res)$fraction <- as.factor(as.character(df.temp$fraction))
      if (anyNA(colData(res)$fraction))
        colData(res)$fraction[is.na(colData(res)$fraction)] <- ""
    }
    colData(res)$replicate <- as.factor(df.temp$rep)
  }
  return(res)
}
