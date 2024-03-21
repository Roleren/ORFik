#' experiment class definition
#'
#' It is an object that simplify and error correct your NGS workflow,
#' creating a single R object that stores and controls all results relevant
#' to a specific experiment.\cr It contains following important parts:
#' \itemize{
#'  \item{filepaths : }{and info for each library in the experiment
#'  (for multiple files formats: bam, bed, wig, ofst, ..)}
#'  \item{genome : }{annotation files of the experiment
#'  (fasta genome, index, gtf, txdb)}
#'  \item{organism : }{name (for automatic GO, sequence analysis..)}
#'  \item{description : }{and author information (list.experiments(),
#'  show all experiments you have made with ORFik, easy to find and load them later)}
#'  \item{API : }{ORFik supports a rich API for using the experiment,
#'  like outputLibs(experiment, type = "wig") will load all libraries
#'  converted to wig format into R,  loadTxdb(experiment) will load the
#'  txdb (gtf) of experiment, transcriptWindow() will automatically plot
#'  metacoverage of all libraries in the experiment, countTable(experiment) will load count tables, etc..)}
#'  \item{Safety : }{It is also a safety in that it verifies your experiments contain
#'   no duplicate, empty or non-accessible files.}
#' }
#' Act as a way of extension of \code{\link{SummarizedExperiment}} by allowing
#' more ease to find not only counts, but rather
#' information about libraries, and annotation, so that more tasks are
#' possible. Like coverage per position in some transcript etc.\cr\cr
#' ## Constructor:\cr
#' Simplest way to make is to call:\cr
#' create.experiment(dir)\cr
#' On some folder with NGS libraries (usually bam files) and see what you get.
#' Some of the fields
#' might be needed to fill in manually. Each resulting row must be unique
#' (not including filepath, they are always unique), that means
#' if it has replicates then that must be said explicit. And all
#' filepaths must be unique and have files with size > 0.\cr\cr
#' Here all the columns in the experiment will be described:
#' name (column info): examples\cr
#' \describe{
#'      \item{libtype}{library type: rna-seq, ribo-seq, CAGE etc}
#'      \item{stage}{stage or tissue: 64cell, Shield, HEK293}
#'      \item{rep}{replicate: 1,2,3 etc}
#'      \item{condition}{treatment or condition: :
#'      WT (wild-type), control, target, mzdicer, starved}
#'      \item{fraction}{fraction of total: 18, 19 (TCP / RCP fractions),
#'      or other ways to split library.\cr}
#'      \item{filepath}{Full filepath to file}
#'      \item{reverse}{optional: 2nd filepath or info, only used if paired files}
#' }
#'
#' @details
#' Special rules:\cr
#' Supported:\cr
#' Single/paired end bam, bed, wig, ofst + compressions of these\cr
#' The reverse column of the experiments says "paired-end" if bam file.
#' If a pair of wig files, forward and reverse strand, reverse is filepath
#' to '-' strand wig file.
#' Paired forward / reverse wig files, must have same name except
#'  _forward / _reverse in name\cr
#' Paired end bam, when creating experiment, set pairedEndBam = c(T, T, T, F).
#' For 3 paired end libraries, then one single end.\cr
#' Naming:
#' Will try to guess naming for tissues / stages, replicates etc.
#' If it finds more than one hit for one file, it will not guess.
#' Always check that it guessed correctly.
#' @importFrom methods new
#' @examples
#' ## To see an internal ORFik example
#' df <- ORFik.template.experiment()
#' ## See libraries in experiment
#' df
#' ## See organism of experiment
#' organism(df)
#' ## See file paths in experiment
#' filepath(df, "default")
#' ## Output NGS libraries in R, to .GlobalEnv
#' #outputLibs(df)
#' ## Output cds of experiment annotation
#' #loadRegion(df, "cds")
#'
#' ## This is how to make it:
#' \dontrun{
#' library(ORFik)
#'
#' # 1. Update path to experiment data  directory (bam, bed, wig files etc)
#' exp_dir = "/data/processed_data/RNA-seq/Lee_zebrafish_2013/aligned/"
#'
#' # 2. Set a short character name for experiment, (Lee et al 2013 -> Lee13, etc)
#' exper_name = "Lee13"
#'
#' # 3. Create a template experiment (gtf and fasta genome)
#' temp <- create.experiment(exp_dir, exper_name, saveDir = NULL,
#'  txdb = "/data/references/Zv9_zebrafish/Danio_rerio.Zv9.79.gtf",
#'  fa = "/data/references/Zv9_zebrafish/Danio_rerio.Zv9.fa",
#'  organism = "Homo sapiens")
#'
#' # 4. Make sure each row(sample) is unique and correct
#' # You will get a view open now, check the data.frame that it is correct:
#' # library type (RNA-seq, Ribo-seq), stage, rep, condition, fraction.
#' # Let say it did not figure out it is RNA-seq, then we do:"
#'
#' temp[5:6, 1] <- "RNA" # [row 5 and 6, col 1] are library types
#'
#' # You can also do this in your spread sheet program (excel, libre office)
#' # Now save new version, if you did not use spread sheet.
#' saveName <- paste0("/data/processed_data/experiment_tables_for_R/",
#'  exper_name,".csv")
#' save.experiment(temp, saveName)
#'
#' # 5. Load experiment, this will validate that you actually made it correct
#' df <- read.experiment(saveName)
#'
#' # Set experiment name not to be assigned in R variable names
#' df@expInVarName <- FALSE
#' df
#' }
#' @return a ORFik experiment
#' @export
#' @family ORFik_experiment
experiment <- setClass("experiment",
                       slots=list(experiment = "character",
                                  txdb = "character",
                                  fafile = "character",
                                  organism = "character",
                                  assembly = "character",
                                  author = "character",
                                  expInVarName = "logical",
                                  envir = "environment",
                                  resultFolder = "character"),
                       contains = "DFrame")

#' experiment show definition
#'
#' Show a simplified version of the experiment.
#' The show function simplifies the view so that any
#' column of data (like replicate or stage) is not shown, if all
#' values are identical in that column. Filepaths are also never shown.
#' @param object an ORFik \code{\link{experiment}}
#' @importFrom withr local_options
#' @export
#' @return print state of experiment
setMethod("show",
          "experiment",
          function(object) {
            type <- ifelse(length(unique(object@listData$libtype)) == 1,
                           "type", "types")
            cat("experiment:", object@experiment, "with",
                length(unique(object@listData$libtype)), "library", type, "and",
                length(object@listData$libtype), "runs","\n")
            if (object@author != "") cat(object@author, "et al. \n")

            obj <- as.data.table(as(object@listData, Class = "DataFrame"))
            withr::local_options(list(datatable.print.class = FALSE))
            if (nrow(obj) > 0) {
              obj <- obj[,-c("filepath", "index")]
              if ("Run" %in% colnames(obj)) obj <- obj[,-c("Run")]
              if (!is.null(obj$reverse)) obj <- obj[,-"reverse"]
              skip <- c()
              for (i in 2:ncol(obj)) {
                if (nrow(unique(obj[,i, with = FALSE])) == 1)
                  skip <- c(skip, i)
              }
              if (length(skip) > 0) {
                show(obj[,-skip, with = FALSE])
              } else show(obj)
            }
          }
)

#' Internal nrow function for ORFik experiment
#' Number of runs in experiment
#' @param x an ORFik \code{\link{experiment}}
#' @return number of rows in experiment (integer)
setMethod("nrow",
          "experiment",
          function(x) {
            length(x$libtype)
          }
)

#' Get name of ORFik experiment
#' @param x an ORFik \code{\link{experiment}}
#' @return character, name of experiment
#' @export
setGeneric("name", function(x) standardGeneric("name"))

#' @inherit name
setMethod("name",
          "experiment",
          function(x) {
            x@experiment
          }
)

#' Get ORFik experiment environment
#'
#' More correctly, get the pointer reference, default is .GlobalEnv
#' @param x an ORFik \code{\link{experiment}}
#' @return environment pointer, name of environment: pointer
#' @export
setGeneric("envExp", function(x) standardGeneric("envExp"))

#' Set ORFik experiment environment
#'
#' More correctly, set the pointer reference, default is .GlobalEnv
#' @param x an ORFik \code{\link{experiment}}
#' @param value environment pointer to assign to experiment
#' @return an ORFik \code{\link{experiment}} with updated environment
#' @export
setGeneric("envExp<-", function(x, value) standardGeneric("envExp<-"))

#' @inherit envExp
setMethod("envExp",
          "experiment",
          function(x) {
            x@envir
          }
)

#' @inherit envExp<-
setMethod("envExp<-",
          "experiment",
          function(x, value) {
            stopifnot(is.environment(value))
            x@envir <- value
            return(x)
          }
)

#' Get ORFik experiment organism
#'
#' If not defined directly, checks the txdb / gtf organism information, if existing.
#' @param object an ORFik \code{\link{experiment}}
#' @return character, name of organism
#' @family ORFik_experiment
#' @export
#' @examples
#' # if you have set organism in txdb of ORFik experiment:
#' df <- ORFik.template.experiment()
#' organism(df)
#'
#' #' If you have not set the organism you can do:
#' #gtf <- "pat/to/gff_or_gff"
#' #txdb_path <- paste0(gtf, ".db") # This file is created in next step
#' #txdb <- makeTxdbFromGenome(gtf, genome, organism = "Homo sapiens",
#' # optimize = TRUE, return = TRUE)
#' # then use this txdb in you ORFik experiment and load:
#' # create.experiment(exper = "new_experiment",
#' #   txdb = txdb_path) ...
#' # organism(read.experiment("new-experiment))
setMethod("organism",
          "experiment",
          function(object) {
            if (object@organism != "") return(object@organism)
            org <- BiocGenerics::organism(loadTxdb(object))
            if (is.na(org))
              message("Organism not set in either of experiment and txdb of gtf")
            return(org)
          }
)

#' Get ORFik experiment main output directory
#'
#' @param x an ORFik \code{\link{experiment}}
#' @return a character path
#' @export
setGeneric("resFolder", function(x) standardGeneric("resFolder"))

#' @inherit resFolder
setMethod("resFolder",
          "experiment",
          function(x) {
            if (is(try(x@resultFolder, silent = TRUE), "try-error"))
              return(libFolder(x))
            if (x@resultFolder != "") {
              return(x@resultFolder)
            } else return(libFolder(x))
          }
)


#' Get ORFik experiment QC folder path
#'
#' @param x an ORFik \code{\link{experiment}}
#' @return a character path
#' @export
setGeneric("QCfolder", function(x) standardGeneric("QCfolder"))

#' @inherit QCfolder
setMethod("QCfolder",
          "experiment",
          function(x) {
            file.path(resFolder(x), "QC_STATS/")
          }
)

#' Get ORFik experiment library folder
#'
#' @param x an ORFik \code{\link{experiment}}
#' @param mode character, default "first". Alternatives: "unique", "all".
#' @return a character path
#' @export
setGeneric("libFolder", function(x, mode = "first") standardGeneric("libFolder"))

#' @inherit libFolder
setMethod("libFolder",
          "experiment",
          function(x, mode = "first") {
            path <-
            if (mode == "first") {
              dirname(x$filepath[1])
            } else if (mode == "unique") {
              unique(dirname(x$filepath))
            } else if (mode == "all") {
              dirname(x$filepath)
            } else stop("argument 'mode', must be either first, unique or all")
            return(path)
          }
)

#' Get SRR/DRR/ERR run ids from ORFik experiment
#'
#' @param x an ORFik \code{\link{experiment}}
#' @return a character vector of runIDs, "" if not existing.
#' @export
setGeneric("runIDs", function(x) standardGeneric("runIDs"))

#' @inherit runIDs
setMethod("runIDs",
          "experiment",
          function(x) {
            if ("Run" %in% colnames(x)) return(x$Run)
            return(rep("", nrow(x)))
          }
)

#' Seqlevels ORFik experiment
#' Extracted from fasta genome index
#' @param x an ORFik \code{\link{experiment}}
#' @return integer vector with names
#' @export
setMethod("seqlevels",
          "experiment",
          function(x) {
            seqlevels(findFa(df@fafile))
          }
)

#' Seqinfo ORFik experiment
#' Extracted from fasta genome index
#' @param x an ORFik \code{\link{experiment}}
#' @return integer vector with names
#' @export
setMethod("seqinfo",
          "experiment",
          function(x) {
            seqinfo(findFa(x@fafile))
          }
)

#' Get ORFik experiment QC folder path
#'
#' @param x an ORFik \code{\link{experiment}}
#' @return a data.table with gene id, gene symbols and tx ids (3 columns)
#' @export
setGeneric("symbols", function(x) standardGeneric("symbols"))

#' @inherit QCfolder
setMethod("symbols",
          "experiment",
          function(x) {
            cand_path <- file.path(dirname(x@txdb), "gene_symbol_tx_table.fst")
            if (file.exists(cand_path)) {
              return(as.data.table(read_fst(cand_path)))
            } else {
              message("Gene symbols not created, run ",
              "ORFik:::makeTxdbFromGenome(gene_symbols = TRUE)")
              return(data.table())
            }
          }
)

#' Get experimental design
#' Find the column/columns that create a separation between samples,
#' by default skips replicate and choose first that is
#' from either: libtype, condition, stage and fraction.
#' @param object an ORFik \code{\link{experiment}}
#' @param batch.correction.design logical, default FALSE. If true,
#' add replicate as a second design factor (only if >= 2 replicates exists).
#' @param as.formula logical, default FALSE. If TRUE, return as formula
#' @param multi.factor logical, default TRUE If FALSE, return first factor only
#' (+ rep, if batch.correction.design is true). Order of picking is:
#' libtype, if not then: stage, if not then: condition, if not then: fraction.
#' @return a character (name of column) or a formula
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' design(df) # The 2 columns that decides the design here
#' # If we subset it changes
#' design(df[df$libtype == "RFP",])
#' # Only single factor design, it picks first
#' design(df, multi.factor = FALSE)
setMethod("design",
          "experiment",
          function(object, batch.correction.design = FALSE,
                   as.formula = FALSE, multi.factor = TRUE) {
            formula <- colnames(object)
            formula <- formula[formula %in% c("libtype", "stage", "rep", "condition", "fraction")]
            if ("rep" %in% formula & !batch.correction.design)
              formula <- formula[-grep("rep", formula)]
            dt <- as.data.table(object)[, ..formula]

            dt <- dt[, apply(dt, FUN = function(i)
                      return(if (length(unique(i)) == 1) {FALSE} else TRUE), MARGIN = 2),
                     with = FALSE]
            if (nrow(dt) == 0)
              stop("Malformed experiment, you need a column that seperates the libraries (> 1 unique value")
            factors <- colnames(dt)

            if (!multi.factor) {
              factors <- c(factors[!(factors %in% "rep")][1], factors[factors %in% "rep"][1])
              factors <- factors[!is.na(factors)]
            }
            if (as.formula) {
              return(as.formula(paste(c("~", paste(factors,
                                                   collapse = " + ")),
                                      collapse = " ")))
            } else return(factors)
          }
)

#' Get experiment design model matrix
#'
#' The function extends stats::model.matrix.
#' @param object an ORFik \code{\link{experiment}}
#' @param design_formula the experiment design, as formula, subset columns, to
#'  change the model.matrix, default: \code{design(object, as.formula = TRUE)}
#' @return a matrix with design and level attributes
#' @export
#' @importFrom stats model.matrix
#' @examples
#' df <- ORFik.template.experiment()
#' model.matrix(df)
setMethod("model.matrix",
          "experiment",
          function(object, design_formula = design(object, as.formula = TRUE)) {
            stats::model.matrix(design_formula, data = object)
})

