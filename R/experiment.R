#' experiment class definition
#'
#' It is an object to massivly simplify your coding, by having a
#' table of all libraries of an experiment. That contains
#' filepaths and info for each library in the experiment. It also tries
#' to guess grouping / types / pairs by the file names.
#' \cr
#' Act as a way of extension of \code{\link{SummarizedExperiment}} by allowing
#' more ease to find not only counts, but rather
#' information about libraries, and annotation, so that more tasks are
#' possible. Like coverage per position in some transcript etc.
#'
#'
#' ## Constructor:
#'
#' Simplest way to make is to call:
#' create.experiment(dir)
#'
#' On some folder with libraries and see what you get. Some of the fields
#' might be needed to fill in manually. Each resulting row must be unique
#' (excluding filepath), that means
#' if it has replicates then that must be said explicit. And all
#' filepaths must be unique and have files with size > 0.
#' Syntax (columns):\cr
#' libtype (library type): rna-seq, ribo-seq, CAGE etc.\cr
#' rep (replicate): 1,2,3 etc\cr
#' condition: WT (wild-type), control, target, mzdicer, starved etc.\cr
#' fraction: 18, 19 (fractinations), or other ways to split library.\cr
#' filepath: Full filepath to file
#' @details
#' Special rules:\cr
#' Supported:\cr
#' Single end bam, bed, wig + compressions of these\cr
#' Paired forward / reverse wig files, must have same name except
#'  _forward / _reverse etc
#'
#' Not supported yet:\cr
#' Paired end bam files not supported yet!
#' @importFrom methods new
#' @examples
#' \dontrun{
#' library(ORFik)
#'
#' # 1. Update path to experiment data  directory (bam, bed, wig files etc)
#' exp_dir = "/data/processed_data/RNA-seq/Lee_zebrafish_2013/aligned/"
#'
#' # 2. Set a 5 character name for experiment, (Lee 2013 -> Lee13, etc)
#' exper_name = "Lee13"
#'
#' # 3. Create a template experiment
#' temp <- create.experiment(exp_dir, exper_name,
#'  txdb = "/data/references/Zv9_zebrafish/Danio_rerio.Zv9.79.gtf",
#'  fa = "/data/references/Zv9_zebrafish/Danio_rerio.Zv9.fa")
#'
#' # 4. Make sure each row(sample) is unique and correct
#' # You will get a view open now, check the data.frame that it is correct:
#' # library type (RNA-seq, Ribo-seq), stage, rep, condition, fraction.
#' # Let say it did not figure out it is RNA-seq, then we do:"
#'
#' temp[5:6, 1] <- "RNA" # [row 5 and 6, col 1] are library types
#'
#' # You can also do this in your spread sheet program (excel, libre..)
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
#' @export
#' @family ORFik_experiment
experiment <- setClass("experiment",
                       slots=list(experiment = "character",
                                  txdb = "character",
                                  fafile = "character",
                                  expInVarName = "logical"),
                       contains = "DataFrame")

#' experiment show definition
#'
#' Show a simplified version of experiment.
#' The show function simplifies the view so that any
#' column of data (like replicate or stage) is not shown, if all
#' values are identical in that column. Filepath is also never shown.
#' @param object an ORFik \code{\link{experiment}}
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

            obj <- as.data.table(as(object@listData, Class = "DataFrame"))
            if (nrow(obj) > 0) {
              obj <- obj[,-"filepath"]
              if (!is.null(obj$reverse)) obj <- obj[,-"reverse"]
              skip <- c()
              for (i in 2:ncol(obj)) {
                if (nrow(unique(obj[,i, with = F])) == 1)
                  skip <- c(skip, i)
              }
              if (length(skip) > 0) {
                show(obj[,-skip, with = F])
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
            nrow(as.data.table(as(x@listData, Class = "DataFrame")))
          }
)

#' Read ORFik \code{\link{experiment}}
#'
#' Read in runs / samples from an experiment as a single R object.
#' To read an ORFik experiment, you must of course make one first.
#' See \code{\link{create.experiment}}
#' The file must be csv and be a valid ORFik experiment
#' @param file a .csv file following ORFik experiment style ("," as seperator)
#' , or a template data.frame from \code{\link{create.experiment}}
#' @return an ORFik \code{\link{experiment}}
#' @export
#' @examples
#' # From file
#' \dontrun{
#' # Read from file
#' df <- read.experiment(filepath) # <- valid .csv file
#' }
#' # Read from (create.experiment() template)
#' # 1. Pick directory
#' dir <- system.file("extdata", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               viewTemplate = FALSE)
#' template$X5[6] <- "heart" # <- fix non unique row
#' # read experiment from template
#' df <- read.experiment(template)
#'
#' # To save it, do:
#' # save.experiment(df, file = "path/to/save/experiment.csv")
#' @family ORFik_experiment
read.experiment <-  function(file) {
  if (is(file, "character")) {
    if (file_ext(file) == "") file <- paste0(file, ".csv")
    info <- read.table(file, sep = ",", nrows = 3, stringsAsFactors = FALSE)
    listData <- read.csv2(file, skip = 3, header = TRUE, sep = ",",
                          stringsAsFactors = FALSE)
  } else if(is(file, "data.frame")) {
    info <- file[seq(3),]
    listData <- file[-seq(4),]
    colnames(listData) <- file[4,]
  } else stop("file must be either character or data.frame template")


  exper <- info[1, 2]
  txdb <- ifelse(is.na(info[2, 2]),  "", info[2, 2])
  fa <- ifelse(is.na(info[3, 2]),  "", info[3, 2])

  df <- experiment(experiment = exper, txdb = txdb, fafile = fa,
                   listData = listData, expInVarName = TRUE)

  validateExperiments(df)
  return(df)
}

#' Create a template for new ORFik \code{\link{experiment}}
#'
#' Create information on runs / samples from an experiment as a single R object.
#' By using files in a folder / folders. It will try to make an experiment table
#' with information per sample. There will be several columns you can fill in,
#' most of there it will try to auto-detect. Like if it is RNA-seq or Ribo-seq,
#' Wild type or mutant etc.
#' You will have to fill in the details that were not autodetected.
#' Easiest way to fill in the blanks are in a csv editor like libre Office
#' or excel. Remember that each row (sample) must have a unique combination
#' of values.
#' An extra column called "reverse" is made if there are paired data,
#' like +/- strand wig files.
#' @param dir Which directory / directories to create experiment from
#' @param exper Short name of experiment, max 5 characters long
#' @param saveDir Directory to save experiment csv file (NULL)
#' @param types Default (bam, bed, wig), which types of libraries to allow
#' @param txdb A path to gff/gtf file used for libraries
#' @param fa A path to fasta genome/sequences used for libraries
#' @param viewTemplate run View() on template when finished, default (TRUE)
#' @return a data.frame, NOTE: this is not a ORFik experiment,
#'  only a template for it!
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
#' # Save with: save.experiment(df, file = "path/to/save/experiment.csv")
#' @family ORFik_experiment
create.experiment <- function(dir, exper, saveDir = NULL,
                              types = c("bam", "bed", "wig"), txdb = "",
                              fa = "", viewTemplate = TRUE) {
  notDir <- !all(dir.exists(dir))
  if (notDir) stop(paste(dir[!dir.exists(dir)], "is not a valid directory!"))
  file_dt <- findLibrariesInFolder(dir, types)
  if (is(file_dt, "data.table")) { # If paired data
    files <- file_dt$forward
    df <- data.frame(matrix(ncol = 7, nrow = length(files) + 4))
    df[4,] <- c("libtype", "stage", "rep", "condition", "fraction","filepath",
                "reverse")
    # set lib column names
    df[5:(5+length(files)-1), 7] <- file_dt$reverse
  } else { # only single libraries
    files <- file_dt
    df <- data.frame(matrix(ncol = 6, nrow = length(files) + 4))
    # set lib column names
    df[4,] <- c("libtype", "stage", "rep", "condition", "fraction","filepath")
  }

  # set file paths
  df[5:(5+length(files)-1), 6] <- files
  # Set library type (RNA-seq etc)
  df[5:(5+length(files)-1), 1] <- findFromPath(files)
  # set stage
  stages <- c("unfertalized", "_fertalized",
              "2-4cell", "2to4cell", "64cell", "256cell", "512cell","1Kcell",
              "2-4Cell", "2to4Cell","64Cell", "256Cell", "512Cell","1KCell",
              "2-4_cell", "2to4_cell", "64_cell", "256_cell", "512_cell",
              "sphere", "shield", "dome", "oblong", "bud",
              "Sphere", "Shield", "Dome", "Oblong", "Bud",
              "_2h", "_4h", "_6h", "_8h", "_12h", "_24h", "_28h", "_48h",
              "_02h", "_04h", "_06h", "_08h",
              "2hpf", "4hpf", "6hpf", "8hpf", "12hpf", "24hpf", "28hpf",
              "48hpf",
              "1dpf", "2dpf", "3dpf", "4dpf", "5dpf", "6dpf", "10dpf",
              "21dpf", "24dpf")
  cell_lines <- c("HEK293", "HeLa", "THP-1", "PC3")
  tissues <- c("adipose", "brain", "bladder", "blood", "breast", "colon",
               "cortex", "eye", "fibroblast", "frontal lobe", "heart",
               "kidney", "liver", "lung", "muscle","ovary", "prostate",
               "rectum", "testis","urunary", "vagina", "skin", "tongue")
  stages <- c(stages, tissues, cell_lines)
  df[5:(5+length(files)-1), 2] <- findFromPath(files, stages)
  # set rep
  df[5:(5+length(files)-1), 3] <- findFromPath(files, c("rep1", "rep2", "rep3",
                                                       "Rep1", "Rep2", "Rep3",
                                                       "run1", "run2", "run3",
                                                       "_r1_", "_r2_", "_r3_"))
  # Set condition
  conditions <- c("WT", "control", "MZ", "dicer", "4Ei", "4ei", "silvesterol",
                  "Silvesterol", "mutant", "Mutant", "cas9", "Cas9")
  df[5:(5+length(files)-1), 4] <- findFromPath(files, conditions)

  df[1, seq(2)] <- c("name", exper)
  df[2, seq(2)] <- c("gff", txdb)
  df[3, seq(2)] <- c("fasta", fa)
  df[is.na(df)] <- ""
  if (!is.null(saveDir))
    save.experiment(df, pasteDir(saveDir, exper,".csv"))
  if (viewTemplate) View(df)
  return(df)
}

#' Save \code{\link{experiment}} to disc
#'
#'
#' @param df an ORFik \code{\link{experiment}}
#' @param file name of file to save df as
#' @export
#' @return NULL (experiment save only)
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
#' # Save with: save.experiment(df, file = "path/to/save/experiment.csv")
#' @family ORFik_experiment
save.experiment <- function(df, file) {
  write.table(x = df, file = file, sep = ",",
              row.names = FALSE, col.names = FALSE)
  return(invisible(NULL))
}

#' Find all candidate library types filenames
#'
#' From the given \code{\link{experiment}}
#' @param filepaths path to all files
#' @param candidates Possible names to search for.
#' @return a candidate library types (character vector)
findFromPath <- function(filepaths, candidates = c("RNA", "rna-seq",
                                                   "Rna-seq", "RNA-seq",
                                                   "RFP", "RPF", "ribo-seq",
                                                   "Ribo-seq", "mrna",
                                                   "mrna-seq", "mRNA-seq",
                                                   "CAGE", "cage", "LSU",
                                                   "SSU", "ATAC", "tRNA",
                                                   "SHAPE", "PRPF")) {
  types <- c()
  for (path in filepaths) {
    hit <- unlist(sapply(candidates, grep, x = path))
    hitRel <- unlist(sapply(candidates, grep, x = gsub(".*/", "", path)))
    type <- if(length(hit) == 1 & length(hitRel) == 0) names(hit)
    over <- hit[names(hit) %in% names(hitRel)]
    type <- ifelse(length(over) == 1, names(over),
                   ifelse(is.null(type), "", type))
    types <- c(types, gsub(pattern = "_", "", type))
  }
  return(types)
}


#' Which type of library type in \code{\link{experiment}}?
#' @param df an ORFik \code{\link{experiment}}
#' @return library types (character vector)
#' @family ORFik_experiment
libraryTypes <- function(df) {
  if (is(df, "experiment")) {
    return(unique(df$libtype))
  } else if (is(df, "character") | is(df, "factor")) {
    return(gsub("_.*", x = df, replacement = ""))
  } else stop("library types must be data.frame or character vector!")
}

#' Validate ORFik \code{\link{experiment}}
#'
#' Check for valid existing, non-empty and all unique.
#' A good way to see if your experiment is valid.
#' @param df an ORFik \code{\link{experiment}}
#' @return NULL (Stops if failed)
#' @family ORFik_experiment
validateExperiments <- function(df) {
  libTypes <- libraryTypes(df)
  if (!is(df, "experiment")) stop("df must be experiment!")
  if (!all((c("stage", "libtype") %in% colnames(df))))
    stop("stage, libtype and experiment must be colnames in df!")
  if (length(libTypes) == 0) stop("df have no valid sequencing libraries!")
  if (nrow(df) == 0) stop("df must have at least 1 row!")
  if (!is.null(df$reverse)) {
    files <- c(df$filepath, df$reverse)
  } else files <- df$filepath

  emptyFiles <- c()
  for (i in files) {
    emptyFiles <- c(emptyFiles, as.numeric(sapply(as.character(i),
                                                  file.size)) == 0)
  }
  if (any(is.na(emptyFiles)))
    stop(paste("File is not existing:\n", files[is.na(emptyFiles)]))
  if (any(emptyFiles)) {
    print(files[emptyFiles])
    stop("Empty files in list, see above for which")
  }
  if (length(bamVarName(df)) != length(unique(bamVarName(df))))
    stop("experiment table has non-unique rows!")
  if (length(files) != length(unique(files)))
    stop("Duplicated filepaths in experiment!")
}

#' Get library variable names from ORFik \code{\link{experiment}}
#'
#' What will each sample be called given the columns of the experiment?
#' @param df an ORFik \code{\link{experiment}}
#' @param skip.replicate a logical (FALSE), don't include replicate
#' in variable name.
#' @param skip.condition a logical (FALSE), don't include condition
#' in variable name.
#' @param skip.stage a logical (FALSE), don't include stage
#' in variable name.
#' @param skip.fraction a logical (FALSE), don't include fraction
#' @param skip.experiment a logical (FALSE), don't include experiment
#' @param skip.libtype a logical (FALSE), don't include libtype
#' @return variable names of libraries (character vector)
#' @export
#' @family ORFik_experiment
bamVarName <- function(df, skip.replicate = length(unique(df$rep)) == 1,
                       skip.condition = length(unique(df$condition)) == 1,
                       skip.stage = length(unique(df$stage)) == 1,
                       skip.fraction = length(unique(df$fraction)) == 1,
                       skip.experiment = !df@expInVarName,
                       skip.libtype = FALSE) {

  varName <- c()
  for (i in 1:nrow(df)) {
    varName <- c(varName, bamVarNamePicker(df[i,], skip.replicate,
                                           skip.condition, skip.stage,
                                           skip.fraction, skip.experiment,
                                           skip.libtype))
  }
  return(varName)
}

#' Get variable name per filepath in experiment
#'
#' @param df an ORFik \code{\link{experiment}}
#' @param skip.replicate a logical (FALSE), don't include replicate
#' in variable name.
#' @param skip.condition a logical (FALSE), don't include condition
#' in variable name.
#' @param skip.stage a logical (FALSE), don't include stage
#' in variable name.
#' @param skip.fraction a logical (FALSE), don't include fraction
#' @param skip.experiment a logical (FALSE), don't include experiment
#' @param skip.libtype a logical (FALSE), don't include libtype
#' @return variable name of library (character vector)
bamVarNamePicker <- function(df, skip.replicate = FALSE,
                             skip.condition = FALSE,
                             skip.stage = FALSE, skip.fraction = FALSE,
                             skip.experiment = FALSE, skip.libtype = FALSE) {
  if(nrow(df) != 1) stop("experiment must only input 1 row")
  lib <- df$libtype
  stage <- df$stage
  cond <- df$condition
  rep <- df$rep
  frac <- df$fraction
  current <- ""
  # Add only underscore if x is not ""
  spaste <- function(x, y, reverse = FALSE) {
    if (reverse)
      return(paste(x, y, sep = ifelse(y == "", "", "_")))
    return(paste(x, y, sep = ifelse(x == "", "", "_")))
  }
  if (!skip.libtype)
    current <- lib
  if(!skip.condition)
    current <- spaste(current, cond)
  if (!skip.stage)
    current <- spaste(current, stage)
  if (!(skip.fraction | is.null(frac) | is.na(frac))) {
    if (frac != "")
      current <- spaste(current, paste0("f", frac))
  }

  if (!(skip.replicate | is.null(rep)))
    current <- spaste(current, paste0("r", rep))
  if (! (skip.experiment | is.null(df@experiment)))
    current <- spaste(df@experiment, current, TRUE)

  current <- gsub(pattern = "__", "_", current)
  return(gsub("_$", "", current))
}

#' Output bam/bed/bedo/wig files to R as variables
#'
#' Variable names defined by df (ORFik experiment DataFrame)
#' Uses multiple cores to load, defined by multicoreParam
#' @param df an ORFik \code{\link{experiment}}
#' @inheritParams matchSeqStyle
#' @param type a character(default: "defualt"), load files in experiment
#' or some precomputed variant, either "bedo" or "pshifted".
#' These are made with ORFik:::simpleLibs()
#' @param envir environment to save to, default (.GlobalEnv)
#' @param BPPARAM how many cores? default: bpparam()
#' @return NULL (libraries set by envir assignment)
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
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
#' # Output to .GlobalEnv with:
#' # outputLibs(df)
#' @family ORFik_experiment
outputLibs <- function(df, chrStyle = NULL, type = "default",
                       envir = .GlobalEnv, BPPARAM = bpparam()) {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)

  for (df in dfl) {
    validateExperiments(df)
    libTypes <- libraryTypes(df)
    varNames <- bamVarName(df)
    loaded <- c()
    for (i in 1:nrow(df)) { # For each stage
      if (exists(x = varNames[i], envir = envir, inherits = FALSE)) {
        loaded <- append(loaded, TRUE)
      } else loaded <- append(loaded, FALSE)
    }
    # Par apply
    if (!all(loaded)) {
      message(paste0("Ouputing libraries from: ",df@experiment))
      libs <- bplapply(df$filepath,
          function(x, chrStyle, df, type) {
            i <- which(df$filepath == x)
            varNames <- bamVarName(df)
            message(paste(i, ": ", varNames[i]))
            if (type == "bedo") {
              out.dir <- paste0(dirname(df$filepath[1]), "/bedo/")
              if (dir.exists(out.dir)) {
                input <- paste0(out.dir,
                                remove.file_ext(x,
                                                basename = TRUE)
                                , ".bedo")
              } else type <- "default"
            } else if (type == "pshifted") {
              out.dir <- paste0(dirname(df$filepath[1]), "/pshifted/")
              if (dir.exists(out.dir)) {
                input <- paste0(out.dir,
                                remove.file_ext(x,
                                                basename = TRUE)
                                , "_pshifted.bedo")
                if (!file.exists(input)) {
                  input <- paste0(out.dir,
                                  remove.file_ext(x,
                                                  basename = TRUE)
                                  , "_pshifted.bed")
                  if (!file.exists(input))
                    type <- "default"
                  }
              } else type <- "default"
            }
            if (type == "default") {
              if (!is.null(df$reverse)) {
                if (df[i,]$reverse != "")
                input <- c(x, df[i,]$reverse)
              } else input <- x
            }
            return(fimport(input, chrStyle))
          },
          BPPARAM = BPPARAM, chrStyle = chrStyle,
          df = df, type = type
      )

      # assign to environment
      for (i in 1:nrow(df)) { # For each stage
        assign(varNames[i], libs[[i]], envir = envir)
      }
    }
  }
  return(invisible(NULL))
}

#' Will make a simplified version of NGS libraries
#'
#' Export files as .bedo files. It is a bed file with 2 score columns.
#' Gives a massive speedup when cigar strings are not needed.
#'
#' See \code{\link{export.bedo}} for information on file format
#' @param df an ORFik \code{\link{experiment}}
#' @param out.dir optional output directory, default:
#' paste0(dirname(df$filepath[1]), "/bedo/"),
#' if it is NULL, it will just reassign R objects to simplified libraries.
#' @inheritParams convertToOneBasedRanges
#' @param must.overlap default (NULL), else a GRanges / GRangesList object, so
#' only reads that overlap (must.overlap) are kept. This is useful when you
#' only need the reads over transcript annotation or subset etc.
#' @return NULL (saves files to disc or R .GlobalEnv)
simpleLibs <- function(df,
                       out.dir = paste0(dirname(df$filepath[1]), "/bedo/"),
                       addScoreColumn = TRUE, addSizeColumn = TRUE,
                       must.overlap = NULL) {
  validateExperiments(df)
  if (!is.null(must.overlap) & !is.gr_or_grl(must.overlap))
    stop("must.overlap must be GRanges or GRangesList object!")
  if (!is.null(out.dir)) {
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(out.dir)) stop("could not create directory!")
    message(paste("Saving .bedo files to:", out.dir))
  }

  outputLibs(df, chrStyle = must.overlap)

  varNames <- bamVarName(df)
  i <- 1
  for (f in varNames) {
    message(f)
    gr <- convertToOneBasedRanges(gr = get(f),
                                  addScoreColumn = addScoreColumn,
                                  addSizeColumn = addSizeColumn,
                                  method = "None")
    if (!is.null(must.overlap)) gr <- optimizeReads(must.overlap, gr)

    if (!is.null(out.dir)) {
      output <- paste0(out.dir,
                       remove.file_ext(df$filepath[i], basename = TRUE),
                       ".bedo")
      export.bedo(gr, output)
    } else {
      assign(x = f, value = gr, envir = .GlobalEnv)
    }
    i <- i + 1
  }
}


#' Remove bam/bed/wig files load in R as variables
#'
#' Variable names defined by df, in envir defined
#' @param df an ORFik \code{\link{experiment}}
#' @param envir environment to save to, default (.GlobalEnv)
#' @return NULL (objects removed from envir specified)
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
#' # Output to .GlobalEnv with:
#' # outputLibs(df)
#' # Then remove them with: remove.experiments(df)
remove.experiments <- function(df, envir = .GlobalEnv) {
  rm(list =  bamVarName(df), envir = envir)
  message(paste0("Removed loaded libraries from experiment:",
                 df@experiment))
  return(invisible(NULL))
}

#' Get all library files in folder/folders of given types
#'
#' Will try to guess paired / unpaired wig, bed, bam files.
#' @param dir The directory/directories to find bam, bed, wig files.
#' @param types All accepted types of bam, bed, wig files..
#' @importFrom tools file_ext
#' @return (data.table) All files found from types parameter.
#' With 2 extra column (logical), is it wig pairs, and paired bam files.
findLibrariesInFolder <- function(dir, types) {
  regex <- paste("\\.", types, collapse = "|", sep = "")
  # Find files in multiple dirs in correct order
  files <- unlist(lapply(dir,
                         FUN = function(d)
                           grep(pattern = regex,
                                x = list.files(d, full.names = TRUE),
                           value = TRUE)))
  files <- pasteDir(files)
  # Remove .bai bam index files etc
  fext <- file_ext(files)
  if (!(all(fext %in% types))) {
    files <- files[fext != "bai"]
    fext <- fext[fext != "bai"]

    compressed = fext %in% c("gzip", "gz", "bgz", "zip")
    if (any(compressed)) {
      fext[compressed] <-file_ext(file_path_sans_ext(files[compressed],
                                                      compression = FALSE))
    }
    files <- files[fext %in% types]
  }
  filesOld <- files
  # Wig pairs
  wig_files <- findNGSPairs(files[fext == "wig"])
  # Bed pairs
  bed_pairs <- findNGSPairs(filesOld[fext == "bed"], format = "bed")
  # Paired end bam
  bam_pairs <- findNGSPairs(filesOld[fext == "bam"],
                                    f = c("_R1_00", "_F", "_Forward", "_forward"),
                                    r = c("_R2_00", "_R", "_Reverse", "_reverse"),
                                    format = "bam")
  pairs <- data.table()
  if (is(wig_files, "data.table")) pairs <- rbind(pairs, wig_files)
  if (is(bed_pairs, "data.table")) pairs <- rbind(pairs, bed_pairs)
  if (is(bam_pairs, "data.table")) pairs <- rbind(pairs, bam_pairs)

  if (nrow(pairs) > 0) {
    if (nrow(pairs)*2 != length(files)) { # if more than just matched pairs
      others <- files[!(files %in% c(pairs$forward, pairs$reverse))]
      file_dt <- data.table(forward = others, reverse = "", match = FALSE)
      files <- rbind(pairs, file_dt)
    } else files <- pairs
    if (nrow(files) == 0) stop("Found no valid files in folder")
  } else {
    if (length(files) == 0) stop("Found no valid files in folder")
  }
  return(files)
}
