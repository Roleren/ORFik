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
#'      \item{fraction}{fraction of total: 18, 19 (TCP / RCP fractinations),
#'      or other ways to split library.\cr}
#'      \item{filepath}{Full filepath to file}
#'      \item{reverse}{optional: only used if paired files,
#'      "paired-end" if bam file or filepath to 2nd file if wig files etc}
#' }
#'
#' @details
#' Special rules:\cr
#' Supported:\cr
#' Single/paired end bam, bed, wig, ofst + compressions of these\cr
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
#' organism.df(df)
#' ## See file paths in experiment
#' filepath(df, "default")
#' ## Output objects in R, to .GlobalEnv
#' #outputLibs(df)
#'
#' ## This is how to make it:
#' \dontrun{
#' library(ORFik)
#'
#' # 1. Update path to experiment data  directory (bam, bed, wig files etc)
#' exp_dir = "/data/processed_data/RNA-seq/Lee_zebrafish_2013/aligned/"
#'
#' # 2. Set a 5 character name for experiment, (Lee et al 2013 -> Lee13, etc)
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
#' @export
#' @family ORFik_experiment
experiment <- setClass("experiment",
                       slots=list(experiment = "character",
                                  txdb = "character",
                                  fafile = "character",
                                  organism = "character",
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
#' @param file relative path to a ORFik experiment. That is a .csv file following
#' ORFik experiment style ("," as seperator).
#' , or a template data.frame from \code{\link{create.experiment}}. Can
#' also be full path to file, then in.dir argument is ignored.
#' @param in.dir Directory to load experiment csv file from, default:
#' "~/Bio_data/ORFik_experiments/" \cr Set to NULL if you don't want to save
#' it to disc. Does not apply if file is not a path, but a data.frame. Also
#' does not apply if file was given as full path.
#' @return an ORFik \code{\link{experiment}}
#' @importFrom utils read.table read.csv2
#' @export
#' @examples
#' # From file
#' \dontrun{
#' # Read from file
#' df <- read.experiment(filepath) # <- valid ORFik .csv file
#' }
#' ## Read from (create.experiment() template)
#' df <- ORFik.template.experiment()
#'
#' ## To save it, do:
#' # save.experiment(df, file = "path/to/save/experiment")
#' ## You can then do:
#' # read.experiment("path/to/save/experiment")
#' # or (identical):
#' # read.experiment("experiment", in.dir = "path/to/save/")
#' @family ORFik_experiment
read.experiment <-  function(file, in.dir = "~/Bio_data/ORFik_experiments/") {
  if (is(file, "character")) {
    if (file_ext(file) == "") file <- paste0(file, ".csv")
    if (!file.exists(file)) file <- pasteDir(in.dir, file)
    if (!file.exists(file)) { # This will only trigger on CBU server @ UIB
      cbu.path <- "/export/valenfs/data/processed_data/experiment_tables_for_R/"
      if (file.exists(pasteDir(cbu.path, basename(file))))
        file <- pasteDir(cbu.path, basename(file))
    }


    info <- read.table(file, sep = ",", nrows = 3, stringsAsFactors = FALSE)
    listData <- read.csv2(file, skip = 3, header = TRUE, sep = ",",
                          stringsAsFactors = FALSE)
  } else if(is(file, "data.frame")) {
    info <- file[seq(3),]
    listData <- file[-seq(4),]
    colnames(listData) <- file[4,]
  } else stop("file must be either character or data.frame template")
  org <- ifelse(info[2,5] == "organism" & !is.na(info[2,6]),
                info[2,6], "")
  exper <- info[1, 2]
  txdb <- ifelse(is.na(info[2, 2]),  "", info[2, 2])
  fa <- ifelse(is.na(info[3, 2]),  "", info[3, 2])

  df <- experiment(experiment = exper, txdb = txdb, fafile = fa,
                   organism = org, listData = listData,
                   expInVarName = TRUE)

  validateExperiments(df)
  return(df)
}

#' Create a ORFik \code{\link{experiment}}
#'
#' Create information on runs / samples from an experiment as a single R object.
#' By using files in a folder / folders. It will try to make an experiment table
#' with information per sample. There will be several columns you can fill in,
#' most of there it will try to auto-detect. Like if it is RNA-seq or Ribo-seq,
#' Wild type or mutant etc.
#' You will have to fill in the details that were not auto detected.
#' Easiest way to fill in the blanks are in a csv editor like libre Office
#' or excel. Remember that each row (sample) must have a unique combination
#' of values.
#' An extra column called "reverse" is made if there are paired data,
#' like +/- strand wig files.
#' @param dir Which directory / directories to create experiment from
#' @param exper Short name of experiment, max 5 characters long
#' @param saveDir Directory to save experiment csv file, default:
#' "~/Bio_data/ORFik_experiments/" \cr Set to NULL if you don't want to save
#' it to disc.
#' @param types Default (bam, bed, wig), which types of libraries to allow
#' @param txdb A path to gff/gtf file used for libraries
#' @param fa A path to fasta genome/sequences used for libraries, remember the
#' file must have a fasta index too.
#' @param viewTemplate run View() on template when finished, default (TRUE)
#' @param organism character, default: "" (no organism set), scientific name
#' of organism. Homo sapiens, Danio rerio, Rattus norvegicus etc.
#' If you have a SRA metadata csv file, you can set this argument to
#' study$ScientificName[1], where study is the SRA metadata for all files
#' that was aligned.
#' @param pairedEndBam logical FALSE, else TRUE, or a logical list of
#' TRUE/FALSE per library you see will be included (run first without and check
#' what order the files will come in) 1 paired end file, then two single will
#' be c(T, F, F). If you have a SRA metadata csv file, you can set this argument to
#' study$LibraryLayout == "PAIRED", where study is the SRA metadata for all files
#' that was aligned.
#' @return a data.frame, NOTE: this is not a ORFik experiment,
#'  only a template for it!
#' @importFrom utils View
#' @export
#' @examples
#' # 1. Pick directory
#' dir <- system.file("extdata", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' # 4. Pick fasta genome of organism
#' fa <- system.file("extdata", "genome.fasta", package = "ORFik")
#' # 5. Set organism (optional)
#' org <- "Homo sapiens"
#'
#' # Create temple not saved on disc yet:
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               saveDir = NULL,
#'                               fa = fa, organism = org,
#'                               viewTemplate = FALSE)
#' ## Now fix non-unique rows: either is libre office, microsoft excel, or in R
#' template$X5[6] <- "heart"
#' # read experiment (if you set correctly)
#' df <- read.experiment(template)
#' # Save with: save.experiment(df, file = "path/to/save/experiment.csv")
#'
#' ## Create and save experiment directly:
#' ## Default location: "~/Bio_data/ORFik_experiments/"
#' #template <- create.experiment(dir = dir, exper, txdb = txdb,
#' #                               fa = fa, organism = org,
#' #                               viewTemplate = FALSE)
#' ## Custom location
#' #template <- create.experiment(dir = dir, exper, txdb = txdb,
#' #                               saveDir = "~/MY/CUSTOME/LOCATION",
#' #                               fa = fa, organism = org,
#' #                               viewTemplate = FALSE)
#' @family ORFik_experiment
create.experiment <- function(dir, exper, saveDir = "~/Bio_data/ORFik_experiments/",
                              txdb = "", fa = "", organism = "",
                              pairedEndBam = FALSE,
                              viewTemplate = TRUE,
                              types = c("bam", "bed", "wig")) {
  notDir <- !all(dir.exists(dir))
  if (notDir) stop(paste(dir[!dir.exists(dir)], "is not a existing directory!"))
  file_dt <- findLibrariesInFolder(dir, types, pairedEndBam)

  if (is(file_dt, "data.table")) { # If paired data
    files <- file_dt$forward
    df <- data.frame(matrix(ncol = 7, nrow = length(files) + 4))
    # set lib column names
    df[4,] <- c("libtype", "stage", "rep", "condition", "fraction","filepath",
                "reverse")
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
  df[5:(5+length(files)-1), 1] <- findFromPath(files, libNames())
  # set stage (sphere, shield etc)
  stages <- rbind(stageNames(), tissueNames(), cellLineNames())
  df[5:(5+length(files)-1), 2] <- findFromPath(files, stages)
  # set rep (1, 2, 3 etc)
  df[5:(5+length(files)-1), 3] <- findFromPath(files, repNames())
  # Set condition (WT, control, mutant etc)
  df[5:(5+length(files)-1), 4] <- findFromPath(files, conditionNames())

  df[1, seq(2)] <- c("name", exper)
  df[2, seq(2)] <- c("gff", txdb)
  df[3, seq(2)] <- c("fasta", fa)
  if (organism != "") df[2, seq(5, 6)] <- c("organism", organism)
  df[is.na(df)] <- ""
  if (!is.null(saveDir)) {
    cbu.path <- "/export/valenfs/data/processed_data/experiment_tables_for_R/"
    if (dir.exists(cbu.path)) { # This will only trigger on CBU server @ UIB
      message("This is on internal CBU server, saving to preset directory:")
      message(cbu.path)
      saveDir <- cbu.path
    }
    save.experiment(df, pasteDir(saveDir, exper,".csv"))
  }

  if (viewTemplate) View(df)
  return(df)
}

#' Save \code{\link{experiment}} to disc
#'
#' @param df an ORFik \code{\link{experiment}}
#' @param file name of file to save df as
#' @export
#' @importFrom utils write.table
#' @return NULL (experiment save only)
#' @examples
#' df <- ORFik.template.experiment()
#' ## Save with:
#' #save.experiment(df, file = "path/to/save/experiment.csv")
#' ## Identical (.csv not needed, can be added):
#' #save.experiment(df, file = "path/to/save/experiment")
#' @family ORFik_experiment
save.experiment <- function(df, file) {
  if (file_ext(file) == "") file <- paste0(file, ".csv")
  if (!dir.exists(dirname(file)))
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)

  write.table(x = df, file = file, sep = ",",
              row.names = FALSE, col.names = FALSE)
  return(invisible(NULL))
}

#' Find all candidate library types filenames
#'
#' From the given \code{\link{experiment}}
#' @param filepaths path to all files
#' @param candidates a data.table with 2 columns,
#' Possible names to search for, see experiment_naming family for candidates.
#' @return a candidate library types (character vector)
findFromPath <- function(filepaths, candidates) {
  dt <- candidates
  candidates <- unlist(dt$allNames)
  types <- c()
  for (path in filepaths) {
    hit <- unlist(sapply(candidates, grep, x = path))
    hitRel <- unlist(sapply(candidates, grep, x = gsub(".*/", "", path)))
    type <- if(length(hit) == 1 & length(hitRel) == 0) names(hit)
    over <- hit[names(hit) %in% names(hitRel)]
    type <- ifelse(length(over) == 1, names(over),
                   ifelse(is.null(type), "", type))
    types <- c(types, type)
  }
  types <- mainNames(types, dt)

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
    stop("stage and libtype must be colnames in df!")
  if (length(libTypes) == 0) stop("df have no valid sequencing libraries!")
  if (nrow(df) == 0) stop("df must have at least 1 row!")
  files <- df$filepath
  if (!is.null(df$reverse)) {
    reversePaths <- df$reverse[!(df$reverse %in% c("", "paired-end"))]
    files <- c(files, reversePaths)
  }

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
#' @examples
#' df <- ORFik.template.experiment()
#' bamVarName(df)
#'
#' ## without libtype
#' bamVarName(df, skip.libtype = TRUE)
#' ## Without experiment name
#' bamVarName(df, skip.experiment = TRUE)
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
      return(paste(x, y, sep = ifelse(y %in% "", "", "_")))
    return(paste(x, y, sep = ifelse(x %in% "", "", "_")))
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
  # TODO: FIX _NA for replicates
  return(gsub("_$", "", current))
}

#' Get filepaths to ORFik experiment
#'
#' If other type than "default" is given and that type is not found,
#' it will return you default filepaths without warning. \cr
#'
#' For pshifted libraries, it will load ".bedo"
#' prioritized over ".bed", if there exists both file types for the same file.
#' @inheritParams outputLibs
#' @param basename logical, default (FALSE).
#' Get relative paths instead of full. Only use for inspection!
#' @return a character vector of paths, or a list of character with 2 paths per,
#' if paired libraries exists
#' @importFrom BiocGenerics basename
#' @export
#' @family ORFik_experiment
#' @examples
#' df <- ORFik.template.experiment()
#' filepath(df, "default")
#' # If you have bedo files, see simpleLibs():
#' # filepath(df, "bedo")
#' # If you have pshifted files, see shiftFootprintsByExperiment():
#' # filepath(df, "pshifted")
filepath <- function(df, type, basename = FALSE) {
  if (!is(df, "experiment")) stop("df must be ORFik experiment!")

  paths <- lapply(df$filepath, function(x, df, type) {
    i <- which(df$filepath == x)
    input <- NULL
    if (type %in% c("bedoc", "bedo", "bed", "ofst")) {
      out.dir <- paste0(dirname(df$filepath[1]), "/",type,"/")
      if (dir.exists(out.dir)) {
        input <- paste0(out.dir,
                        remove.file_ext(x,
                                        basename = TRUE)
                        , ".", type)
        if (!file.exists(input)) type <- "default"
      } else type <- "default"
    } else if (type == "pshifted") {
      out.dir <- paste0(dirname(df$filepath[1]), "/pshifted/")
      if (dir.exists(out.dir)) {
        input <- paste0(out.dir,
                        remove.file_ext(x,
                                        basename = TRUE)
                        , "_pshifted.ofst")
        if (!file.exists(input)) { # if not ofst
          input <- paste0(out.dir,
                          remove.file_ext(x,
                                          basename = TRUE)
                          , "_pshifted.bedo")
          if (!file.exists(input)) { # if not bedo
            input <- paste0(out.dir,
                            remove.file_ext(x,
                                            basename = TRUE)
                            , "_pshifted.bed")
            if (!file.exists(input))
              type <- "default"
          }
        }
      } else type <- "default"
    }
    if (type == "default") {
      input <- x
      if (!is.null(df$reverse)) { # If reverse exists
        if (df[i,]$reverse != "")
          input <- c(x, df[i,]$reverse)
      }
    }
    if (is.null(input)) stop("filepath type not valid!")
    if (basename) input <- basename(input)
    return(input)
  }, df = df, type = type)
  if (all(lengths(paths) == 1)) {
    paths <- unlist(paths)
  }
  return(paths)
}

#' Output bam/bed/bedo/bedoc/ofst/wig files to R as variables
#'
#' Variable names defined by df (ORFik experiment DataFrame)
#' Uses multiple cores to load, defined by multicoreParam
#' @param df an ORFik \code{\link{experiment}}
#' @inheritParams matchSeqStyle
#' @param type a character(default: "default"), load files in experiment
#' or some precomputed variant, either "bedo", "bedoc", "ofst or "pshifted".
#' These are made with ORFik:::simpleLibs(), shiftFootprintsByExperiment()..
#' @param envir environment to save to, default (.GlobalEnv)
#' @param BPPARAM how many cores/threads to use? default: bpparam().
#' To see number of threads used, do \code{bpparam()$workers}
#' @return NULL (libraries set by envir assignment)
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @export
#' @examples
#' ## Load a template ORFik experiment
#' df <- ORFik.template.experiment()
#' ## Default library type load, usually bam files
#' # outputLibs(df, type = "default")
#' ## .ofst file load, if ofst files does not exists
#' ## it will load default
#' # outputLibs(df, type = "ofst")
#' ## .wig file load, if wiggle files does not exists
#' ## it will load default
#' # outputLibs(df, type = "wig")
#' @family ORFik_experiment
outputLibs <- function(df, chrStyle = NULL, type = "default",
                       envir = .GlobalEnv, BPPARAM = bpparam()) {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)

  for (df in dfl) {
    validateExperiments(df)
    varNames <- bamVarName(df)
    loaded <- c()
    for (i in 1:nrow(df)) { # For each stage
      if (exists(x = varNames[i], envir = envir, inherits = FALSE,
                 mode = "S4")) {
        loaded <- append(loaded, TRUE)
      } else loaded <- append(loaded, FALSE)
    }
    # Par apply
    if (!all(loaded)) {
      message(paste0("Outputting libraries from: ", df@experiment))
      paths <- filepath(df, type)
      libs <- bplapply(seq_along(paths),
                       function(i, paths, df, chrStyle) {
        varNames <- bamVarName(df)
        message(paste(i, ": ", varNames[i]))
        fimport(paths[i], chrStyle)
      }, BPPARAM = BPPARAM, paths = paths, chrStyle = chrStyle, df = df)

      # assign to environment
      for (i in 1:nrow(df)) { # For each stage
        assign(varNames[i], libs[[i]], envir = envir)
      }
    }
  }
  return(invisible(NULL))
}

#' Converted format of NGS libraries
#'
#' Export as either .ofst, .bedo or .bedoc files.\cr
#' Export files as .bedo files: It is a bed file with 2 score columns.
#' Gives a massive speedup when cigar string and bam flags are not needed.\cr
#' Export files as .bedoc files: If cigar is needed, gives you replicates
#' and cigar, so a fast way to load a GAlignment object, other bam flags
#' are lost. If type is bedoc addSizeColumn and method will be ignored.
#'
#' See \code{\link{export.bedo}} and \code{\link{export.bedoc}}
#' for information on file formats
#' @param df an ORFik \code{\link{experiment}}
#' @param out.dir optional output directory, default:
#' dirname(df$filepath[1]),
#' if it is NULL, it will just reassign R objects to simplified libraries.
#' @param addScoreColumn logical, default TRUE, if FALSE will not add
#' replicate numbers as score column, see ORFik::convertToOneBasedRanges.
#' @param addSizeColumn logical, default TRUE, if FALSE will not add
#' size (width) as size column, see ORFik::convertToOneBasedRanges.
#' Does not apply for .ofst or .bedoc.
#' @param must.overlap default (NULL), else a GRanges / GRangesList object, so
#' only reads that overlap (must.overlap) are kept. This is useful when you
#' only need the reads over transcript annotation or subset etc.
#' @param method character, default "None", the method to reduce ranges,
#' for more info see \code{\link{convertToOneBasedRanges}}
#' @param type a character of format, default "ofst".
#' Alternatives: "ofst", "wig","bedo" or "bedoc". Which format you want.
#' Will make a folder within out.dir with this name containing the files.
#' @return NULL (saves files to disc or R .GlobalEnv)
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' #convertLibs(df)
#' # Keep only 5' ends of reads
#' #convertLibs(df, method = "5prime")
convertLibs <- function(df,
                       out.dir = dirname(df$filepath[1]),
                       addScoreColumn = TRUE, addSizeColumn = TRUE,
                       must.overlap = NULL, method = "None",
                       type = "ofst") {
  if (!(type %in% c("ofst", "bedo", "bedoc", "wig")))
    stop("type must be either ofst, bedo or bedoc")

  validateExperiments(df)
  if (!is.null(must.overlap) & !is.gr_or_grl(must.overlap))
    stop("must.overlap must be GRanges or GRangesList object!")
  if (!is.null(out.dir)) {
    out.dir <- paste0(out.dir, "/", type, "/")
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(out.dir)) stop("could not create directory!")
    message(paste("Saving,", type, "files to:", out.dir))
  }

  outputLibs(df, type = "ofst", chrStyle = must.overlap)

  varNames <- bamVarName(df)
  i <- 1
  message("--------------------------")
  message("Converting to new format:")
  for (f in varNames) {
    message(f)
    if (type == "bedo") { # bedo
    gr <- convertToOneBasedRanges(gr = get(f),
                                  addScoreColumn = addScoreColumn,
                                  addSizeColumn = addSizeColumn,
                                  method = method)
    } else if (type %in% c("bedoc", "ofst")) {
      gr <- collapseDuplicatedReads(x = get(f),
                                    addScoreColumn = addScoreColumn)
    }

    if (!is.null(must.overlap)) gr <- optimizeReads(must.overlap, gr)

    if (!is.null(out.dir)) {
      output <- paste0(out.dir,
                       remove.file_ext(df$filepath[i], basename = TRUE),
                       ".", type)
      if (type == "bedo"){
        export.bedo(gr, output)
      } else if (type == "ofst"){
        export.ofst(gr, file = output)
      } else if (type == "wig"){
        export.wiggle(gr, file = output)
      } else # Must be bedoc, check done
        export.bedoc(gr, output)

    } else {
      assign(x = f, value = gr, envir = .GlobalEnv)
    }
    i <- i + 1
  }
}

# Keep for legacy purpose for now
#' @inherit convertLibs
#' @export
simpleLibs <- convertLibs


#' Remove bam/bed/wig files load in R as variables
#'
#' Variable names defined by df, in envir defined
#' @param df an ORFik \code{\link{experiment}}
#' @param envir environment to save to, default (.GlobalEnv)
#' @return NULL (objects removed from envir specified)
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' # Output to .GlobalEnv with:
#' # outputLibs(df)
#' # Then remove them with:
#' # remove.experiments(df)
remove.experiments <- function(df, envir = .GlobalEnv) {
  rm(list =  bamVarName(df), envir = envir)
  message(paste0("Removed loaded libraries from experiment:",
                 df@experiment))
  return(invisible(NULL))
}

#' Get all library files in folder/folders of given types
#'
#' Will try to guess paired / unpaired wig, bed, bam files.
#'
#' Set pairedEndBam if you have paired end reads as a single bam file.
#' @inheritParams create.experiment
#' @importFrom tools file_ext
#' @return (data.table) All files found from types parameter.
#' With 2 extra column (logical), is it wig pairs, and paired bam files.
findLibrariesInFolder <- function(dir, types, pairedEndBam = FALSE) {
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
      if (any(pairedEndBam)) {
        if (length(pairedEndBam) == 1) {
          files[files$reverse == "",] <- "paired-end"
        } else if (length(pairedEndBam) == nrow(files)) {
          files[(files$reverse == "") & pairedEndBam,] <- "paired-end"
        } else stop("pairedEndBam must be either length 1 or
                    number of files in experiment")
      }
    } else files <- pairs
    if (nrow(files) == 0) stop("Found no valid files in folder")
  } else if (length(files) > 0) {
    if (any(pairedEndBam)) {
      files <- data.table(forward = files, reverse = "", match = TRUE)
      files[pairedEndBam == TRUE,]$reverse <- "paired-end"
    }
  } else stop("Found no valid files in folder")
  return(files)
}

#' Get organism of the ORFik experiment
#'
#' Uses the txdb / gtf organism information, if existing.
#' @param df an ORFik \code{\link{experiment}}
#' @return organism (character vector), if no organism set: NA
#' @family ORFik_experiment
#' @importFrom BiocGenerics organism
#' @export
#' @examples
#' # if you have set organism in txdb of
#' # ORFik experiment:
#' df <- ORFik.template.experiment()
#' #organism.df(df)
#'
#' #' If you have not set the organism you can do:
#' #txdb <- GenomicFeatures::makeTxDbFromGFF("pat/to/gff_or_gff")
#' #BiocGenerics::organism(txdb) <- "Homo sapiens"
#' #saveDb(txdb, paste0("pat/to/gff_or_gff", ".db"))
#' # then use this txdb in you ORFik experiment and load:
#' # create.experiment(exper = "new_experiment",
#' #   txdb = paste0("pat/to/gff_or_gff", ".db")) ...
#' # organism.df(read.experiment("new-experiment))
organism.df <- function(df) {
  if (df@organism != "") return(df@organism)
  org <- BiocGenerics::organism(loadTxdb(df))
  if (is.na(org))
    message("Organism not set in either of experiment and txdb of gtf")
  return(org)
}

#' List current experiment available
#'
#' Will only search .csv extension, also exclude any experiment
#' with the word template.
#' @param dir directory for ORFik experiments: default:
#' "~/Bio_data/ORFik_experiments/"
#' @param pattern allowed patterns in experiment file name:
#' default ("*", all experiments)
#' @param libtypeExclusive search for experiments with exclusivly this
#' libtype, default (NULL, all)
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a data.table, 1 row per experiment with columns experiment (name), libtypes
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @export
#' @examples
#' ## Make your experiments
#' df <- ORFik.template.experiment(TRUE)
#' df2 <- df[1:6,] # Only first 2 libs
#' ## Save them
#' # save.experiment(df, "~/Bio_data/ORFik_experiments/exp1.csv")
#' # save.experiment(df2, "~/Bio_data/ORFik_experiments/exp1_subset.csv")
#' ## List all experiment you have:
#' ## Path above is default path, so no dir argument needed
#' #list.experiments()
#' #list.experiments(pattern = "subset")
#' ## For non default directory experiments
#' #list.experiments(dir = "MY/CUSTOM/PATH)
list.experiments <- function(dir =  "~/Bio_data/ORFik_experiments/",
                             pattern = "*", libtypeExclusive = NULL,
                             BPPARAM = bpparam()) {
  experiments <- list.files(path = dir, pattern = "\\.csv")
  if (length(experiments) == 0) { # This will only trigger on CBU server @ UIB
    cbu.path <- "/export/valenfs/data/processed_data/experiment_tables_for_R/"
    if (dir.exists(cbu.path)) { # If on UIB SERVER
      dir <- cbu.path
      experiments <- list.files(path = dir, pattern = "\\.csv")
    }
  }

  experiments <- grep(experiments, pattern = pattern, value = TRUE)
  experiments <- experiments[grep(experiments, pattern = "template", value = FALSE, invert = TRUE)]
  if (length(experiments) == 0) {
    message(paste("Searching for experiments in dir:", dir))
    stop("No experiments found, have you made any ?")
  }

  info <- bplapply(experiments, function(x, dir) { # Open each experiment in parallell
    e <- read.experiment(x, dir)
    list(libtype = unique(e$libtype), runs = length(e$libtype), organism = e@organism)
  }, dir = dir)

  info <- unlist(info, recursive = FALSE)
  libtypes <- info[grep("libtype", names(info))]
  samples <- unlist(info[names(info) == "runs"])
  organism <- unlist(info[names(info) == "organism"])

  dt <- data.table(name = gsub(".csv", "", experiments), libtypes, samples, organism)
  dt <- dt[order(organism, name),]
  if (!is.null(libtypeExclusive)) {
    message(paste("subset on libtype:", libtypeExclusive))
    match <- lapply(dt$libtypes, function(i) any(libtypeExclusive %in% i))
    match <- unlist(match)

    dt <- dt[match,]
  }
  return(dt)
}

#' An ORFik experiment to see how it looks
#'
#' NOTE! This experiment should only be used for testing, since
#' it is just sampled data internal in ORFik.
#' @param as.temp logical, default FALSE, load as ORFik experiment.
#' If TRUE, loads as data.frame template of the experiment.
#' @return an ORFik \code{\link{experiment}}
#' @export
#' @family ORFik_experiment
#' @examples
#' ORFik.template.experiment()
ORFik.template.experiment <- function(as.temp = FALSE) {
  dir <- system.file("extdata", "", package = "ORFik")
  # 2. Pick an experiment name
  exper <- "ORFik"
  # 3. Pick .gff/.gtf location
  txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
  fa <- system.file("extdata", "genome.fasta", package = "ORFik")
  template <- create.experiment(dir = dir, saveDir = NULL,
                                exper, txdb = txdb, fa = fa,
                                organism = "Homo sapiens",
                                viewTemplate = FALSE)
  # read experiment
  if (as.temp) return(template)
  return(read.experiment(template))
}
