

#' Find all candidate library types filenames
#'
#' From the given \code{\link{experiment}}
#' @param filepaths path to all files
#' @param candidates a data.table with 2 columns,
#' Possible names to search for, see experiment_naming family for candidates.
#' @param slot character, default "auto". If auto, use auto guessing of slot,
#' else must be a character vector of length 1 or equal length as filepaths.
#' @return a candidate library types (character vector)
#' @keywords internal
findFromPath <- function(filepaths, candidates, slot = "auto") {
  if (all(slot != "auto", na.rm = TRUE)) { # If not auto guess
    if(length(slot) != 1 & length(slot) != length(filepaths)) {
      stop("When experiment slot is not auto, length must be 1 or length(files)!")
    } else return(slot)
  }
  dt <- candidates
  candidates <- unlist(dt$allNames)
  types <- c()
  for (path in filepaths) {
    hit <- names(unlist(sapply(candidates, grep, x = path)))
    if (length(hit) > 0) # Remove multiple hits from same group:
      hit <- unique(mainNames(hit, dt))
    hitRel <- names(unlist(sapply(candidates, grep, x = gsub(".*/", "", path))))
    if (length(hitRel) > 0) # Remove multiple hits from same group:
      hitRel <- unique(mainNames(hitRel, dt))

    # Assign the one (if existing) with 1 match
    type <- if(length(hit) == 1 & length(hitRel) == 0) {
      hit
    } else if(length(hit) == 0 & length(hitRel) == 1) {
      hitRel
    }
    over <- hit[hit %in% hitRel] #overlaps
    # If exactly one overlap, set to that one
    type <- ifelse(length(over) == 1, over,
                   ifelse(is.null(type), "", type))
    types <- c(types, type)
  }

  return(types)
}


#' Which type of library type in \code{\link{experiment}}?
#' @param df an ORFik \code{\link{experiment}}
#' @param uniqueTypes logical, default TRUE. Only return unique lib types.
#' @return library types (character vector)
#' @family ORFik_experiment
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' libraryTypes(df)
#' libraryTypes(df, uniqueTypes = FALSE)
libraryTypes <- function(df, uniqueTypes = TRUE) {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  libtypes <- c()
  for (df in dfl) {
    if (is(df, "experiment")) {
      libtypes <- if (uniqueTypes) {
        c(libtypes, unique(df$libtype))
      } else c(libtypes, df$libtype)

    } else if (is(df, "character") | is(df, "factor")) {
      libtypes <- if (uniqueTypes) {
        c(libtypes, unique(gsub("_.*", x = df, replacement = "")))
      } else c(libtypes, gsub("_.*", x = df, replacement = ""))
    } else stop("library types must be data.frame or character vector!")
  }
  return(libtypes)
}

#' Validate ORFik \code{\link{experiment}}
#'
#' Check for valid existing, non-empty and all unique.
#' A good way to see if your experiment is valid.
#' @inheritParams outputLibs
#' @return NULL (Stops if failed)
#' @family ORFik_experiment
#' @keywords internal
validateExperiments <- function(df, library.names = bamVarName(df)) {
  libTypes <- libraryTypes(df)
  if (!is(df, "experiment")) stop("df must be experiment!")
  if (!all((c("stage", "libtype") %in% colnames(df))))
    stop("stage and libtype must be colnames in df!")
  if (length(libTypes) == 0) stop("df have no valid sequencing libraries!")
  if (nrow(df) == 0) stop("df must have at least 1 row!")
  files <- df$filepath
  if (length(df$filepath) == 0) stop("df have no filepaths!")
  if (!is.null(df$reverse)) {
    reversePaths <- df$reverse[!(df$reverse %in% c("", "paired-end"))]
    files <- c(files, reversePaths)
  }

  emptyFiles <- file.size(files) == 0
  if (any(is.na(emptyFiles))) {
    message("Error in experiment:", name(df))
    stop(paste("File does not exist:\n", files[is.na(emptyFiles)]))
  }

  if (any(emptyFiles)) {
    print(files[emptyFiles])
    stop("Empty files in list, see above for which")
  }
  stopifnot(is.character(library.names) && (length(library.names) == nrow(df)))
  names <- library.names
  if (length(names) != length(unique(names))) {
    message("Duplicated rows: ", paste(names[duplicated(names)],
                                      collapse = " ; "))
    stop("Experiment table has non-unique rows!",
         " Update either replicate, stage, condition or fraction,",
         " to get unique rows!")
  }

  if (length(files) != length(unique(files))) {
    message("Duplicated filepaths: ", paste(files[duplicated(files)],
                                            collapse = " ; "))
    stop("Duplicated filepaths in experiment!")
  }
}

#' Get library variable names from ORFik \code{\link{experiment}}
#'
#' What will each sample be called given the columns of the experiment?
#' A column is included if more than 1 unique element value exist in that column.
#' @inheritParams bamVarNamePicker
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
                       skip.libtype = FALSE,
                       fraction_prepend_f = TRUE) {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  varName <- character()
  for (df in dfl) {
    res <- vapply(seq(1, nrow(df), length.out = nrow(df)),
                  function(i) bamVarNamePicker(df[i,], skip.replicate,
                                               skip.condition, skip.stage,
                                               skip.fraction,
                                               skip.experiment, skip.libtype,
                                               fraction_prepend_f), character(1))
    varName <- c(varName, res)
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
#' @param fraction_prepend_f a logical (TRUE), include "f" in front of
#' fraction, useful for knowing what fraction is.
#' @return variable name of library (character vector)
#' @keywords internal
bamVarNamePicker <- function(df, skip.replicate = FALSE,
                             skip.condition = FALSE,
                             skip.stage = FALSE, skip.fraction = FALSE,
                             skip.experiment = FALSE, skip.libtype = FALSE,
                             fraction_prepend_f = TRUE) {
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
  if(!(skip.condition | is.na(cond)))
    current <- spaste(current, cond)
  if (!(skip.stage | is.na(stage)))
    current <- spaste(current, stage)
  if (!(skip.fraction | is.null(frac) | is.na(frac))) {
    if (frac != "") {
      if (fraction_prepend_f) {
        current <- spaste(current, paste0("f", frac))
      } else current <- spaste(current, frac)
    }

  }
  # TODO: FIX _NA for replicates
  if (!(skip.replicate | is.null(rep) | is.na(rep) | (rep == "")))
    current <- spaste(current, paste0("r", rep))
  if (! (skip.experiment | is.null(df@experiment)))
    current <- spaste(df@experiment, current, TRUE)

  current <- gsub(pattern = "__", "_", current)

  return(gsub("_$", "", current))
}

#' Get filepaths to ORFik experiment
#'
#' If other type than "default" is given and that type is not found
#' (and 'fallback' is TRUE), it will return you ofst files, if they do not exist,
#' then default filepaths without warning. \cr
#'
#' For pshifted libraries, if "pshifted" is specified as type: if
#'  if multiple formats exist it will use a priority:
#'  ofst -> bigwig -> wig -> bed. For formats outside default, all files
#' must be stored in the directory of the first file:
#' \code{base_folder <- libFolder(df)}
#' @inheritParams outputLibs
#' @param basename logical, default (FALSE).
#' Get relative paths instead of full. Only use for inspection!
#' @param fallback logical, default: type %in% c("pshifted", "bed", "ofst", "bedoc", "bedo").
#'  If TRUE, will use type fallback, see above for info.
#' @param suffix_stem character, default "AUTO". Which is "" for all except
#' type = "pshifted". Then it is "_pshifted" appended to end of names before
#' format. Can be vector, then it searches suffixes in priority: so if you insert
#'  c("_pshifted", ""), it will look for suffix _pshifted, then the empty suffix.
#' @param base_folders character vector, default libFolder(df),
#' path to base folder to search for library variant directories.
#' If single path (length == 1), it will apply to all libraries in df.
#' If df is a collection, an experiment where libraries are put in different
#' folders and library variants like pshifted are put inside those respective
#' folders, set base_folders = libFolder(df, mode = "all")
#' @return a character vector of paths, or a list of character with 2 paths per,
#' if paired libraries exists
#' @importFrom BiocGenerics basename
#' @export
#' @family ORFik_experiment
#' @examples
#' df <- ORFik.template.experiment()
#' filepath(df, "default")
#' # Subset
#' filepath(df[9,], "default")
#' # Other format path
#' filepath(df[9,], "ofst")
#' ## If you have pshifted files, see shiftFootprintsByExperiment()
#' filepath(df[9,], "pshifted") # <- falls back to ofst
filepath <- function(df, type, basename = FALSE,
                     fallback = type %in% c("pshifted", "bed", "ofst", "bedoc", "bedo"),
                     suffix_stem = "AUTO",
                     base_folders = libFolder(df)) {
  if (!is(df, "experiment")) stop("df must be ORFik experiment!")
  stopifnot(length(type) == 1)
  rel_folder <- c(cov = "cov_RLE/", covl = "cov_RLE_List/", bigwig = "bigwig/",
                  pshifted = "pshifted/")
  if (length(base_folders) == 1) base_folders <- rep(base_folders, nrow(df))

  paths <- lapply(df$filepath, function(x, df, type) {
    i <- which(df$filepath == x)
    base_folder <- base_folders[i]
    name_stem <- remove.file_ext(x, basename = TRUE)
    found_valid_file <- FALSE
    input <- NULL
    if (type == "pshifted") {
      if (all(suffix_stem == "AUTO")) suffix_stem <- "_pshifted"
      out.dir.type <- file.path(base_folder, rel_folder["pshifted"])
      type_dir_exists <- dir.exists(out.dir.type)
      if (type_dir_exists) {
        formats <- paste0(".", c("ofst", "bigWig", "wig", "bed"))
        for (suf_stem in suffix_stem) {
          input_stem <- paste0(out.dir.type, name_stem, suf_stem)
          for (format in formats) {
            input <-  paste0(input_stem, format)
            found_valid_file <- file.exists(input)
            if (found_valid_file) break
          }
        }
      }
      # If not hit, fall back to default ofst files
      if (!found_valid_file) {
        if (fallback) {
          type <- "ofst"
        } else stop(filepath_errors(type))
    }} else if (all(suffix_stem == "AUTO")) suffix_stem <- ""


    ext <- c(cov = ".covrds", covl = ".covrds", bigwig = ".bigWig")
    paired_files <- list(bigwig = c("_forward", "_reverse"))
    if (type %in% c("cov", "covl", "bigwig")) {
      t <- type
      out.dir.type <- file.path(base_folder, rel_folder[t])
      paired_file <- if(is.null(paired_files[[t]])) {""} else paired_files[[t]]
      for (suf_stem in suffix_stem) {
        input <- paste0(out.dir.type, name_stem, suf_stem, paired_file, ext[t])
        found_valid_file <- all(file.exists(input))
        if (found_valid_file) break
      }

      if (!found_valid_file) {
        if (fallback) {
          type <- "ofst"
        } else stop(filepath_errors(t))
      }
    }

    if (type %in% c("bed", "ofst", "bedoc", "bedo")) {
      out.dir <- paste0(base_folder, "/",type,"/")
      if (dir.exists(out.dir)) {
        input <- paste0(out.dir, remove.file_ext(x,basename = TRUE), ".", type)
        if (!file.exists(input)) type <- "default"
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

filepath_errors <- function(format) {
  stopifnot(is(format, "character"))
  stem <- "File did not exist,"
  candidates <- c(cov = paste(stem, "did you create covRle yet?"),
                  covl = paste(stem, "did you create covRleList yet?"),
                  bigwig = paste(stem, "did you create bigwig yet?",
                    " (only supports naming _forward.bigWig etc for now)"),
                  pshifted = paste(stem, "did you create pshifted files yet?"))
  stopifnot(format %in% names(candidates))
  return(candidates[format])
}

#' Output NGS libraries to R as variables
#'
#' By default loads the original files of the experiment into
#' the global environment, named by the rows of the experiment
#' required to make all libraries have unique names.\cr
#' Uses multiple cores to load, defined by multicoreParam
#'
#' The functions checks if the total set of libraries have already been loaded:
#' i.e. Check if all names from 'library.names' exists as S4 objects in
#' environment of experiment.
#' @param df an ORFik \code{\link{experiment}}
#' @inheritParams fimport
#' @param type a character(default: "default"), load files in experiment
#' or some precomputed variant, like "ofst" or "pshifted".
#' These are made with ORFik:::convertLibs(),
#' shiftFootprintsByExperiment(), etc.
#' Can also be custom user made folders inside the experiments bam folder.
#' It acts in a recursive manner with priority: If you state "pshifted",
#' but it does not exist, it checks "ofst". If no .ofst files, it uses
#' "default", which always must exists.\cr Presets are (folder is relative
#' to default lib folder, some types fall back to other formats if folder does not exist):\cr
#' - "default": load the original files for experiment, usually bam.\cr
#' - "ofst": loads ofst files from the ofst folder, relative to lib folder (falls back to default)\cr
#' - "pshifted": loads ofst, wig or bigwig from pshifted folder (falls back to ofst, then default)\cr
#' - "cov": Load covRle objects from cov_RLE folder (fail if not found)\cr
#' - "covl": Load covRleList objects, from cov_RLE_List folder (fail if not found)\cr
#' - "bed": Load bed files, from bed folder (falls back to default)\cr
#' - Other formats must be loaded directly with fimport
#' @param paths character vector, the filpaths to use,
#' default \code{filepath(df, type)}. Change type argument if not correct.
#' If that is not enough, then you can also update this argument.
#' But be careful about using this directly.
#' @param naming a character (default: "minimum"). Name files as minimum
#' information needed to make all files unique. Set to "full" to get full
#' names. Set to "fullexp", to get full name with experiment name as prefix,
#' the last one guarantees uniqueness.
#' @param library.names character vector, names of libraries, default:
#' name_decider(df, naming)
#' @param output.mode character, default "envir". Output libraries to environment.
#' Alternative: "list", return as list. "envirlist", output to envir and return
#' as list. If output is list format, the list elements are named from:
#' \code{bamVarName(df.rfp)} (Full or minimum naming based on 'naming' argument)
#' @param envir environment to save to, default
#' \code{envExp(df)}, which defaults to .GlobalEnv, but can be set with
#' \code{envExp(df) <- new.env()} etc.
#' @param verbose logical, default TRUE, message about library output status.
#' @param force logical, default TRUE If TRUE, reload library files even if
#' matching named variables are found in environment used by experiment
#'  (see \code{\link{envExp}}) A simple way to make
#' sure correct libraries are always loaded. FALSE is faster if data
#' is loaded correctly already.
#' @param BPPARAM how many cores/threads to use? default: bpparam().
#' To see number of threads used, do \code{bpparam()$workers}.
#' You can also add a time remaining bar, for a more detailed pipeline.
#' @return NULL (libraries set by envir assignment), unless output.mode is
#' "list" or "envirlist": Then you get a list of the libraries.
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
#' ## Load as list
#' outputLibs(df, output.mode = "list")
#' ## Load libs to new environment (called ORFik in Global)
#' # outputLibs(df, envir = assign(name(df), new.env(parent = .GlobalEnv)))
#' ## Load to hidden environment given by experiment
#' # envExp(df) <- new.env()
#' # outputLibs(df)
#'
#' @family ORFik_experiment
outputLibs <- function(df, type = "default", paths = filepath(df, type),
                       param = NULL, strandMode = 0, naming = "minimum",
                       library.names = name_decider(df, naming),
                       output.mode = "envir", chrStyle = NULL,
                       envir = envExp(df), verbose = TRUE, force = TRUE,
                       BPPARAM = bpparam()) {
  stopifnot(output.mode %in% c("envir", "list", "envirlist"))
  stopifnot(is.character(type))

  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  all_libs <- NULL
  for (df in dfl) {
    validateExperiments(df, library.names)
    loaded <- libs_are_loaded(library.names, envir)
    varNames <- library.names
    import_is_needed <- !all(loaded) | force
    if (import_is_needed) {
      if (verbose) message(paste0("Outputting libraries from: ", name(df)))
      if (is(BPPARAM, "SerialParam")) {
        libs <- lapply(seq_along(paths),
                       function(i, paths, df, chrStyle, param, strandMode, varNames, verbose) {
                           if (verbose) message(paste(i, ": ", varNames[i]))
                           fimport(paths[i], chrStyle, param, strandMode)
                         }, paths = paths, chrStyle = chrStyle, df = df,
                         param = param, strandMode = strandMode, varNames = varNames,
                         verbose = verbose)
      } else {
        libs <- bplapply(seq_along(paths),
                         function(i, paths, df, chrStyle, param, strandMode, varNames, verbose) {
                           if (verbose) message(paste(i, ": ", varNames[i]))
                           fimport(paths[i], chrStyle, param, strandMode)
                         }, paths = paths, chrStyle = chrStyle, df = df,
                         param = param, strandMode = strandMode, varNames = varNames,
                         verbose = verbose, BPPARAM = BPPARAM)
      }

      # assign to environment
      if (output.mode %in% c("envir", "envirlist")) {
        for (i in 1:nrow(df)) { # For each stage
          assign(varNames[i], libs[[i]], envir = envir)
        }
      }
      if (output.mode %in% c("list", "envirlist")) {
        names(libs) <- varNames
        all_libs <- c(all_libs, libs)
      }

    } else if (output.mode %in% c("list", "envirlist")) {
      libs <- lapply(varNames, function(x) get(x, envir = envir))
      names(libs) <- varNames
      all_libs <- c(all_libs, libs)
    }
  }
  if (!is.null(all_libs)) return(all_libs)
  return(invisible(NULL))
}

name_decider <- function(df, naming) {
  stopifnot(naming %in% c("minimum", "full", "fullexp"))
  varNames <-
    if (naming == "minimum") {
       bamVarName(df)
    } else if (naming == "full") {
      bamVarName(df, FALSE, FALSE, FALSE, FALSE)
    } else bamVarName(df, FALSE, FALSE, FALSE, FALSE, FALSE)

  return(varNames)
}

libs_are_loaded <- function(lib_names, envir) {
  vapply(lib_names, function(lib) exists(x = lib, envir = envir, inherits = FALSE,
                                           mode = "S4"), logical(1))
}

#' Converted format of NGS libraries
#'
#' Export as either .ofst, .wig, .bigWig,.bedo (legacy format) or .bedoc (legacy format) files:\cr
#' Export files as .ofst for fastest load speed into R.\cr
#' Export files as .wig / bigWig for use in IGV or other genome browsers.\cr
#' The input files are checked if they exist from: \code{envExp(df)}.\cr
#'
#' We advice you to not use this directly, as other function are more safe
#' for library type conversions. See family description below. This is
#' mostly used internally in ORFik. It is only adviced to use if large bam files
#' are already loaded in R and conversions are wanted from those.
#'
#' See \code{\link{export.ofst}}, \code{\link{export.wiggle}},
#' \code{\link{export.bedo}} and \code{\link{export.bedoc}}
#' for information on file formats.\cr
#' If libraries of the experiment are
#' already loaded into environment (default: .globalEnv) is will export
#' using those files as templates. If they are not in environment the
#' .ofst files from the bam files are loaded (unless you are converting
#' to .ofst then the .bam files are loaded).
#' @inheritParams outputLibs
#' @param out.dir optional output directory, default: libFolder(df),
#' if it is NULL, it will just reassign R objects to simplified libraries.
#' Will then create a final folder specfied as: paste0(out.dir, "/", type, "/").
#' Here the files will be saved in format given by the type argument.
#' @param addScoreColumn logical, default TRUE, if FALSE will not add
#' replicate numbers as score column, see ORFik::convertToOneBasedRanges.
#' @param addSizeColumn logical, default TRUE, if FALSE will not add
#' size (width) as size column, see ORFik::convertToOneBasedRanges.
#' Does not apply for (GAlignment version of.ofst) or .bedoc. Since they
#' contain the original cigar.
#' @param must.overlap default (NULL), else a GRanges / GRangesList object, so
#' only reads that overlap (must.overlap) are kept. This is useful when you
#' only need the reads over transcript annotation or subset etc.
#' @param method character, default "None", the method to reduce ranges,
#' for more info see \code{\link{convertToOneBasedRanges}}
#' @param type character, output format, default "ofst".
#' Alternatives: "ofst", "bigWig", "wig","bedo" or "bedoc". Which format you want.
#' Will make a folder within out.dir with this name containing the files.
#' @param input.type character, input type "ofst". Remember this function
#' uses the loaded libraries if existing, so this argument is usually ignored.
#' Only used if files do not already exist.
#' @param libs list, output of outputLibs as list of
#'  GRanges/GAlignments/GAlignmentPairs objects. Set input.type and force arguments to define parameters.
#' @param reassign.when.saving logical, default FALSE. If TRUE, will reassign
#' library to converted form after saving. Ignored when out.dir = NULL.
#' @return invisible NULL (saves files to disc or R .GlobalEnv)
#' @export
#' @family lib_converters
#' @examples
#' df <- ORFik.template.experiment()
#' #convertLibs(df, out.dir = NULL)
#' # Keep only 5' ends of reads
#' #convertLibs(df, out.dir = NULL, method = "5prime")
convertLibs <- function(df,
                        out.dir = libFolder(df),
                        addScoreColumn = TRUE, addSizeColumn = TRUE,
                        must.overlap = NULL, method = "None",
                        type = "ofst", input.type = "ofst",
                        reassign.when.saving = FALSE,
                        envir = envExp(df),
                        force = TRUE,
                        library.names = bamVarName(df),
                        libs = outputLibs(df, type = input.type, chrStyle = must.overlap, library.names = library.names,
                                          output.mode = "list", force = force, BPPARAM = BPPARAM),
                        BPPARAM = bpparam()) {
  if (!(type %in% c("ofst", "bedo", "bedoc", "wig", "bigWig")))
    stop("type must be either ofst, bedo or bedoc")
  validateExperiments(df, library.names)
  stopifnot(length(libs) == nrow(df))
  if (!is.null(must.overlap) & !is.gr_or_grl(must.overlap))
    stop("must.overlap must be GRanges or GRangesList object!")
  if (!is.null(out.dir)) {
    out.dir <- paste0(out.dir, "/", type, "/")
    dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(out.dir)) stop("could not create directory!")
    message(paste("Saving,", type, "files to:", out.dir))
  }


  message("--------------------------")
  message("Converting libraries to new format: ", type)
  lapply(seq_along(libs), function(i) {
    f <- library.names[i]
    message(f)
    if (type %in% c("bedo", "wig")) { # bedo, wig
      gr <- convertToOneBasedRanges(gr = libs[[i]],
                                    addScoreColumn = addScoreColumn,
                                    addSizeColumn = addSizeColumn,
                                    method = method)
    } else if (type %in% c("bedoc", "ofst")) {
      gr <- collapseDuplicatedReads(x = libs[[i]],
                                    addScoreColumn = addScoreColumn)
    }

    if (!is.null(must.overlap)) gr <- optimizeReads(must.overlap, gr)

    if (!is.null(out.dir)) {
      output <- paste0(out.dir,
                       remove.file_ext(df$filepath[i], basename = TRUE),
                       ".", type)
      if (type == "bedo") {
        export.bedo(gr, output)
      } else if (type == "ofst"){
        export.ofst(gr, file = output)
      } else if (type == "wig"){
        export.wiggle(gr, file = output)
      } else if (type == "bigWig"){
        export.bigWig(gr, file = output)
      } else # Must be bedoc, check done
        export.bedoc(gr, output)
      if (reassign.when.saving)
        assign(x = f, value = gr, envir = envir)
    } else {
      assign(x = f, value = gr, envir = envir)
    }
    return(invisible(NULL))
  })
  return(invisible(NULL))
}

# Keep for legacy purpose for now
#' @inherit convertLibs
#' @export
simpleLibs <- convertLibs

#' Merge and save libraries of experiment
#'
#' Aggregate count of reads (from the "score" column) by making a merged library.
#' Only allowed for .ofst files!
#' @inheritParams outputLibs
#' @inheritParams ofst_merge
#' @param out_dir Ouput directory, default \code{file.path(dirname(df$filepath[1]), "ofst_merged")},
#' saved as "all.ofst" in this folder if mode is "all". Use a folder called pshifted_merged, for
#' default Ribo-seq ofst files.
#' @param mode character, default "all". Merge all or "rep" for collapsing replicates only, or
#' "lib" for collapsing all per library type.
#' @return NULL, files saved to disc. A data.table with a score column that now contains the sum
#' of scores per merge setting.
#' @export
#' @examples
#' df2 <- ORFik.template.experiment()
#' df2 <- df2[df2$libtype == "RFP",]
#' # Merge all
#' #mergeLibs(df2, tempdir(), mode = "all", type = "default")
#' # Read as GRanges with mcols
#' #fimport(file.path(tempdir(), "all.ofst"))
#' # Read as direct fst data.table
#' #read_fst(file.path(tempdir(), "all.ofst"))
#' # Collapse replicates
#' #mergeLibs(df2, tempdir(), mode = "rep", type = "default")
#' # Collapse by lib types
#' #mergeLibs(df2, tempdir(), mode = "lib", type = "default")
mergeLibs <- function(df, out_dir = file.path(libFolder(df), "ofst_merged"), mode = "all",
                      type = "ofst", keep_all_scores = TRUE, paths = filepath(df, type)) {
  stopifnot(mode %in% c("all", "rep", "lib"))
  stopifnot(nrow(df) == length(paths))
  filepaths <- paths
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lib_names_full <- bamVarName(df, skip.libtype = FALSE)

  if (mode == "rep") {
    lib_names <- bamVarName(df, skip.libtype = FALSE, skip.replicate = TRUE)
    libs <- lapply(unique(lib_names), function(x) grep(x, lib_names))
    names(libs) <- unique(lib_names)
  } else if (mode == "lib") {
    lib_names <- bamVarName(df, TRUE, TRUE, TRUE, TRUE, TRUE)
    libs <- lapply(unique(lib_names), function(x) grep(x, lib_names))
    names(libs) <- unique(lib_names)
  } else {
    libs <- list(all = seq(nrow(df)))
  }

  for (name in names(libs)) {
    specific_paths <- filepaths[libs[[name]]]
    specific_names <- lib_names_full[libs[[name]]]
    save_path <- file.path(out_dir, paste0(name, ".ofst"))
    write_fst(ofst_merge(specific_paths, specific_names, keep_all_scores), save_path)
  }
  return(invisible(NULL))
}

#' Remove ORFik experiment libraries load in R
#'
#' Variable names defined by df, in envir defined
#' @inheritParams outputLibs
#' @return NULL (objects removed from envir specified)
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' # Output to .GlobalEnv with:
#' # outputLibs(df)
#' # Then remove them with:
#' # remove.experiments(df)
remove.experiments <- function(df, envir = envExp(df)) {
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
#' @keywords internal
findLibrariesInFolder <- function(dir, types, pairedEndBam = FALSE) {
  notDir <- !all(dir.exists(dir))
  if (notDir) stop(paste(dir[!dir.exists(dir)], "is not a existing directory!"))

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

#' List current experiment available
#'
#' Will only search .csv extension, also exclude any experiment
#' with the word template.
#' @param dir directory for ORFik experiments: default:
#' ORFik::config()["exp"], which by default is:
#' "~/Bio_data/ORFik_experiments/"
#' @inheritParams read.experiment
#' @param pattern allowed patterns in experiment file name:
#' default ("*", all experiments)
#' @param libtypeExclusive search for experiments with exclusivly this
#' libtype, default (NULL, all)
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a data.table, 1 row per experiment with columns:\cr
#'  - experiment (name),\cr
#'  - organism\cr
#'  - author \cr
#'  - libtypes\cr
#'  - number of samples
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
list.experiments <- function(dir =  ORFik::config()["exp"],
                             pattern = "*", libtypeExclusive = NULL,
                             validate = TRUE,
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

  info <- bplapply(experiments, function(x, dir, validate) { # Open each experiment in parallell
    e <- read.experiment(x, dir, validate)
    list(libtype = unique(e$libtype), runs = length(e$libtype), organism = e@organism,
         author = e@author)
  }, dir = dir, validate = validate, BPPARAM = BPPARAM)

  info <- unlist(info, recursive = FALSE)
  libtypes <- info[grep("libtype", names(info))]
  samples <- unlist(info[names(info) == "runs"])
  organism <- unlist(info[names(info) == "organism"])
  author <- unlist(info[names(info) == "author"])

  dt <- data.table(name = gsub(".csv", "", experiments), organism, author,libtypes, samples)
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
#' Toy-data created to resemble Zebrafish genes:\cr
#' Number of genes: 150\cr
#' Ribo-seq: 1 library
#' @param as.temp logical, default FALSE, load as ORFik experiment.
#' If TRUE, loads as data.frame template of the experiment.
#' @return an ORFik \code{\link{experiment}}
#' @export
#' @family ORFik_experiment
#' @examples
#' ORFik.template.experiment.zf()
ORFik.template.experiment.zf <- function(as.temp = FALSE) {
  sample_dir <- system.file("extdata/Danio_rerio_sample", "", package = "ORFik")
  # 2. Pick an experiment name
  exper <- "ORFik"
  # 3. Pick .gff/.gtf location
  txdb <- system.file("extdata/references/danio_rerio",
                      "annotations.gtf", package = "ORFik")
  fa <- system.file("extdata/references/danio_rerio",
                    "genome_dummy.fasta", package = "ORFik")
  template <- create.experiment(dir = sample_dir, saveDir = NULL,
                                exper, txdb = txdb, fa = fa,
                                organism = "Danio rerio",
                                author = "Tjeldnes",
                                viewTemplate = FALSE)
  # read experiment
  if (as.temp) return(template)
  return(read.experiment(template))
}

#' An ORFik experiment to see how it looks
#'
#' Toy-data created to resemble human genes:\cr
#' Number of genes: 6\cr
#' Genome size: 1161nt x 6 chromosomes = 6966 nt\cr
#' Experimental design (2 replicates, Wild type vs Mutant)\cr
#' CAGE: 4 libraries\cr
#' PAS (poly-A): 4 libraries\cr
#' Ribo-seq: 4 libraries\cr
#' RNA-seq: 4 libraries\cr
#' @param as.temp logical, default FALSE, load as ORFik experiment.
#' If TRUE, loads as data.frame template of the experiment.
#' @return an ORFik \code{\link{experiment}}
#' @export
#' @family ORFik_experiment
#' @examples
#' ORFik.template.experiment()
ORFik.template.experiment <- function(as.temp = FALSE) {
  sample_dir <- system.file("extdata/Homo_sapiens_sample", "", package = "ORFik")
  # 2. Pick an experiment name
  exper <- "ORFik"
  # 3. Pick .gff/.gtf location
  txdb <- system.file("extdata/references/homo_sapiens",
                      "Homo_sapiens_dummy.gtf.db", package = "ORFik")
  fa <- system.file("extdata/references/homo_sapiens",
                    "Homo_sapiens_dummy.fasta", package = "ORFik")
  template <- create.experiment(dir = sample_dir, saveDir = NULL,
                                exper, txdb = txdb, fa = fa,
                                organism = "Homo sapiens",
                                author = "Tjeldnes",
                                viewTemplate = FALSE, types = "ofst")
  # read experiment
  if (as.temp) return(template)
  return(read.experiment(template))
}
