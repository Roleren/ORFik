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
#' \code{ORFik::config()["exp"]}, which has default "~/Bio_data/ORFik_experiments/"\cr
#' Set to NULL if you don't want to save it to disc.
#' Does not apply if file argument is not a path (can also be a data.frame).
#' Also does not apply if file argument was given as full path.
#' @param validate logical, default TRUE. Abort if any library files does not exist.
#' Do not set this to FALSE, unless you know what you are doing!
#' @param output.env an environment, default .GlobalEnv. Which environment
#' should ORFik output libraries to (if this is done),
#' can be updated later with \code{envExp(df) <- new.env()}.
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
read.experiment <-  function(file, in.dir = ORFik::config()["exp"],
                             validate = TRUE, output.env = .GlobalEnv) {
  stopifnot(is.environment(output.env))
  stopifnot(is.logical(validate))
  # Split up metadata and library data
  parse_list <- experiment_parse_list_info(file, in.dir)

  # Parse metadata info from list
  info <- parse_list$info
  assembly <- ifelse(info[1,5] == "assembly" & !is.na(info[1,6]),
                     info[1,6], "")
  org <- ifelse(info[2,5] == "organism" & !is.na(info[2,6]),
                info[2,6], "")
  author <- ifelse(info[3,5] == "author" & !is.na(info[3,6]),
                   info[3,6], "")
  resultFolder <- ifelse(info[1, 3] == "results" & !is.na(info[1,4]),
                          info[1,4], "")
  exper <- info[1, 2]
  txdb <- ifelse(is.na(info[2, 2]),  "", info[2, 2])
  fa <- ifelse(is.na(info[3, 2]),  "", info[3, 2])

  df <- experiment(experiment = exper, txdb = txdb, fafile = fa,
                   organism = org, author = author,
                   assembly = assembly,
                   listData = parse_list$listData,
                   expInVarName = FALSE,
                   resultFolder = resultFolder,
                   envir = output.env)

  if (validate) validateExperiments(df)
  return(df)
}

#' Create an ORFik \code{\link{experiment}}
#'
#' Create a single R object that stores and controls all results relevant to
#' a specific Next generation sequencing experiment.
#' Click the experiment link above in the title if you are not sure what an
#' ORFik experiment is.\cr\cr
#' By using files in a folder / folders. It will make an experiment table
#' with information per sample, this object allows you to use the extensive API in
#' ORFik that works on experiments. \cr\cr
#' Information Auto-detection:\cr
#' There will be several columns you can fill in, when creating the object,
#' if the files have logical names like (RNA-seq_WT_rep1.bam) it will try to auto-detect
#' the most likely values for the columns. Like if it is RNA-seq or Ribo-seq,
#' Wild type or mutant, is this replicate 1 or 2 etc.\cr
#' You will have to fill in the details that were not auto detected.
#' Easiest way to fill in the blanks are in a csv editor like libre Office
#' or excel. You can also remake the experiment and specify the
#' specific column manually.
#' Remember that each row (sample) must have a unique combination
#' of values.
#' An extra column called "reverse" is made if there are paired data,
#' like +/- strand wig files.
#' @param dir Which directory / directories to create experiment from,
#' must be a directory with NGS data from your experiment. Will include
#' all files of file type specified by "types" argument. So do not mix
#' files from other experiments in the same folder!
#' @param exper Short name of experiment. Will be name used to load
#' experiment, and name shown when running \code{\link{list.experiments}}
#' @param saveDir Directory to save experiment csv file, default:
#' \code{ORFik::config()["exp"]}, which has default: "~/Bio_data/ORFik_experiments/".
#' Set to NULL if you don't want to save it to disc.
#' @param types Default \code{c("bam", "bed", "wig", "bigWig","ofst")},
#' which types of libraries to allow as NGS data.
#' @param txdb A path to TxDb (prefered) or gff/gtf (not adviced, slower)
#' file with transcriptome annotation for the organism.
#' @param fa A path to fasta genome/sequences used for libraries, remember the
#' file must have a fasta index too.
#' @param viewTemplate run View() on template when finished, default (FALSE).
#' Usually gives you a better view of result than using print().
#' @param organism character, default: "" (no organism set), scientific name
#' of organism. Homo sapiens, Danio rerio, Rattus norvegicus etc.
#' If you have a SRA metadata csv file, you can set this argument to
#' study$ScientificName[1], where study is the SRA metadata for all files
#' that was aligned.
#' @param assembly character, default: "" (no assembly set).
#' The genome assembly name, like GRCh38 etc. Useful to add if you want
#' detailed metadata of experiment analysis.
#' @param pairedEndBam logical FALSE, else TRUE, or a logical list of
#' TRUE/FALSE per library you see will be included (run first without and check
#' what order the files will come in) 1 paired end file, then two single will
#' be c(T, F, F). If you have a SRA metadata csv file, you can set this argument to
#' study$LibraryLayout == "PAIRED", where study is the SRA metadata for all files
#' that was aligned.
#' @param libtype character, default "auto". Library types,
#' must be length 1 or equal length of number of libraries.
#' "auto" means ORFik will try to guess from file names.
#' Example: RFP (Ribo-seq), RNA (RNA-seq), CAGE, SSU (TCP-seq 40S),
#' LSU (TCP-seq 80S).
#' @param stage character, default "auto". Developmental stage, tissue or
#' cell line, must be length 1 or equal length of number of libraries.
#' "auto" means ORFik will try to guess from file names.
#' Example: HEK293 (Cell line), Sphere (zebrafish stage), ovary (Tissue).
#' @param rep character, default "auto". Replicate numbering,
#' must be length 1 or equal length of number of libraries.
#' "auto" means ORFik will try to guess from file names.
#' Example: 1 (rep 1), 2 rep(2). Insert only numbers here!
#' @param condition character, default "auto". Library conditions,
#' must be length 1 or equal length of number of libraries.
#' "auto" means ORFik will try to guess from file names.
#' Example: WT (wild type), mutant, etc.
#' @param fraction character, default "auto". Fractionation of library,
#' must be length 1 or equal length of number of libraries.
#' "auto" means ORFik will try to guess from file names. This columns
#' is used to make experiment unique, if the other columns are not sufficient.
#' Example: cyto (cytosolic fraction), dmso (dmso treated fraction), etc.
#' @param author character, default "". Main author of experiment,
#' usually last name is enough. When printing will state "author et al" in info.
#' @param files character vector or data.table of library paths in dir.
#' Default: \code{findLibrariesInFolder(dir, types, pairedEndBam)}.
#' Do not touch unless you want to do some subsetting, it will automatically
#' remove files that are not of file format defined by 'type' argument.
#' Note that sorting on number that: 10 is before 2, so 1, 2, 10, is sorted as:
#' 1, 10, 2. If you want to fix this, you could update this argument with:
#' ORFik:::findLibrariesInFolder()[1,3,2] to get order back to 1,2,10 etc.
#' @param result_folder character, default NULL. The folder to output analysis
#' results like QC, count tables etc. By default the libFolder(df) folder is used,
#' the folder of first library in experiment. If you are making a new experiment
#' which is a collection of other experiments, set this to a new folder,
#' to not contaminate your other experiment directories.
#' @param runIDs character ids, usually SRR, ERR, or DRR identifiers, default is to search for any of these 3 in the filename by:
#'  \code{extract_run_id(files)}. They are optional.
#' @return a data.frame, NOTE: this is not a ORFik experiment,
#'  only a template for it!
#' @importFrom utils View
#' @export
#' @examples
#' # 1. Pick directory
#' dir <- system.file("extdata/Homo_sapiens_sample", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata/references/homo_sapiens",
#'                     "Homo_sapiens_dummy.gtf.db", package = "ORFik")
#' # 4. Pick fasta genome of organism
#' fa <- system.file("extdata/references/homo_sapiens",
#'                   "Homo_sapiens_dummy.fasta", package = "ORFik")
#' # 5. Set organism (optional)
#' org <- "Homo sapiens"
#'
#' # Create temple not saved on disc yet:
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               saveDir = NULL,
#'                               fa = fa, organism = org,
#'                               viewTemplate = FALSE)
#' ## Now fix non-unique rows: either is libre office, microsoft excel, or in R
#' template$X5[6] <- "heart" # here a dummy example, even though data is correct
#' # read experiment (if you set correctly)
#' df <- read.experiment(template)
#'
#' ## Default location of experiments is ORFik::config()["exp"]
#' # default_experiments_path <- ORFik::config()["exp"]
#' # exp_path <- file.path(default_experiments_path, paste0("exper", ".csv"))
#' # Save with: save.experiment(df, file = exp_path)
#' # Then you can simply load with read.experiment(exper),
#' # since you saved in the default directory
#'
#' ## Custom location (If you work in a team, use a shared folder)
#' # Remember to update ORFik::config() to ripple the effect through whole
#' # of ORFik if you want to use this as default
#' new_dir <- tempdir() # Here we just use tempdir
#' create.experiment(dir = dir, exper, txdb = txdb,
#'                   saveDir = new_dir, fa = fa, organism = org)
#' template_loaded <- read.experiment(exper,  new_dir)
#' # The csv template paths (from index 5) is equal to file paths of loaded exp
#' identical(template$X6[-seq(4)], filepath(template_loaded, "default"))
#'
#' @family ORFik_experiment
create.experiment <- function(dir, exper, saveDir = ORFik::config()["exp"],
                              txdb = "", fa = "", organism = "", assembly = "",
                              pairedEndBam = FALSE,
                              viewTemplate = FALSE,
                              types = c("bam", "bed", "wig", "bigWig","ofst"),
                              libtype = "auto", stage = "auto", rep = "auto",
                              condition = "auto", fraction = "auto",
                              author = "",
                              files = findLibrariesInFolder(dir, types, pairedEndBam),
                              result_folder = NULL,
                              runIDs = extract_run_id(files)) {
  if (!(is(files, "character") | is(files, "data.table")))
    stop("'files' must be of class character or data.table")
  stopifnot(is(exper, "character") & length(exper) == 1)
  stopifnot(is(runIDs, "character"))
  stopifnot(is(txdb, "character") & length(txdb) == 1)
  stopifnot(is(fa, "character") & length(fa) == 1)
  stopifnot(is(organism, "character") & length(organism) == 1)
  stopifnot(is(assembly, "character") & length(assembly) == 1)


  is_paired_end <- is(files, "data.table")
  runID_exists <- ifelse(all(runIDs != ""), TRUE, FALSE)
  exper <- sub("\\.csv$", "", exper)

  # Define data.frame size
  df <- data.frame(matrix(ncol = 6 + runID_exists + is_paired_end,
                          nrow = length(files) + 4))
  # set file paths
  df <- add_file_paths(df, files, is_paired_end, runIDs, runID_exists)
  ## Specify library information columns
  lib_metadata_columns <- c(c("libtype", "stage", "rep", "condition",
                              "fraction","filepath"),
                            c("reverse", "Run")[c(is_paired_end, runID_exists)])
  df <- guess_metadata_from_filepaths(df, files, libtype, stage, rep,
                                      condition, fraction, lib_metadata_columns)
  df <- add_reference_info(df, exper, txdb, fa, assembly, organism)

  # Additional information
  if (author != "") df[3, seq(5, 6)] <- c("author", author)
  if (!is.null(result_folder)) {
    df[1, seq(3,4)] <- c("results", result_folder)
  }
  df[is.na(df)] <- ""

  if (!is.null(saveDir)) {
    cbu.path <- "/export/valenfs/data/processed_data/experiment_tables_for_R/"
    if (dir.exists(cbu.path)) { # This will only trigger on CBU server @ UIB
      message("This is on internal CBU server, saving to preset directory:")
      message(cbu.path)
      saveDir <- cbu.path
    }
    save.experiment(df, pasteDir(saveDir, exper, ".csv"))
  }
  if (viewTemplate) View(df)
  return(df)
}

add_reference_info <- function(df, exper, txdb, fa, assembly, organism) {
  ## Reference assembly information
  df[1, seq(2)] <- c("name", exper)
  df[2, seq(2)] <- c("gff", txdb)
  df[3, seq(2)] <- c("fasta", fa)
  if (assembly != "") df[1, seq(5, 6)] <- c("assembly", assembly)
  if (organism != "") df[2, seq(5, 6)] <- c("organism", organism)
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

experiment_parse_list_info <- function(file, in.dir) {
  if (is(file, "character")) {
    if (length(file) != 1) stop("Experiment name must be single string!")
    if (file == "") stop("Experiment name is empty: ''")
    if (file_ext(file) == "") file <- paste0(file, ".csv")
    if (!file.exists(file)) {
      stopifnot(is.character(in.dir))
      file <- pasteDir(in.dir, file)
    }
    if (!file.exists(file)) { # This will only trigger on CBU server @ UIB
      cbu.path <- "/export/valenfs/data/processed_data/experiment_tables_for_R"
      if (file.exists(file.path(cbu.path, basename(file))))
        file <- file.path(cbu.path, basename(file))
    }

    info <- read.table(file, sep = ",", nrows = 3, stringsAsFactors = FALSE)
    listData <- read.csv2(file, skip = 3, header = TRUE, sep = ",",
                          stringsAsFactors = FALSE)
  } else if(is(file, "data.frame")) {
    if (nrow(file) < 5) stop("For data.frame input, file must have > 4 rows")
    info <- file[seq(3),]
    listData <- file[-seq(4),]
    colnames(listData) <- file[4,]
  } else stop("file must be either character or data.frame template")
  listData <- cbind(listData, index = as.integer(seq.int(nrow(listData))))

  return(list(listData = listData, info = info))
}

add_file_paths <- function(df, files, is_paired_end, runIDs,
                           runID_exists) {
  file_dt <- files
  files <- if(is_paired_end) { file_dt$forward } else {file_dt}
  df[5:(5+length(files)-1), 6] <- files
  if (is_paired_end) df[5:(5+length(files)-1), 7] <- file_dt$reverse
  if (runID_exists) df[5:(5+length(files)-1), ncol(df)] <- runIDs
  return(df)
}
