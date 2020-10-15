#' Download sra toolkit
#' @param folder default folder, "~/bin"
#' @param version a string, default "2.10.8"
#' @return path to fastq-dump in sratoolkit
#' @importFrom utils untar
#' @references https://ncbi.github.io/sra-tools/fastq-dump.html
#' @export
#' @examples
#' # install.sratoolkit()
#' ## Custom version
#' # install.sratoolkit(version = "2.10.7")
#'
install.sratoolkit <- function(folder = "~/bin", version = "2.10.8") {
  if (.Platform$OS.type != "unix")
    stop("fastp does not work on Windows, try RSubread")
  folder <- path.expand(folder)
  is_linux <- Sys.info()[1] == "Linux" # else it is mac


  path.final <- ifelse(is_linux,
                       paste0(folder, "/sratoolkit.", version, "-centos_linux64"),
                       paste0(folder, "/sratoolkit.", version, "-mac64"))
  path.final <- paste0(path.final, "/bin/fastq-dump")
  if (file.exists(path.final)) {
    message(paste("Using fastq-dump at location:",
                  path.final))
    return(path.final)
  }
  message("Downloading and configuring SRA-toolkit for you,
          this is done only once!")

  url <- paste0("https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/", version, "/")
  url <- ifelse(is_linux,
                paste0(url, "sratoolkit.", version, "-centos_linux64.tar.gz"),
                paste0(url, "sratoolkit.", version, "-mac64.tar.gz"))
  path <- paste0(folder, "/sratoolkit.tar.gz")

  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  utils::download.file(url, destfile = path)
  untar(path, exdir = folder)

  # Update access rights
  system(paste("chmod a+x", path.final))
  # Make config file, will give ignorable seqmentation faul warning
  message("Ignore the following config warning: SIGNAL - Segmentation fault ")
  conf <- suppressWarnings(system(paste0(dirname(path.final), "/vdb-config -i"),
                                  intern = TRUE))

  return(path.final)
}

#' Download SRR libraries from SRA
#'
#' Multicore version download, see documentation for SRA toolkit for more information.
#' @param info character vector of only SRR numbers or
#' a data.frame with SRA metadata information including the SRR numbers in a column called
#' "Run" or "SRR".
#' If only SRR numbers can not rename, since no additional information is given.
#' @param outdir a string, default: cbu server
#' @param rename logical, default TRUE. Rename files according to info, only if info has a column called
#' `Experiment Title`
#' @param fastq.dump.path path to fastq-dump binary, default: path returned
#' from install.sratoolkit()
#' @param settings a string of arguments for fastq-dump,
#' default: paste("--gzip", "--skip-technical", "--split-files")
#' @param subset a numeric or NULL, default NULL (no subset). If defined as
#' a numeric will download only the first n reads specified by subset.
#' @param BPPARAM how many cores/threads to use? default: bpparam().
#' To see number of threads used, do \code{bpparam()$workers}
#' @return a character vector of download files filepaths
#' @references https://ncbi.github.io/sra-tools/fastq-dump.html
#' @export
#' @examples
#' \dontrun{
#' ## Simple single SRR run of YEAST
#' SRR <- c("SRR453566") # Can be more than one
#' outdir <- tempdir() # Specify output directory
#' # Download, get 5 first reads
#' download.SRA(SRR, outdir, subset = 5)
#'
#' ## Using metadata column to get SRR numbers and to be able to rename samples
#' info <- download.SRA.metadata("SRP226389") # By study id
#' outdir <- tempdir() # Specify output directory
#' # Download, 5 first reads of each library and rename
#' download.SRA(info, outdir, subset = 5)
#' }
download.SRA <- function(info, outdir, rename = TRUE,
                         fastq.dump.path = install.sratoolkit(),
                         settings =  paste("--gzip", "--skip-technical", "--split-files"),
                         subset = NULL,
                         BPPARAM = bpparam()) {
  # If character presume SRR, if not check for column Run or SRR
  SRR <- if (is.character(info)) { # if character
    info
  } else { # else metadata
    if (is.null(info$Run)) { # If not called Run
      info$SRR
    } else  { # If called Run
      info$Run
    }
  }

  if (is.null(SRR)) stop("Could not find SRR numbers in 'info'")



  if (is.null(SRR)) stop("No SRR column found in info!")

  fastq.dump <- fastq.dump.path
  settings <- paste("--outdir", outdir, settings)
  if (!is.null(subset)) {
    if(!is.numeric(subset)) stop("subset must be numeric if not NULL")
    settings <- paste(settings, "-X", subset)
  }
  message("Starting download of SRA runs:")
  BiocParallel::bplapply(SRR, function(i, fastq.dump, settings) {
    message(i)
    system(command = paste(fastq.dump, i, settings),
           wait = TRUE)
  }, fastq.dump = fastq.dump, settings = settings, BPPARAM = BPPARAM)


  files <- unlist(lapply(SRR, function(S)
    dir(outdir, paste0(S, ".*", "\\.fastq\\.gz"), full.names = TRUE))
  )

  rename <- ifelse(!is.character(info) & rename, TRUE, FALSE)
  if (rename) {
    files <- rename.SRA.files(files, info)
  }
  return(files)
}

#' Downloads metadata from SRA
#'
#' @param SRP a string, a study ID as either the SRP, ERP or PRJ of the study,
#' examples would be "SRP226389" or "ERP116106".
#' @param outdir directory to save file,
#' The file will be called "SraRunInfo_SRP.csv", where SRP is
#' the SRP argument.
#' The directory will be created if not existing.
#' @return a data.table of the opened file
#' @importFrom utils download.file
#' @importFrom data.table fread
#' @export
#' @examples
#' ## Originally on SRA
#' # outdir <- tempdir() # Specify output directory
#' # download.SRA.metadata("SRP226389", outdir)
#' ## ORiginally on ENA
#' # download.SRA.metadata("ERP116106", outdir)
download.SRA.metadata <- function(SRP, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  destfile <- paste0(outdir, "/SraRunInfo_", SRP, ".csv")
  if (file.exists(destfile)) {
    message(paste("Existing metadata file found in dir:", outdir, "will not download"))
  } else {
    url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term="
    url <- paste0(url, SRP)
    download.file(url, destfile = destfile)
  }
  file <- fread(destfile)
  if (nrow(file) == 0) {
    warning(paste("No runs found from experiment:", SRP))
  } else {
    file <- file[, -c("ReleaseDate", "LoadDate", "download_path", "RunHash", "ReadHash")]
  }
  return(file)
}

#' Rename SRA files from metadata
#'
#' @param files a character vector, all the files
#' @param new_names a character vector of new names or
#' a data.table with metadata to use to rename
#' @return a character vector of new file names
rename.SRA.files <- function(files, new_names) {
  if (is.character(new_names)) {
    # Do nothing
  } else { # Find from meta data
    info <- new_names
    if (!is.null(info$`Experiment Title`)) {
      message("Renaming files:")
      new_names <- info$`Experiment Title`
      new_names <- gsub(".*: ", "", new_names)
      new_names <- gsub(";.*", "", new_names)
      new_names <- paste0(dirname(files), "/",new_names, ".fastq.gz")
    } else if (!is.null(info$LibraryName)) {
      if (!any(is.na(info$LibraryName))) {
        new_names <- info$LibraryName
        message("Renaming files:")
      } else new_names <- NULL
    } else { # No valid naming, set to NULL
      new_names <- NULL
    }
  }

  if (any(duplicated(new_names))) {
    new_names <- make.unique(new_names, sep = "_rep")
  }

  if (!is.null(new_names)) {
    if (length(new_names) != length(files))
      stop("Length of files and new_names to rename by is not equal!")

    new_names <- gsub(" |\\(|\\)", "_", new_names)
    new_names <- gsub("__", "_", new_names)
    for (i in seq(length(files))) {
      file.rename(files[i], new_names[i])
    }
  } else {
    warning("Did not find a way for valid renaming, returning without renaming!")
    return(files)
  }
  return(new_names)
}
