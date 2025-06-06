#' Download sra toolkit
#'
#' Currently supported for Linux (64 bit centos and ubunutu is tested to work)
#' and Mac-OS(64 bit). If other linux distro, centos binaries will be used.
#' @param folder default folder, "~/bin"
#' @param version a string, default "2.11.3"
#' @return path to fastq-dump in sratoolkit
#' @importFrom utils untar
#' @references https://ncbi.github.io/sra-tools/fastq-dump.html
#' @family sra
#' @export
#' @examples
#' # install.sratoolkit()
#' ## Custom folder and version (not adviced)
#' folder <- "/I/WANT/IT/HERE/"
#' # install.sratoolkit(folder, version = "2.10.9")
#'
install.sratoolkit <- function(folder = "~/bin", version = "2.11.3") {
  if (.Platform$OS.type != "unix")
    stop("sratoolkit is not currently supported for windows by ORFik,
         download manually or use WSL (windows subsystem linux)")
  folder <- path.expand(folder)
  is_linux <- Sys.info()[1] == "Linux" # else it is mac
  os <- if(!is_linux) {
    "mac64"
    } else {
     is_ubuntu <- grepl("^Ubuntu", osVersion, ignore.case = TRUE)
     is_centos <- grepl("^cent", osVersion, ignore.case = TRUE)
     found_exact_match <- sum(c(is_ubuntu, is_centos)) == 1
     if (!found_exact_match | is_centos) {
       "centos_linux64"
     } else "ubuntu64"
    }
  # TODO; Check if ubuntu compliation is needed for safer download ->
  #length(grep("Ubuntu", system("cat /etc/*release", intern = TRUE)[1])) == 1

  path.final <- paste0(folder, "/sratoolkit.", version, "-", os)
  path.final <- paste0(path.final, "/bin/fastq-dump")
  if (file.exists(path.final)) {
    message(paste("Using fastq-dump at location:",
                  path.final))
    return(path.final)
  }
  message("Downloading and configuring SRA-toolkit for you,
          this is done only once!")

  url <- paste0("https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/", version, "/")
  url <- paste0(url, "sratoolkit.", version, "-", os, ".tar.gz")

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

#' Download read libraries from SRA
#'
#' Multicore version download, see documentation for SRA toolkit for more information.
#' @param info character vector of only SRR numbers or
#' a data.frame with SRA metadata information including the SRR numbers in a column called
#' "Run" or "SRR". Can be SRR, ERR or DRR numbers.
#' If only SRR numbers can not rename, since no additional information is given.
#' @param outdir directory to store runs,
#' files are named by default (rename = TRUE) by information from SRA metadata table,
#' if (rename = FALSE) named according to SRR numbers.
#' @param rename logical or character, default TRUE (Auto guess new names). False: Skip
#' renaming. A character vector of equal size as files wanted can also be given.
#' Priority of renaming from
#' the metadata is to check for unique names in the LibraryName column,
#' then the sample_title column if no valid names in LibraryName.
#' If new names found and still duplicates, will
#' add "_rep1", "_rep2" to make them unique. If no valid names, will not
#' rename, that is keep the SRR numbers, you then can manually rename files
#' to something more meaningful.
#' @param fastq.dump.path path to fastq-dump binary, default: path returned
#' from install.sratoolkit()
#' @param settings a string of arguments for fastq-dump,
#' default: paste("--gzip", "--skip-technical", "--split-files")
#' @param subset an integer or NULL, default NULL (no subset). If defined as
#' a integer will download only the first n reads specified by subset. If subset is
#' defined, will force to use fastq-dump which is slower than ebi download.
#' @param compress logical, default TRUE. Download compressed files ".gz".
#' @param use.ebi.ftp logical, default: is.null(subset). Use ORFiks much faster download
#' function that only works when subset is null,
#' if subset is defined, it uses fastqdump, it is slower but supports subsetting.
#' Force it to use fastqdump by setting this to FALSE.
#' @param ebiDLMethod character, default "auto". Which download protocol
#' to use in download.file when using ebi ftp download. Sometimes "curl"
#' is might not work (the default auto usually), in those cases use wget.
#' See "method" argument of ?download.file, for more info.
#' @param timeout 5000, how many seconds before killing download if
#' still active? Will overwrite global option until R session is closed.
#' Increase value if you are on a very slow connection or downloading a large dataset.
#' @param BPPARAM how many cores/threads to use? default: bpparam().
#' To see number of threads used, do \code{bpparam()$workers}
#' @return a character vector of download files filepaths
#' @references https://ncbi.github.io/sra-tools/fastq-dump.html
#' @family sra
#' @export
#' @examples
#' SRR <- c("SRR453566") # Can be more than one
#' \donttest{
#' ## Simple single SRR run of YEAST
#' outdir <- tempdir() # Specify output directory
#' # Download, get 5 first reads
#' #download.SRA(SRR, outdir, rename = FALSE, subset = 5)
#'
#' ## Using metadata column to get SRR numbers and to be able to rename samples
#' outdir <- tempdir() # Specify output directory
#' info <- download.SRA.metadata("SRP226389", outdir) # By study id
#' ## Download, 5 first reads of each library and rename
#' #files <- download.SRA(info, outdir, subset = 5)
#' #Biostrings::readDNAStringSet(files[1], format = "fastq")
#'
#' ## Download full libraries of experiment
#' ## (note, this will take some time to download!)
#' #download.SRA(info, outdir)
#' }
download.SRA <- function(info, outdir, rename = TRUE,
                         fastq.dump.path = install.sratoolkit(),
                         settings =  paste("--skip-technical", "--split-files"),
                         subset = NULL,
                         compress = TRUE,
                         use.ebi.ftp = is.null(subset),
                         ebiDLMethod = "auto",
                         timeout = 5000,
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
  if (is.null(SRR) | (length(SRR) == 0))
    stop("Could not find SRR numbers in 'info'")

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(outdir)) stop("Could not create output directory, is the disc not mounted?")

  settings <- paste("--outdir", outdir, settings)
  if (!is.null(subset)) {
    if(!is.numeric(subset)) stop("subset must be numeric if not NULL")
    subset <- as.integer(subset)
    settings <- paste(settings, "-X", subset)
  } else if (use.ebi.ftp){
    files <- download.ebi(info, outdir, rename, ebiDLMethod, timeout, BPPARAM)
    if (length(files) > 0) return(files)
    message("Checking for fastq files using fastq-dump")
  }
  if (compress) {
    settings <- paste(settings, "--gzip")
  }
  fastq.dump <- fastq.dump.path
  message("Starting download of SRA runs:")
  BiocParallel::bplapply(SRR, function(i, fastq.dump, settings) {
    message(i)
    system(command = paste(fastq.dump, i, settings),
           wait = TRUE)
  }, fastq.dump = fastq.dump, settings = settings, BPPARAM = BPPARAM)

  search_it <- ifelse(compress, "\\.fastq\\.gz$", "\\.fastq$")
  files <- unlist(lapply(SRR, function(S)
    dir(outdir, paste0(S, ".*", search_it), full.names = TRUE))
  )

  valid <- TRUE
  if (length(files) == 0) valid <- FALSE
  any.paired <- length(grep("_[2]\\.fastq\\.gz", files))
  # TODO: validate that this will work in download of mixed
  paired <- ifelse(any.paired,
                   length(grep("_[1-2]\\.fastq\\.gz", files)),
                   0)

  if (length(SRR) != (paired/2 + length(files) - paired))
    valid <- FALSE
  if (!valid) {
    warning("Some of the files specified was not downloaded,",
            " are you behind a strict firewall?")
    message("If only few files remaining, subset to those SRR numbers and run again")
  }

  if (is.logical(rename)) { # Renaming
    # Set to false if no metadata
    if (is.character(info) & rename) {
      rename <- FALSE
      warning("rename = TRUE, but no metadata given. Can not rename!")
    } else if (rename) files <- rename.SRA.files(files, info)
  } else { # else names were assigned manually
    files <- rename.SRA.files(files, rename)
  }
  return(files)
}

#' Rename SRA files from metadata
#'
#' @param files a character vector, with full path to all the files
#' @param new_names a character vector of new names or
#' a data.table with metadata to use to rename (usually from SRA metadata).
#' Priority of renaming from
#' the metadata is to check for unique names in the LibraryName column,
#' then the sample_title column if no valid names in LibraryName.
#' If found and still duplicates, will
#' add "_rep1", "_rep2" to make them unique. Paired end data will get a extension
#' of _p1 and _p2. If no valid names, will not
#' rename, that is keep the SRR numbers, you then can manually rename files
#' to something more meaningful.
#' @return a character vector of new file names
#' @family sra
#' @keywords internal
rename.SRA.files <- function(files, new_names) {
  info <- NULL # Set to default
  if (!is.character(new_names)) { # Then auto-guess from meta data
    message("Auto-guessing new names from metadata, check that they are valid")
    info <- new_names
    new_names <- NULL

    valid_libraryName_column <- !is.null(info$LibraryName) &
      !any(is.na(info$LibraryName)) & !any("" %in% info$LibraryName)
    if (valid_libraryName_column) {
        new_names <- info$LibraryName
    }
    not_defined_yet <- is.null(new_names)
    valid_sample_column <- !is.null(info$sample_title) &
      !any(is.na(info$sample_title)) & !any("" %in% new_names)
    if (not_defined_yet & valid_sample_column) {
      new_names <- info$sample_title
      new_names <- gsub(".*: ", "", new_names)
      new_names <- gsub(";.*", "", new_names)
    }
    libStrat <- info$LibraryStrategy
    libSelect <- info$LibrarySelection
    libStrat_usable <- !is.null(libStrat) &
      !any(is.na(libStrat)) & !all(c("") %in% libStrat) &
      !all(c("OTHER") %in% libStrat) & !all(c("other") %in% libStrat) &
      !all(c("unspecified") %in% libStrat)

    libSelect_usable <- !is.null(libSelect) &
      !any(is.na(libSelect)) & !all(c("") %in% libSelect) &
      !all(c("OTHER") %in% libSelect) & !all(c("other") %in% libSelect) &
      !all(c("unspecified") %in% libSelect)

    if (!is.null(new_names)) {
      new_names <- paste0(toupper(substr(new_names, 1, 1)),
                           substr(new_names, 2, nchar(new_names)))
    }
  }

  if (any(duplicated(new_names))) {
    new_names <- make.unique(new_names, sep = "_rep")
  }

  if (!is.null(new_names)) {
    message("Renaming files:")
    if (!is.null(info)) { # If metadata given, update if paired end
      if (any("PAIRED" %in% info$LibraryLayout)) {
        new_names <- lapply(seq_along(info$LibraryLayout),
                            function(x) if(info$LibraryLayout[x] == "PAIRED") {
                              c(paste0(new_names[x], "_p1"),
                                paste0(new_names[x],"_p2"))
                            } else new_names[x])
        new_names <- unlist(new_names)
      }
    }

    if (length(new_names) != length(files))
      stop("Length of files and new_names to rename by is not equal!",
           " If manual assign of paired end name, repeat each element twice!")
    new_names <- gsub(",", "_", new_names)
    new_names <- gsub(" |\\(|\\)", "_", new_names)
    new_names <- gsub("__", "_", new_names)
    new_names <- gsub("/", "", new_names)
    is_gzipped <- grep("\\.fastq\\.gz", files)

    new_names <- paste0(dirname(files), "/", basename(new_names), ".fastq")

    new_names[is_gzipped] <- paste0(new_names, ".gz")
    for (i in seq(length(files))) {
      file.rename(files[i], new_names[i])
    }
  } else {
    warning("Did not find a way for valid renaming, returning without renaming!")
    return(files)
  }
  names(new_names) <- basename(files)
  return(new_names)
}

#' Faster download of fastq files
#'
#' Uses ftp download from vol1 drive on EBI ftp server,
#'  for faster download of ERR, SRR or DRR files.
#' But does not support subsetting or custom settings of files!
#' @inheritParams download.SRA
#' @return character, full filepath of downloaded  files
#' @family sra
#' @keywords internal
download.ebi <- function(info, outdir, rename = TRUE,
                         ebiDLMethod = "auto", timeout = 5000,
                         BPPARAM = bpparam()) {

  study <- NULL
  # If character presume SRR, if not check for column Run or SRR
  SRR <- if (is.character(info)) { # if character
    info
  } else { # else metadata
    # Check if study is specified
    if (length(unique(info$BioProject)) == 1)
      study <- info$BioProject[1]
    if (is.null(info$Run)) { # If not called Run
      info$SRR
    } else  { # If called Run
      info$Run
    }
  }
  if (is.null(SRR) | (length(SRR) == 0))
    stop("Could not find SRR numbers in 'info'")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  urls <- find_url_ebi(SRR, study = study)
  if (length(urls) == 0) {
    message("None of the Fastq files specified found on ebi")
    return(NULL)
  } else if (length(urls) < length(SRR)) {
    message("Not all fastq files specified found on ebi")
    return(NULL)
  }

  files <- file.path(outdir, basename(urls))
  message("Starting download of EBI runs:")
  method <- ebiDLMethod
  withr::local_options(timeout = timeout)
  BiocParallel::bplapply(urls, function(i, outdir, method) {
    message(i)
    download.file(i, destfile = file.path(outdir, basename(i)),
                  method = method, quiet = TRUE)
  }, outdir = outdir, method = method, BPPARAM = BPPARAM)

  if (is.logical(rename)) {
    # Set to false if not metadata
    if (is.character(info) & rename) {
      rename <- FALSE
      warning("rename = TRUE, but no metadata given. Can not rename!")
    } else if (rename) files <- rename.SRA.files(files, info)
  } else { # else manual assign names
    files <- rename.SRA.files(files, rename)
  }
  return(files)
}

#' Locates and check if fastq files exists in ebi
#'
#' Look for files in ebi file servers,
#' Paired end and single end fastq files.\cr
#' Fastq ftp url: \code{ftp://ftp.sra.ebi.ac.uk/vol1/fastq}\cr
#' SRA   ftp url: \code{ftp://ftp.sra.ebi.ac.uk/vol1/srr}\cr
#' Fastq ASCP url: \code{era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq}\cr
#' SRA   ASCP url: \code{era-fasp@fasp.sra.ebi.ac.uk:vol1/srr}\cr
#' @param SRR character, SRR, ERR or DRR numbers.
#' @param stop.on.error logical FALSE, if TRUE will stop
#'  if all files are not found. If FALSE returns empty character vector if error
#'  is catched.
#' @param study default NULL, optional PRJ (study id) to speed up search
#' for URLs.
#' @param ebi_file_format character, format of run download, default is fastq (ftp):
#'  \code{c("fastq_ftp", "sra_ftp")[1]}
#' @param convert_to_ascp logical, default FALSE. If TRUE use server:
#' \code{era-fasp@fasp.sra.ebi.ac.uk:}
#' @return full url to fastq files, same length as input
#' (2 urls for paired end data). Returns empty character() if all
#' files not found.
#' @export
#' @examples
#' # Test the 3 ways to get fastq files from EBI
#' # Both single end and paired end data
#'
#' # Most common: SRR(3 first)/0(2 last)/whole
#' # Single
#' ORFik:::find_url_ebi("SRR10503056")
#' # Paired
#' ORFik:::find_url_ebi("SRR10500056")
#'
#' # less common: SRR(3 first)/00(1 last)/whole
#' # Single
#' #ORFik:::find_url_ebi("SRR1562873")
#' # Paired
#' #ORFik:::find_url_ebi("SRR1560083")
#' # least common SRR(3 first)/whole
#' # Single
#' #ORFik:::find_url_ebi("SRR105687")
#' # Paired
#' #ORFik:::find_url_ebi("SRR105788")
find_url_ebi <- function(SRR, stop.on.error = FALSE, study = NULL,
                         ebi_file_format = c("fastq_ftp", "sra_ftp")[1],
                         convert_to_ascp = FALSE) {
  message("Finding optimal download urls from ebi...")
  if (!is.null(study)) {
    SRR <- study
  }
  return(find_url_ebi_safe(SRR, stop.on.error = stop.on.error,
                           ebi_file_format = ebi_file_format,
                           convert_to_ascp = convert_to_ascp))
}

#' Find URL for EBI fastq files
#'
#' Safer version
#' @inheritParams find_url_ebi
#' @param accession character: (PRJ, SRP, ERP, DRP, SRX, SRR, ERR,..). For studies or samples,
#' it returns all runs per study or sample.
#' @param SRR character, which SRR numbers to subset by (can also be ERR or DRR numbers)
#' @return character (1 element per SRR number)
#' @keywords internal
find_url_ebi_safe <- function(accession, SRR = NULL, stop.on.error = FALSE,
                              ebi_file_format = c("fastq_ftp", "sra_ftp")[1],
                              convert_to_ascp = FALSE) {
  stopifnot(ebi_file_format %in% c("fastq_ftp", "sra_ftp"))
  a <- data.table()
  for (i in accession) {
    search_url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                         i, "&result=read_run&fields=run_accession,", ebi_file_format)
    temp <- suppressWarnings(temp <- fread(search_url, header = TRUE))
    a <- rbindlist(list(a, temp))
  }
  if (!is.null(SRR)) {
    if (!all(SRR %in% a$run_accession)) {
      if (stop.on.error) stop("Study does not contain some of the SRR numbers given!")
      return(character())
    }
    a <- a[run_accession %in% SRR,]
  }
  paths <- a[, colnames(a) == ebi_file_format, with = FALSE][[1]]
  paths <- unlist(strsplit(paths, ";"))
  if (convert_to_ascp) {
    paths <- sub("ftp.sra.ebi.ac.uk/", "era-fasp@fasp.sra.ebi.ac.uk:", paths)
  }
  return(paths)
}

#' Extract SRR/ERR/DRR run IDs from string
#'
#' @param x character vector to search through.
#' @param search the regex search, default: \code{"(SRR[0-9]+|DRR[0-9]+|ERR[0-9]+)"}
#' @param only_valid logical, default FALSE. If TRUE, return only the hits.
#' @return a character vector of run accepted run ids according to search,
#' if only_valid named character vector for which indices are returned
#' @examples
#' search <- c("SRR1230123_absdb", "SRR1241204124_asdasd", "asd_ERR1231230213",
#'  "DRR12412412_asdqwe", "ASDASD_ASDASD", "SRRASDASD")
#' ORFik:::extract_run_id(search)
#' ORFik:::extract_run_id(search, only_valid = TRUE)
extract_run_id <- function(x, search = "(SRR[0-9]+|DRR[0-9]+|ERR[0-9]+)", only_valid = FALSE) {
  hits <- gsub(paste0(".*", search), "\\1", gsub(paste0(search, ".*"), "\\1", x))
  match <- grep(search, hits)
  if (only_valid) {
    hits <- hits[match]
    names(hits) <- match
  } else hits[!(seq_along(hits) %in% match)] <- ""
  return(hits)
}
