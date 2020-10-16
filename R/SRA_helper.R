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
#' @param rename logical, default TRUE. Priority of renaming from
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
#' outdir <- tempdir() # Specify output directory
#' info <- download.SRA.metadata("SRP226389", outdir) # By study id
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
#' @param remove.invalid logical, default TRUE. Remove Runs with 0 reads (spots)
#' @return a data.table of the opened file
#' @importFrom utils download.file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom xml2 read_xml
#' @importFrom xml2 as_list
#' @export
#' @examples
#' ## Originally on SRA
#' # outdir <- tempdir() # Specify output directory
#' # download.SRA.metadata("SRP226389", outdir)
#' ## ORiginally on ENA
#' # download.SRA.metadata("ERP116106", outdir)
download.SRA.metadata <- function(SRP, outdir, remove.invalid = TRUE) {
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

  msg <- paste("Found Runs with 0 reads (spots) in metadata, will not be able
              to download the run/s:", file[spots == 0,]$Run)
  if (any(file$spots == 0)) {
    warning(msg)
    if (remove.invalid) {
      warning("Removing invalid Runs from final metadata list")
      file <- file[spots > 0,]
    }
  }


  if (nrow(file) == 0) {
    warning(paste("No valid runs found from experiment:", SRP))
    return(file)
  } else {
    if ("sample_title" %in% colnames(file)) return(file)

    file <- file[, -c("ReleaseDate", "LoadDate", "download_path", "RunHash", "ReadHash", "Consent")]
    # Download xml and add more data
    url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=xml&term="
    url <- paste0(url, SRP)
    destfile_xml <- paste0(outdir, "/SraRunInfo_", SRP, ".xml")
    download.file(url, destfile = destfile_xml)
    a <- xml2::read_xml(destfile_xml)
    a <- xml2::as_list(a)

    dt <- data.table()
    for(i in seq_along(a$EXPERIMENT_PACKAGE_SET)) {
      dt <- rbind(dt,
                  cbind(unlist(a$EXPERIMENT_PACKAGE_SET[i]$EXPERIMENT_PACKAGE$SAMPLE$TITLE),
                        unlist(a$EXPERIMENT_PACKAGE_SET[i]$EXPERIMENT_PACKAGE$RUN_SET$RUN$IDENTIFIERS$PRIMARY_ID))
      )
    }
    colnames(dt) <- c("sample_title", "Run")
    dt <- dt[Run %in% file$Run]
    if (length(dt) > 0) {
      file <- data.table::merge.data.table(file, dt, by = "Run")
    }
    # Remove xml and keep runinfo
    file.remove(destfile_xml)
    fwrite(file, destfile)
  }
  return(file)
}

#' Rename SRA files from metadata
#'
#' @param files a character vector, all the files
#' @param new_names a character vector of new names or
#' a data.table with metadata to use to rename. Priority of renaming from
#' the metadata is to check for unique names in the LibraryName column,
#' then the sample_title column if no valid names in LibraryName.
#' If found and still duplicates, will
#' add "_rep1", "_rep2" to make them unique. If no valid names, will not
#' rename, that is keep the SRR numbers, you then can manually rename files
#' to something more meaningful.
#' @return a character vector of new file names
rename.SRA.files <- function(files, new_names) {
  if (is.character(new_names)) {
    # Do nothing
  } else { # Find from meta data
    info <- new_names
    if (!is.null(info$LibraryName) & !any(is.na(info$LibraryName))) {
        new_names <- info$LibraryName
        message("Renaming files:")
    } else if (!is.null(info$sample_title) & !any(is.na(info$sample_title))) {
       new_names <- info$sample_title
       new_names <- gsub(".*: ", "", new_names)
       new_names <- gsub(";.*", "", new_names)
       message("Renaming files:")
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
  return(new_names)
}
