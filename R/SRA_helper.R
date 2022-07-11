#' Download sra toolkit
#'
#' Currently supported for Linux (64 bit centos and ubunutu is tested to work)
#' and Mac-OS(64 bit)
#' @param folder default folder, "~/bin"
#' @param version a string, default "2.10.9"
#' @return path to fastq-dump in sratoolkit
#' @importFrom utils untar
#' @references https://ncbi.github.io/sra-tools/fastq-dump.html
#' @family sra
#' @export
#' @examples
#' # install.sratoolkit()
#' ## Custom folder and version
#' folder <- "/I/WANT/IT/HERE/"
#' # install.sratoolkit(folder, version = "2.10.7")
#'
install.sratoolkit <- function(folder = "~/bin", version = "2.10.9") {
  if (.Platform$OS.type != "unix")
    stop("sratoolkit is not currently supported for windows by ORFik, download manually")
  folder <- path.expand(folder)
  is_linux <- Sys.info()[1] == "Linux" # else it is mac
  # TODO; Check if ubuntu compliation is needed for safer download ->
  #length(grep("Ubuntu", system("cat /etc/*release", intern = TRUE)[1])) == 1

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
  url <- paste0(url, "sratoolkit.", version)
  url <- ifelse(is_linux,
                paste0(url, "-centos_linux64.tar.gz"),
                paste0(url, "-mac64.tar.gz"))
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
#' @param outdir a string, default: cbu server
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
#' #download.SRA(SRR, outdir, subset = 5)
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

  settings <- paste("--outdir", outdir, settings)
  if (!is.null(subset)) {
    if(!is.numeric(subset)) stop("subset must be numeric if not NULL")
    subset <- as.integer(subset)
    settings <- paste(settings, "-X", subset)
  } else if (use.ebi.ftp){
    files <- download.ebi(info, outdir, rename, ebiDLMethod, BPPARAM)
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

#' Downloads metadata from SRA
#'
#' Given a experiment identifier, query information from different locations of SRA
#' to get a complete metadata table of the experiment. It first finds Runinfo for each
#' library, then sample info,
#' if pubmed id is not found searches for that and searches for author through pubmed.
#'
#' A common problem is that the project is not linked to an article, you will then not
#' get a pubmed id.
#' @param SRP a string, a study ID as either the PRJ, SRP, ERP, DRPor GSE of the study,
#' examples would be "SRP226389" or "ERP116106". If GSE it will try to convert to the SRP
#' to find the files. The call works as long the runs are registered on the efetch server,
#' as their is a linked SRP link from bioproject or GSE. Example which fails is "PRJNA449388",
#' which does not have a linking like this.
#' @param outdir directory to save file, default: tempdir().
#' The file will be called "SraRunInfo_SRP.csv", where SRP is
#' the SRP argument. We advice to use bioproject IDs "PRJNA...".
#' The directory will be created if not existing.
#' @param remove.invalid logical, default TRUE. Remove Runs with 0 reads (spots)
#' @param auto.detect logical, default FALSE. If TRUE, ORFik will add additional columns:\cr
#' LIBRARYTYPE: (is this Ribo-seq or mRNA-seq, CAGE etc), \cr
#' REPLICATE: (is this replicate 1, 2 etc),\cr
#' STAGE: (Which time point, cell line or tissue is this, HEK293, TCP-1, 24hpf etc),\cr
#' CONDITION: (is this Wild type control or a mutant etc).\cr
#' These values are only qualified guesses from the metadata, so always double check!
#' @param abstract character, default "printsave". If abstract for project exists,
#' print and save it (save the file to same directory as runinfo).
#' Alternatives: "print", Only print first time downloaded,
#' will not be able to print later.\cr
#' save" save it, no print\cr
#' "no" skip download of abstract
#' @return a data.table of the metadata, 1 row per sample,
#'  SRR run number defined in Run column.
#' @importFrom utils download.file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom xml2 read_xml
#' @importFrom xml2 as_list
#' @references doi: 10.1093/nar/gkq1019
#' @family sra
#' @export
#' @examples
#' ## Originally on SRA
#' download.SRA.metadata("SRP226389")
#' ## Now try with auto detection (guessing additional library info)
#' ## Need to specify output dir as tempfile() to re-download
#' #download.SRA.metadata("SRP226389", tempfile(), auto.detect = TRUE)
#' ## Originally on ENA (RCP-seq data)
#' # download.SRA.metadata("ERP116106")
#' ## Originally on GEO (GSE) (save to directory to keep info with fastq files)
#' # download.SRA.metadata("GSE61011", "/path/to/fastq.folder/")
download.SRA.metadata <- function(SRP, outdir = tempdir(), remove.invalid = TRUE,
                                  auto.detect = FALSE, abstract = "printsave") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  destfile <- paste0(outdir, "/SraRunInfo_", SRP, ".csv")
  abstract_destfile <- paste0(outdir, "/abstract_", SRP, ".csv")

  if (file.exists(destfile)) {
    message(paste("Existing metadata file found in dir:", outdir, " ,will not download"))
    if (abstract %in% c("print", "printsave")) { # Print abstract if wanted
      if (file.exists(abstract_destfile)) {
        print("Study abstract:")
        print(read.table(abstract_destfile, header = TRUE)$abstract)
      }
    }
    return(fread(destfile))
  }
  is_GSE <- length(grep("GSE", x = SRP)) == 1
  if (is_GSE) { # Find SRP from GSE
    message("GSE inserted, trying to find SRP from the GSE")
    url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
    url <- paste0(url, SRP, "&form=text")
    # Get GSE info
    a <- fread(url, sep = "!", header = FALSE, col.names = c("row.info", "info"))
    # Now find SRP from GSE
    SRP_line <- grepl("Series_relation = ", a$info)
    if (length(SRP_line) == 0) stop("GSE does not have a recorded Bioproject or SRP; check that it is correct!")
    b <- a[SRP_line,]$info
    d <- b[grepl("=SRP", b)]
    if (length(d) == 0) {
      d <- b[grepl("= BioProject:", b)]
      SRP <- gsub(".*/", replacement = "", d)
    } else SRP <- gsub(".*term=", replacement = "", d)
    if (length(d) == 0) stop("GSE does not have a recorded Bioproject or SRP; check that it is correct!")
    SRP <- unique(SRP)
    if (length(SRP) > 1) stop("Found multiple non identical SRPs for GSE:", SRP)
    message(paste("Found SRP, will continue using:", SRP))
  }
  # Download runinfo from Trace / SRA server
  url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term="
  url <- paste0(url, SRP)
  download.file(url, destfile)
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
  }
  # Remove unwanted columns, and search for more info
  file[, MONTH := substr(ReleaseDate, 6, 7)]
  file[, YEAR := gsub("-.*", "", ReleaseDate)]
  file <- file[, -c("AssemblyName", "ReleaseDate", "LoadDate",
                    "download_path", "RunHash", "ReadHash", "Consent")]
  # Download xml and add more data
  url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/exp?term="
  url <- paste0(url, SRP)
  a <- xml2::read_xml(url)
  a <- xml2::as_list(a)

  dt <- data.table()
  for(i in seq_along(a$EXPERIMENT_PACKAGE_SET)) { # Per sample in study
    EXP_SAMPLE <- a$EXPERIMENT_PACKAGE_SET[i]$EXPERIMENT_PACKAGE
    # Get Sample title
    xml.TITLE <- unlist(EXP_SAMPLE$SAMPLE$TITLE)
    xml.AUTHOR <- unlist(EXP_SAMPLE$Organization$Contact$Name$Last)
    # For each run in sample
    xml.RUN <- c()
    for (j in seq_along(EXP_SAMPLE$RUN_SET)) {
      xml.RUN <- c(xml.RUN, unlist(EXP_SAMPLE$RUN_SET[j]$RUN$IDENTIFIERS$PRIMARY_ID))
    }
    # Check if sample has alternative name (SOURCE)
    xml.SOURCE <- c()
    for (j in seq_along(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES)) {
      tag <- unlist(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES[[j]]$TAG)
      if (!is.null(tag)) {
        if (tag == "source_name") {
          xml.SOURCE <- c(xml.SOURCE, unlist(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES[[j]]$VALUE))
        }
      }
    }
    xml.TITLE <- ifelse(is.null(xml.TITLE), "", xml.TITLE)
    xml.AUTHOR <- ifelse(is.null(xml.AUTHOR), "",
                         ifelse(xml.AUTHOR %in% c("Curators", "GEO"), "", xml.AUTHOR))
    xml.SOURCE <- ifelse(is.null(xml.SOURCE), "", xml.SOURCE)
    if (length(xml.RUN) == 0) xml.RUN <- ""
    dt <- rbind(dt, cbind(xml.AUTHOR, xml.SOURCE, xml.TITLE, xml.RUN))
  }

  colnames(dt) <- c("AUTHOR", "sample_source", "sample_title", "Run")
  dt <- dt[Run %in% file$Run]
  if (length(dt) > 0) {
    file <- data.table::merge.data.table(file, dt, by = "Run")
  }
  if (abstract %in% c("print", "save", "printsave")) {
    if (!is.null(EXP_SAMPLE) && !is.null(EXP_SAMPLE$STUDY$DESCRIPTOR$STUDY_ABSTRACT[[1]])) {
      abstract_text <- EXP_SAMPLE$STUDY$DESCRIPTOR$STUDY_ABSTRACT[[1]]
      if (abstract %in% c("print", "printsave")) {
        print("Study abstract:")
        print(abstract_text)
      }
      if (abstract %in% c("save", "printsave")) {
        fwrite(data.table(abstract = abstract_text), abstract_destfile)
      }
    } else message("Could not find abstract for project")
  }
  pubmed.id <- unlist(EXP_SAMPLE$STUDY$STUDY_LINKS$STUDY_LINK$XREF_LINK$DB[[1]])
  if (!is.null(pubmed.id)) {
    if (pubmed.id == "pubmed") {
      file$Study_Pubmed_id <- as.integer(unlist(EXP_SAMPLE$STUDY$STUDY_LINKS$STUDY_LINK$XREF_LINK$ID[[1]]))
    }
  }
  # Find Study: example: 24476825
  if (!is.null(file$Study_Pubmed_id)) {
    if (!is.na(file$Study_Pubmed_id[1]) & (file$Study_Pubmed_id[1] != 3)) { # 3 is error code
      url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id="
      url <- paste0(url, file$Study_Pubmed_id[1])
      a <- xml2::as_list(xml2::read_xml(url))
      authors <- a$eSummaryResult$DocSum[[5]]
      if (!is.null(unlist(authors)[1])) {
        file$AUTHOR <- gsub(" .*", "", unlist(authors)[1])
        #TODO: Decide if I add Last author:
        # file$LAST_AUTHOR <- gsub(" .*", "", tail(unlist(authors), 1))
      }
    }
  }


  # Create ORFik guess columns from metadata:
  if (auto.detect) {
    file$LIBRARYTYPE <- findFromPath(file$sample_title,
                                     libNames(), "auto")
    if (any(file$LIBRARYTYPE %in% c(""))){ # Check if valid library strategy
      file[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LIBRARYTYPE <-
        findFromPath(file[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LibraryStrategy,
                     libNames(), "auto")
    }
    if (any(file$LIBRARYTYPE %in% c(""))){ # Check if valid library name
      file[LIBRARYTYPE == "",]$LIBRARYTYPE <-
        findFromPath(file[LIBRARYTYPE == "",]$LibraryName,
                     libNames(), "auto")
    }

    file$REPLICATE <- findFromPath(file$sample_title,
                                   repNames(), "auto")
    stages <- rbind(stageNames(), tissueNames(), cellLineNames())
    file$STAGE <- findFromPath(file$sample_title, stages, "auto")
    file$CONDITION <- findFromPath(file$sample_title,
                                   conditionNames(), "auto")
    file$INHIBITOR <- findFromPath(file$sample_title,
                                   inhibitorNames(), "auto")
  }
  fwrite(file, destfile)
  return(file)
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
                         ebiDLMethod = "auto", BPPARAM = bpparam()) {

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

  urls <- ORFik:::find_url_ebi(SRR, study = study)
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
#' Look for files in ebi following url: ftp://ftp.sra.ebi.ac.uk/vol1/fastq
#' Paired end and single end fastq files.\cr
#' EBI uses 3 ways to organize data inside vol1/fastq:\cr
#' - 1: Most common: SRR(3 first)/0(2 last)/whole\cr
#' - 2: less common: SRR(3 first)/00(1 last)/whole\cr
#' - 3: least common SRR(3 first)/whole
#' @param SRR character, SRR, ERR or DRR numbers.
#' @param stop.on.error logical FALSE, if TRUE will stop
#'  if all files are not found. If FALSE returns empty character vector if error
#'  is catched.
#' @param study default NULL, optional PRJ (study id) to speed up search
#' for URLs.
#' @return full url to fastq files, same length as input
#' (2 urls for paired end data). Returns empty character() if all
#' files not found.
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
find_url_ebi <- function(SRR, stop.on.error = FALSE, study = NULL) {
  message("Finding optimal download urls from ebi...")
  final.path <- if (!is.null(study)) {
    find_url_ebi_safe(study, SRR, stop.on.error = stop.on.error)
  } else find_url_ebi_safe(SRR, stop.on.error = stop.on.error)
  return(final.path)
  # TODO: remove when not needed
  ebi_server <- "ftp://ftp.sra.ebi.ac.uk"
  # Check that we can connect to ebi
  exists.ftp.dir.fast(ebi_server, report.error = TRUE)

  # Create candidate directories
  SRR_first_3 <- substring(SRR, 1, 6)
  SRR_last_3 <- paste0("0", reverse(substring(reverse(SRR), 1, 2)))
  SRR_last_1 <- paste0("00", reverse(substring(reverse(SRR), 1, 1)))
  SRR_default <- file.path(ebi_server, "vol1/fastq", SRR_first_3)

  SRR_fastq <- paste0(SRR, ".fastq.gz")
  SRR_fastq_paired <- c(paste0(SRR, c("_1"), ".fastq.gz"),
                        paste0(SRR, c("_2"), ".fastq.gz"))
  # method 1
  SRR_paths <- file.path(SRR_default, SRR_last_3, SRR, SRR_fastq)
  SRR_paths_paired <- file.path(SRR_default, SRR_last_3, SRR, SRR_fastq_paired)
  # method 2:
  SRR_paths_spec2 <- file.path(SRR_default, SRR_last_1, SRR, SRR_fastq)
  SRR_paths_paired_spec2 <- file.path(SRR_default, SRR_last_1, SRR, SRR_fastq_paired)
  # Method 3: Special location
  SRR_paths_spec <- file.path(SRR_default, SRR, SRR_fastq)
  SRR_paths_spec_paired <- file.path(SRR_default, SRR, SRR_fastq_paired)
  # Check what format the files are found in (3 types: 2 each)
  lib_counter <- 0
  SRR_possibilities <- list(SRR_paths, SRR_paths_paired, SRR_paths_spec2,
                            SRR_paths_paired_spec2, SRR_paths_spec,
                            SRR_paths_spec_paired)
  url.exists <- c()
  for (i in seq_along(SRR_possibilities)) { # For each url area
    if (lib_counter == length(SRR)) break
    url.temp <-  sapply(SRR_possibilities[[i]], function(x)
      exists.ftp.file.fast(x))
    url.temp <- url.temp[url.temp]
    if (i %% 2 == 1) {
      lib_counter <- lib_counter + length(url.temp)
    } else {
      lib_counter <- lib_counter + (length(url.temp)/2)
    }
    url.exists <- c(url.exists, url.temp)
  }

  final.path.temp <- names(url.exists[url.exists])
  # Sort them correctly as input
  final.path <- unlist(sapply(c(SRR, "asdasd"), function(x, final.path.temp) {
    sort(final.path.temp[grepl(x, final.path.temp)])
  }, final.path.temp = final.path.temp), use.names = FALSE)

  valid <- TRUE
  if (lib_counter != length(SRR)) valid <- FALSE
  paired <- length(grep("_[1-2]\\.fastq\\.gz",final.path))
  if (length(SRR) != (paired/2 + length(final.path) - paired))
    valid <- FALSE
  if (!valid & stop.on.error) stop("Did not find all fastq files on ENA",
                                   "set use.ebi.ftp to FALSE,
                                   to use fastq-dump instead")
  if (!valid) final.path <- character()

  return(final.path)
}

#' Find URL for EBI fastq files
#'
#' Safer version
#' @param accession character: (PRJ, SRP, ERP, DRP, SRX, SRR, ERR,..). For studies or samples,
#' it returns all runs per study or sample.
#' @param SRR character, which SRR numbers to subset by (can also be ERR or DRR numbers)
#' @param stop.on.error logical FALSE, if TRUE will stop
#'  if all files are not found. If FALSE returns empty character vector if error
#'  is catched.
#' @return character (1 element per SRR number)
#' @keywords internal
find_url_ebi_safe <- function(accession, SRR = NULL, stop.on.error = FALSE) {
  a <- data.table()
  for (i in accession) {
    search_url <- paste0("http://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                         i, "&result=read_run&fields=run_accession,fastq_ftp")
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

  return(unlist(strsplit(a$fastq_ftp, ";")))
}
