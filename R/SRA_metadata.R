#' Query eutils for bioproject IDs
#' @param term character, default "Ribosome Profiling human"
#' @param RetMax integer, default 10000. How many IDs to return maximum
#' @return character vector of IDs
#' @export
#' @references https://www.ncbi.nlm.nih.gov/books/NBK25501/
#' @family sra
#' @examples
#' term <- "Ribosome Profiling Saccharomyces cerevisiae"
#' # get_bioproject_candidates(term)
get_bioproject_candidates <- function(term = "Ribosome Profiling human", RetMax = 10000) {
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  esearch_base_url <- "esearch.fcgi?db=bioproject"
  fixed_term <- gsub(" ", "+", term)
  search <- paste0("&retmode=xml&RetMax=", RetMax, "&term=", fixed_term)
  url <- paste0(base_url, esearch_base_url, search)
  hits <- xml2::read_xml(url)
  hits <- unlist(xml2::as_list(hits)$eSearchResult$IdList, use.names = FALSE)
  return(hits)
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
#' @param force logical, default FALSE. If TRUE, will redownload
#' all files needed even though they exists. Useuful if you wanted
#' auto.detection, but already downloaded without it.
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
#' # download.SRA.metadata("GSE61011")
download.SRA.metadata <- function(SRP, outdir = tempdir(), remove.invalid = TRUE,
                                  auto.detect = FALSE, abstract = "printsave",
                                  force = FALSE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  destfile <- paste0(outdir, "/SraRunInfo_", SRP, ".csv")
  abstract_destfile <- paste0(outdir, "/abstract_", SRP, ".csv")

  if (file.exists(destfile) & !force) {
    message(paste("Existing metadata file found in dir:", outdir, ",will not download"))
    if (abstract %in% c("print", "printsave")) { # Print abstract if wanted
      if (file.exists(abstract_destfile)) {
        cat("Study Abstract:\n")
        cat(readLines(abstract_destfile)[-1], " \n")
        cat("------------------\n")
      }
    }
    return(fread(destfile))
  }
  is_GSE <- length(grep("GSE", x = SRP)) == 1
  # Convert GSE to SRP
  if (is_GSE) SRP <- SRP_from_GSE(SRP)
  # Download runinfo from Trace / SRA server
  file <- study_runinfo_download(SRP, destfile)
  if (nrow(file) == 0) {
    warning("Experiment not found on SRA, are you sure it is public?")
    return(file)
  }

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
  # Remove unwanted columns
  file[, MONTH := substr(ReleaseDate, 6, 7)]
  file[, YEAR := gsub("-.*", "", ReleaseDate)]
  file <- file[, -c("AssemblyName", "ReleaseDate", "LoadDate",
                    "download_path", "RunHash", "ReadHash", "Consent")]

  file <- sample_info_append_SRA(file, SRP, abstract_destfile, abstract)

  # Create ORFik guess columns from metadata:
  if (auto.detect) {
    file <- metadata.autnaming(file)
  }
  fwrite(file, destfile)
  return(file)
}

SRP_from_GSE <- function(SRP) {
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
  return(SRP)
}

sample_info_append_SRA <- function(file, SRP, abstract_destfile, abstract) {
  # Download xml and add more data
  a <- sample_info_download(SRP)

  dt <- data.table()
  is_EXP_SET <- !is.null(a$EXPERIMENT_PACKAGE_SET)
  EXP <- if(is_EXP_SET) {a$EXPERIMENT_PACKAGE_SET} else a

  for(i in seq_along(EXP)) { # Per sample in study
    EXP_SAMPLE <- EXP[i]$EXPERIMENT_PACKAGE
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
    # Add title and author information
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
        cat("Study Abstract:\n")
        cat(abstract_text, " \n")
        cat("------------------\n")
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
  return(file)
}

study_runinfo_download <- function(SRP, destfile) {
  url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term="
  url <- paste0(url, SRP)
  download.file(url, destfile)
  return(fread(destfile))
}

sample_info_download <- function(SRP) {
  # Trace db is dead, now we have to get ids first, then get samples
  url_prefix <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term="
  url_suffix <- "&retmax=2000&retmode=xml"
  url <- paste0(url_prefix, SRP, url_suffix)
  ids_xml <- xml2::as_list(xml2::read_xml(url))
  ids <- unlist(ids_xml$eSearchResult$IdList, use.names = FALSE)
  url_prefix <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="
  url_suffix <- "&rettype=fasta&retmode=xml"
  url <- paste0(url_prefix, paste(ids, collapse = ","), url_suffix)
  return(xml2::as_list(xml2::read_xml(url)))
}

#' Guess SRA metadata columns
#' @param file a data.table of SRA metadata
#' @return a data.table of SRA metadata with additional columns:
#'       LIBRARYTYPE, REPLICATE, STAGE, CONDITION, INHIBITOR
metadata.autnaming <- function(file) {
  ## Library type
  # First test sample_title
  file$LIBRARYTYPE <- findFromPath(file$sample_title,
                                   libNames(), "auto")
  # Test LIBRARYTYPE
  if (any(file$LIBRARYTYPE %in% c(""))){ # Check if valid library strategy
    file[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LIBRARYTYPE <-
      findFromPath(file[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LibraryStrategy,
                   libNames(), "auto")
  }
  # Test LIBRARYNAME
  if (any(file$LIBRARYTYPE %in% c(""))){ # Check if valid library name
    file[LIBRARYTYPE == "",]$LIBRARYTYPE <-
      findFromPath(file[LIBRARYTYPE == "",]$LibraryName,
                   libNames(), "auto")
  }
  # The other columns
  file$REPLICATE <- findFromPath(file$sample_title,
                                 repNames(), "auto")
  stages <- rbind(stageNames(), tissueNames(), cellLineNames())
  file$STAGE <- findFromPath(file$sample_title, stages, "auto")
  file$CONDITION <- findFromPath(file$sample_title,
                                 conditionNames(), "auto")
  file$INHIBITOR <- findFromPath(file$sample_title,
                                 inhibitorNames(), "auto")
  return(file)
}
