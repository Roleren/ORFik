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

fetch_xml_attributes <- function(xml, xml_path) {
  dt <- XML::getNodeSet(xml, xml_path)
  dt <-  lapply(dt, XML::xmlToDataFrame)
  dt <-  lapply(dt, function(x) as.data.table(t(x)))
  dt <- lapply(dt, function(x) {colnames(x) <- as.character(x[1,]); return(x[2,])})
  dt <- rbindlist(dt, fill = TRUE)
  to_keep <- sapply(dt, function(x) !all(x == "NA"))
  dt <- dt[,..to_keep]
  return(dt)
}

sample_info_append_SRA <- function (SRP, destfile, abstract_destfile, abstract, remove.invalid,
                                    rich.format = FALSE, fetch_GSE = FALSE)
{
  xml <- sample_info_download(SRP)
  if (fetch_GSE) {
    all_ext <- xml2::xml_text(xml2::xml_find_all(xml, ".//EXTERNAL_ID"))
    gse <- all_ext[grep("GSE",all_ext)]
    gse <- unique(gse)
    if (length(gse) > 1) warning("Multiple GSE identifiers found, using first")
    if (length(gse) == 1) gse <- gse[1]
    if (length(gse) == 0) gse <- NA_character_
  }
  if (rich.format) {
    sample_xml <- XML::xmlParse(xml)
    sample_dt <- fetch_xml_attributes(sample_xml, "//SAMPLE/SAMPLE_ATTRIBUTES")
    exp_attr_dt <- fetch_xml_attributes(sample_xml, "//EXPERIMENT/EXPERIMENT_ATTRIBUTES")
  }
  xml <- xml2::as_list(xml)
  dt <- data.table()
  is_EXP_SET <- !is.null(xml$EXPERIMENT_PACKAGE_SET)
  EXP <- if (is_EXP_SET) {
    xml$EXPERIMENT_PACKAGE_SET
  }  else xml
  for (i in seq_along(EXP)) {
    EXP_SAMPLE <- EXP[i]$EXPERIMENT_PACKAGE
    dt_single <- sample_info_single(EXP_SAMPLE)
    dt <- rbind(dt, dt_single)
  }
  no_linker <- !is.null(EXP$eFetchResult)
  if (no_linker) {
    warning("Could not find SRA linker, falling back to deprecated search",
            "This might not find all the samples!")
    dt <- study_runinfo_download_deprecated(SRP, destfile)
  }
  else {
    dt <- add_pubmed_id(EXP_SAMPLE, dt)
  }
  dt <- add_author(dt)
  if (fetch_GSE) dt[,GEO := gse]
  if (rich.format)
    dt <- cbind(dt, sample_dt, exp_attr_dt)
  abstract_save(EXP_SAMPLE, abstract, abstract_destfile)
  return(filter_empty_runs(dt, remove.invalid, SRP))
}

sample_info_single <- function(EXP_SAMPLE) {

  # For each run in sample
  dt_run_all <- data.table()
  for (j in seq_along(EXP_SAMPLE$RUN_SET)) {
    RUN <- EXP_SAMPLE$RUN_SET[j]$RUN
    xml.RUN <- unlist(RUN$IDENTIFIERS$PRIMARY_ID)

    spots <- as.integer(attr(RUN, "total_spots"))
    bases <- as.numeric(attr(RUN, "total_bases"))
    # spots_with_mates <- 0
    avgLength <- as.integer(attr(RUN$Statistics$Read, "average"))
    size_MB <- floor(as.numeric(attr(RUN, "size"))/1024^2)
    Experiment <- EXP_SAMPLE$EXPERIMENT$IDENTIFIERS$PRIMARY_ID[[1]]
    # if (length(xml.RUN) == 0) xml.RUN <- ""
    dt_run <- data.table(Run = xml.RUN, spots, bases,
                         avgLength, size_MB, Experiment)
    dt_run_all <- rbind(dt_run_all, dt_run)
  }
  # Sample info (only 1 row per sample)
  # Library
  lib_design <- EXP_SAMPLE$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR
  if (length(lib_design$LIBRARY_NAME) == 0) {
    LibraryName <- NA
  } else LibraryName <- lib_design$LIBRARY_NAME[[1]]

  LibraryStrategy <- lib_design$LIBRARY_STRATEGY[[1]]
  LibrarySelection <- lib_design$LIBRARY_SELECTION[[1]]
  LibrarySource <- lib_design$LIBRARY_SOURCE
  LibraryLayout <- names(lib_design$LIBRARY_LAYOUT)
  dt_library <- data.table(LibraryName, LibraryStrategy,
                           LibrarySelection, LibrarySource, LibraryLayout)
  # Insert
  InsertSize <- 0
  InsertDev <- 0
  dt_insert <- data.table(InsertSize, InsertDev)
  # Instrument
  Platform_info <- EXP_SAMPLE$EXPERIMENT$PLATFORM
  Platform <- names(Platform_info)
  Model <- unlist(Platform_info[[1]]$INSTRUMENT_MODEL)
  dt_platform <- data.table(Platform, Model)
  # Sample
  SRAStudy <- EXP_SAMPLE$STUDY$IDENTIFIERS$PRIMARY_ID[[1]]
  BioProject <- attr(EXP_SAMPLE$STUDY$IDENTIFIERS$EXTERNAL_ID, "namespace")
  BioProject <- EXP_SAMPLE$STUDY$IDENTIFIERS$EXTERNAL_ID[[1]]
  if (is.null(BioProject)) {
    BioProject <- attr(EXP_SAMPLE$EXPERIMENT$STUDY_REF$IDENTIFIERS$EXTERNAL_ID, "namespace")
    BioProject <- EXP_SAMPLE$EXPERIMENT$STUDY_REF$IDENTIFIERS$EXTERNAL_ID[[1]]
  }
  Study_Pubmed_id <- as.numeric(NA)
  ProjectID <- ""
  Sample_info <- EXP_SAMPLE$SAMPLE
  Sample <- Sample_info$IDENTIFIERS$PRIMARY_ID[[1]]
  BioSample <- Sample_info$IDENTIFIERS$EXTERNAL_ID[[1]]
  SampleType <- "simple"
  dt_sample <- data.table(SRAStudy, BioProject, Study_Pubmed_id, ProjectID,
                          Sample, BioSample, SampleType)
  # Organism
  TaxID <- Sample_info$SAMPLE_NAME$TAXON_ID[[1]]
  ScientificName <- Sample_info$SAMPLE_NAME$SCIENTIFIC_NAME[[1]]
  dt_organism <- data.table(TaxID, ScientificName)
  # Submission
  SampleName <- attr(Sample_info, "alias")
  CenterName <- attr(EXP_SAMPLE$SUBMISSION, "center_name")
  Submission <- attr(EXP_SAMPLE$SUBMISSION, "accession")
  ReleaseDate <- attr(RUN, "published")
  MONTH <- substr(ReleaseDate, 6, 7)
  YEAR <- gsub("-.*", "", ReleaseDate)
  AUTHOR <- unlist(EXP_SAMPLE$Organization$Contact$Name$Last)
  AUTHOR <- ifelse(is.null(AUTHOR), "",
                   ifelse(AUTHOR %in% c("Curators", "GEO"), "", AUTHOR))

  dt_submission <- data.table(SampleName, CenterName, Submission,
                              MONTH, YEAR, AUTHOR)
  # Source and title
  # Check if sample has alternative name (SOURCE)
  sample_source <- c()
  for (j in seq_along(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES)) {
    tag <- unlist(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES[[j]]$TAG)
    if (!is.null(tag)) {
      if (tag == "source_name") {
        sample_source <- c(sample_source, unlist(EXP_SAMPLE$SAMPLE$SAMPLE_ATTRIBUTES[[j]]$VALUE))
      }
    }
  }
  sample_source <- ifelse(is.null(sample_source), "", sample_source)
  # Get Sample title
  sample_title <- unlist(EXP_SAMPLE$SAMPLE$TITLE)
  sample_title <- ifelse(is.null(sample_title), "", sample_title)
  dt_additional <- data.table(sample_source, sample_title)


  dt_single <- data.table(cbind(dt_run_all, dt_library, dt_insert, dt_platform,
                                dt_sample, dt_organism, dt_submission,
                                dt_additional))
}

study_runinfo_download_deprecated <- function(SRP, destfile) {
  url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term="
  url <- paste0(url, SRP)
  download.file(url, destfile)
  return(fread(destfile))
}

sample_info_download <- function(SRP, max_samples = 2000) {
  # Trace db is dead, now we have to get ids first, then get samples
  url_prefix <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term="
  url_suffix <- paste0("&retmax=", max_samples, "&retmode=xml&tool=ORFik")
  url <- paste0(url_prefix, SRP, url_suffix)
  ids_xml <- xml2::as_list(xml2::read_xml(url))
  ids <- sort(unlist(ids_xml$eSearchResult$IdList, use.names = FALSE))
  if (is.null(ids)) {
    stop("SRA server returned the following error: \n",
         ids_xml$eSearchResult$WarningList$OutputMessage[[1]])
  }
  url_prefix <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="
  url_suffix <- "&rettype=fasta&retmode=xml&tool=ORFik"
  url <- paste0(url_prefix, paste(ids, collapse = ","), url_suffix)
  message("Downloading metadata from:")
  message(url)
  return(xml2::read_xml(url))
}

abstract_save <- function(EXP_SAMPLE, abstract, abstract_destfile) {
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
    } else message("Note: Could not find abstract for project")
  }
}

add_pubmed_id <- function(EXP_SAMPLE, file) {
  pubmed.id <- unlist(EXP_SAMPLE$STUDY$STUDY_LINKS$STUDY_LINK$XREF_LINK$DB[[1]])
  if (!is.null(pubmed.id)) {
    if (pubmed.id == "pubmed") {
      file$Study_Pubmed_id <- as.integer(unlist(EXP_SAMPLE$STUDY$STUDY_LINKS$STUDY_LINK$XREF_LINK$ID[[1]]))
    }
  } else {
    file$Study_Pubmed_id <- pmid_from_bioproject(file$BioProject[1])
  }
  return(file)
}

pmid_from_bioproject <- function(id) {
  if (is.null(id)) return(as.numeric(NA))
  id_short <- gsub("[^0-9]", "", id)
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=bioproject&db=pubmed&id=", id_short,
                "&retmode=xml&tool=ORFik")
  xml <- xml2::as_list(xml2::read_xml(url))
  pmid <- xml$eLinkResult$LinkSet$LinkSetDb$Link$Id[[1]]
  if (is.null(pmid)) {
    pmid <- as.numeric(NA)
  }
  names(pmid) <- id
  return(as.numeric(pmid))
}

add_author <- function(file) {
  if (!is.null(file$Study_Pubmed_id)) {
    # 3 is error code
    valid_pubmed_id <- !is.na(file$Study_Pubmed_id[1]) & (file$Study_Pubmed_id[1] != 3)
    if (valid_pubmed_id) {
      url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id="
      url <- paste0(url, file$Study_Pubmed_id[1], "&tool=ORFik")
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

filter_empty_runs <- function(dt, remove.invalid = TRUE, SRP) {
  # Filter
  if (nrow(dt) == 0) {
    warning("Experiment not found on SRA, are you sure it is public?")
    return(dt)
  }

  msg <- paste("Found Runs with 0 reads (spots) in metadata, will not be able
              to download the run/s:", dt[spots == 0,]$Run)
  dt[is.na(spots), spots := 0]
  if (any(dt$spots == 0)) {
    warning(msg)
    if (remove.invalid) {
      message("Removing invalid Runs from final metadata list")
      dt <- dt[spots > 0,]
    }
  }

  if (nrow(dt) == 0) {
    warning(paste("No valid runs found from experiment:", SRP))
    return(dt)
  }
  return(dt)
}
