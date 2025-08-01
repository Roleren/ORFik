#' Downloads metadata from SRA
#'
#' Given a experiment identifier, query information from different locations of SRA
#' to get a complete metadata table of the experiment. It first finds Runinfo for each
#' library, then sample info,
#' if pubmed id is not found searches for that and searches for author through pubmed.
#'
#' A common problem is that the project is not linked to an article, you will then not
#' get a pubmed id.
#'
#' The algorithm works like this:\cr
#' If GEO identifier, find the SRP.\cr
#' Then search Entrez for project and get sample identifier.\cr
#' From that extract the run information and collect into a final table.\cr
#'
#' @param SRP character string, a study ID as either the PRJ, SRP, ERP, DRP, GSE or SRA of the study,
#' examples would be "SRP226389" or "ERP116106". If GSE it will try to convert to the SRP
#' to find the files. The call works as long the runs are registered on the efetch server,
#' as their is a linked SRP link from bioproject or GSE. Example which fails is "PRJNA449388",
#' which does not have a linking like this.
#' @param outdir character string, directory to save file, default: tempdir().
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
#' @param rich.format logical, default FALSE. If TRUE, will fetch all Experiment and Sample attributes.
#' It means, that different studies can have different set of columns if set to TRUE.
#' @param fetch_GSE logical, default FALSE. Search for GSE, if exists, appends a column
#' called GEO. Will be included even though this study is not from GEO, then it
#' sets all to NA.
#' @return a data.table of the metadata, 1 row per sample,
#'  SRR run number defined in 'Run' column.
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
#' ## Bioproject ID
#' # download.SRA.metadata("PRJNA231536")
download.SRA.metadata <- function(SRP, outdir = tempdir(), remove.invalid = TRUE,
                                  auto.detect = FALSE, abstract = "printsave",
                                  force = FALSE, rich.format = FALSE,
                                  fetch_GSE = FALSE) {
  stopifnot(length(SRP) == 1)
  stopifnot(length(outdir) == 1)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  destfile <- paste0(outdir, "/SraRunInfo_", SRP, ".csv")
  abstract_destfile <- paste0(outdir, "/abstract_", SRP, ".csv")

  if (file.exists(destfile) & !force) {
    message(paste("Existing metadata file found in dir:", outdir, ", will not download"))
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


  file <- sample_info_append_SRA(SRP, destfile, abstract_destfile, abstract,
                                 remove.invalid, rich.format = rich.format,
                                 fetch_GSE = fetch_GSE)

  # Create ORFik guess columns from metadata:
  if (auto.detect) {
    file <- metadata.autnaming(file)
  }
  fwrite(file, destfile)
  return(file)
}

#' Query eutils for bioproject IDs
#'
#' The default query of Ribosome Profiling human, will result in internal
#' entrez search of:
#' Ribosome[All Fields] AND Profiling[All Fields] AND ("Homo sapiens"[Organism]
#' OR human[All Fields])
#' @param term character, default "Ribosome Profiling human".
#' A space is translated into AND, that means "Ribosome AND Profiling AND human",
#' will give same as above. To do OR operation, do:
#' "Ribosome OR profiling OR human".
#' @param as_accession logical, default TRUE. Get bioproject accessions:
#' PRJNA, PRJEB, PRJDB values, or IDs (FALSE), numbers only. Accessions
#' are usually the thing needed for most tools.
#' @param add_study_title logical, default FALSE. If TRUE, return as data table
#' with 2 columns: id: ID or accessions. title: The title of the study.
#' @param RetMax integer, default 10000. How many IDs to return maximum
#' @return character vector of Accessions or IDs. If add_study_title is TRUE,
#' returns a data.table.
#' @export
#' @references https://www.ncbi.nlm.nih.gov/books/NBK25501/
#' @family sra
#' @examples
#' term <- "Ribosome Profiling Saccharomyces cerevisiae"
#' # get_bioproject_candidates(term)
get_bioproject_candidates <- function(term = "Ribosome Profiling human",
                                      as_accession = TRUE,
                                      add_study_title = FALSE,
                                      RetMax = 10000) {
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  esearch_base_url <- "esearch.fcgi?db=bioproject"
  fixed_term <- gsub(" ", "+", term)
  search <- paste0("&retmode=xml&RetMax=", RetMax, "&term=", fixed_term)
  url <- paste0(base_url, esearch_base_url, search)
  message("- Retrieving IDs")
  hits <- xml2::read_xml(url)
  hits <- unlist(xml2::as_list(hits)$eSearchResult$IdList, use.names = FALSE)
  message("done")

  needs_esummary <- (as_accession | add_study_title) & length(hits) > 0
  if (needs_esummary) {
    message("- Retrieving Accessions from IDs")
    esummary_base_url <- "esummary.fcgi?db=bioproject&id="
    # Max 500 ids per call
    ids_list <- split(hits, ceiling(seq_along(hits)/ 500))
    hits <- if(add_study_title) {data.table()} else c()
    for (ids_subset in ids_list) {
      ids <- paste(ids_subset, collapse = ",")
      url <- paste0(base_url, esummary_base_url, ids)
      temp_hits <- xml2::read_xml(url)
      xml_list <- xml2::as_list(temp_hits)

      if (as_accession) {
        temp_hits <- unlist(lapply(xml_list$eSummaryResult$DocumentSummarySet,
                              function(x) x$Project_Acc[[1]]), use.names = FALSE)
      }
      if (add_study_title) {
        titles <- unlist(lapply(xml_list$eSummaryResult$DocumentSummarySet,
                                function(x) x$Project_Title[[1]]), use.names = FALSE)
        temp_hits <- data.table(id = temp_hits, title = titles)
        hits <- rbindlist(list(hits, temp_hits))
      } else hits <- c(hits, temp_hits)
    } # Max ids per time is 500

    message("done")
  }
  return(hits)
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

#' Open SRA in browser for specific bioproject
#' @param x character, bioproject ID.
#' @inheritParams utils::browseURL
#' @return invisible(NULL), opens webpage only
#' @family sra
#' @export
#' @examples
#' #browseSRA("PRJNA336542")
#'
#' #' # For windows make sure a valid browser is defined:
#' browser <- getOption("browser")
#' #browseSRA("PRJNA336542", browser)
browseSRA <- function(x, browser = getOption("browser")) {
  browseURL(paste0("https://www.ncbi.nlm.nih.gov/Traces/study/?acc=", x, "&o=acc_s%3Aa"),
            browser = browser)
}

#' Download summary information of a gene
#'
#' Uses ncbi gene database summary from RefSeq
#' @param gene character, gene name (symbol)
#' @param organism, default NULL. Scientific name (e.g. Homo sapiens)
#' @param by character, default symbol (search by gene symbol name).
#' If "ensembl id", it seraches as it is ensembl gene id ENSG.. etc.
#' @return character, summary text for gene from the database.
#' @export
#' @importFrom jsonlite read_json
#' @examples
#' # Wrap in 'try' to avoid wrong bioc test error
#' try(download_gene_info(gene = "CCND1"))
#' try(download_gene_info("ENSG00000110092", by = "ensembl_id")) # By ensembl id
#' try(download_gene_info(gene = "CCND1", organism = "Mus musculus"))
download_gene_info <- function(gene = "CCND1", organism = "Homo sapiens", by = "symbol") {
  stopifnot(by %in% c("symbol", "ensembl_id"))
  format <- ifelse(by == "symbol", "[Gene%20Name]+", "[Source%20ID]+")
  esearch_gene_api_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene"
  if (!is.null(organism)) organism <- gsub(" ", "%20", organism)
  url_search <- paste0(esearch_gene_api_url, "&term=", gene, format, organism, "[Organism]&retmode=json&tool=ORFik")
  json <- jsonlite::read_json(url_search)
  hits <- json$esearchresult$count[[1]]
  if (hits == 0) {
    message("No id found for this gene, returning empty string")
    return("")
  }
  id <- json$esearchresult$idlist[[1]]
  esummary_gene_api_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene"
  url_summary <- paste0(esummary_gene_api_url, "&id=", id,"&retmode=json&tool=ORFik")
  json <- jsonlite::read_json(url_summary)
  summary <- json$result[[id]]$summary
  if (length(summary) == 0 || summary == "") {
    summary <- paste(json$result[[id]]$description, "| Chromosome:", json$result[[id]]$chromosome)
  }
  return(summary)
}

#' Download homologue information of a gene
#'
#' Uses ncbi gene database for vertebrates
#' @param gene_id character, gene name (ensembl gene id, not symbol!)
#' @param organism, default NULL. Scientific name (e.g. Homo sapiens)
#' @return character, summary text for gene from the database.
#' @export
download_gene_homologues <- function(gene_id = "ENSG00000110092", organism = "Homo sapiens") {
  organism <- tolower(gsub(" ", "_", organism))
  rest_api_url <- "https://rest.ensembl.org/homology/id/"
  organism_url <- paste0(organism, "/", gene_id)
  settings_url <- "?content-type=application/json&sequence=none&type=orthologues&cigar_line=0&format=condensed"
  url_summary <- paste0(rest_api_url, organism_url, settings_url)
  json <- jsonlite::read_json(url_summary)
  homologues <- json$data[[1]]$homologies
  dt <- rbindlist(lapply(homologues, function(species) {
    data.table::data.table(species = species$species, gene_id = species$id)
    }))[!duplicated(species)]
  return(dt)
}
