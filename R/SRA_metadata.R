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


  file <- sample_info_append_SRA(SRP, destfile, abstract_destfile, abstract)

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
