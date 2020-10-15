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
#' @param info a data.frame with information and SRR numbers or character vector of only SRR numbers.
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
#' ## With additional info for renaming, info download manually from SRA
#' library(data.table)
#' # The SRR numbers as a list in file downloaded manually
#' SRR.ids <- "/export/valenfs/data/raw_data/TCP-seq/sel-TCPSeq_Wagner/SRR.ids"
#' # The sra results as a file downloaded manually
#' info <- fread("/export/valenfs/data/raw_data/TCP-seq/sel-TCPSeq_Wagner/sra_result.csv")
#' SRR <- unlist(read.table(SRR.ids,
#'                          stringsAsFactors = F), use.names = FALSE)
#' info <- cbind(info, SRR)
#' info <- info[1,]
#' # Download, 5 first reads
#' download.SRA(info, outdir, subset = 5)
#' }
download.SRA <- function(info, outdir, rename = TRUE,
                         fastq.dump.path = install.sratoolkit(),
                         settings =  paste("--gzip", "--skip-technical", "--split-files"),
                         subset = NULL,
                         BPPARAM = bpparam()) {
  SRR <- ifelse(is.character(info), info, info$SRR)
  rename <- ifelse(!is.character(info) & rename, TRUE, FALSE)

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

  if (rename) {
    if (is.null(info$`Experiment Title`)) {
      warning("Did not find the column Experiment Title in info, returning without renaming!")
      return(files)
    }
    message("Renaming files:")
    new_names <- info$`Experiment Title`
    new_names <- gsub(".*: ", "", new_names); new_names <- gsub(";.*", "", new_names)
    new_names <- gsub(" ", "_", new_names)
    new_names <- paste0(dirname(files), "/",new_names, ".fastq.gz")
    for (i in seq(length(files))) {
      file.rename(files[i], new_names[i])
    }
    files <- new_names
  }
  return(files)
}
