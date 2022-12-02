#' Merge groups of Fastq /Fasta files
#'
#' Will use multithreading to speed up process.
#' Only works for Unix OS (Linux and Mac)
#' @export
#' @param in_files \code{character} specify the full path to the individual fastq.gz files.
#' Seperated by space per file in group: For 2 output files from 4 input files:
#' in_files <- c("file1.fastq file2.fastq". "file3.fastq file4.fastq")
#' @param out_files \code{character} specify the path to the FASTQ directory
#' For 2 output files: out_files <- c("/merged/file1&2.fastq", "/merged/file3&4.fastq")
#' @inheritParams download.SRA
#' @return invisible(NULL).
#' @export
#' @examples
#' fastq.folder <- tempdir() # <- Your fastq files
#' infiles <- dir(fastq.folder, "*.fastq", full.names = TRUE)
#' \dontrun{
#' # Seperate files into groups (here it is 4 output files from 12 input files)
#' in_files <- c(paste0(grep(infiles, pattern = paste0("ribopool-",
#'                seq(11, 14), collapse = "|"), value = TRUE), collapse = " "),
#'               paste0(grep(infiles, pattern = paste0("ribopool-",
#'                seq(18, 19), collapse = "|"), value = TRUE), collapse = " "),
#'               paste0(grep(infiles, pattern = paste0("C11-",
#'                seq(11, 14), collapse = "|"), value = TRUE), collapse = " "),
#'               paste0(grep(infiles, pattern = paste0("C11-",
#'                seq(18, 19), collapse = "|"), value = TRUE), collapse = " "))
#'
#' out_files <- paste0(c("SSU_ribopool", "LSU_ribopool", "SSU_WT", "LSU_WT"), ".fastq.gz")
#' merged.fastq.folder <- file.path(fastq.folder, "merged/")
#' out_files <- file.path(merged.fastq.folder, out_files)
#'
#' mergeFastq(in_files, out_files)
#' }
mergeFastq <- function(in_files, out_files, BPPARAM = bpparam()) {
  # TODO: Make it work on windows
  if (.Platform$OS.type != "unix") stop("Merge does not work on windows OS")
  if (length(in_files) != length(out_files)) stop("Not equal length of in_files and out_files!")

  bplapply(seq_along(in_files), function(x, in_files, out_files) {
    dir.create(dirname(out_files[x]), showWarnings = FALSE, recursive = TRUE)
    system(paste("cat", in_files[x], ">", out_files[x]))
  }, in_files = in_files, out_files = out_files, BPPARAM = BPPARAM)

  message("Merge Done")
  return(invisible(NULL))
}

#' Very fast fastq/fasta collapser
#'
#' For each unique read in the file, collapse into 1 and state in the fasta header
#' how many reads existed of that type. This is done after trimming usually, works
#' best for reads < 50 read length. Not so effective for 150 bp length mRNA-seq etc.
#' @param files paths to fasta / fastq files to collapse. I tries to detect format per file,
#' if file does not have .fastq, .fastq.gz, .fq or fq.gz extensions, it will be treated
#' as a .fasta file format.
#' @param outdir outdir to save files, default:
#' \code{file.path(dirname(files[1]), "collapsed")}.
#' Inside same folder as input files, then create subfolder "collapsed", and add a prefix
#' of "collapsed_" to the output names in that folder.
#' @param header.out.format character, default "ribotoolkit", else must be "fastx".
#' How the read header of the output fasta should be formated: ribotoolkit: ">seq1_x55",
#' sequence 1 has 55 duplicated reads collapsed.
#' fastx: ">1-55", sequence 1 has 55 duplicated reads collapsed
#' @param prefix character, default "collapsed_"
#' Prefix to name of output file.
#' @param compress logical, default FALSE
#' @return invisible(NULL), files saved to disc in fasta format.
#' @export
#' @examples
#'
#' fastq.folder <- tempdir() # <- Your fastq files
#' infiles <- dir(fastq.folder, "*.fastq", full.names = TRUE)
#' # collapse.fastq(infiles)
#'
collapse.fastq <- function(files, outdir = file.path(dirname(files[1]), "collapsed"),
                           header.out.format = "ribotoolkit", compress = FALSE,
                           prefix = "collapsed_") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }

  format <- rep("fasta", length(files))
  format[grep("\\.fastq|\\.fq", files)] <- "fastq"

  for (f in seq_along(files)) {
    file <- files[f]
    message("File ", f, "/", length(files), ":  ", file)
    fasta_name <- gsub(pattern = "\\.fastq$", replacement = ".fasta",
                       basename(file))
    new_file_name <- paste0(prefix, fasta_name,
                            ifelse(compress, ".gz", ""))
    writeXStringSet(collapse.fastq.internal(
              readDNAStringSet(file, format = format[f], use.names = FALSE), header.out.format),
                    file.path(outdir, new_file_name),
                    compress = compress, format = "fasta")
  }
  return(invisible(NULL))
}

collapse.fastq.internal <- function(seqs, header.out.format = "ribotoolkit") {
  # Fast collapser using data.table
  replicates <- data.table(seqs = as.character(seqs))
  # Much faster with 1 core actually, strange..
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(1)
  replicates <- replicates[, .N, by = seqs][order(N, decreasing = TRUE),]
  data.table::setDTthreads(old_threads)
  if (header.out.format == "fastx") {
    headers <- paste0(seq.int(nrow(replicates)), "-", replicates$N)
  } else if (header.out.format == "ribotoolkit") {
    headers <- paste0("seq", seq.int(nrow(replicates)), "_x", replicates$N)
  } else stop("format must be 'fastx' or 'ribotoolkit'")
  new_seqs <- replicates$seqs
  names(new_seqs) <- headers
  return(DNAStringSet(new_seqs, use.names = TRUE))
}
