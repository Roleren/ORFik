#' Merge groups of Fastq /Fasta files
#'
#' Will use multithreading to speed up process.
#' Only works for Unix OS (Linux and Mac)
#' @export
#' @param in_files_by_out_file_list \code{list of character vectors},
#' Per element of list are  full path to the individual fastq.gz files to that
#' output file.
#' @inheritParams download.SRA
#' @importFrom BiocParallel bpmapply
#' @return invisible(NULL).
#' @export
#' @examples
#'
#' # Make small example fastq files
#' fastq.folder <- tempdir() # <- Your fastq files
#' # Seperate files into groups (here it is 4 output files from 12 input files)
#' in_files <- paste0(LETTERS[1:16], ".fastq.gz")
#' in_files <- file.path(fastq.folder, in_files)
#' samples_dna_letters <- vapply(seq_along(in_files), function(x)
#'   paste(sample(DNA_ALPHABET[1:4], 12, replace = TRUE), collapse = ""), character(1))
#' # Write example input files to temp
#' lapply(seq_along(in_files), function(i) {
#'  seq <- DNAStringSet(samples_dna_letters[i])
#'  names(seq) <- basename(in_files[i])
#'  writeXStringSet(seq, in_files[i])
#' })
#'
#' out_files <- paste0(c("SSU_ribopool", "LSU_ribopool", "SSU_WT", "LSU_WT"), ".fastq.gz")
#' merged.fastq.folder <- file.path(fastq.folder, "merged/")
#' out_files <- file.path(merged.fastq.folder, out_files)
#'
#' in_files_by_out_file_list <- split(in_files, rep(out_files, each = 4))
#' mergeFastq(in_files_by_out_file_list, BiocParallel::SerialParam())
#' lapply(out_files, readDNAStringSet)
mergeFastq <- function(in_files_by_out_file_list, BPPARAM = bpparam()) {
  # TODO: Make it work on windows
  if (.Platform$OS.type != "unix") stop("Merge does not work on windows OS")
  all_input_files_exist <- all(file.exists(unlist(in_files_by_out_file_list)))
  stopifnot(is(in_files_by_out_file_list, "list"))
  stopifnot(length(in_files_by_out_file_list) > 0)
  stopifnot(all_input_files_exist)
  if (any(lengths(in_files_by_out_file_list) == 0))
    stop("An output file group had 0 input files!")

  equal_compression_per_group <- max(sapply(in_files_by_out_file_list, function(files) length(unique(tools::file_ext(files)))))
  if (any(equal_compression_per_group != 1))
    stop("You cant mix compressed and uncompressed input files per output file!")


  BiocParallel::bpmapply(function(in_files, out_file) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    status <- system(paste("cat", paste(shQuote(in_files), collapse = " "),
                           ">", shQuote(out_file)))
    if (status != 0)
      stop("Failed to merge fastq files into: ", out_file)
  }, in_files = in_files_by_out_file_list,
     out_file = names(in_files_by_out_file_list), BPPARAM = BPPARAM)

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
#' seqs <- rep(c("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGG"), times = c(2,3))
#' fastq.folder <- tempfile()
#' dir.create(fastq.folder)
#' fastq.file <- file.path(fastq.folder, "test.fasta")
#' writeXStringSet(DNAStringSet(seqs), filepath = fastq.file)
#' infiles <- dir(fastq.folder, "fasta|fastq", full.names = TRUE)
#' collapse.fastq(infiles)
#' readDNAStringSet(file.path(fastq.folder, "collapsed", "collapsed_test.fasta"))
#' # You see names says x3 of read 1 (GGGG...) and x2 of read 2 (AAAA....)
collapse.fastq <- function(files, outdir = file.path(dirname(files[1]), "collapsed"),
                           header.out.format = "ribotoolkit", compress = FALSE,
                           prefix = "collapsed_") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }
  format <- get_fast_format(files)


  for (f in seq_along(files)) {
    file <- files[f]
    message("File ", f, "/", length(files), ":  ", file)
    fasta_name <- gsub(pattern = fastq_regex_pattern(), replacement = ".fasta",
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

collapse.fastq.internal <- function(seqs, header.out.format = "ribotoolkit",
                                    asDNAStringSet = TRUE) {
  # Fast collapser using data.table
  replicates <- data.table(seqs = as.character(seqs))
  local_DTthreads(1)
  replicates <- replicates[, .N, by = seqs][order(N, decreasing = TRUE),]

  new_seqs <- replicates$seqs
  names(new_seqs) <- collapse_header_set(replicates$N, seq.int(nrow(replicates)),
                                         header.out.format)
  if (!asDNAStringSet) return(new_seqs)
  return(DNAStringSet(new_seqs, use.names = TRUE))
}

collapse_header_set <- function(N, indices = seq_along(N),
                                header.out.format = "ribotoolkit") {
  if (header.out.format == "fastx") {
    headers <- sprintf("%d-%d", indices, N)
  } else if (header.out.format == "ribotoolkit") {
    headers <- sprintf("seq%d_x%d", indices, N)
  } else stop("format must be 'fastx' or 'ribotoolkit'")
  return(headers)
}

get_fast_format <- function(files, allowed_compressions = "gz") {
  format <- rep("fasta", length(files))
  format[file_is_fastq(files, allowed_compressions)] <- "fastq"
  return(format)
}

file_is_fastq <- function(files, allowed_compressions = "gz") {
  grep(fastq_regex_pattern(allowed_compressions), files)
}

fastq_regex_pattern <- function(allowed_compressions = "gz") {
  paste0("\\.(fastq|fq)(", paste("\\.", allowed_compressions, sep = "", collapse = "|"), ")?$")
}
