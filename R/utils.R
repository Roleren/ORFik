#' Get the transcripts that have minimum lengths of leaders, cds and trailer.
#'
#' Filter transcripts to those who have 5' UTR, CDS, 3' UTR of some lengths,
#' pick the longest per gene.
#'
#' If a transcript does not have a 3' UTR, then the length is 0,
#' so they will be filtered out. So only transcripts with leaders, cds and
#' 3' UTRs will be returned. You can set the integer to 0, that will return all
#' within that group.
#' @param txdb a TxDb object from gtf
#' @param minFiveUTR (integer) minimum bp for 5' UTR during filtering for the
#' transcripts
#' @param minCDS (integer) minimum bp for CDS during filtering for the
#' transcripts
#' @param minThreeUTR (integer) minimum bp for 3' UTR during filtering for the
#' transcripts
#' @param stopOnEmpty logical TRUE, stop if no valid names are found ?
#' @return a character vector of valid tramscript names
#' @export
#' @examples
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
#' txNames <- filterTranscripts(txdb)
#'
filterTranscripts <- function(txdb, minFiveUTR = 30L, minCDS = 150L,
                              minThreeUTR = 30L, stopOnEmpty = TRUE) {
  if(!is(txdb, "TxDb")) stop("txdb must be a TxDb object")

  tx <- data.table::setDT(
    GenomicFeatures::transcriptLengths(
      txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
  tx <- tx[tx$utr5_len >= minFiveUTR & tx$cds_len >= minCDS &
             tx$utr3_len >= minThreeUTR, ]
  gene_id <- cds_len <- NULL
  data.table::setorder(tx, gene_id, -cds_len)
  tx <- tx[!duplicated(tx$gene_id), ]
  tx <- tx[!is.na(tx$gene_id)]
  if (stopOnEmpty & length(tx$tx_name) == 0)
    stop("No transcript has leaders and trailers of specified minFiveUTR",
         "minCDS, minThreeUTR")

  return(tx$tx_name)
}

#' Helper function to find overlaping seqlevels
#'
#' Useful to avoid warnings in bioC
#' @param grl a GRangesList or GRanges object
#' @param reads a GRanges or GAlignment object
#' @return a character vector of valid seqlevels
validSeqlevels <- function(grl, reads) {
  readNames <- unique(seqnames(reads))
  seqMatch <- readNames %in%
    unique(seqnamesPerGroup(grl, FALSE))
  return(readNames[seqMatch])
}

#' Helper function to check for GRangesList
#' @param class the class you want to check if is GRL,
#' either a character from class or the object itself.
#' @return a boolean
#' @family utils
#'
is.grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return((class == "GRangesList" || class == "CompressedGRangesList"))
}


#' Helper function to check for GRangesList or GRanges class
#' @param class the class you want to check if is GRL or GR,
#'  either a character from class or the object itself.
#' @return a boolean
#' @family utils
#'
is.gr_or_grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return(is.grl(class) || class == "GRanges")
}

#' Check if all requirements for an ORFik ORF is accepted.
#' @param grl a GRangesList or GRanges to check
#' @return a logical (TRUE/FALSE)
is.ORF <- function(grl){
  if (is.gr_or_grl(class(grl))){
     if (is.grl(grl)) {
       names <- unlist(grl[1], use.names = FALSE)$names
     } else {
         names <-grl[1]$names
     }
    return(any(grep("_", names)))
  }
  return(FALSE)
}

#' Helper Function to check valid GRangesList input
#' @param class as character vector the given class of
#'  supposed GRangesList object
#' @param type a character vector, is it gtf, cds, 5', 3', for messages.
#' @param checkNULL should NULL classes be checked and return indeces of these?
#' @return either NULL or indices (checkNULL == TRUE)
#' @family utils
#'
validGRL <- function(class, type = "grl", checkNULL = FALSE) {
  if(length(class) != length(type)) stop("not equal length of classes",
                                         " and types, see validGRL")
  if (checkNULL) {
    indeces <- "NULL" == class
    class <- class[!indeces]
    if (length(class) == 0) return(rep(TRUE, length(type)))
    type <- type[!indeces]
  }
  for (classI in seq_along(class)) {
    if (!is.grl(class[classI])) {
      messageI <- paste(type[classI], "must be given and be type GRangesList")
      stop(messageI)
    }
  }
  if (checkNULL) {
    return(indeces)
  }
}

#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An data.frame from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a GRanges object from bed
#' @family utils
#'
bedToGR <- function(x, bed6 = TRUE){

  if (!bed6) {
    gr <- GRanges(x[, 1L], IRanges(x[, 2L] + 1L, x[, 3L]))
    return(gr)
  }
  starts <- x[, 2L] + 1L
  ends <- x[, 3L]
  gr <- GRanges(x[, 1L], IRanges(starts, ends),
                strand = x[, 6L])
  mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = x[, 5L])
  if (ncol(x) > 6L) mcols(gr) <- x[, 7L:ncol(x)]
  return(gr)
}


#' Load bed file as GRanges.
#'
#' Wraps around rtracklayer::import.bed and tries to speed up loading with the
#' use of data.table. Supports gzip, gz, bgz and bed formats.
#' @param filePath The location of the bed file
#' @importFrom data.table fread setDF
#' @importFrom tools file_ext
#' @importFrom rtracklayer import.bed
#' @return a GRanges object
#' @export
#' @family utils
#' @examples
#' # path to example CageSeq data from hg19 heart sample
#' cageData <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'                         package = "ORFik")
#' fread.bed(cageData)
#'
fread.bed <- function(filePath) {

  if (.Platform$OS.type == "unix") {
    if (file.exists(filePath)) {
      if (any(file_ext(filePath) == c("gzip", "gz", "bgz"))) {
        bed <- bedToGR(setDF(
          fread(cmd = paste("gunzip -c", filePath), sep = "\t")))
      } else if (file_ext(filePath) == "bed"){
        bed <- bedToGR(setDF(fread(filePath, sep = "\t")))
      } else {
        bed <- import.bed(con =  filePath)
      }
    } else {stop("Filepath specified does not name existing file.")}
  } else {
    ## NB: Windows user will have slower loading
    bed <- import.bed(con =  filePath)
  }
  return(bed)
}

#' Convert a GRanges Object to 1 width reads
#'
#' There are 4 ways of doing this
#' 1. Take 5' ends, reduce away rest (5prime)
#' 2. Take 3' ends, reduce away rest (3prime)
#' 3. Tile and include all (tileAll)
#' 4. Take middle point per GRanges (middle)
#'
#' Many other ways to do this have their own functions, like startCodons and
#' stopCodons.
#' @param gr GRanges, GAlignment Object to reduce
#' @param method the method to reduce, see info. (5prime defualt)
#' @param addScoreColumn logical (FALSE), if TRUE, add a score column that
#'  sums up the hits per position.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read.
#' @return  Converted GRanges object
#' @family utils
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE){
  if (addSizeColumn) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr), size = readWidths(gr))
  }
  gr <- GRanges(gr)
  if (method == "5prime") {
    gr <- resize(gr, width = 1, fix = "start")
  } else if(method == "3prime") {
    gr <- resize(gr, width = 1, fix = "end")
  } else if(method == "tileAll") {
    gr <- unlist(tile(gr, width = 1), use.names = FALSE)

  } else if (method == "middle") {
    ranges(gr) <- IRanges(start(gr) + ceiling((end(gr) - start(gr)) / 2),
                          width = 1)
  } else stop("method not defined")

  if (addScoreColumn) {
    pos <- strandBool(gr)
    posGr <- gr[pos]
    dt <- as.data.table(posGr)[, .N, .(seqnames, start)]
    posGr <- GRanges(dt$seqnames, IRanges(dt$start, width = 1), "+")
    score <- dt$N
    negGr <- gr[!pos]
    dt <- as.data.table(negGr)[, .N, .(seqnames, end)]
    negGr <- GRanges(dt$seqnames, IRanges(dt$end, width = 1), "-")
    score <- as.integer(c(score, dt$N))

    gr <- c(posGr, negGr)
    gr$score <- NULL
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = score)
  }
  return(gr)
}

#' Convenience wrapper for Rsamtools FaFile
#' @param faFile a character path or FaFile
#' @importFrom Rsamtools FaFile
#' @return a FaFile or BSgenome
#' @family utils
#' @importFrom methods is
#'
findFa <- function(faFile) {
  if (is.character(faFile)) {
    if (file.exists(faFile)) {
      return(FaFile(faFile))
    } else {
      stop("faFile does not name a valid fasta/index file")
    }
  } else if (is(faFile, "FaFile") || is(faFile, "BSgenome")) {
    return(faFile)
  }
  stop("faFile must be FaFile, BSgenome or valid filePath")
}
