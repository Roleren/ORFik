#' Get the transcripts with accepted lengths of leaders, cds and trailer.
#'
#' Filter transcripts to those who have leaders, CDS, trailers of some lengths,
#' you can also pick the longest per gene.
#'
#' If a transcript does not have a trailer, then the length is 0,
#' so they will be filtered out. So only transcripts with leaders, cds and
#' trailers will be returned. You can set the integer to 0, that will return
#' all within that group.
#'
#' If your annotation does not have leaders or trailers, set them to NULL.
#' @inheritParams loadTxdb
#' @param minFiveUTR (integer) minimum bp for 5' UTR during filtering for the
#' transcripts. Set to NULL if no 5' UTRs exists for annotation.
#' @param minCDS (integer) minimum bp for CDS during filtering for the
#' transcripts
#' @param minThreeUTR (integer) minimum bp for 3' UTR during filtering for the
#' transcripts. Set to NULL if no 3' UTRs exists for annotation.
#' @param longestPerGene logical (TRUE), return only longest valid transcript
#' per gene.
#' @param stopOnEmpty logical TRUE, stop if no valid names are found ?
#' @return a character vector of valid tramscript names
#' @export
#' @examples
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
#' txNames <- filterTranscripts(txdb)
#'
filterTranscripts <- function(txdb, minFiveUTR = 30L, minCDS = 150L,
                              minThreeUTR = 30L, longestPerGene = TRUE,
                              stopOnEmpty = TRUE) {
  txdb <- loadTxdb(txdb)
  five <- !is.null(minFiveUTR)
  three <- !is.null(minThreeUTR)

  tx <- data.table::setDT(
    GenomicFeatures::transcriptLengths(
      txdb, with.cds_len = TRUE, with.utr5_len = five, with.utr3_len = three))
  five <- rep(five, nrow(tx))
  three <- rep(three, nrow(tx))

  tx <- tx[ifelse(five, utr5_len >= minFiveUTR, TRUE) & cds_len >= minCDS &
             ifelse(three, utr3_len >= minThreeUTR, TRUE), ]

  gene_id <- cds_len <- NULL
  data.table::setorder(tx, gene_id, -cds_len)
  if (longestPerGene) {
    tx <- tx[!duplicated(tx$gene_id), ]
  }
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

#' Find optimized subset of valid reads
#' @inheritParams validSeqlevels
#' @return the reads as GRanges or GAlignment
optimizeReads <- function(grl, reads) {
  seqMatch <- validSeqlevels(grl, reads)
  reads <- keepSeqlevels(reads, seqMatch, pruning.mode = "coarse")
  if (length(reads) > 1e6) { # speedup on big riboseq libraries
    reads <- reads[countOverlaps(reads, grl, type = "within") > 0]
    reads <- sort(reads)
  }
  return(reads)
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
#' Custome bam reader
#'
#' Safer version that handles the most important error done.
#' In the future will use a faster .bam loader for big .bam files in R.
#' @param path a character path to .bam file
#' @param chrStyle a GRanges object, or a characterPath (Default: NULL) to
#' get seqlevelsStyle from. Is chromosome 1 called chr1 or 1, is mitocondrial
#' chromosome called MT or chrM etc.
#' @return a GAlignment object of bam file
readBam <- function(path, chrStyle = NULL) {
  reads <- readGAlignments(path)
  if (!is.null(chrStyle)) {
    if (is.character(chrStyle)) {
      seqlevelsStyle(reads) <- chrStyle
    } else if (is.gr_or_grl(chrStyle)) {
      seqlevelsStyle(reads) <- seqlevelsStyle(chrStyle)
    } else stop("chrStyle must be valid GRanges object, or a valid chr style!")
  }
}

#' Convert a GRanges Object to 1 width reads
#'
#' There are 4 ways of doing this
#' 1. Take 5' ends, reduce away rest (5prime)
#' 2. Take 3' ends, reduce away rest (3prime)
#' 3. Tile to 1-mers and include all (tileAll)
#' 4. Take middle point per GRanges (middle)
#'
#'
#' Many other ways to do this have their own functions, like startSites and
#' stopSites etc.
#' @param gr GRanges, GAlignment Object to reduce
#' @param method the method to reduce, see info. (5prime defualt)
#' @param addScoreColumn logical (FALSE), if TRUE, add a score column that
#'  sums up the hits per unique range This will make each read unique, so
#'  that each read is 1 time, and score column gives the number of hits.
#'  A useful compression.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read. Useful if you need
#'  original read lengths. This takes care of soft clips etc.
#' @return  Converted GRanges object
#' @family utils
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE){
  if (addSizeColumn) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr), size = readWidths(gr))
  }
  if (addScoreColumn) {
    if (!is.null(cigar(reads))) {
      dt <- as.data.table(gr)[, .N, .(seqnames, start, end, strand, cigar)]
    } else {
      dt <- as.data.table(gr)[, .N, .(seqnames, start, end, strand)]
    }
    gr$score <- NULL
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = as.integer(dt$N))
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

  return(gr)
}

#' Convenience wrapper for Rsamtools FaFile
#' @param faFile FaFile, BSgenome or fasta/index file path used to find the
#'  transcripts
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
