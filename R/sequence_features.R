#' Get distances between ORF Start and TSS of its transcript
#'
#' Matching is done by transcript names.
#' This is applicable practically to the upstream (fiveUTRs) ORFs.
#' If ORF is not within specified search space in fiveUTRs, this function
#' will crash.
#' @references doi: 10.1074/jbc.R116.733899
#' @param ORFs orfs as \code{\link{GRangesList}},
#' names of orfs must be txname_[rank]
#' @param fiveUTRs 5' leaders as \code{\link{GRangesList}}.
#' @return an integer vector, 1 means on TSS, 2 means second base of Tx.
#' @family features
#' @export
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(5, 10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(2, 20), "+"))
#' distToTSS(grl, fiveUTRs)
#'
distToTSS <- function(ORFs, fiveUTRs){
  validGRL(class(ORFs), "ORFs")

  startSites <-  startSites(ORFs, asGR = TRUE, keep.names = TRUE,
                            is.sorted = TRUE)
  return(start(asTX(startSites, fiveUTRs)))
}

#' Get distances between ORF ends and starts of their transcripts cds'.
#'
#' Will calculate distance between each ORF end and begining of the
#' corresponding cds. Matching is done by transcript names.
#' This is applicable practically to the upstream (fiveUTRs) ORFs.
#' The cds start site, will be presumed to be on + 1 of end of fiveUTRs.
#' @references doi: 10.1074/jbc.R116.733899
#' @param ORFs orfs as \code{\link{GRangesList}},
#' names of orfs must be transcript names
#' @param fiveUTRs fiveUTRs as \code{\link{GRangesList}},
#' remember to use CAGE version of 5' if you did CAGE reassignment!
#' @param cds cds' as \code{\link{GRangesList}},
#' only add if you have ORFs going into CDS.
#' @return an integer vector, +1 means one base upstream of cds, -1 means
#' 2nd base in cds, 0 means orf stops at cds start.
#' @family features
#' @export
#' @examples
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1, 10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1, 20), "+"))
#' distToCds(grl, fiveUTRs)
#'
distToCds <- function(ORFs, fiveUTRs, cds = NULL){
  validGRL(class(ORFs), "ORFs")

  cdsStarts <- widthPerGroup(fiveUTRs[
    txNames(ORFs)], FALSE) + 1

  lastExons <- lastExonPerGroup(ORFs)
  if (is.grl(cds)) {
    fiveUTRs <- addCdsOnLeaderEnds(fiveUTRs, cds)
  }
  orfsTx <- asTX(lastExons, fiveUTRs)

  # this is ok, since it is tx not genomic ->
  orfEnds <- lastExonEndPerGroup(orfsTx, FALSE)

  return(cdsStarts - orfEnds)
}


#' Make a score for each ORFs start region by proximity to Kozak
#'
#' The closer the sequence is to the Kozak sequence
#' the higher the score, based on the experimental pwms
#' from article referenced.
#' Minimum score is 0 (worst correlation), max is 1 (the best
#' base per column was chosen).
#' @references doi: https://doi.org/10.1371/journal.pone.0108475
#' @param grl a \code{\link{GRangesList}} grouped by ORF
#' @param faFile a FaFile from the fasta file, see ?FaFile.
#'  Can also be path to fastaFile with fai file in same dir.
#' @param species ("human"), which species to use,
#' currently supports human, zebrafish and mouse (m. musculus).
#' You can also specify a pfm for your own species.
#' Syntax of pfm is an rectangular integer matrix,
#' where all columns must sum to the same value, normally 100.
#' See example for more information.
#' Rows are in order: c("A", "C", "G", "T")
#' @param include.N logical (F), if TRUE, allow N bases to be counted as hits,
#' score will be average of the other bases. If True, N bases will be
#' added to pfm, automaticly, so dont include them if you make your own pfm.
#' @return a numeric vector with values between 0 and 1
#' @return an integer vector, one score per orf
#' @family features
#' @importFrom Biostrings PWM
#' @export
#' @examples
#' # Usually the ORFs are found in orfik, which makes names for you etc.
#' # Here we make an example from scratch
#' seqName <- "Chromosome"
#' ORF1 <- GRanges(seqnames = seqName,
#'                    ranges = IRanges(c(1007, 1096), width = 60),
#'                    strand = c("+", "+"))
#' ORF2 <- GRanges(seqnames = seqName,
#'                     ranges = IRanges(c(400, 100), width = 30),
#'                     strand = c("-", "-"))
#' ORFs <- GRangesList(tx1 = ORF1, tx2 = ORF2)
#' ORFs <- makeORFNames(ORFs) # need ORF names
#' # get faFile for sequences
#' faFile <- FaFile(system.file("extdata", "genome.fasta", package = "ORFik"))
#' kozakSequenceScore(ORFs, faFile)
#' # For more details see vignettes.
kozakSequenceScore <- function(grl, faFile, species = "human",
                               include.N = FALSE) {
  faFile <- findFa(faFile)
  firstExons <- firstExonPerGroup(grl)
  kozakLocation <- promoters(firstExons, upstream = 9, downstream = 6)

  sequences <- as.character(txSeqsFromFa(kozakLocation,
                                         faFile, is.sorted = TRUE))
  if (!all(nchar(sequences) == 15)) {
    stop("not all ranges had valid kozak sequences length, not 15")
  }

  if(is(species, "matrix")){
    # self defined pfm
    pfm <- species
  } else if (species == "human") {
    # human pfm, see article reference
    pfm <- t(matrix(as.integer(c(20,20,21,21,19,24,46,29,19,22,28,16,
                                 27,33,32,23,32,38,10,38,45,15,39,26,
                                 35,29,28,39,30,26,37,20,28,49,18,37,
                                 18,18,19,17,19,12,7,13,8,14,15,21)),
                    ncol = 4))
  } else if (species == "mouse") {
    # zebrafish pfm, see article reference
    pfm <- t(matrix(as.integer(c(20,19,21,20,18,25,49,28,17,23,28,15,
                                 27,34,31,23,32,38,9,39,47,14,40,26,
                                 34,28,27,39,29,25,36,20,28,49,18,37,
                                 19,19,21,18,21,12,6,13,8,14,14,22)),
                    ncol = 4))
  } else if (species == "zebrafish") {
    # zebrafish pfm, see article reference
    pfm <- t(matrix(as.integer(c(29,26,28,26,22,35,62,39,28,24,27,17,
                                 21,26,24,16,28,32,5,23,35,12,42,21,
                                 25,24,22,33,22,19,28,17,27,47,16,34,
                                 25,24,26,25,28,14,5,21,10,17,15,28)),
                    ncol = 4))
  } else  {
    stop("Either input species as a matrix
         or name of presupported pfm organism")
  }

  bases <- c("A", "C", "G", "T")
  rownames(pfm) <- bases
  pwm <- PWM(pfm)

  if (include.N) {
    bases <- c(bases, "N")
    pwm <- rbind(pwm, colMeans(pwm))
    rownames(pwm) <- bases
  }

  # exclude start codon
  s <- paste0(substr(x = sequences, 1, 9), substr(x = sequences, 13, 15))
  # split strings and relist as letters of 9 rows
  subSplit <- strsplit(s, split = "")
  # this will not when ATG is on start of chr
  mat <- t(matrix(unlist(subSplit, use.names = FALSE), ncol = length(s)))

  scores <- rep(0., length(s))
  for (i in seq(ncol(mat))) {
    for (n in seq_along(bases)) {
      match <- mat[, i] == bases[n]
      scores[match] <- scores[match] + pwm[n, i]
    }
  }
  return(scores)
}

#' Find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function
#' [distToCds()]
#'
#' possible outputs:
#' 0: orf is in frame with cds
#' 1: 1 shifted from cds
#' 2: 2 shifted from cds
#'
#' @references doi: 10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @return a logical vector
#' @family features
#' @examples
#' # simple example
#' isInFrame(c(3,6,8,11,15))
#'
#' # GRangesList example
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs)
#' isInFrame <- isInFrame(dist)
#' @export
#'
isInFrame <- function(dists){
  return((dists - 1) %% 3)
}


#' Find frame for each orf relative to cds
#'
#' Input of this function, is the output of the function
#' [distToCds()]
#' @references doi: 10.1074/jbc.R116.733899
#' @param dists a vector of distances between ORF and cds
#' @family features
#' @return a logical vector
#' @export
#' @examples
#' # simple example
#' isOverlapping(c(-3,-6,8,11,15))
#'
#' # GRangesList example
#' grl <- GRangesList(tx1_1 = GRanges("1", IRanges(1,10), "+"))
#' fiveUTRs <- GRangesList(tx1 = GRanges("1", IRanges(1,20), "+"))
#' dist <- distToCds(grl, fiveUTRs)
#' isOverlapping <- isOverlapping(dist)
isOverlapping <- function(dists) {
  return(dists < 0)
}


#' ORF rank in transcripts
#'
#' @references doi: 10.1074/jbc.R116.733899
#' @description ig. second orf _2 -> 2
#' @param grl a \code{\link{GRangesList}} object with ORFs
#' @return a numeric vector of integers
#' @family features
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' grl <- ORFik:::makeORFNames(grl)
#' rankOrder(grl)
rankOrder <- function(grl) {
  gr <- unlist(grl, use.names = FALSE)

  if (is.null(names(grl))) {
    if (is.null(gr$names)) {
      if (is.null(names(gr))) {
        stop("no valid names to find ranks")
      } else {
        orfName <- names(gr)
        if (length(orfName) > length(grl)) {
          orfName <- names(groupGRangesBy(gr, names(gr)))
        }
      }
    } else {
      orfName <- gr$names
      if (length(orfName) > length(grl)) {
        orfName <- names(groupGRangesBy(gr, gr$names))
      }
    }
  } else {
    orfName <- names(grl)
    if (suppressWarnings(anyNA(as.integer(sub(".*_", "", orfName,
                                              perl = TRUE))))) {
      if (!is.null(gr$names)) {
        orfName <- names(groupGRangesBy(gr, gr$names))
      }
    }
  }
  if (length(orfName) > length(grl)) {
    stop("did not find a valid column to find ranks, easiest way to fix is",
         " set grl to: ORFik:::groupGRangesBy(grl, names), ",
         "where names are the orf names with _* in them-")
  }

  if (is.null(orfName)) stop("grl must have column called names")
  ranks <- as.integer(sub(".*_", "", orfName, perl = TRUE))
  if (anyNA(ranks)) {
    stop("no valid names to find ranks, check for orf _* names eg.",
         "tx_1, tx_2.")
  }
  return(ranks)
}

#' Fraction Length
#' @description Fraction Length is defined as
#' \preformatted{(lengths of grl)/(length of tx_len)}
#' so that each group in
#' the grl is divided by the corresponding transcript.
#' @references doi: 10.1242/dev.098343
#' @param grl a \code{\link{GRangesList}} object
#' with usually either leaders,
#' cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param tx_len the transcript lengths of the transcripts,
#' a named (tx names) vector of integers.
#' If you have the transcripts as GRangesList,
#' call `ORFik:::widthPerGroup(tx, TRUE)`.
#'
#' If you used CageSeq to reannotate leaders, then the tss for the the leaders
#' have changed, therefore the tx lengths have changed. To account for that
#' call: `tx_len <- widthPerGroup(extendLeaders(tx, cageFiveUTRs))`
#' and calculate graction length using `fractionLength(grl, tx_len)`.
#' @return a numeric vector of ratios
#' @family features
#' @export
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' # grl must have same names as cds + _1 etc, so that they can be matched.
#' tx <-  GRangesList(tx1 = GRanges("1", IRanges(1, 50), "+"))
#' fractionLength(grl, ORFik:::widthPerGroup(tx, keep.names = TRUE))
#'
fractionLength <- function(grl, tx_len){
  grl_len <- widthPerGroup(grl, FALSE)
  tx_len <- tx_len[txNames(grl)]
  names(tx_len) <- NULL
  return(grl_len / tx_len)
}