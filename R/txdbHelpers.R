#' Get new exon ids
#' @param txList a list, call of as.list(txdb)
#' @return a new valid ordered list of exon ids (integer)
remakeTxdbExonIds <- function(txList) {
  # remake exon ids
  DT <- data.table(start = txList$splicings$exon_start,
                   end = txList$splicings$exon_end,
                   chr = txList$splicings$exon_chr,
                   strand = txList$splicings$exon_strand,
                   id = seq.int(1,length(txList$splicings$exon_start)))
  setkeyv(DT, c("start", "end", "chr", "strand"))
  d <- duplicated(DT[,.(start, end, chr, strand)])
  c <- rep.int(1L, length(d))
  # find unique exon ids
  for(x in seq.int(2, length(d))){
    if (d[x]) {
      c[x] <- c[x-1L]
    } else {
      c[x] <- c[x-1L] + 1L
    }
  }
  return(as.integer(c[order(DT$id)]))
}

#' Update exon ranks of exon data.frame
#'
#' @param exons a data.frame, call of as.list(txdb)$splicings
#' @return a data.frame, modified call of as.list(txdb)
updateTxdbRanks <- function(exons) {

  exons$exon_rank <-  unlist(lapply(runLength(Rle(exons$tx_id)),
                         function(x) seq.int(1L, x)), use.names = FALSE)
  return(exons)

}

#' Remove exons in txList that are not in fiveUTRs
#'
#' @param txList a list, call of as.list(txdb)
#' @param fiveUTRs a GRangesList of 5' leaders
#' @return a list, modified call of as.list(txdb)
#' @importFrom data.table setDT
removeTxdbExons <- function(txList, fiveUTRs) {
  # remove old "dead" exons
  # get fiveUTR exons
  gr <- unlistGrl(fiveUTRs)
  dt <- as.data.table(gr)
  dt$names <- names(gr)
  # find the ones that does not start on 1
  d <-  dt[, .I[which.min(exon_rank)], by = names]
  d$ranks <- dt$exon_rank[d$V1]
  d <- d[ranks > 1L,]
  if(nrow(d) == 0) return(txList)
  # delete these in txList
  rankHits <- (txList$transcripts$tx_name[txList$splicings$tx_id] %in% d$names)
  e <- setDT(txList$splicings[rankHits,])
  f <- data.table(rank = e$exon_rank,
                  names = txList$transcripts$tx_name[e$tx_id])
  f$minRank <- rep.int(d$ranks, runLength(Rle(f$names)))
  f$remove <- f$minRank > f$rank

  txList$splicings <- txList$splicings[-which(rankHits)[f$remove],]
  txList$splicings <- updateTxdbRanks(txList$splicings)

  return(txList)
}

#' Remove specific transcripts in txdb List
#'
#' Remove all transcripts, except the ones in fiveUTRs.
#' @inheritParams updateTxdbStartSites
#' @return a txList
removeTxdbTranscripts <- function(txList, fiveUTRs) {
  # Transcripts
  match <- txList$transcripts$tx_name %in% names(fiveUTRs)
  ids <- which(match)
  txList$transcripts <- txList$transcripts[match, ]
  # Splicing
  match <- txList$splicings$tx_id %in% ids
  txList$splicings <- txList$splicings[match, ]
  # Genes
  match <- txList$genes$tx_id %in% ids
  txList$genes <- txList$genes[match, ]
  return(txList)
}

#' Update start sites of leaders
#'
#' @param txList a list, call of as.list(txdb)
#' @param fiveUTRs a GRangesList of 5' leaders
#' @param removeUnused logical (FALSE), remove leaders that did not have any
#' cage support. (standard is to set them to original annotation)
#' @return a list, modified call of as.list(txdb)
updateTxdbStartSites <- function(txList, fiveUTRs, removeUnused) {
  if (removeUnused) {
    txList <- removeTxdbTranscripts(txList, fiveUTRs)
  }
  txList <- removeTxdbExons(txList, fiveUTRs)

  starts <- startSites(fiveUTRs, keep.names = TRUE)

  # find all transcripts with 5' UTRs
  hitsTx <- txList$transcripts$tx_name %in% names(fiveUTRs)
  if(sum(hitsTx) != length(fiveUTRs)) {
    stop("Number of transcripts and leaders does not match, check naming!")
  }
  idTx <- which(hitsTx)
  # reassign starts for positive strand
  txList$transcripts$tx_start[hitsTx][strandBool(fiveUTRs)] <-
    starts[strandBool(fiveUTRs)]
  txList$splicings$exon_start[(txList$splicings$tx_id %in% idTx) &
                                (txList$splicings$exon_rank == 1) &
                                (txList$splicings$exon_strand == "+")] <-
    starts[strandBool(fiveUTRs)]
  # reassign stops for negative strand
  txList$transcripts$tx_end[hitsTx][!strandBool(fiveUTRs)] <-
    starts[!strandBool(fiveUTRs)]
  txList$splicings$exon_end[(txList$splicings$tx_id %in% idTx) &
                              (txList$splicings$exon_rank == 1) &
                              (txList$splicings$exon_strand == "-")] <-
    starts[!strandBool(fiveUTRs)]


  return(txList)
}

#' General loader for txdb
#'
#' Useful to allow fast TxDb loader like .db
#' @param txdb a TxDb file or a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite)
#' @return a TxDb object
#' @importFrom AnnotationDbi loadDb
#' @export
#' @examples
#' library(GenomicFeatures)
#' # Get the gtf txdb file
#' txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                         package = "GenomicFeatures")
#' txdb <- loadDb(txdbFile)
#'
loadTxdb <- function(txdb) {
  #TODO: Check that is is an open connection!
  if (is(txdb, "character")) {
    f <- file_ext(txdb)
    if (f == "gff" | f == "gff2" | f == "gff3" | f == "gtf") {
      txdb <- GenomicFeatures::makeTxDbFromGFF(txdb)
    } else if(f == "db" | f == "sqlite") {
      txdb <- loadDb(txdb)
    } else stop("when txdb is path, must be one of .gff, .gtf and .db")

  } else if(!is(txdb, "TxDb")) stop("txdb must be path or TxDb")
  return(txdb)
}

#' Load transcript region
#'
#' Load GRangesList if input is not already GRangesList.
#' @param txdb a TxDb file or a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite), if it is a GRangesList,
#'  it will return it self.
#' @param part a character, one of: tx, leader, cds, trailer, intron, mrna
#' NOTE: difference between tx and mrna is that tx are all transcripts, while
#' mrna are all transcripts with a cds
#' @return a GrangesList of region
#' @export
#' @examples
#' gtf <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' loadRegion(gtf, "intron")
loadRegion <- function(txdb, part = "tx") {
  if (is.grl(txdb)) return(txdb)
  txdb <- loadTxdb(txdb)
  if (part == "tx") {
    return(exonsBy(txdb, by = "tx", use.names = TRUE))
  } else if (part %in% c("leader", "leaders")) {
    return(fiveUTRsByTranscript(txdb, use.names = TRUE))
  } else if (part == "cds") {
    return(cdsBy(txdb, by = "tx", use.names = TRUE))
  } else if (part %in% c("trailer", "trailers")) {
    return(threeUTRsByTranscript(txdb, use.names = TRUE))
  } else if (part %in% c("intron", "introns")) {
    return(intronsByTranscript(txdb, use.names = TRUE))
  }  else if (part %in% c("mrna", "mRNA")) {
    txNames <- filterTranscripts(txdb, 0, 1, 0, FALSE)
    return(exonsBy(txdb, by = "tx", use.names = TRUE)[txNames])
  } else stop("invalid: must be tx, leader, cds, trailer, introns or mrna")
}

#' Load transcripts of given biotype
#'
#' Like rRNA, snoRNA etc.
#' NOTE: Only works on gtf/gff, not .db object for now.
#' Also note that these anotations are not perfect, some rRNA annotations
#' only contain 5S rRNA etc. If your gtf does not contain evertyhing you need,
#' use a resource like repeatmasker and download a gtf:
#' https://genome.ucsc.edu/cgi-bin/hgTables
#' @references doi: 10.1002/0471250953.bi0410s25
#' @param path path to gtf/gff
#' @param part a character, default rRNA. Can also be:
#' snoRNA, tRNA etc. As long as that biotype is defined in the gtf.
#' @param tx a GRangesList of transcripts (Optional, default NULL),
#'  add to save run time.
#' @return a GRangesList of transcript of that type
loadTranscriptType <- function(path, part = "rRNA", tx = NULL) {
  if (!is.character(path)) stop("path must be a file path to gtf/gff")
  type <- import(path)

  valids <- type[grep(x = type$transcript_biotype, pattern = part)]
  if (length(valids) == 0) stop("found no valid transcript of type", part)
  if (is.null(tx)) tx <- loadRegion(path)

  return(tx[unique(valids$transcript_id)])
}

#' Convert transcript names to gene names
#'
#' Works for ensembl etc.
#' @param txNames character vector, the transcript names to convert.
#' @param txdb the transcript database to use or gtf/gff path to it.
#' @return character vector of gene names
txNamesToGeneNames <- function(txNames, txdb) {
  txdb <- loadTxdb(txdb)
  g <- mcols(transcripts(txdb, columns = c("tx_name", "gene_id")))
  match <- g[g$tx_name %in% txNames,]
  return(as.character(match$gene_id))
}

#' Filter transcripts by lengths
#'
#' Filter transcripts to those who have leaders, CDS, trailers of some lengths,
#' you can also pick the longest per gene.
#'
#' If a transcript does not have a trailer, then the length is 0,
#' so they will be filtered out if you minThreeUTR to 1. So only transcripts
#' with leaders, cds and trailers will be returned. You can set the integer
#' to 0, that will return all within that group.
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
#' @param stopOnEmpty logical TRUE, stop if no valid transcripts are found ?
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
