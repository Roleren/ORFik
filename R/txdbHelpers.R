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
  if (is(txdb, "character")) {
    f <- file_ext(txdb)
    if ( f == "gff" | f == "gff2" | f == "gff3" | f == "gtf") {
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
#' @param part a character, one of: tx, leader, cds, trailer, introns
#' @return a GrangesList of region
loadRegion <- function(txdb, part = "tx") {
  if (is.grl(txdb)) return(txdb)
  txdb <- loadTxdb(txdb)
  if (part == "tx") {
    return(exonsBy(txdb, by = "tx", use.names = TRUE))
  } else if (part == "leader") {
    return(fiveUTRsByTranscript(txdb, use.names = TRUE))
  } else if(part == "cds") {
    return(cdsBy(txdb, by = "tx", use.names = TRUE))
  } else if(part == "trailer") {
    return(threeUTRsByTranscript(txdb, use.names = TRUE))
  } else if(part == "introns") {
    return(intronsByTranscript(txdb)(txdb, use.names = TRUE))
  } else stop("invalid part, must be tx, leader, cds or trailer")
}

#' Load type of transcript
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
#' snoRNA, tRNA etc. As long as that type is defined in the gtf.
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
