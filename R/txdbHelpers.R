#' Get new exon ids after update of txdb
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

#' Update exon ranks of exon data.frame inside txdb object
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
#' @param txdb a TxDb file, a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite)
#'  or an ORFik experiment
#' @inheritParams matchSeqStyle
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
loadTxdb <- function(txdb, chrStyle = NULL) {
  if (is(txdb, "experiment")) {
    txdb <- txdb@txdb
  }
  #TODO: Check that is is an open connection!
  if (is(txdb, "character")) {
    f <- file_ext(txdb)
    if (f == "gff" | f == "gff2" | f == "gff3" | f == "gtf") {
      txdb <- GenomicFeatures::makeTxDbFromGFF(txdb)
    } else if(f == "db" | f == "sqlite") {
      txdb <- loadDb(txdb)
    } else stop("when txdb is path, must be one of .gff, .gtf and .db")

  } else if(!is(txdb, "TxDb")) stop("txdb must be path or TxDb")
  return(matchSeqStyle(txdb, chrStyle))
}

#' Load transcript region
#'
#' Usefull to simplify loading of standard regions, like cds' and leaders.
#'
#' Load as GRangesList if input is not already GRangesList.
#' @param txdb a TxDb file or a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite), if it is a GRangesList,
#'  it will return it self.
#' @param part a character, one of: tx, leader, cds, trailer, intron, mrna
#' NOTE: difference between tx and mrna is that tx are all transcripts, while
#' mrna are all transcripts with a cds
#' @param names.keep a character vector of subset of names to keep. Example:
#' loadRegions(txdb, names = ENST1000005), will return only that transcript.
#' Remember if you set by to "gene", then this list must be with gene names.
#' @param by a character, default "tx" Either "tx" or "gene". What names to
#' output region by, the transcript name "tx" or gene names "gene"
#' @return a GrangesList of region
#' @export
#' @examples
#' gtf <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' loadRegion(gtf, "cds")
#' loadRegion(gtf, "intron")
loadRegion <- function(txdb, part = "tx", names.keep = NULL, by = "tx") {
  if (is.grl(txdb)) return(txdb)
  if (length(part) != 1) stop("argument: (path) must be length 1")
  txdb <- loadTxdb(txdb)
  region <-
    if (part %in% c("tx", "transcript", "transcripts")) {
      exonsBy(txdb, by = "tx", use.names = TRUE)
    } else if (part %in% c("leader", "leaders", "5'", "5", "5utr",
                           "fiveUTRs", "5pUTR")) {
      fiveUTRsByTranscript(txdb, use.names = TRUE)
    } else if (part %in% c("cds", "CDS", "mORF")) {
      cdsBy(txdb, by = "tx", use.names = TRUE)
    } else if (part %in% c("trailer", "trailers", "3'", "3", "3utr",
                           "threeUTRs", "3pUTR")) {
      threeUTRsByTranscript(txdb, use.names = TRUE)
    } else if (part %in% c("intron", "introns")) {
      intronsByTranscript(txdb, use.names = TRUE)
    }  else if (part %in% c("mrna", "mrnas", "mRNA", "mRNAs")) {
      loadRegion(txdb, "tx")[names(cdsBy(txdb, use.names = TRUE))]
    } else stop("invalid: must be tx, leader, cds, trailer, introns or mrna")

  if (by == "gene") {
    names(region) <- txNamesToGeneNames(names(region), txdb)
  }
  if (!is.null(names.keep)) { # If subset
    subset <- names(region) %in% names.keep
    if (length(subset) == 0)
      stop(paste("Found no transcripts kepts, for region:", part))
    region <- region[subset]
  }
  return(region)
}

#' Get all regions of transcripts specified to environment
#'
#' By default loads all parts to .GlobalEnv (global environemnt)
#' Useful to not spend time on finding the functions to load regions.
#' @inheritParams loadTxdb
#' @inheritParams loadRegion
#' @param parts the transcript parts you want, default:
#' c("mrna", "leaders", "cds", "trailers").\cr See ?loadRegion for more info
#' on this argument.
#' @param extension What to add on the name after leader, like: B -> leadersB
#' @param envir Which environment to save to, default: .GlobalEnv
#' @return invisible(NULL) (regions saved in envir)
#' @export
#' @examples
#' # Load all mrna regions to Global environment
#' gtf <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' loadRegions(gtf, parts = c("mrna", "leaders", "cds", "trailers"))
loadRegions <- function(txdb, parts = c("mrna", "leaders", "cds", "trailers"),
                        extension = "", names.keep = NULL,
                        by = "tx",
                        envir = .GlobalEnv) {
  txdb <- loadTxdb(txdb)
  for (i in parts) {
    assign(x = paste0(i, extension),
           value = loadRegion(txdb, i,names.keep, by),
           envir = envir)
  }
  return(invisible(NULL))
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
#' @param object a TxDb, ORFik experiment or path to gtf/gff,
#' @param part a character, default rRNA. Can also be:
#' snoRNA, tRNA etc. As long as that biotype is defined in the gtf.
#' @param tx a GRangesList of transcripts (Optional, default NULL,
#' all transcript of that type), else it must be names a list to subset on.
#' @return a GRangesList of transcript of that type
loadTranscriptType <- function(object, part = "rRNA", tx = NULL) {
  type <- importGtfFromTxdb(object)

  valids <- type[grep(x = type$transcript_biotype, pattern = part)]
  if (length(valids) == 0) stop("found no valid transcript of type", part)
  if (is.null(tx)) tx <- loadRegion(object)

  return(tx[unique(valids$transcript_id)])
}

#' Convert transcript names to gene names
#'
#' Works for ensembl, UCSC and other standard annotations.
#' @param txNames character vector, the transcript names to convert. Can also be
#' a named object with tx names (like a GRangesList), will then extract names.
#' @param txdb the transcript database to use or gtf/gff path to it.
#' @return character vector of gene names
#' @export
#' @examples
#' gtf <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- loadTxdb(gtf)
#' loadRegions(txdb, "cds") # using tx names
#' txNamesToGeneNames(cds, txdb)
#' # Identical to:
#' loadRegions(txdb, "cds", by = "gene")
txNamesToGeneNames <- function(txNames, txdb) {
  if (!is(txNames, "character")) txNames <- names(txNames)

  txdb <- loadTxdb(txdb)
  g <- mcols(transcripts(txdb, columns = c("tx_name", "gene_id")))
  match <- chmatch(txNames, g$tx_name)
  if (anyNA(match)) stop("Not all txNames exists in txdb, check for spelling errors!")
  return(as.character(g$gene_id)[match])
}

#' Import the GTF / GFF that made the txdb
#'
#' @param txdb a TxDb, path to txdb / gff or ORFik experiment object
#' @return data.frame, the gtf/gff object imported with rtracklayer::import
importGtfFromTxdb <- function(txdb) {
  if (is(txdb, "experiment")) txdb <- txdb@txdb
  if(is(txdb, "character")) {
    if (!(file_ext(txdb) %in% c("gtf", "gff", "gff3", "gff2"))) {
      # It means it shoud be a TxDb path
      txdb <- loadTxdb(txdb)
    }
  }
  if (is(txdb, "TxDb")) {
    if (!(file_ext(metadata(txdb)[3,2]) %in%
          c("gtf", "gff", "gff3", "gff2"))) {
      message("This is error txdb ->")
      message("It should be found by: metadata(txdb)[3,2]")
      print(txdb)
      stop("Could not find valid gtf / gff file, only data base object!")
    }
    txdb <- metadata(txdb)[3,2]
  }
  return(import(txdb))
}


#' Filter transcripts by lengths
#'
#' Filter transcripts to those who have leaders, CDS, trailers of some lengths,
#' you can also pick the longest per gene.
#'
#' If a transcript does not have a trailer, then the length is 0,
#' so they will be filtered out if you set minThreeUTR to 1.
#' So only transcripts with leaders, cds and trailers will be returned.
#' You can set the integer to 0, that will return all within that group.
#'
#' If your annotation does not have leaders or trailers, set them to NULL,
#' since 0 does mean there must exist a column called utr3_len etc.
#' Genes with gene_id = NA will be be removed.
#' @inheritParams loadRegion
#' @param minFiveUTR (integer) minimum bp for 5' UTR during filtering for the
#' transcripts. Set to NULL if no 5' UTRs exists for annotation.
#' @param minCDS (integer) minimum bp for CDS during filtering for the
#' transcripts
#' @param minThreeUTR (integer) minimum bp for 3' UTR during filtering for the
#' transcripts. Set to NULL if no 3' UTRs exists for annotation.
#' @param longestPerGene logical (TRUE), return only longest valid transcript
#' per gene.
#' @param stopOnEmpty logical TRUE, stop if no valid transcripts are found ?
#' @return a character vector of valid transcript names
#' @export
#' @examples
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
#' txNames <- filterTranscripts(txdb, minFiveUTR = 1, minCDS = 30,
#'                              minThreeUTR = 1)
#' loadRegion(txdb, "mrna")[txNames]
#' loadRegion(txdb, "5utr")[txNames]
#'
filterTranscripts <- function(txdb, minFiveUTR = 30L, minCDS = 150L,
                              minThreeUTR = 30L, longestPerGene = TRUE,
                              stopOnEmpty = TRUE,
                              by = "tx") {
  if (!(by %in% c("tx", "gene"))) stop("by must be either tx or gene!")
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

  if (by == "gene") return(tx$gene_id)
  return(tx$tx_name)
}
