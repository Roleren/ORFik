
#' Make txdb from genome
#'
#' Make a Txdb with defined seqlevels and
#' seqlevelsstyle from the fasta genome.
#' This makes it more fail safe than standard Txdb creation.
#' Example is that you can not create a coverage window outside
#' the chromosome boundary, this is only possible if you have
#' set the seqlengths.
#' @param gtf path to gtf file
#' @param genome character, default NULL. Path to fasta genome
#' corresponding to the gtf. If NULL, can not set seqlevels.
#' If value is NULL or FALSE, it will be ignored.
#' @param organism Scientific name of organism, first letter
#' must be capital! Example: Homo sapiens. Will force first letter
#' to capital and convert any "_" (underscore) to " " (space)
#' @param optimize logical, default FALSE. Create a folder
#' within the folder of the gtf, that includes optimized objects
#' to speed up loading of annotation regions from up to 15 seconds
#' on human genome down to 0.1 second. ORFik will then load these optimized
#' objects instead. Currently optimizes filterTranscript() function and
#' loadRegion() function for 5' UTRs, 3' UTRs, CDS,
#'  mRNA (all transcript with CDS) and tx (all transcripts).
#' @return NULL, Txdb saved to disc named paste0(gtf, ".db")
#' @export
#' @examples
#' gtf <- "/path/to/local/annotation.gtf"
#' genome <- "/path/to/local/genome.fasta"
#' #makeTxdbFromGenome(gtf, genome, organism = "Saccharomyces cerevisiae")
makeTxdbFromGenome <- function(gtf, genome = NULL, organism,
                               optimize = FALSE) {
  message("Making txdb of GTF")
  organismCapital <- paste0(toupper(substr(organism, 1, 1)),
                            substr(organism, 2, nchar(organism)))
  organismCapital <- gsub("_", " ", organismCapital)

  if (!is.logical(genome) & !is.null(genome)) {
    fa <- FaFile(genome)
    fa.seqinfo <- seqinfo(fa)
    if ("MT" %in% names(fa.seqinfo)) { # Ensembl
      isCircular(fa.seqinfo)["MT"] <- TRUE
    } else if ("chrM" %in% names(fa.seqinfo)) { # UCSC
      isCircular(fa.seqinfo)["chrM"] <- TRUE
    } else if ("Mito" %in% names(fa.seqinfo)) { # UCSC
      isCircular(fa.seqinfo)["Mito"] <- TRUE
    }

    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, organism = organismCapital,
                                             chrominfo = fa.seqinfo)
    if (seqlevelsStyle(txdb)[1] != seqlevelsStyle(fa)[1]) {
      seqlevelsStyle(txdb) <- seqlevelsStyle(fa)[1]
    }

  } else {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, organism = organismCapital)
  }

  txdb_file <- paste0(gtf, ".db")
  AnnotationDbi::saveDb(txdb, txdb_file)
  message("Txdb stored at: ", txdb_file)

  if (optimize) {
    base_path <- optimized_txdb_path(txdb, create.dir = TRUE)
    message("--------------------------")
    message("Optimizing annotation, saving to: ", dirname(base_path))
    # Save all transcript, cds and UTR lengths as .fst
    optimizedTranscriptLengths(txdb, create.fst.version = TRUE)
    # Save RDS version of all transcript regions
    message("Creating rds speedup files for transcript regions")
    parts <- c("tx", "mrna", "leaders", "cds", "trailers")
    for (i in parts) {
      saveRDS(loadRegion(txdb, i, by = "tx"),
              file = paste0(base_path, "_", i, ".rds"))
    }
  }
  return(invisible(NULL))
}

#' Get new exon ids after update of txdb
#' @param txList a list, call of as.list(txdb)
#' @return a new valid ordered list of exon ids (integer)
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' Adds another safety in that seqlevels will be set
#'
#' Load as GRangesList if input is not already GRangesList.
#' @param txdb a TxDb file or a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite), if it is a GRangesList,
#'  it will return it self.
#' @param part a character, one of: tx, leader, cds, trailer, intron, mrna
#' NOTE: difference between tx and mrna is that tx are all transcripts, while
#' mrna are all transcripts with a cds
#' @param names.keep a character vector of subset of names to keep. Example:
#' loadRegions(txdb, names = "ENST1000005"), will return only that transcript.
#' Remember if you set by to "gene", then this list must be with gene names.
#' @param by a character, default "tx" Either "tx" or "gene". What names to
#' output region by, the transcript name "tx" or gene names "gene".
#' NOTE: this is note the same as cdsBy(txdb, by = "gene"), cdsBy would then
#' only give 1 cds per Gene, loadRegion gives all isoforms, but with gene names.
#' @param skip.optimized logical, default FALSE. If TRUE, will not search for optimized
#' rds files to load created from ORFik::makeTxdbFromGenome(..., optimize = TRUE). The
#' optimized files are ~ 100x faster to load for human genome.
#' @return a GrangesList of region
#' @export
#' @examples
#' gtf <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' loadRegion(gtf, "cds")
#' loadRegion(gtf, "intron")
loadRegion <- function(txdb, part = "tx", names.keep = NULL, by = "tx",
                       skip.optimized = FALSE) {
  if (is.grl(txdb)) return(txdb)
  if (length(part) != 1) stop("argument: (part) must be length 1,
                              use loadRegions for multi region loading!")
  txdb <- loadTxdb(txdb)
  # Check for optimized paths
  optimized_path <- optimized_txdb_path(txdb, stop.error = FALSE)
  optimized <- !is.null(optimized_path) & !skip.optimized
  region <-
    if (part %in% c("tx", "transcript", "transcripts")) {
      optimized.rds <- paste0(optimized_path, "_", "tx", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else exonsBy(txdb, by = "tx", use.names = TRUE)
    } else if (part %in% c("leader", "leaders", "5'", "5", "5utr",
                           "fiveUTRs", "5pUTR")) {
      optimized.rds <- paste0(optimized_path, "_", "leaders", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else fiveUTRsByTranscript(txdb, use.names = TRUE)
    } else if (part %in% c("cds", "CDS", "mORF")) {
      optimized.rds <- paste0(optimized_path, "_", "cds", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else cdsBy(txdb, by = "tx", use.names = TRUE)
    } else if (part %in% c("trailer", "trailers", "3'", "3", "3utr",
                           "threeUTRs", "3pUTR")) {
      optimized.rds <- paste0(optimized_path, "_", "trailers", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else threeUTRsByTranscript(txdb, use.names = TRUE)
    } else if (part %in% c("intron", "introns")) {
      intronsByTranscript(txdb, use.names = TRUE)
    }  else if (part %in% c("mrna", "mrnas", "mRNA", "mRNAs")) {
      optimized.rds <- paste0(optimized_path, "_", "mrna", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else exonsBy(txdb, by = "tx", use.names = TRUE)[names(cdsBy(txdb, use.names = TRUE))]
    } else stop("invalid: must be tx, leader, cds, trailer, introns or mrna")

  if (by == "gene") {
    names(region) <- txNamesToGeneNames(names(region), txdb)
  }
  if (!is.null(names.keep)) { # If subset
    subset <- names(region) %in% names.keep
    if (length(subset) == 0)
      stop(paste("Found no kept transcripts, for region:", part))
    region <- region[subset]
  }
  if (all(seqlevels(region) %in% seqlevels(txdb))) { # Avoid warnings
    seqlevels(region) <- seqlevels(txdb)
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
                        by = "tx", skip.optimized = FALSE,
                        envir = .GlobalEnv) {
  txdb <- loadTxdb(txdb)
  for (i in parts) {
    assign(x = paste0(i, extension),
           value = loadRegion(txdb, i, names.keep, by, skip.optimized),
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
#' @export
#' @examples
#' gtf <- "path/to.gtf"
#' #loadTranscriptType(gtf, part = "rRNA")
#' #loadTranscriptType(gtf, part = "miRNA")
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
#' @inheritParams getGtfPathFromTxdb
#' @param txdb a TxDb, path to txdb / gff or ORFik experiment object
#' @return data.frame, the gtf/gff object imported with rtracklayer::import.
#' Or NULL, if stop.on.error is FALSE, and no GTF file found.
importGtfFromTxdb <- function(txdb, stop.on.error = TRUE) {
  if (is(txdb, "experiment")) txdb <- txdb@txdb
  if(is(txdb, "character")) {
    if (!(file_ext(txdb) %in% c("gtf", "gff", "gff3", "gff2"))) {
      # It means it shoud be a TxDb path
      txdb <- loadTxdb(txdb)
    }
  }
  if (is(txdb, "TxDb")) {
    txdb <- getGtfPathFromTxdb(txdb, stop.error = stop.error)
    if (!stop.error) return(txdb)
  }
  if (!file.exists(txdb)) {
    message <- paste("Could not open gtf, did you rename folder of gtf?")
  }
  return(import(txdb))
}

#' Get path of GTF that created txdb
#'
#' Will crash and report proper error if no gtf is found
#' @param txdb a loaded TxDb object
#' @param stop.error logical TRUE, stop if Txdb does not have a gtf.
#' If FALSE, return NULL.
#' @return a character file path, returns NULL if not valid
#' and stop.error is FALSE.
#' @keywords internal
getGtfPathFromTxdb <- function(txdb, stop.error = TRUE) {
  genome <- metadata(txdb)[metadata(txdb)[,1] == "Data source", 2]
  valid <- TRUE
  if (length(genome) == 0) valid <- FALSE
  if (valid) {
    if (!(file_ext(genome) %in%
          c("gtf", "gff", "gff3", "gff2"))) {
      valid <- FALSE
    }
  }

  if (!valid) {
    if (stop.error) {
      message("This is error txdb ->")
      message("It should be found by: metadata(txdb)[metadata(txdb)[,1] == 'Data source',]")
      print(txdb)
      stop("Your Txdb does not point to a valid gtf/gff (no valid Data source defined)")
    } else return(NULL)
  }
  return(genome)
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
#' since 0 means there must exist a column called utr3_len etc.
#' Genes with gene_id = NA will be be removed.
#' @inheritParams loadRegion
#' @param minFiveUTR (integer) minimum bp for 5' UTR during filtering for the
#' transcripts. Set to NULL if no 5' UTRs exists for annotation.
#' @param minCDS (integer) minimum bp for CDS during filtering for the
#' transcripts
#' @param minThreeUTR (integer) minimum bp for 3' UTR during filtering for the
#' transcripts. Set to NULL if no 3' UTRs exists for annotation.
#' @param longestPerGene logical (TRUE), return only longest valid transcript
#' per gene. NOTE: This is by priority longest cds isoform, if equal then pick
#' longest total transcript. So if transcript is shorter but cds is longer,
#'  it will still be the one returned.
#' @param stopOnEmpty logical TRUE, stop if no valid transcripts are found ?
#' @param create.fst.version logical, FALSE. If TRUE, creates a .fst version
#' of the transcript length table (if it not already exists),
#' reducing load time from ~ 15 seconds to
#' ~ 0.01 second next time you run filterTranscripts with this txdb object.
#' The file is stored in the
#' same folder as the genome this txdb is created from, with the name:\cr
#' \code{paste0(ORFik:::remove.file_ext(metadata(txdb)[3,2]), "_",
#'        gsub(" \\(.*| |:", "", metadata(txdb)[metadata(txdb)[,1] ==
#'         "Creation time",2]), "_txLengths.fst")}\cr
#' Some error checks are done to see this is a valid location, if the txdb
#' data source is a repository like UCSC and not a local folder, it will not
#' be made.
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
                              by = "tx", create.fst.version = FALSE) {
  if (!(by %in% c("tx", "gene"))) stop("by must be either tx or gene!")
  txdb <- loadTxdb(txdb)
  five <- !is.null(minFiveUTR)
  three <- !is.null(minThreeUTR)

  tx <- optimizedTranscriptLengths(txdb, five, three, create.fst.version)
  five <- rep(five, nrow(tx))
  three <- rep(three, nrow(tx))

  tx <- tx[ifelse(five, utr5_len >= minFiveUTR, TRUE) & cds_len >= minCDS &
             ifelse(three, utr3_len >= minThreeUTR, TRUE), ]

  gene_id <- cds_len <- NULL
  tx <- data.frame(tx)
  tx <- tx[order(tx$gene_id, -rank(tx$cds_len), -rank(tx$tx_len)), ]
  # can't be used due to crashes of R, no errors reported...
  # data.table::setorder(tx, gene_id, -cds_len, -tx_len)
  if (longestPerGene) {
    tx <- tx[!duplicated(tx$gene_id), ]
  }
  tx <- tx[!is.na(tx$gene_id), ]

  if (stopOnEmpty & length(tx$tx_name) == 0)
    stop("No transcript has leaders and trailers of specified minFiveUTR",
         " minCDS, minThreeUTR, create pseudo-UTRs, see annotation vignette!")
  if (by == "gene")
    return(tx$gene_id)
  return(tx$tx_name)
}

#' Get path for optimization files for txdb
#' @param txdb a loaded TxDb object
#' @param create.dir logical FALSE, if TRUE create the
#' optimization directory, this should only be called first time used.
#' @param stop.error logical TRUE
#' @return a character file path, returns NULL if not valid
#' and stop.error is FALSE.
#' @keywords internal
optimized_txdb_path <- function(txdb, create.dir = FALSE, stop.error = TRUE) {
  genome <- getGtfPathFromTxdb(txdb, stop.error = stop.error)
  if (is.null(genome)) {
    return(NULL)
  }
  base_dir <- file.path(dirname(genome), "ORFik_optimized")
  base_path <- file.path(base_dir, basename(remove.file_ext(genome)))
  creation.time <- gsub(" \\(.*| |:", "",
                        metadata(txdb)[metadata(txdb)[,1] == "Creation time",2])
  if (create.dir) {
    dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(paste0(base_path, "_", creation.time))
}

#' Load length of all transcripts
#'
#' A speedup wrapper around \code{\link{transcriptLengths}},
#' default load time of lengths is ~ 15 seconds, if ORFik fst
#' optimized lengths object has been made, load that file instead:
#' load time reduced to ~ 0.1 second.
#' @inheritParams filterTranscripts
#' @param with.utr5_len logical TRUE, include length of 5' UTRs,
#'  ignored if .fst exists
#' @param with.utr3_len logical TRUE, include length of 3' UTRs,
#'  ignored if .fst exists
#' @return a data.table of loaded lengths 8 columns,
#' 1 row per transcript isoform.
#' @importFrom data.table setDT
#' @importFrom fst read_fst
#' @importFrom fst write_fst
#' @keywords internal
optimizedTranscriptLengths <- function(txdb, with.utr5_len = TRUE,
                                       with.utr3_len = TRUE,
                                       create.fst.version = FALSE) {
  txdb <- loadTxdb(txdb)

  optimized_path <- optimized_txdb_path(txdb, stop.error = FALSE)
  found_gtf <- !is.null(optimized_path)
  possible_fst <- paste0(optimized_path, "_txLengths.fst")
  if (file.exists(possible_fst) & found_gtf) { # If fst exists
    return(setDT(fst::read_fst(possible_fst)))
  } else if (create.fst.version & found_gtf) { # If make fst
    tx <- data.table::setDT(
      GenomicFeatures::transcriptLengths(
        txdb, with.cds_len = TRUE,
        with.utr5_len = with.utr5_len,
        with.utr3_len = with.utr3_len))
    # Validation save location
    if (dir.exists(dirname(possible_fst))) {
        message("Creating fst speedup file for transcript lengths, at location:")
        message(possible_fst)
        fst::write_fst(tx, possible_fst)
    } else {
        stop("Optimization directory does not exist: ", dirname(possible_fst))
      }

    return(tx)
  }
  # Else load normally
  return(data.table::setDT(
    GenomicFeatures::transcriptLengths(
      txdb, with.cds_len = TRUE,
      with.utr5_len = with.utr5_len,
      with.utr3_len = with.utr3_len)))
}
