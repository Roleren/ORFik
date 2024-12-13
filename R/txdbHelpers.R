
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
#' within the output folder (defined by txdb_file_out_path),
#' that includes optimized objects
#' to speed up loading of annotation regions from up to 15 seconds
#' on human genome down to 0.1 second. ORFik will then load these optimized
#' objects instead. Currently optimizes filterTranscript() function and
#' loadRegion() function for 5' UTRs, 3' UTRs, CDS,
#'  mRNA (all transcript with CDS) and tx (all transcripts).
#' @param gene_symbols logical default FALSE. If TRUE, will download
#' and store all gene symbols for all transcripts (coding and noncoding)-
#' In a file called: "gene_symbol_tx_table.fst" in same folder as txdb.
#' hgcn for human, mouse symbols for mouse and rat, more to be added.
#' @param uniprot_id logical default FALSE.  If TRUE, will download
#' and store all uniprot id for all transcripts (coding and noncoding)-
#' In a file called: "gene_symbol_tx_table.fst" in same folder as txdb.
#' @param return logical, default FALSE. If TRUE, return TXDB object,
#' else invisible(NULL).
#' @param txdb_file_out_path character path, default paste0(gtf, ".db").
#' Set to NULL to not write file to disc.
#' @param symbols_file_out_path character path, default
#' file.path(dirname(gtf), "gene_symbol_tx_table.fst").
#' Must be defined as character if "gene_symbols" is TRUE. Ignored if
#' "gene_symbols" is FALSE.
#' @inheritParams add_pseudo_5utrs_txdb_if_needed
#' @return logical, default is.null(txdb_file_out_path),
#' Txdb saved to disc named default paste0(gtf, ".db").
#' Set 'return' argument to TRUE, to also get txdb back as an object.
#' @export
#' @examples
#' gtf <- "/path/to/local/annotation.gtf"
#' genome <- "/path/to/local/genome.fasta"
#' #makeTxdbFromGenome(gtf, genome, organism = "Saccharomyces cerevisiae")
#' # Runnable full example
#' df <- ORFik.template.experiment()
#' gtf <- sub("\\.db$", "", df@txdb)
#' genome <- df@fafile
#' txdb <- makeTxdbFromGenome(gtf, genome, organism = "Saccharomyces cerevisiae",
#'   txdb_file_out_path = NULL)
#' ## Add pseudo UTRs if needed (< 30% of cds have a defined 5'UTR)
makeTxdbFromGenome <- function(gtf, genome = NULL, organism,
                               optimize = FALSE, gene_symbols = FALSE,
                               uniprot_id = FALSE,
                               pseudo_5UTRS_if_needed = NULL,
                               minimum_5UTR_percentage = 30,
                               return = is.null(txdb_file_out_path),
                               txdb_file_out_path = paste0(gtf, ".db"),
                               symbols_file_out_path = file.path(dirname(gtf), "gene_symbol_tx_table.fst")) {

  txdb <- makeTxdbTemplate(gtf, genome, organism)
  txdb <- add_pseudo_5utrs_txdb_if_needed(txdb, pseudo_5UTRS_if_needed,
                                          minimum_5UTR_percentage)
  if (!is.null(txdb_file_out_path)) {
    AnnotationDbi::saveDb(txdb, txdb_file_out_path)
    message("--------------------------")
    message("Txdb stored at: ", txdb_file_out_path)
  }

  if (optimize) {
    message("--------------------------")
    optimized_path <- optimized_txdb_path(txdb, create.dir = TRUE,
                                          stop.error = TRUE,
                                          txdb_file_out_path)
    message("Optimizing annotation, saving to: ", dirname(optimized_path))
    # Save all transcript, cds and UTR lengths as .fst
    optimizedTranscriptLengths(txdb, create.fst.version = TRUE,
                               optimized_path = optimized_path)
    # Save RDS version of all transcript regions
    optimizeTranscriptRegions(txdb, optimized_path)
  }
  if (gene_symbols) {
    stopifnot(is(symbols_file_out_path, "character") &
                length(symbols_file_out_path) == 1)
    symbols <- geneToSymbol(txdb, include_tx_ids = TRUE,
                            uniprot_id = uniprot_id)

    if (nrow(symbols) > 0) fst::write_fst(symbols, symbols_file_out_path)
  }
  if (return) return(txdb)
  return(invisible(NULL))
}

makeTxdbTemplate <- function(gtf, genome = NULL, organism) {

  message("Making optimized txdb using ", toupper(tools::file_ext(gtf)))
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

    txdb <- loadTxdb(gtf, organism = organismCapital, chrominfo = fa.seqinfo)

    supported_species <- try(seqlevelsStyle(txdb)[1], silent = TRUE)
    if (!is(supported_species, "try-error")) {
      if (seqlevelsStyle(txdb)[1] != seqlevelsStyle(fa)[1]) {
        seqlevelsStyle(txdb) <- seqlevelsStyle(fa)[1]
      }
    } else message("Species not supported by seqlevelsStyle,",
                   "skipping safety test!")


  } else {
    txdb <- loadTxdb(gtf, organism = organismCapital)
  }
  return(txdb)
}

#' add_pseudo_5utrs_txdb_if_needed
#' @param txdb a TxDb object
#' @param pseudo_5UTRS_if_needed integer, default NULL. If defined > 0,
#' will add pseudo 5' UTRs of maximum this length if 'minimum_5UTR_percentage" (default 30%) of
#' mRNAs (coding transcripts) do not have a leader. (NULL and 0 are both the ignore command)
#' @param minimum_5UTR_percentage numeric, default 30. What minimum percentage
#' of mRNAs most have a 5' UTRs (leaders), to not do the pseudo_UTR addition.
#' If percentage is higher, addition is ignored, set to 101 to always do it.
#' @return txdb (new txdb if it was done, old if not)
add_pseudo_5utrs_txdb_if_needed <- function(txdb, pseudo_5UTRS_if_needed = NULL, minimum_5UTR_percentage = 30) {
  pseudo5 <- pseudo_5UTRS_if_needed
  if (!is.null(pseudo5)) {
    if (pseudo5 > 0 & is.numeric(pseudo5)) {
      message("---- Adding Pseudo 5' UTRs")
      leaders <- loadRegion(txdb, "leaders")

      message("- Total leaders: ", length(leaders))
      if (length(leaders) > 0) {
        message("- Leader length statistics:")
        print(summary(widthPerGroup(leaders, FALSE)))
        cds <- loadRegion(txdb, "cds")
        if (length(cds) == 0)
          stop("Can not add pseudo 5' UTRs to a genome without Coding sequences!")
        percentage <- round((length(leaders) / length(cds))*100, 1)
        message("- Percentage of CDS' with leaders: ", percentage, "%")
      } else {
        message("-- This genome has no predefined leaders")
        percentage <- 0
      }
      if (percentage < minimum_5UTR_percentage) {
        message("- Adding pseudo")
        txdb <- assignTSSByCage(txdb, cage = NULL, pseudoLength = pseudo5)
      }
    } else warning("Malformed 'pseudo_5UTRS_if_needed', ignoring it")
  }
  return(txdb)
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
#' @param organism character, default NA. Scientific name of organism.
#'  Only used if input is path to gff.
#' @param chrominfo Seqinfo object, default NULL.
#'  Only used if input is path to gff.
#' @return a TxDb object
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures makeTxDbFromGFF makeTxDb
#' @export
#' @examples
#' library(GenomicFeatures)
#' # Get the gtf txdb file
#' txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                         package = "GenomicFeatures")
#' txdb <- loadTxdb(txdbFile)
#'
loadTxdb <- function(txdb, chrStyle = NULL, organism = NA,
                     chrominfo = NULL) {
  if (is(txdb, "experiment")) {
    txdb <- txdb@txdb
  }
  #TODO: Check that is is an open connection!
  if (is(txdb, "character")) {
    f <- file_ext(txdb)
    if (f %in% c("gff", "gff2", "gff3", "gtf")) {
      txdb <- makeTxDbFromGFF(txdb, organism = organism, chrominfo = chrominfo)
    } else if(f %in% c("db", "sqlite")) {
      txdb <- loadDb(txdb)
    } else {
      stop("when txdb is path, must be one of:
           gff - (.gff, .gff2, .gff3, .gtf)
           txdb - (.db, .sqlite")
    }
  } else if (is.list(txdb)) {
    must_have_cols <- c("transcripts", "splicings", "genes", "chrominfo")
    if (!all(must_have_cols %in% names(txdb))) {
      stop("When txdb is list input, must have names:",
           paste(must_have_cols, collapse = ", "))
    }
    return(do.call(makeTxDb, txdb))
  } else if(!is(txdb, "TxDb")) stop("txdb must be path, list or TxDb")
  return(matchSeqStyle(txdb, chrStyle))
}

#' Load transcript region
#'
#' Usefull to simplify loading of standard regions, like cds' and leaders.
#' Adds another safety in that seqlevels will be set
#'
#' Load as GRangesList if input is not already GRangesList.
#' @param txdb a TxDb object, ORFik experiment object or a path to one of:
#'  (.gtf ,.gff, .gff2, .gff2, .db or .sqlite),
#'  Only in the loadRegion function: if it is a GRangesList, it will return it self.
#' @param part a character, one of: tx, ncRNA, mrna, leader, cds, trailer, intron,
#' NOTE: difference between tx and mrna is that tx are all transcripts, while
#' mrna are all transcripts with a cds, respectivly ncRNA are all tx without a cds.
#' @param names.keep a character vector of subset of names to keep. Example:
#' loadRegions(txdb, names = "ENST1000005"), will return only that transcript.
#' Remember if you set by to "gene", then this list must be with gene names.
#' @param by a character, default "tx" Either "tx" or "gene". What names to
#' output region by, the transcript name "tx" or gene names "gene".
#' NOTE: this is not the same as cdsBy(txdb, by = "gene"), cdsBy would then
#' only give 1 cds per Gene, loadRegion gives all isoforms, but with gene names.
#' @param skip.optimized logical, default FALSE. If TRUE, will not search for optimized
#' rds files to load created from ORFik::makeTxdbFromGenome(..., optimize = TRUE). The
#' optimized files are ~ 100x faster to load for human genome.
#' @return a GRangesList of region
#' @importFrom GenomicFeatures exonsBy cdsBy intronsByTranscript
#' fiveUTRsByTranscript threeUTRsByTranscript
#' @export
#' @examples
#' # GTF file is slow, but possible to use
#' gtf <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                         package = "GenomicFeatures")
#' txdb <- loadTxdb(gtf)
#' loadRegion(txdb, "cds")
#' loadRegion(txdb, "intron")
#' # Use txdb from experiment
#' df <- ORFik.template.experiment()
#' txdb <- loadTxdb(df)
#' loadRegion(txdb, "leaders")
#' # Use ORFik experiment directly
#' loadRegion(df, "mrna")
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
    }  else if (part %in% c("ncRNA", "ncrna")) {
      optimized.rds <- paste0(optimized_path, "_", "ncRNA", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else {
        tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
        tx[!(names(tx) %in% names(cdsBy(txdb, use.names = TRUE)))]
      }
    } else if (part %in% c("mrna", "mrnas", "mRNA", "mRNAs")) {
      optimized.rds <- paste0(optimized_path, "_", "mrna", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else exonsBy(txdb, by = "tx", use.names = TRUE)[names(cdsBy(txdb, use.names = TRUE))]
    } else if (part %in% c("uorf", "uorfs")) {
      optimized.rds <- paste0(optimized_path, "_", "uorfs", ".rds")
      if (optimized & file.exists(optimized.rds)) {
        readRDS(optimized.rds)
      } else stop("uORFs must always exist before calling this function, see ?findUORFs_exp")
    } else stop("invalid: must be tx, leader, cds, trailer, introns, ncRNA, uorf or mrna")

  if (by == "gene") {
    names(region) <- txNamesToGeneNames(names(region), txdb)
  }
  if (!is.null(names.keep)) { # If subset
    subset <-
    if (part %in% c("uorf", "uorfs")) {
      txNames(region) %in% names.keep
    } else names(region) %in% names.keep
    if (length(subset) == 0)
      stop("Found no kept transcript with given subset of part:", part)
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
#' gtf <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                         package = "GenomicFeatures")
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
#' @importFrom GenomicFeatures transcripts
#' @examples
#' df <- ORFik.template.experiment()
#' txdb <- loadTxdb(df)
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

#' Get gene symbols from Ensembl gene ids
#'
#' If your organism is not in this list of supported
#' organisms, manually assign the input arguments.
#' There are 2 main fetch modes:\cr
#' By gene ids (Single accession per gene)\cr
#' By tx ids (Multiple accessions per gene)\cr
#' Run the mode you need depending on your required attributes.
#' \cr\cr
#' Will check for already existing table of all genes, and use that instead
#' of re-downloading every time (If you input valid experiment or txdb
#' and have run \code{\link{makeTxdbFromGenome}}
#' with symbols = TRUE, you have a file called gene_symbol_tx_table.fst) will
#' load instantly. If df = NULL, it can still search cache to load a bit slower.
#' @param df an ORFik \code{\link{experiment}} or TxDb object with defined organism slot.
#' If set will look for file at path of txdb / experiment reference path named:
#' 'gene_symbol_tx_table.fst' relative to the txdb/genome directory.
#' Can be set to NULL if gene_ids and organism is defined manually.
#' @param organism_name default, \code{organism(df)}.
#' Scientific name of organism, like ("Homo sapiens"),
#' remember capital letter for first name only!
#' @param gene_ids default, \code{filterTranscripts(df, by = "gene", 0, 0, 0)}.
#'  Ensembl gene IDs to search for (default all transcripts coding and noncoding)
#'  To only get coding do: filterTranscripts(df, by = "gene", 0, 1, 0)
#' @param org.dataset default, \code{paste0(tolower(substr(organism_name, 1, 1)), gsub(".* ", replacement = "", organism_name), "_gene_ensembl")}
#' the ensembl dataset to use. For Homo sapiens, this converts to default as: hsapiens_gene_ensembl
#' @param ensembl default, \code{useEnsembl("ensembl",dataset=org.dataset)} .The mart connection.
#' @param attribute default, "external_gene_name", the biomaRt column / columns
#' default(primary gene symbol names). These are always from specific database, like
#' hgnc symbol for human, and mgi symbol for mouse and rat, sgd for yeast etc.
#' @param include_tx_ids logical, default FALSE, also match tx ids, which then returns as the 3rd column.
#' Only allowed when 'df' is defined. If
#' @param uniprot_id logical, default FALSE. Include uniprotsptrembl and/or uniprotswissprot.
#' If include_tx_ids you will get per isoform if available, else you get canonical uniprot id
#' per gene. If both uniprotsptrembl and uniprotswissprot exists, it will make a merged
#' uniprot id column with rule: if id exists in uniprotswissprot, keep.
#' If not, use uniprotsptrembl column id.
#' @param force logical FALSE, if TRUE will not look for existing file made through \code{\link{makeTxdbFromGenome}}
#' corresponding to this txdb / ORFik experiment stored with name "gene_symbol_tx_table.fst".
#' @param verbose logical TRUE, if FALSE, do not output messages.
#' @return data.table with 2, 3 or 4 columns: gene_id, gene_symbol, tx_id and uniprot_id named after
#' attribute, sorted in order of gene_ids input.
#' (example: returns 3 columns if include_tx_ids is TRUE),
#' and more if additional columns are specified in 'attribute' argument.
#' @export
#' @importFrom biomaRt useEnsembl getBM searchAttributes
#' @examples
#' ## Without ORFik experiment input
#' gene_id_ATF4 <- "ENSG00000128272"
#' #geneToSymbol(NULL, organism_name = "Homo sapiens", gene_ids = gene_id_ATF4)
#' # With uniprot canonical isoform id:
#' #geneToSymbol(NULL, organism_name = "Homo sapiens", gene_ids = gene_id_ATF4, uniprot_id = TRUE)
#'
#' ## All genes from Organism using ORFik experiment
#' # df <- read.experiment("some_experiment)
#' # geneToSymbol(df)
#'
#' ## Non vertebrate species (the ones not in ensembl, but in ensemblGenomes mart)
#' #txdb_ylipolytica <- loadTxdb("txdb_path")
#' #dt2 <- geneToSymbol(txdb_ylipolytica, include_tx_ids = TRUE,
#' #   ensembl = useEnsemblGenomes(biomart = "fungi_mart", dataset = "ylipolytica_eg_gene"))
#'
geneToSymbol <- function(df, organism_name = organism(df),
                         gene_ids = filterTranscripts(df, by = "gene", 0, 0, 0),
                         org.dataset = paste0(tolower(substr(organism_name, 1, 1)), gsub(".* ", replacement = "", organism_name), "_gene_ensembl"),
                         ensembl = biomaRt::useEnsembl("ensembl",dataset=org.dataset),
                         attribute = "external_gene_name",
                         include_tx_ids = FALSE, uniprot_id = FALSE,
                         force = FALSE, verbose = TRUE) {
  include_uniprot_id <- uniprot_id
  need_to_fetch_tx_ids <- include_tx_ids &
    (!all(attribute == "external_gene_name") | include_uniprot_id)

  need_annotation_file <- need_to_fetch_tx_ids |
    (include_tx_ids & !need_to_fetch_tx_ids)
  if (need_annotation_file) {
    if (is.null(df)) stop("'df' can not be NULL when 'include_tx_ids' is TRUE")
  }
  if (anyNA(attribute)) {
    warning("You inserted a species not supported at the moment,
            skipping download of symbols")
    return(data.table())
  }
  if (!is.null(df)) {
    file <- ifelse(is(df, "experiment"), df@txdb, getGtfPathFromTxdb(df))
    file <- file.path(dirname(file), "gene_symbol_tx_table.fst")
    if (file.exists(file) & !force) {
      if (verbose) message("Loading pre-existing symbols from file")
      return(as.data.table(fst::read_fst(file)))
    }
  }
  if (verbose) {
    message("- Symbols extracted from:")
    message("Organism: ", organism_name)
    message("Dataset: ", attr(ensembl, "dataset"))
    message("Attribute: ", attribute)
  }

  # Find correct uniprot column

  if (include_uniprot_id) {
    if (verbose) message("Uniprot id: ", appendLF = FALSE)
    all_attributes <- biomaRt::searchAttributes(ensembl, "uniprot")
    uni_match <- all_attributes$name[grep("uniprot", all_attributes$name)]
    uni_id <- c()
    if ("uniprotswissprot" %in% uni_match) {
      uni_id <- "uniprotswissprot"
    }
    if ("uniprotsptrembl" %in% uni_match) {
      uni_id <- c(uni_id, "uniprotsptrembl")
    }
    if (is.null(uni_id)) stop("This organism does not have uniprot ids, debug this function
                to see all attributes available!")
    if (verbose) message(paste(uni_id, collapse = " "))
  }

  # Only fetch tx ids if some biomaRt column depends on them, else just load

  value_ids <- gene_ids
  filter_column_id <- 'ensembl_gene_id'
  other_attributes <- "ensembl_gene_id"
  all_attr_to_get <- c(other_attributes, attribute)
  if (need_to_fetch_tx_ids) {
    filter_column_id <- 'ensembl_transcript_id'
    tx_ids <- optimizedTranscriptLengths(df, TRUE, TRUE)[,2:3]
    gene_subset <- tx_ids[gene_id %in% gene_ids,]
    value_ids <- gene_subset$tx_name
    gene_ids <- gene_subset$gene_id
    all_attr_to_get <- c(all_attr_to_get, "ensembl_transcript_id")
  }
  if (include_uniprot_id) all_attr_to_get <- c(all_attr_to_get, uni_id)
  # Call to biomaRt
  if (verbose) message("- Calling biomaRt:")
  gene_id <- as.data.table(biomaRt::getBM(attributes=all_attr_to_get,
                                          filters = filter_column_id,
                                          values = value_ids,
                                          mart = ensembl))
  # Filtering
  gene_id <- gene_id[chmatch(value_ids,
                             gene_id[, filter_column_id, with = FALSE][[1]]),]
  which_na <- is.na(gene_id[, filter_column_id, with = FALSE][[1]])
  gene_id[which_na, ensembl_gene_id := gene_ids[which_na]]
  if (need_to_fetch_tx_ids) {
    unique_symbol <- gene_id[!which_na,  unique(external_gene_name),
                             by = ensembl_gene_id]

    gene_id[which_na, external_gene_name := unique_symbol[chmatch(gene_ids[which_na], ensembl_gene_id)]$V1]
    gene_id[which_na, ensembl_transcript_id := value_ids[which_na]]
    colnames(gene_id)[colnames(gene_id) == "ensembl_transcript_id"] <-
      "ensembl_tx_name"
  }


  if ("external_gene_name" %in% all_attr_to_get) {
    if (anyNA(gene_id[,"external_gene_name"])) {
      if (verbose) message("Symbol detection rate: ",
       (1 - sum(is.na(gene_id[,"external_gene_name"]))/nrow(gene_id))*100, " %")
    } else if (verbose) message("Found symbol for all genes")
  }


  if (include_uniprot_id) {
    colnames(gene_id)[ncol(gene_id)] <- "uniprot_id"
    if (length(uni_id) == 2) {
      gene_id[uniprotswissprot != "", uniprot_id := uniprotswissprot]
      gene_id[, uniprotswissprot := NULL]
    }

    if (anyNA(gene_id[,"uniprot_id"])) {
      if (verbose) message("Uniprot id detection rate: ",
             (1 - sum(is.na(gene_id[,"uniprot_id"]))/nrow(gene_id))*100, " %")
    } else if (verbose) message("Found Uniprot id for all genes")
  }

  if (include_tx_ids & !need_to_fetch_tx_ids) {
    if (!exists("tx_ids")) {
      tx_ids <- optimizedTranscriptLengths(df, TRUE, TRUE)[,2:3]
    }
    colnames(tx_ids) <- paste0("ensembl_", colnames(tx_ids))
    gene_id <- merge.data.table(gene_id, tx_ids, by = filter_column_id, all_x = TRUE)
  }
  gene_id[]
  return(gene_id)
}

#' Import the GTF / GFF that made the txdb
#' @inheritParams getGtfPathFromTxdb
#' @param txdb a TxDb, path to txdb / gff or ORFik experiment object
#' @return data.frame, the gtf/gff object imported with rtracklayer::import.
#' Or NULL, if stop.error is FALSE, and no GTF file found.
importGtfFromTxdb <- function(txdb, stop.error = TRUE) {
  if (is(txdb, "experiment")) txdb <- txdb@txdb
  if(is(txdb, "character")) {
    if (!(file_ext(txdb) %in% c("gtf", "gff", "gff3", "gff2"))) {
      # It means it shoud be a TxDb path
      txdb <- loadTxdb(txdb)
    }
  }
  if (is(txdb, "TxDb")) {
    txdb <- getGtfPathFromTxdb(txdb, stop.error = stop.error)
    if (is.null(txdb)) return(txdb)
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
  txdb <- loadTxdb(txdb)
  genome <- metadata(txdb)[metadata(txdb)[,1] == "Data source", 2]
  valid <- TRUE
  if (length(genome) == 0) valid <- FALSE
  if (valid) {
    if (!file.exists(genome)) {
      genome <- tools::file_path_sans_ext(txdb$conn@dbname)
      if (!file.exists(genome)) valid <- FALSE
    }

    if (!(file_ext(genome) %in%
          c("gtf", "gff", "gff3", "gff2"))) {
      valid <- FALSE
    }
  }

  if (!valid) {
    if (stop.error) {
      message("This is error txdb ->")
      message("It should be found by: metadata(txdb)[metadata(txdb)[,1] == 'Data source',]")
      message("Or here: tools::file_path_sans_ext(txdb$conn@dbname)")
      print(txdb)
      stop("Your Txdb does not point to a valid gtf/gff")
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
#' df <- ORFik.template.experiment.zf()
#' txdb <- loadTxdb(df)
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
  five <- !is.null(minFiveUTR)
  three <- !is.null(minThreeUTR)

  tx <- optimizedTranscriptLengths(txdb, five, three, create.fst.version)
  five <- rep(five, nrow(tx))
  three <- rep(three, nrow(tx))

  tx <- tx[ifelse(five, utr5_len >= minFiveUTR, TRUE) & cds_len >= minCDS &
             ifelse(three, utr3_len >= minThreeUTR, TRUE), ]

  gene_id <- cds_len <- NULL
  data.table::setorderv(tx, c("gene_id", "cds_len", "tx_len"), c(1,-1,-1))
  # can't be used due to crashes of R, no errors reported...
  # data.table::setorder(tx, gene_id, -cds_len, -tx_len)
  if (longestPerGene) {
    tx <- tx[!duplicated(gene_id), ]
  }
  tx <- tx[!is.na(gene_id), ]

  if (stopOnEmpty & length(tx$tx_name) == 0)
    stop("No transcript has leaders and trailers of specified minFiveUTR",
         " minCDS, minThreeUTR, create pseudo-UTRs, see annotation vignette!")
  if (by == "gene")
    return(tx$gene_id)
  return(tx$tx_name)
}

#' Get path for optimization files for txdb
#' @inheritParams loadRegion
#' @param create.dir logical FALSE, if TRUE create the
#' optimization directory, this should only be called first time used.
#' @param stop.error logical TRUE
#' @param gtf_path path to gtf where output should be stored in subfolder
#'  "./ORFik_optimized"
#' @return a character file path, returns NULL if not valid
#' and stop.error is FALSE.
#' @keywords internal
optimized_txdb_path <- function(txdb, create.dir = FALSE, stop.error = TRUE,
                                gtf_path = getGtfPathFromTxdb(txdb, stop.error = stop.error)) {
  if (is.null(gtf_path)) {
    return(NULL)
  }
  base_dir <- file.path(dirname(gtf_path), "ORFik_optimized")
  base_path <- file.path(base_dir, basename(remove.file_ext(gtf_path)))
  creation.time <- gsub(" \\(.*| |:", "",
                        metadata(txdb)[metadata(txdb)[,1] == "Creation time",2])
  if (create.dir) {
    dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  }
  return(paste0(base_path, "_", creation.time))
}

#' Load length and names of all transcripts
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
#' @param optimized_path character, path to optimized txdb objects,
#' default: optimized_txdb_path(txdb, stop.error = FALSE). If no existing file,
#' will be slower and load lengths through \code{\link{transcriptLengths}}.
#' @return a data.table of loaded lengths 8 columns,
#' 1 row per transcript isoform.
#' @importFrom data.table setDT
#' @importFrom fst read_fst
#' @importFrom fst write_fst
#' @export
#' @examples
#' dt <- optimizedTranscriptLengths(ORFik.template.experiment())
#' dt
#' dt[cds_len > 0,] # All mRNA
optimizedTranscriptLengths <- function(txdb, with.utr5_len = TRUE,
                                       with.utr3_len = TRUE,
                                       create.fst.version = FALSE,
                                       optimized_path = optimized_txdb_path(txdb, stop.error = FALSE)) {
  txdb <- loadTxdb(txdb)

  found_gtf <- !is.null(optimized_path)
  possible_fst <- paste0(optimized_path, "_txLengths.fst")
  if (found_gtf && file.exists(possible_fst)) { # If fst exists
    return(setDT(fst::read_fst(possible_fst)))
  } else if (found_gtf && create.fst.version) { # If make fst
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

#' Make optimized GRangesList objects saved to disc
#'
#' Much faster to load
#' @inheritParams filterTranscripts
#' @param base_path Directy and file prefix for files, will append "_region.rds", where region is
#' specific region.
#' @param regions character, default: c("tx", "mrna", "leaders", "cds", "trailers", "ncRNA").
#' Valid options specified by loadRegion.
#' @return invisible(NULL)
optimizeTranscriptRegions <- function(txdb, base_path = optimized_txdb_path(txdb, create.dir = TRUE),
                                       regions = c("tx", "mrna", "leaders", "cds", "trailers", "ncRNA")) {
  message("Creating rds speedup files for transcript regions")
  for (region in regions) {
    saveRDS(loadRegion(txdb, region, by = "tx"),
            file = paste0(base_path, "_", region, ".rds"))
  }
  return(invisible(NULL))
}
