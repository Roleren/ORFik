
#' Shift ribo-seq reads using cigar string
#'
#' @param cigar the cigar of the reads
#' @param shift the shift as numeric
#' @param is_plus_strand logical
#' @return the shifted read
parse_cigar <- function(cigar, shift, is_plus_strand) {
  c_signs <- unlist(explodeCigarOps(cigar))
  c_counts <- unlist(explodeCigarOpLengths(cigar))

  i = ifelse(is_plus_strand, 0, length(c_signs) + 1)
  increment = ifelse(is_plus_strand, 1, -1)
  limit = 0
  plusShift = 0
  while (shift >= limit) {
    i = i + increment
    if (c_signs[i] == "M") {
      limit = limit + c_counts[i]
    } else if (c_signs[i] == "N" || c_signs[i] == "D") {
      plusShift = plusShift + c_counts[i]
    } else if (c_signs[i] == "I") {
      plusShift = plusShift - c_counts[i]
    } else {
      warning(paste0("Not supported sign:", c_signs[i]))
    }
  }
  shift = shift + plusShift
  return(shift)
}

#' Get only unique mapping reads
#'
#' Helper for detectRibosomeShifts
#' @param file a bam or bed file path or GAlignment object
#' @param filterNH1 logical(T), only keep unique reads, require a
#'  meta column called NH
#' @return a GAlignments object
#' @importFrom tools file_ext
#' @importFrom Rsamtools ScanBamParam
#' @export
#' @examples
#' bam <- system.file("extdata", "example.bam",
#'        package = "ORFik") ## location of the bam file
#' readFootprints(bam)
#'
readFootprints <- function(file, filterNH1=TRUE) {

  if (class(file) == "GAlignments" || class(file) == "GRanges") {
    alignments <- file
  } else if(is.character(file)) {
    message("Loading reads into memory...")
    if (file_ext(file) == "bam") {
      param <- ScanBamParam(tag="NH")
      alignments <- readGAlignments(file, param=param)

    } else if (file_ext(file) == "bed") {
      alignments <- fread.bed(file)
    }
    message("Reads loaded.")
  } else {
    stop("could not get alignments from file, check input.")
  }

  if (filterNH1) { ## filter for unique reads
    if(!is.null(mcols(alignments)$NH)) {
      alignments <- alignments[mcols(alignments)$NH == 1]
    }
  }

  return(alignments)
}

#' Get the longest transcripts
#'
#' These settings must be true:
#' utr5 & utr3 > 30 width
#' cds > 150 width
#' If these are accepted, it is returned
#' @param txdb a TxDb object from gtf
#' @return a character vector of valid tramscript names
#' @importFrom data.table setDT setorder
get_longest_tx <- function(txdb) {
  if(class(txdb) != "TxDb") stop("txdb must be a TxDb object")

  txLengths <- setDT(transcriptLengths(txdb, with.cds_len=TRUE,
                                       with.utr5_len=TRUE, with.utr3_len=TRUE))
  ## avoid check warning
  gene_id <- NULL
  cds_len <- NULL
  setorder(txLengths, gene_id, -cds_len)

  longest_tx <- txLengths[!duplicated(txLengths$gene_id),]
  rownames(longest_tx) <- longest_tx$tx_name
  # get CDS longer than 150, UTRs longer than 30
  longest_tx <- longest_tx[longest_tx$utr5_len > 30,]
  longest_tx <- longest_tx[longest_tx$cds_len > 150,]
  longest_tx <- longest_tx[longest_tx$utr3_len > 30,]
  return(longest_tx$tx_name)
}

#' Get cds' for shoelaces periodicity
#'
#' Gives you the first 150 bases per cds
#' @param txdb a TxDb object from gtf
#' @param txNames a character vector of the longest
#'  valid transcript names
#' @return a GRangesList of tiled cds' to 150 length
get_cds150 <- function(txdb, txNames) {
  cds <- cdsBy(txdb, by="tx", use.names=TRUE)[txNames]
  tiledCDS <- tile1(groupGRangesBy(unlist(cds, use.names = TRUE)))

  return(phead(tiledCDS, 150L))
}

#' Get start and stops of cds windows
#'
#' For each cds in longest_tx, get a window around start and stop
#' of each cds. +29 and -30 from cds start
#' @param txdb a txdb object from gtf
#' @param txNames a character vector of the longest
#'  valid transcript names
#' @param start logical (T), include starts
#' @param stop logical (F), include stops
#' @return a list of starts and stops with cds windows
get_start_stop <- function(txdb, txNames, start=TRUE, stop=TRUE) {
  cds <- cdsBy(txdb, by="tx", use.names=TRUE)[txNames]
  cdsTiled <- tile1(groupGRangesBy(unlist(cds, use.names = TRUE)))

  if (start) {
    fiveUTRs <- fiveUTRsByTranscript(txdb, use.names=TRUE)[txNames]
    ##get tiled cds
    cdsHead <- phead(cdsTiled, 30L)

    ##get tiled fiveUTRs
    fiveTiled <- tile1(groupGRangesBy(unlist(fiveUTRs, use.names = TRUE)))
    fiveTail <- ptail(fiveTiled, 30L)

    merged <- unlist(c(fiveTail, cdsHead), use.names = FALSE)
    start <- split(merged, names(merged))
  }

  if (stop) {
    threeUTRs <- threeUTRsByTranscript(txdb, use.names=TRUE)[txNames]
    ##get tiled cds
    cdsTail <- ptail(cdsTiled, 30L)

    ##get tiled threeUTRs
    threeTiled <- tile1(groupGRangesBy(unlist(threeUTRs, use.names = TRUE)))
    threeHead <- phead(threeTiled, 30L)

    merged <- unlist(c(cdsTail, threeHead), use.names = FALSE)
    stop <- split(merged, names(merged))
  }
  ss <- list(start, stop)
  names(ss) <- c("start", "stop")

  return(ss)
}

#' Get the start of each read
#'
#' @param alignments the ribo-seq footprints
#' @param read_length the read length to search for in alignments
#' @return a GRanges object with a hit column
get_5ends <- function(alignments, read_length) {

  if(class(alignments) == "GAlignments"){
    aln <- alignments[qwidth(alignments) == read_length]
  } else {
    aln <- alignments[width(alignments) == read_length]
  }
  plus <- aln[strand(aln) == "+"]
  minus <- aln[strand(aln) == "-"]
  ## get 5'ends separately for plus and minus strand
  ends_plus <- GRanges(seqnames = seqnames(plus),
                       ranges = IRanges(start = start(plus), width = rep(1, length(plus))),
                       strand = strand(plus))
  ends_minus <- GRanges(seqnames = seqnames(minus),
                        ranges = IRanges(start = end(minus), width = rep(1, length(minus))),
                        strand = strand(minus))
  ends_uniq_plus <- unique(ends_plus)
  ends_uniq_plus$hits<- countOverlaps(ends_uniq_plus, ends_plus, type="equal") # get the coverage
  ends_uniq_minus <- unique(ends_minus)
  ends_uniq_minus$hits<- countOverlaps(ends_uniq_minus, ends_minus, type="equal") # get the coverage
  ### 5'ends
  ends_uniq <- c(ends_uniq_plus, ends_uniq_minus)
  return(ends_uniq)
}

#' Get top percentage of transcript riboseq reads
#'
#' @param ribo a list of lists with overlaps of transcripts and ribo-seq reads
#' @param top_tx numeric(10), which top percentage of transcripts to use.
#' @return a list of lists, with accepted lists in ribo
#' @importFrom S4Vectors tail
get_top_tx <- function(ribo, top_tx) {
  idx <- ceiling(length(ribo)/(100/top_tx))
  sums <- sapply(ribo,sum)
  thr <- as.numeric(tail(sort(sums),idx+1)[1])
  ribo <- ribo[sums >= thr]

  return(ribo)
}

#' Find if there are periodicity in sample
#'
#' @param cds150_over a GRanges of tiled cds' to 150 length
#' @param ends_uniq a GRanges object with 5' ends of reads
#' @param top_tx numeric(10), which top percentage of transcripts to use.
#' @return a logical, if it is periodic.
#' @importFrom stats fft spec.pgram
periodicity <- function(cds150_over, ends_uniq, top_tx) {
  cds150_overlaps <- findOverlaps(cds150_over, ends_uniq)
  cds150_over$riboseq <- rep(0, length(cds150_over))
  cds150_over$riboseq[from(cds150_overlaps)] <-
    ends_uniq[to(cds150_overlaps)]$hits
  ## split cds150 per 150nt (seqnames), add up riboseq
  ribo <- split(cds150_over$riboseq, ceiling(seq_along(cds150_over$riboseq)/150))
  #meta <- Reduce("+", ribo) ## all tx
  ribo <- get_top_tx(ribo, top_tx) ## top n% tx
  meta <- Reduce("+", ribo) ## from base
  ## FFT
  amplitudes <- abs(fft(meta)) #!!!!!!!!!!!! WHAT happens here???
  amp <- amplitudes[2:(length(amplitudes)/2+1)]

  periods <- 1/spec.pgram(x = meta, plot = F)$freq
  if ((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1)) {
    periodic <- TRUE
  } else {
    periodic <- FALSE
  }
  return(periodic)
}

#' Get number of ribo reads per tile position of cds start window
#' @param start a GRangesList or logical, of tiled cds
#'  start windows of 60 length
#' @param ends_uniq a GRanges object with 5' ends of reads
#' @param top_tx numeric(10), which top percentage of transcripts to use.
#' @return a list of reads per position
start_ribo <- function(start, ends_uniq, top_tx) {
  start_over <- unlist(start)
  start_overlaps <- findOverlaps(start_over, ends_uniq)
  start_over$riboseq <- rep(0, length(start_over))
  start_over$riboseq[from(start_overlaps)] <-
    ends_uniq[to(start_overlaps)]$hits
  ## split start per 60nt (seqnames), add up riboseq
  start_ribo <- split(start_over$riboseq,
                      ceiling(seq_along(start_over$riboseq)/60))
  #start_meta <- Reduce("+", start_ribo) ## all tx
  start_ribo <- get_top_tx(start_ribo, top_tx) ## top n% tx
  start_meta <- Reduce("+", start_ribo)
  return(start_meta)
}

#' Get number of ribo reads per tile position of cds stop window
#' @param stop a GRangesList or logical, of tiled cds
#'  stop windows of 60 length
#' @param ends_uniq a GRanges object with 5' ends of reads
#' @param top_tx numeric(10), which top percentage of transcripts to use.
#' @return a list of reads per position
stop_ribo <- function(stop, ends_uniq, top_tx) {
  stop_over <- unlist(stop)
  stop_overlaps <- findOverlaps(stop_over, ends_uniq)
  stop_over$riboseq <- rep(0, length(stop_over))
  stop_over$riboseq[from(stop_overlaps)] <- ends_uniq[to(stop_overlaps)]$hits
  ## split start per 60nt (seqnames), add up riboseq
  stop_ribo <- split(stop_over$riboseq,
                     ceiling(seq_along(stop_over$riboseq)/60))
  #start_meta <- Reduce("+", start_ribo) ## all tx
  stop_ribo <- get_top_tx(stop_ribo, top_tx) ## top n% tx
  stop_meta <- Reduce("+", stop_ribo)
  return(stop_meta)
}

#' Get the offset for specific ribo-seq read width
#' @param df a data.frame with points to analyse,
#'  from function get_start_stop
#' @param feature a character, (start) or stop
#' @return a single numeric offset
change_point_analysis <- function(df, feature="start") {
  if (feature == "start") {
    meta <- df$count[1:40]
    means <- c()
    for (j in 15:35) {
      m <- mean(meta[j:40]) - mean(meta[1:(j-1)])
      means <- c(means, m)
    }
    shift <- which.max(abs(means)) + 14
    offset <- df$codon[shift]
  } else if (feature == "stop") {
    meta <- df$count[1:40]
    shift <- which.max(meta)
    offset <- df$codon[shift] + 3
  }
  return(offset)
}
