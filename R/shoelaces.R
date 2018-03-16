### shoelaces for ORFik
library(GenomicFeatures)
library(GenomicAlignments)
library(plyr)
library(GeneCycle)
library(ggplot2)

### read in BAM and GTF files
gtf <- "/Users/kasia/Documents/shoelaces/Data/example.gtf"
bam <- "/Users/kasia/Documents/shoelaces/Data/example.bam"
#noise <- #
shifts <- shoelaces(bam, gtf, start=T, stop=F, offset_plots=T, top_tx=10)


read_bam <- function(bam) {
  ########## BAM (get only uniquely mapping reads) ##########
  param <- ScanBamParam(tag="NH")
  alignments <- readGAlignments(bam, param=param)
  alignments <- alignments[mcols(alignments)$NH == 1]
  return(alignments)
}

get_longest_tx <- function(txdb) {
  ## tx with longest CDS
  txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
  rownames(txLengths) <- txLengths$tx_name
  longest_tx <- arrange(txLengths, gene_id, desc(cds_len))
  longest_tx <- longest_tx[!duplicated(longest_tx$gene_id),]
  rownames(longest_tx) <- longest_tx$tx_name
  # get CDS longer than 150, UTRs longer than 30
  longest_tx <- longest_tx[longest_tx$utr5_len > 30,]
  longest_tx <- longest_tx[longest_tx$cds_len > 150,]
  longest_tx <- longest_tx[longest_tx$utr3_len > 30,]
  return(longest_tx)
}

get_cds150 <- function(txdb, longest_tx) {
  ########## CDS150 for periodicity ##########
  cds <- cdsBy(txdb, by="tx", use.names=TRUE)
  cds <- cds[order(names(cds))]
  cds <- cds[names(cds) %in% longest_tx$tx_name]
  ucds <- unlist(cds)
  tucds <- tile(ucds, width =  1L)
  names(tucds) <- names(ucds)
  ucds <- unlist(tucds)
  ucdsl <- split(ucds, names(ucds))
  ## reverse minus strand
  ucdsl_plus <- ucdsl[as.logical(runValue(strand(ucdsl) == "+"))]
  ucdsl_minus <- ucdsl[as.logical(runValue(strand(ucdsl) == "-"))]
  ## get first 150nt
  cds150_plus <- GRangesList(sapply(ucdsl_plus, function(x){x[1:150]}))
  cds150_minus <- GRangesList(sapply(ucdsl_minus, function(x){rev(x[(length(x)-149):(length(x))])}))
  cds150 <- c(cds150_plus, cds150_minus)
  cds150@unlistData@ranges@NAMES <- NULL
  return(cds150)
}

get_start_stop <- function(txdb, longest_tx, start=T, stop=T) {
  ########## START/STOP ##########
  exons <- exonsBy(txdb, by="tx", use.names=TRUE)
  exons <- exons[order(names(exons))]
  exons <- exons[names(exons) %in% longest_tx$tx_name]
  uexons <- unlist(exons)
  tuexons <- tile(uexons, width =  1L)
  names(tuexons) <- names(uexons)
  uexons <- unlist(tuexons)
  uexonsl <- split(uexons, names(uexons))
  ## get start and end separately for minus strand
  uexonsl_plus <- uexonsl[as.logical(runValue(strand(uexonsl) == "+"))]
  uexonsl_plus <- uexonsl_plus[order(names(uexonsl_plus))]
  uexonsl_minus <- uexonsl[as.logical(runValue(strand(uexonsl) == "-"))]
  uexonsl_minus <- uexonsl_minus[order(names(uexonsl_minus))]
  ##
  tx_plus <- longest_tx[rownames(longest_tx) %in% names(uexonsl_plus),]
  tx_plus <- tx_plus[order(rownames(tx_plus)),]
  tx_minus <- longest_tx[rownames(longest_tx) %in% names(uexonsl_minus),]
  tx_minus <- tx_minus[order(rownames(tx_minus)),]
  
  ########## START ##########
  if (start) {
    start_plus <- tx_plus$utr5_len + 1
    start_minus <- tx_minus$tx_len - tx_minus$utr5_len
    ## subset
    start_plus <- GRangesList(mapply(function(x,y){x[(y-30):(y+29)]}, uexonsl_plus, start_plus))
    start_minus <- GRangesList(mapply(function(x,y){rev(x[(y-29):(y+30)])}, uexonsl_minus, start_minus))
    start <- c(start_plus, start_minus)
    start@unlistData@ranges@NAMES <- NULL
  }
  
  ########## STOP ##########
  if (stop) {
    stop_plus <- tx_plus$utr5_len + tx_plus$cds_len - 2
    stop_minus <- tx_minus$utr3_len + 3
    ## subset
    stop_plus <- GRangesList(mapply(function(x,y){x[(y-30):(y+29)]}, uexonsl_plus, stop_plus))
    stop_minus <- GRangesList(mapply(function(x,y){rev(x[(y-29):(y+30)])}, uexonsl_minus, stop_minus))
    stop <- c(stop_plus, stop_minus)
    stop@unlistData@ranges@NAMES <- NULL
  }
  
  ss <- list(start, stop)
  names(ss) <- c("start", "stop")
  
  return(ss)
}

get_5ends <- function(alignments, read_length) {
  aln <- alignments[qwidth(alignments) == read_length]
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

get_top_tx <- function(ribo, top_tx) {
  idx <- ceiling(length(ribo)/(100/top_tx))
  sums <- sapply(ribo,sum)
  thr <- as.numeric(tail(sort(sums),idx+1)[1])
  ribo <- ribo[sums > thr]
  return(ribo)
}

periodicity <- function(cds150, ends_uniq, top_tx) {
  cds150_over <- unlist(cds150)
  cds150_overlaps <- findOverlaps(cds150_over, ends_uniq)
  cds150_over$riboseq <- rep(0, length(cds150_over))
  cds150_over$riboseq[from(cds150_overlaps)] <- ends_uniq[to(cds150_overlaps)]$hits
  ribo <- split(cds150_over$riboseq, ceiling(seq_along(cds150_over$riboseq)/150)) ## split cds150 per 150nt (seqnames), add up riboseq
  #meta <- Reduce("+", ribo) ## all tx 
  ribo <- get_top_tx(ribo, top_tx) ## top n% tx
  meta <- Reduce("+", ribo)
  ## FFT
  amplitudes <- abs(fft(meta))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(meta, method = "clone")$freq
  if ((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1)) {
    periodic <- T
  } else {
    periodic <- F
  }
  return(periodic)
}

start_ribo <- function(start, ends_uniq, top_tx) {
  start_over <- unlist(start)
  start_overlaps <- findOverlaps(start_over, ends_uniq)
  start_over$riboseq <- rep(0, length(start_over))
  start_over$riboseq[from(start_overlaps)] <- ends_uniq[to(start_overlaps)]$hits
  start_ribo <- split(start_over$riboseq, ceiling(seq_along(start_over$riboseq)/60)) ## split start per 60nt (seqnames), add up riboseq
  #start_meta <- Reduce("+", start_ribo) ## all tx
  start_ribo <- get_top_tx(start_ribo, top_tx) ## top n% tx
  start_meta <- Reduce("+", start_ribo)
  return(start_meta)
}

stop_ribo <- function(stop, ends_uniq, top_tx) {
  stop_over <- unlist(stop)
  stop_overlaps <- findOverlaps(stop_over, ends_uniq)
  stop_over$riboseq <- rep(0, length(stop_over))
  stop_over$riboseq[from(stop_overlaps)] <- ends_uniq[to(stop_overlaps)]$hits
  stop_ribo <- split(stop_over$riboseq, ceiling(seq_along(stop_over$riboseq)/60)) ## split start per 60nt (seqnames), add up riboseq
  #start_meta <- Reduce("+", start_ribo) ## all tx
  stop_ribo <- get_top_tx(stop_ribo, top_tx) ## top n% tx
  stop_meta <- Reduce("+", stop_ribo)
  return(stop_meta)
}

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
    offset <- df$codon[shift]
  }
  return(offset)
}



################# PLOT ##################
plot_periodic_lengths_start <- function(start_df) {
  colnames(start_df) <- c("start_codon", "count", "flength")
  start_df$fill <- factor(rep(c(0,1,1),10))
  start_df$start_codon <- as.character(start_df$start_codon)
  start_df$start_codon <- factor(start_df$start_codon, levels = start_df$start_codon[1:60])
  p <- ggplot(start_df, aes(x=start_codon, y=count, fill=fill)) + geom_bar(stat="identity") +
    scale_fill_manual(values=c("#E72B3C", "#A9A9A9")) + facet_grid( ~ flength) + theme(legend.position="none")
  rect <- data.frame(xmin=30.5, xmax=33.5, ymin=-Inf, ymax=Inf) # rectangle
  p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                     fill="green", alpha=0.25, inherit.aes = FALSE)
  ll <- levels(start_df$start_codon)[c(T, F, F)] # breaks and labels
  p <- p + scale_x_discrete(breaks=ll, labels = ll)
  return(p)
}

plot_periodic_lengths_stop <- function(stop_df) {
  colnames(stop_df) <- c("stop_codon", "count", "flength")
  stop_df$fill <- factor(rep(c(0,1,1),10))
  stop_df$stop_codon <- as.character(stop_df$stop_codon)
  stop_df$stop_codon <- factor(stop_df$stop_codon, levels = stop_df$stop_codon[1:60])
  p <- ggplot(stop_df, aes(x=stop_codon, y=count, fill=fill)) + geom_bar(stat="identity") +
    scale_fill_manual(values=c("#E72B3C", "#A9A9A9")) + facet_grid( ~ flength) + theme(legend.position="none")
  rect <- data.frame(xmin=30.5, xmax=33.5, ymin=-Inf, ymax=Inf) # rectangle
  p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                     fill="red", alpha=0.25, inherit.aes = FALSE)
  ll <- levels(stop_df$stop_codon)[c(T, F, F)] # breaks and labels
  p <- p + scale_x_discrete(breaks=ll, labels = ll)
  return(p)
}


shoelaces <- function(bam, gtf, start=T, stop=F, offset_plots=T, top_tx=10) {
  ## initialize
  alignments <- read_bam(bam)
  txdb <- makeTxDbFromGFF(gtf, format = "gtf")
  longest_tx <- get_longest_tx(txdb)
  cds150 <- get_cds150(txdb, longest_tx)
  ss <- get_start_stop(txdb, longest_tx, start=start, stop=stop)
  
  ## footprint lengths
  all_lengths <- sort(unique(qwidth(alignments)))
  selected_lengths <- c()
  offsets_start <- c()
  offsets_stop <- c()
  
  if (start) {
    start_df <- data.frame(start_codon=c(), count, flength=c())
  }
  if (stop) {
    stop_df <- data.frame(stop_codon=c(), count, flength=c())
  }
  
  for(i in 1:length(all_lengths)) {
    ### 5'ends
    ends_uniq <- get_5ends(alignments, all_lengths[i])
    periodic <- periodicity(cds150, ends_uniq, top_tx=10)
    if (periodic) {
      selected_lengths <- c(selected_lengths, all_lengths[i])
      if (start) {
        start_meta <- start_ribo(ss[[1]], ends_uniq, top_tx)
        df <- data.frame(codon=c(-30:30)[-31], count=start_meta, flength=rep(all_lengths[i], 60))
        start_df <- rbind(start_df, df)
        ### change point analysis
        offset <- change_point_analysis(df, feature="start")
        offsets_start <- c(offsets_start, offset)
      }
      if (stop) {
        stop_meta <- stop_ribo(ss[[2]], ends_uniq, top_tx)
        df <- data.frame(codon=c(-27:33)[-28], count=stop_meta, flength=rep(all_lengths[i], 60))
        stop_df <- rbind(stop_df, df)
        ### change point analysis
        offset <- change_point_analysis(df, feature="stop")
        offsets_stop <- c(offsets_stop, offset)
      }
    }
  }
  
  ### print plots
  if (offset_plots) {
    if (start) {
      pstart <- plot_periodic_lengths_start(start_df)
      print(pstart)
    }
    if (stop) {
      pstop <- plot_periodic_lengths_stop(stop_df)
      print(pstop)
    }
  }
  
  ## return lengths and offsets
  shifts <- data.frame(fragment_length = selected_lengths)
  if (start) {
    shifts$offsets_start <-  offsets_start
  }
  if (stop) {
    shifts$offsets_stop <-  offsets_stop
  }
  return(shifts)
  
}


