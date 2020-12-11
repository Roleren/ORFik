#' Make sequence region heatmap relative to scoring
#'
#' Given sequences, DNA or RNA.
#' And some score, ribo-seq fpkm, TE etc.
#' Create a heatmap divided per letter in seqs, by how strong the score is.
#'
#' It will create blocks around the highest rate per position
#' @param seqs the sequences (character vector, DNAStringSet)
#' @param rate a scoring vector (equal size to seqs)
#' @param start position in seqs to start at (first is 1), default 1.
#' @param stop position in seqs to stop at (first is 1),
#'  default max(nchar(seqs)), that is the longest sequence length
#' @param center position in seqs to center at (first is 1), center will
#' be +1 in heatmap
#' @param min.observations How many observations per position per letter
#' to accept? numeric or quantile, default (">q1", bigger than quartile 1
#' (25 percentile)). You can do (10), to get all with more than
#' 10 observations.
#' @param skip.startCodon startCodon is defined as after centering
#'  (position 1, 2 and 3).
#' Should they be skipped ? default (FALSE).
#' Not relevant if you are not doing Translation initiation sites (TIS).
#' @param xlab Region you are checking, default (TIS)
#' @param type What type is the rate scoring ? default (ribo-seq)
#' @return a ggplot of the heatmap
#' @importFrom data.table melt
#' @export
#' @examples
#'
#' \dontrun{
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'   txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#'   #Extract sequences of Coding sequences.
#'   cds <- loadRegion(txdbFile, "cds")
#'   tx <- loadRegion(txdbFile, "mrna")
#'
#'   # Get region to check
#'   kozakRegions <- startRegionString(cds, tx, BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#'                                     , upstream = 4, 5)
#'   # Some toy ribo-seq fpkm scores on cds
#'   set.seed(3)
#'   fpkm <- sample(1:115, length(cds), replace = TRUE)
#'   kozakHeatmap(kozakRegions, fpkm, 1, 9, skip.startCodon = F)
#' }
#' }
#'
kozakHeatmap <- function(seqs, rate, start = 1, stop = max(nchar(seqs)),
                         center = ceiling((stop - start + 1)/2),
                         min.observations = ">q1", skip.startCodon = FALSE,
                         xlab = "TIS", type = "ribo-seq") {
  if (length(seqs) != length(rate)) stop("Length of rate and seq must be equal!")
  if (length(seqs) == 0 | length(rate) == 0) stop("Length of rate and seq must be > 0!")
  if (start > stop) stop("Stop must be bigger or equal to start")
  if (is.null(names(seqs))) names(seqs) <- as.character(seq(1, length(seqs)))

  # Update to only full length seqs
  hits <- nchar(seqs) >= stop
  seqs <- seqs[hits]; rate <- rate[hits]
  if (!all(hits)) { # inform that not all seqs are full length
    message("Not all sequences can be used, not full length")
    print(summary(hits))
  } else if (all(!hits)) stop("No seqs long enough for interval start, stop!")

  dt <- data.table(X.gene_id = names(seqs))
  vars <- c()
  for (i in seq(start, stop)) {
    dt[,paste0("seq", i)] <- substring(seqs, i, i)
    vars <- c(vars, paste0("seq", i))
  }
  dt$rate <- rate


  dt.melt <- melt(dt, id.vars = c("X.gene_id", "rate"))
  codon.table <- dt.melt[, .(median_score = median(rate, na.rm = TRUE),
                             count_seq_pos_with_count = .N),
                         by = .(variable, value)]

  uniques <- rev(unique(dt.melt$value))
  codon.table$value <- factor(codon.table$value, uniques)
  xPos <- seq(start, stop + 1) - start + 1 - center
  xPos <- xPos[-center]
  xPos[xPos > 0] <- paste0("+", xPos[xPos > 0])
  codon.table$variable <- factor(codon.table$variable, levels = vars,
                                 labels = xPos)

  if (skip.startCodon) {
    codon.table[codon.table$variable %in% c("+1", "+2", "+3"), ]$median_score <- NA
  }

  #output position where we have more than 1st quartile observations etc
  print("Distribution of observations per position per letter")
  print(summary(codon.table$count_seq_pos_with_count))
  if (min.observations == ">q1") {
    quart <- summary(codon.table$count_seq_pos_with_count)[2]
  } else if (is.numeric(min.observations)) {
    quart <- min.observations
  } else stop("min.observations must be >q1 or a numeric!")

  print(paste0("Picking: >", quart, " observations"))
  codon.table.filtered <- codon.table[count_seq_pos_with_count > quart,]
  if (nrow(codon.table.filtered) == 0) stop("No rows passed min.observations!")

  plot_matrix2_log <- ggplot(data=codon.table.filtered,
                             aes(x=variable, y=value, fill=log2(median_score))) +
    theme(panel.background=element_rect(fill="lightgrey", colour="lightgrey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_tile(color = "lightgrey") +
    scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = 'lightgrey',
                         name = paste0("log2(median ", type,")")) +
    xlab(paste0("Position relative to ", xlab)) +
    ylab("Nucleotide")

  rows <- seq(from = length(unique(codon.table.filtered$value)) - 0.5, to = 0.5)
  names(rows) <- unique(codon.table.filtered$value)
  xmin <- -0.5
  # Per position/column (set the black max box)
  for(col in unique(codon.table.filtered$variable)) {
    mat <- codon.table.filtered[codon.table.filtered$variable == col,]
    # Bug here, if not all letters are remaining

    highest <- rows[names(rows) == mat$value[which.max(mat$median_score)]]
    xmin = xmin + 1;xmax = xmin + 1; ymin = highest;ymax = ymin + 1
    if (length(highest) == 0) next
    input <- paste0("geom_rect(aes(xmin=",xmin,",xmax=",xmax,",ymin=",
                    ymin,",ymax=",ymax,"), color='black', size=0.5, fill=NA)")
    plot_matrix2_log <- plot_matrix2_log +
      eval(parse(text = input))
  }
  return(plot_matrix2_log)
}

#' Rank kozak initiation sequences
#'
#' Defined as region (-4, -1) relative to TIS
#' @param cds_k cds ranges (GRangesList)
#' @param mrna mrna ranges (GRangesList)
#' @param dt.ir data.table with a column called IR, initiation rate
#' @param group.min numeric, default 10. Minimum transcripts per initation
#' group to be included
#' @inheritParams kozakSequenceScore
#' @return a ggplot grid object
#' @export
kozak_IR_ranking <- function(cds_k, mrna, dt.ir, group.min = 10,
                             species = "human") {
  seqs <- startRegionString(cds_k, tx = mrna, faFile = df, upstream = 4, downstream = -1)
  dt.ir$seqs <- seqs

  kozak_scores <- kozakSequenceScore(cds_k, mrna, df, species = "zebrafish")
  dt.ir$upstream_kozak_strength <- kozak_scores

  merged_rates <- dt.ir[, .(mean_IR = mean(IR),
                            upstream_kozak_strength = mean(upstream_kozak_strength),
                            count = .N)
                        , by = seqs]
  print("Distribution of observations group")
  print(summary(merged_rates$count))
  message(paste("Originally", nrow(merged_rates), "groups of sequences"))
  merged_rates <- merged_rates[ count > group.min,]
  message(paste("Filter down to", nrow(merged_rates), "groups of sequences"))

  merged_rates$seqs <- factor(merged_rates$seqs, levels = merged_rates[order(mean_IR),]$seqs)
  scales_kozak <- c(ceiling(min(merged_rates$upstream_kozak_strength * 10)) / 10)
  scales_kozak <- c(scales_kozak, scales_kozak + 0.1)

  mean_IR_20 <- ggplot(data=merged_rates, aes(x=seqs, y=1, fill=log2(mean_IR))) +
    geom_tile() +
    scale_fill_gradientn(colors = c("blue", "white", "red1", "red2", "red3")) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "IR") +
    xlab("Initiation sequence") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.key.size = unit(0.4, "cm"))
  mean_IR_20
  tile_KS_20 <- ggplot(data=merged_rates, aes(x=seqs, y=1, fill=upstream_kozak_strength)) +
    geom_tile() +
    # scale_fill_gradientn(colors = c("#146DA8", "white", "#A87225")) +
    scale_fill_gradientn(colors = c("blue", "white", "red"), breaks = scales_kozak) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "Kozak") +
    xlab("Initiation sequence") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.key.size = unit(0.4, "cm"))
  tile_KS_20
  ranking <- gridExtra::grid.arrange(mean_IR_20, tile_KS_20, ncol = 1)
  return(ranking)
}

#' TOP Motif ecdf plot
#'
#' Given sequences, DNA or RNA.
#' And some score, scanning efficiency (SE), ribo-seq fpkm, TE etc.
#'
#' Top motif defined as a TSS of C and 4 T's or C's (pyrimidins) downstream
#' of TSS C.
#'
#' The right plot groups:
#' C nucleotide, TOP motif (C, then 4 pyrimidines) and
#' OTHER (all other TSS variants).
#' @param seqs the sequences (character vector, DNAStringSet),
#' of 5' UTRs (leaders). See example below for input.
#' @param rate a scoring vector (equal size to seqs)
#' @param start position in seqs to start at (first is 1), default 1.
#' @param stop position in seqs to stop at (first is 1),
#'  default max(nchar(seqs)), that is the longest sequence length
#' @param xlim What interval of rate values you want to show
#' type: numeric or quantile of length 2,
#' 1. default c("q10","q99"). bigger than 10 percentile
#'  and less than 99 percentile.
#' 2. Set to numeric values, like c(5, 1000),
#' 3. Set to NULL if you want all values. Backend uses coord_cartesian.
#' @param type What type is the rate scoring ? default ("Scanning efficiency")
#' @param legend.position.1st adjust left plot label position, default c(0.75, 0.28),
#' ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.position.motif adjust right plot label position, default c(0.75, 0.28),
#' ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @return a ggplot gtable of the TOP motifs in 2 plots
#' @importFrom cowplot plot_grid
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#'
#' \dontrun{
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'   txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package = "GenomicFeatures")
#'   #Extract sequences of Coding sequences.
#'   leaders <- loadRegion(txdbFile, "leaders")
#'
#'   # Should update by CAGE if not already done
#'   cageData <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'                           package = "ORFik")
#'   leadersCage <- reassignTSSbyCage(leaders, cageData)
#'   # Get region to check
#'   seqs <- startRegionString(leadersCage, NULL,
#'         BSgenome.Hsapiens.UCSC.hg19::Hsapiens, 0, 4)
#'   # Some toy ribo-seq fpkm scores on cds
#'   set.seed(3)
#'   fpkm <- sample(1:115, length(leadersCage), replace = TRUE)
#'   # Standard arguments
#'   TOP.Motif.ecdf(seqs, fpkm, type = "ribo-seq FPKM",
#'                  legend.position.1st = "bottom",
#'                  legend.position.motif = "bottom")
#'   # with no zoom on x-axis:
#'   TOP.Motif.ecdf(seqs, fpkm, xlim = NULL,
#'                  legend.position.1st = "bottom",
#'                  legend.position.motif = "bottom")
#' }
#' }
#'
TOP.Motif.ecdf <- function(seqs, rate, start = 1, stop = max(nchar(seqs)),
                         xlim = c("q10","q99"),
                         type = "Scanning efficiency",
                         legend.position.1st = c(0.75, 0.28),
                         legend.position.motif = c(0.75, 0.28)) {
  if (length(seqs) != length(rate)) stop("Length of rate and seq must be equal!")
  if (length(seqs) == 0 | length(rate) == 0) stop("Length of rate and seq must be > 0!")
  if (start > stop) stop("Stop must be bigger or equal to start")

  # Update to only full length seqs
  hits <- nchar(seqs) >= stop
  seqs <- seqs[hits]; rate <- rate[hits]
  if (!all(hits)) { # inform that not all seqs are full length
    message("Not all sequences can be used, here is how many valid (TRUE):")
    print(summary(hits))
  } else if (all(!hits)) stop("No seqs long enough for interval start, stop!")

  dt <- topMotif(seqs, start, stop)
  message("Distribution of TOP motifs:")
  print(table(dt$TOP))

  dt$rate <- rate
  # Define plot settings:
  new_pallet_1 <- c("#E495A5","#ABB065","#39BEB1","#ACA4E2")
  new_pallet_2 <- c("#39BEB1", "#ACA4E2", "#E495A5")

  tl <- theme(legend.position = legend.position.1st, legend.background=element_blank())
  tlb <- theme(legend.position = legend.position.motif, legend.background=element_blank(),
               axis.ticks.y=element_blank(), axis.text.y = element_blank())
  tit <- labs(color = "1st nucleotide")
  titb <- labs(color = "Motif")

  se1 <- ggplot(data=dt, aes((rate), colour = seq1)) +
    stat_ecdf() +
    scale_x_log10() +
    scale_color_manual(values=new_pallet_1) +
    ylab("") +
    xlab("") +
    theme_bw() +
    tl + tit + theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

  se2 <- ggplot(data=dt, aes((rate), colour = TOP)) +
    stat_ecdf() +
    scale_x_log10() +
    scale_color_manual(values=new_pallet_2) +
    ylab("") +
    xlab("") +
    theme_bw() +
    tlb + titb + theme(plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"))

  if (!is.null(xlim)) {
    if (length(xlim) != 2) stop("xlim must be length 2 if defined!")
    if (is.character(xlim)) {
      xlim <- as.numeric(gsub("q", "", xlim)) / 100
      if (any(xlim > 1 | xlim < 0)) stop("when xlim is percentiles,",
                                         "range must be 0 to 100")
      xlim <- quantile(dt$rate, xlim)
    }
    if (!is.numeric(xlim)) stop("xlim must be valid character or numeric!")

    se1 <- se1 + coord_cartesian(xlim = xlim)
    se2 <- se2 + coord_cartesian(xlim = xlim)
  }

  comb <- plot_grid(se1, se2, nrow = 1, align = "v")
  comb <- grid.arrange(comb, bottom = paste0("log10(", type, ")"))
  return(comb)
}

