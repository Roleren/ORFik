#' Make sequence region heatmap relative to scoring
#'
#' Given sequences, DNA or RNA.
#' And some score, ribo-seq fpkm, TE etc.
#' Create a heatmap divided per letter in seqs, by how strong the score is.
#'
#' It will create blocks around best rate per position
#' @param seqs the sequences (character vector, DNAStringSet)
#' @param rate a scoring vector (equal size to seqs)
#' @param start position in seqs to start at (first is 1)
#' @param stop position in seqs to stop at (first is 1)
#' @param center position in seqs to center at (first is 1), center will
#' be +1 in heatmap
#' @param min.observations How many observations per position per letter to accept?
#' numeric or quantile, default (">q1", bigger than quartile 1 (25 percentile)).
#' You can do (10), to get all with more than 10 observations.
#' @param skip.startCodon startCodon is defined as after centering (position 1, 2 and 3).
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
kozakHeatmap <- function(seqs, rate, start, stop, center = ceiling((stop - start + 1)/2),
                         min.observations = ">q1", skip.startCodon = FALSE,
                         xlab = "TIS", type = "ribo-seq") {
  if (length(seqs) != length(rate)) stop("Length of rate and seq must be equal!")
  if (length(seqs) == 0 | length(rate) == 0) stop("Length of rate and seq must be > 0!")
  if (is.null(names(seqs))) names(seqs) <- as.character(seq(1, length(seqs)))

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
  codon.table$variable <- factor(codon.table$variable, levels=vars,
                                 labels=xPos)

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
    theme(panel.background=element_rect(fill="lightgrey", colour="lightgrey")) +
    geom_tile(color = "lightgrey") +
    scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = 'lightgrey',
                         name = paste0("log2(median_", type,")")) +
    xlab(paste0("Position realitive to ", xlab)) +
    ylab("Nucleotide")

  rows <- seq(from = length(unique(codon.table.filtered$value)) - 0.5, to = 0.5)
  names(rows) <- unique(codon.table.filtered$value)
  xmin <- -0.5
  for(col in unique(codon.table.filtered$variable)){
    mat <- codon.table.filtered[codon.table.filtered$variable == col,]
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
