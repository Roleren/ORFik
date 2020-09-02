#' Filter out transcript by a median filter
#'
#' For removing extreme peaks in coverage plots, use high quantiles, like
#' 99. If you want to find general stall sites of Ribo-seq, you should use a
#' pre_filter minimum and
#' @param tx a GRangesList
#' @param reads a GAlignments or GRanges
#' @param multiplier a character or numeric, default "0.99",
#' either a quantile if input is string[0-1],
#' like "0.99", or numeric value if input is numeric.
#' How much bigger than median / mean counts per gene,
#' must a value be to be defined as extreme ?
#' @param min_cutoff a character or numeric, default "0.999",
#' either a quantile if input is string[0-1],
#' like "0.999", or numeric value if input is numeric. Lowest allowed value
#' @param pre_filter_minimum numeric, default 0. If value is x,
#' will remove all positions in all genes with coverage < x,
#' before median filter is applied. Set to 1 to remove all 0 positions.
#' @param average character, default "median". Alternative: "mean".
#' How to scale the multiplier argument, from median or mean of gene coverage.
#' @return GRangesList (filtered)
#' @export
filterExtremePeakGenes <-
  function(tx, reads, upstream = NULL, downstream = NULL,
           multiplier = "0.99", min_cutoff = "0.999",
           pre_filter_minimum = 0, average = "median") {
  if (!is.null(upstream) & !is.null(downstream)) {
    region <- startRegion(tx,
                          extendLeaders(tx, extension = max(0, upstream + 1)),
                          upstream = upstream, downstream = downstream)
  } else region <- tx

  coverage <- coveragePerTiling(region, reads, TRUE, as.data.table = TRUE)
  if (pre_filter_minimum > 0)
    coverage <- coverage[count >= pre_filter_minimum]
  if (average == "median") {
    coverage[, median_per_gene := median(count), by = genes]
  } else if (average == "mean") {
    coverage[, median_per_gene := mean(count), by = genes]
  } else stop("average must be median or mean!")

  if (is(min_cutoff, "character")) {
    if (is.na(as.numeric(min_cutoff)))
      stop("min_cutoff must be numeric or character of numeric value")
    min_cutoff <- quantile(coverage$count, as.numeric(min_cutoff))
  }
  if (is(multiplier, "character")) {
    if (is.na(as.numeric(multiplier)))
      stop("multiplier must be numeric or character of numeric value")
    multiplier <- quantile(coverage$count, as.numeric(multiplier))
  }
  toFilter <- coverage[count > median_per_gene*multiplier + min_cutoff]
  print("Number of 0-score positions and ratio:")
  print(coverage[, .(zeropos_num = sum(count == 0),
                     ratio = sum(count == 0) / .N)])
  print(paste("Using", average, "multiplier and minimum cutoff of:"))
  print(c(multiplier, min_cutoff))
  print("Read count distribution:");print(summary(coverage$count))
  print("Hit positions in genes that will be filtered:");print(toFilter)
  print("Transcripts filtered out:");print(names(tx[unique(toFilter$genes)]))
  if (length(toFilter$genes) == 0) {
    message("No genes filtered out, reduce filters if you want hits")
    return(tx)
  }
  return(tx[-unique(toFilter$genes)])
  }

#' Find max peak per gene
#'
#' For finding the highest peak per gene. With some default filters.
#' @param tx a GRangesList
#' @param reads a GAlignments or GRanges
#' @param top_tx numeric, default 0.50 (50\% top tx by reads).
#' @param min_reads_per_tx numeric, default 10. Gene must have at least
#'  10 reads
#' @return a data.table
#' @export
findMaxPeakPerGenes <- function(tx, reads, top_tx = 0.50,
                                min_reads_per_tx = 10) {
    # coverage track
    coverage <- coveragePerTiling(tx, reads, TRUE, as.data.table = TRUE)
    coverage[, sum_per_gene := sum(count), by = genes]
    coverage <- coverage[sum_per_gene >= max(quantile(sum_per_gene, tx_to_use),
                                             min.reads.per.tx),]
    coverage[, mean_per_gene := mean(count), by = genes]
    coverage[, sd_per_gene := sd(count), by = genes]
    coverage[, zscore := (count - mean_per_gene) / sd_per_gene, by = genes]
    summary(coverage$zscore)
    coverage[, gene_id := names(tx)[genes]]
    coverage <- coverage[coverage[, .I[zscore == max(zscore)], by = genes]$V1]
    return(coverage)
}
