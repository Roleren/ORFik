seq_usage <- function(dt, seqs, genes, input.dt.length = 1, output.seq.length = 3,
                      seqs.order.table = NULL, dispersion_method = "MME") {
  stopifnot(is(dt, "data.table"))
  if ("genes" %in% colnames(dt)) {
    stop("dt can not contain column called 'genes'",
         "it should be given seperatly in argument 'genes'!")
  }
  if (length(genes) != nrow(dt)) stop("Length of dt and genes must match!")
  if (length(seqs) != (nrow(dt) / output.seq.length))
    stop("Size of 'seqs' must be nrow(dt) / output.seq.length")
  codon_sums <- copy(dt)
  n.samples <- ncol(dt)
  if (n.samples > 1) {
    codon_sums <- suppressWarnings(melt.data.table(codon_sums,
                                                   value.name = "score"))
  } else {
    stopifnot("score" %in% colnames(codon_sums))
    codon_sums[, variable := "library 1"]
  }
  codon_sums[, genes := rep(genes, length.out = .N)]
  we_must_collapse <- input.dt.length != output.seq.length
  if (we_must_collapse) {
    codon_sums[, codon_sum := frollsum(x = score, n = output.seq.length, align = "left")]
    # Keep only first per seq group
    codon_sums <- codon_sums[rep(c(T, rep(FALSE, output.seq.length - 1)), length.out = .N),]
  }

  codon_sums[, seqs := rep(seqs, length.out = .N)]

  # Normalize
  codon_sums[, `:=`(gene_sum = sum(score, na.rm = TRUE)), by = .(variable, genes)]
  codon_sums[, `:=`(N_AA_of_type_per_gene = .N), by = .(variable, genes, seqs)]
  codon_sums[, `:=`(as_prob_normalized = score / gene_sum / N_AA_of_type_per_gene )]
  codon_sums[, `:=`(as_prob_normalized = as_prob_normalized / sum(as_prob_normalized)),
             by = .(variable, genes)]
  codon_sums <- codon_sums[, .(score = sum(score), as_prob_normalized = sum(as_prob_normalized),
                               N_total = .N), by = .(variable, genes, seqs)]
  seq_scores <- codon_sums[, .(sum = sum(score, na.rm = TRUE),
                               sum_txNorm = sum(as_prob_normalized, na.rm = TRUE),
                               var = var(score, na.rm = TRUE),
                               var_txNorm = var(as_prob_normalized, na.rm = TRUE),
                               N = .N, N.total = sum(N_total)), by = .(variable, seqs)]

  seq_scores[, mean := sum / N]
  seq_scores[!is.finite(mean), mean := 0]
  seq_scores[, mean_txNorm := sum_txNorm / N]
  seq_scores[!is.finite(mean_txNorm), mean_txNorm := 0]
  if (dispersion_method == "MLE") {
    stop("MLE dispersion not implemented yet!")
  } else if (dispersion_method == "MME") {
    seq_scores[, dispersion := (mean^2) / (var - mean)]
    seq_scores[!is.finite(dispersion) | dispersion < 0, dispersion := 0.001]
    seq_scores[, dispersion_txNorm := (mean_txNorm^2) / (var_txNorm - mean_txNorm)]
    seq_scores[!is.finite(dispersion_txNorm) | dispersion_txNorm < 0, dispersion_txNorm := 0.001]
  }
  seq_scores[, mean_percentage := (mean / sum(mean))*100, by = variable]
  seq_scores[, mean_txNorm_prob := (mean_txNorm / sum(mean_txNorm)), by = variable]
  seq_scores[, mean_txNorm_percentage := mean_txNorm_prob*100]
  seq_scores[, alpha := dirichlet_params(mean_txNorm_prob, sqrt(var_txNorm)),
             by = variable]
  seq_scores[, relative_to_max_score := (mean_txNorm_percentage /
                  max(mean_txNorm_percentage)*100), by = variable]
  if (is.null(seqs.order.table)) {
    setorderv(seq_scores, c("variable","mean_txNorm_percentage"), order = c(1,-1))
  } else {
    seq_scores[, seqs := factor(seqs, levels = seqs.order.table, ordered = TRUE)]
    setorderv(seq_scores, c("variable","seqs"), order = c(1,1))
  }

  seq_scores[]
  return(seq_scores)
}

# dirichlet MOM estimater
dirichlet_params <- function(p.mean, sigma) {
  n.params <- length(p.mean)
  if(n.params != length(sigma)){
    stop("Length of mean different from length of sigma")
  }
  # Compute second moment
  p.2 <- sigma^2 + p.mean^2
  # Initialize alpa vector
  alpha <- numeric(n.params)
  for (i in 1:(n.params-1)){
    alpha[i] <- (p.mean[1] - p.2[1])*p.mean[i]/(p.2[1] - p.mean[1]^2)
  }
  alpha[n.params] <- (p.mean[1] - p.2[1])*(1-sum(p.mean[-n.params]))/(p.2[1] - p.mean[1]^2)
  return(alpha)
}

codon_usage <- function(reads, cds,
                        mrna, faFile,
                        filter_table, filter_cds_mod3 = TRUE,
                        min_counts_cds_filter = 1e4,
                        with_A_sites = TRUE, code = GENETIC_CODE) {
  stopifnot(is(filter_table, "matrix"))
  stopifnot(length(cds) == nrow(filter_table))
  stopifnot(with_A_sites == TRUE)
  stopifnot(filter_cds_mod3 == TRUE)
  message("-- Codon usage analysis")
  message("- Starting with ", length(cds), " CDS sequences")
  # Filter CDSs
  if (filter_cds_mod3) {
    is_mod3 <- widthPerGroup(cds) %% 3 == 0
    cds_filtered <- cds[is_mod3]
    message("- ", length(cds_filtered), " CDSs are multiple of 3")
  }
  cds_filtered <- cds_filtered[filter_table[,1][names(cds_filtered)] > min_counts_cds_filter]
  message("- ", length(cds_filtered), " after filtering by counts filter")
  if (length(cds_filtered) == 0) stop("Filter is too strict, set a lower filter!")
  # Create merged A and P site ranges
  cds_with_A_site <- windowPerGroup(startSites(cds_filtered, T, T, T), mrna,
                                    downstream = widthPerGroup(cds_filtered, F) - 1, upstream = 3)
  # Get coverage
  if (is.list(reads)) {
    dt_samp <- lapply(reads, function(lib) {
      coveragePerTiling(cds_with_A_site, reads = lib,
                        is.sorted = TRUE, as.data.table = T)[,1][[1]]
    })
    dt_samp <- setDT(dt_samp)
  } else {
    dt_samp <- coveragePerTiling(cds_with_A_site, reads = reads,
                                 is.sorted = TRUE, as.data.table = T)[,1]
    colnames(dt_samp) <- "score"
  }
  ## Make P and A site tables
  cds_lengths <- widthPerGroup(cds_filtered, FALSE)
  cds_with_A_lengths <- cds_lengths + 3
  stopifnot(sum(cds_with_A_lengths) == nrow(dt_samp))
  gene_split_sites_general <- cumsum(cds_with_A_lengths)
  # P sites
  gene_split_sites <- c(0, gene_split_sites_general[-length(gene_split_sites_general)])
  a_start_codons <- sort(unlist(lapply(c(1,2,3), function(x) x + gene_split_sites)))
  dt_all_real <- dt_samp[-a_start_codons,]
  stopifnot(nrow(dt_all_real) == sum(cds_lengths))
  # A sites
  gene_split_sites <- gene_split_sites_general
  p_stop_codons <- sort(unlist(lapply(c(1,2,3), function(x) gene_split_sites - x)))
  dt_all_real_A <- dt_samp[-p_stop_codons,]
  stopifnot(nrow(dt_all_real_A) == sum(cds_lengths))

  # Get sequences and translate
  #browser()
  seqs <- orf_coding_table(cds_filtered, faFile, code)
  seqs[, merged := paste0(AA, ":", codon)]
  # Calc Codon usage
  genes_pos_index <- rep.int(seq_along(cds_filtered), times = cds_lengths)
  ifelse(filter_cds_mod3, "fast", "exact")
  codon_Psite <- seq_usage(dt_all_real, seqs$merged, genes_pos_index)
  codon_Asite <- seq_usage(dt_all_real_A, seqs$merged, genes_pos_index,
                seqs.order.table = unique(as.character(codon_Psite$seqs)))
  return(rbindlist(list(codon_Psite[, type := as.factor("P")],
                        codon_Asite[, type := as.factor("A")])))
}

GENETIC_CODE_ORFik <- function(code = Biostrings::GENETIC_CODE,
                               as.dt = FALSE, with.charge = FALSE,
                               init_as_hash = TRUE) {
  x <- Biostrings::GENETIC_CODE
  stops <- names(x)[x == "*"]
  x_stop <- x[(names(x) %in% stops)]
  start_codons <- if(init_as_hash) {c("###" = "#")} else NULL
  x <- c(start_codons, x[!(names(x) %in% stops)], x_stop)
  if (as.dt) x <- data.table(AA = x, codon = names(x))
  return(x)
}

orf_coding_table <- function(grl, faFile, code = GENETIC_CODE, as.factors = TRUE) {
  seqs_nt <- txSeqsFromFa(grl, faFile, is.sorted = TRUE)
  seqs_AA <- translate(seqs_nt, genetic.code = code)
  subseq(seqs_AA, 1, 1) <- "#"
  seqs_AA <- unlist(strsplit(as.character(unlist(seqs_AA, use.names = FALSE)),
                  split = ""))
  seqs_nt <- unlist(seqs_nt)
  seqs_codon <- substring(seqs_nt,
                          seq(1, nchar(seqs_nt)-2, by = 3),
                          seq(3, nchar(seqs_nt), by = 3))
  return(data.table(AA = as.factor(seqs_AA), codon = as.factor(seqs_codon)))
}

#' Codon analysis for ORFik experiment
#'
#' Per AA / codon, analyse the coverage, get a multitude of features.
#' For both A sites and P-sites (Input reads must be P-sites for now)
#' This function takes inspiration from the codonDT package, and among
#' others returns the negative binomial estimtes, but in addition many
#' other features.
#' @inheritParams outputLibs
#' @inheritParams txSeqsFromFa
#' @param reads either a single library (GRanges, GAlignment, GAlignmentPairs),
#' or a list of libraries returned from \code{outputLibs(df, )} with p-sites.
#' @param cds a GRangesList, the coding sequences, default:
#' \code{loadRegion(df, "cds", filterTranscripts(df))}, longest isoform
#' per gene.
#' @param mrna a GRangesList, the full mRNA sequences (matching by names
#' the cds sequences), default:
#' \code{loadRegion(df, "mrna", names(cds))}.
#' @param filter_cds_mod3 logical, default TRUE. Remove all ORFs that are not
#' mod3, this speeds up the computation a lot, and usually removes malformed
#' ORFs you would not want anyway.
#' @param filter_table an numeric(integer) matrix, where rownames are the
#' names of the full set of mRNA transcripts. This will be subsetted to the
#' cds subset you use. Then CDSs are filtered from this table by the
#' 'min_counts_cds_filter' argument.
#' @param min_counts_cds_filter numeric, default 1e4. Minimum number of counts
#' from the 'filter_table'  argument.
#' @param with_A_sites logical, default TRUE. Not used yet, will also return
#' A site scores.
#' @return a data.table of rows per codon / AA. The columns are:
#' seq: name of codon / AA
#' And many features returned.
#' @importFrom data.table frollsum setorderv
#' @family codon
#' @references https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7196831/
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
#' ## For single library
#' res <- codon_usage_exp(df, fimport(filepath(df[1,], "pshifted")),
#'                  min_counts_cds_filter = 10)
#' ## For multiple libs
#' res2 <- codon_usage_exp(df, outputLibs(df, type = "pshifted", output.mode = "list"),
#'                  min_counts_cds_filter = 10)
codon_usage_exp <- function(df, reads, cds = loadRegion(df, "cds", filterTranscripts(df)),
                            mrna = loadRegion(df, "mrna", names(cds)),
                            filter_cds_mod3 = TRUE, filter_table = assay(countTable(df, type = "summarized")[names(cds)]),
                            faFile = df@fafile,
                            min_counts_cds_filter = max(min(quantile(filter_table, 0.50), 100), 100),
                            with_A_sites = TRUE, code = GENETIC_CODE) {
  codon_usage(reads = reads, cds = cds,
              mrna = mrna,
              filter_cds_mod3 = filter_cds_mod3, filter_table = filter_table,
              faFile = faFile, min_counts_cds_filter = min_counts_cds_filter,
              with_A_sites = with_A_sites)
}

#' Plot codon_usage
#'
#' @param res a data.table of output from a codon_usage function
#' @param score_column numeric, default: res$relative_to_max_score.
#'  Which parameter to use as score column.
#' @param ylab character vector, names for libraries to show on Y axis
#' @param legend.position character, default "none", do not display legend.
#' @param coord_flip logical, default TRUE. Flip so that seqs are y coordinates
#' @return a ggplot object
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
#' ## For multiple libs
#' res2 <- codon_usage_exp(df, outputLibs(df, type = "pshifted", output.mode = "list"),
#'                  min_counts_cds_filter = 10)
#' codon_usage_plot(res2)
codon_usage_plot <- function(res, score_column = res$relative_to_max_score,
                             ylab = "Ribo-seq library",
                             legend.position = "none", coord_flip = TRUE) {
  if (is.null(score_column)) stop("score_column can not be NULL!")
  plot <- ggplot(res, aes(seqs, type, fill = score_column)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "orange", mid = "white",
                         midpoint = 0, limit = c(0,100), space = "Lab",
                         name="Rel. Cov.") +
    theme_classic()+ # minimal theme
    theme(legend.position=legend.position) +
    ylab(ylab)
  if (length(unique(res$variable)) > 1)
    plot <- plot + facet_wrap(res$variable, ncol = 4)
  if (coord_flip) plot + coord_flip()
}

# library(ORFik); library(data.table); library(ggplot2)
# df <- read.experiment("human_all_merged_l50")
# reads <- fimport(filepath(df, "cov"))
# res <- codon_usage_exp(df, reads)
#
