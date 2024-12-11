seq_usage <- function(dt, seqs, genes, input.dt.length = 1, output.seq.length = 3,
                      seqs.order.table = NULL, dispersion_method = "MME", aligned_position = "left") {
  stopifnot(is(dt, "data.table"))
  if ("genes" %in% colnames(dt)) {
    stop("dt can not contain column called 'genes'",
         "it should be given seperatly in argument 'genes'!")
  }
  if (length(genes) != (nrow(dt))) stop("Length of dt and genes must match!")
  if (length(seqs) != ((nrow(dt)) / output.seq.length))
    stop("Size of 'seqs' must be nrow(dt) / output.seq.length")
  codon_sums <- copy(dt)
  n.samples <- ncol(dt)
  has_names <- !all(colnames(codon_sums) %in% c("score", "count"))
  if (n.samples > 1 | has_names) {
    codon_sums <- suppressWarnings(melt.data.table(codon_sums,
                                                   value.name = "score"))
  } else {
    stopifnot(is(codon_sums[[1]], "numeric") | is(codon_sums[[1]], "integer"))
    if (colnames(codon_sums) == "count") colnames(codon_sums) <- "score"
    codon_sums[, variable := "library 1"]
  }
  codon_sums[, genes := rep(genes, length.out = .N)]
  we_must_collapse <- input.dt.length != output.seq.length
  if (we_must_collapse) {
    codon_sums[, codon_sum := frollsum(x = score, n = output.seq.length, align = aligned_position)]
    # Keep only second (aligned == "center") or first (aligned == "left") per seq group
    keep <- rep(FALSE, output.seq.length)

    keep[ifelse(aligned_position == "left",1, ceiling(output.seq.length/2))] <- TRUE
    codon_sums <- codon_sums[rep(keep,
                                 length.out = .N),]
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

coverage_A_and_P <- function(cds_filtered, mrna, reads,
                             cds_lengths = widthPerGroup(cds_filtered, FALSE), aligned_position = "center") {
  # Create merged A and P site ranges
  A_site_extension <- ifelse(aligned_position == "center", 4,3)
  cds_with_A_site <- startRegion(cds_filtered, mrna,
                                 downstream = cds_lengths - A_site_extension + 2,
                                 upstream = A_site_extension)

  # Get coverage
  if (is.list(reads)) {
    dt_samp <- lapply(reads, function(lib) {
      coveragePerTiling(cds_with_A_site, reads = lib,
                        is.sorted = TRUE, as.data.table = TRUE)[,1][[1]]
    })
    dt_samp <- setDT(dt_samp)
  } else {
    dt_samp <- coveragePerTiling(cds_with_A_site, reads = reads,
                                 is.sorted = TRUE, as.data.table = TRUE)[,1]
    colnames(dt_samp) <- "score"
  }
  ## Make P and A site tables
  cds_with_A_lengths <- cds_lengths + 3
  stopifnot(sum(cds_with_A_lengths) == nrow(dt_samp))
  gene_split_sites_general <- cumsum(cds_with_A_lengths)
  # P sites
  gene_split_sites <- c(0, gene_split_sites_general[-length(gene_split_sites_general)])
  a_start_codons <- sort(unlist(lapply(c(1,2,3), function(x) x + gene_split_sites)))
  dt_all_real <- dt_samp[-a_start_codons,]
  if (nrow(dt_all_real) != sum(cds_lengths))
    stop("nrow(dt_all_real) != sum(cds_lengths), did you include 5' UTRs?")
  # A sites
  gene_split_sites <- gene_split_sites_general
  p_stop_codons <- sort(unlist(lapply(c(1,2,3), function(x) gene_split_sites - x)))
  dt_all_real_A <- dt_samp[-p_stop_codons,]
  if (nrow(dt_all_real_A) != sum(cds_lengths))
    stop("nrow(dt_all_real_A) != sum(cds_lengths), did you include 3' UTRs?")
  return(list(A = dt_all_real_A, P = dt_all_real))
}

# dirichlet MOM estimator
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

filter_CDS_by_counts <- function(cds, filter_table,
                                 min_counts_cds_filter = 1000,
                                 filter_cds_mod3 = TRUE,
                                 minimum_5p_flank = FALSE,
                                 mrna = NULL) {
  message("- Starting with ", length(cds), " CDS sequences")
  # Filter CDSs
  if (filter_cds_mod3) {
    is_mod3 <- widthPerGroup(cds) %% 3 == 0
    cds_filtered <- cds[is_mod3]
    message("- ", length(cds_filtered), " CDSs are multiple of 3")
  }
  if (minimum_5p_flank & !is.null(mrna)) {
    cds_with_flank <- startRegion(cds_filtered, mrna,
                                   upstream = minimum_5p_flank)
    flank_size <- widthPerGroup(cds_with_flank, FALSE)
    has_valid_flank <- flank_size == max(flank_size)
    cds_filtered <- cds_filtered[has_valid_flank]
    message("- ", length(cds_filtered), " after filtering 5' flank of size ", minimum_5p_flank)
  }
  average_cds_counts <- rowMeans(as.matrix(filter_table[names(cds_filtered), , drop = FALSE]))

  cds_filtered <- cds_filtered[average_cds_counts > min_counts_cds_filter]
  message("- ", length(cds_filtered), " after filtering by counts filter")
  if (length(cds_filtered) == 0) stop("Filter is too strict, set a lower filter!")
  return(cds_filtered)
}

#' @inherit codon_usage_exp
#' @title Codon usage
#' @param cds a GRangesList
#' @param mrna a GRangesList
#' @param faFile a FaFile from genome
#' @param filter_table a matrix / vector of length equal to cds
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
#'
#' ## For single library
#' reads <- fimport(filepath(df[1,], "pshifted"))
#' cds <- loadRegion(df, "cds", filterTranscripts(df))
#' mrna <- loadRegion(df, "mrna", names(cds))
#' filter_table <- assay(countTable(df, type = "summarized")[names(cds)])
#' faFile <- findFa(df)
#' res <- codon_usage(reads, cds, mrna, faFile = faFile,
#'              filter_table = filter_table, min_counts_cds_filter = 10)
codon_usage <- function(reads, cds,
                        mrna, faFile,
                        filter_table, filter_cds_mod3 = TRUE,
                        min_counts_cds_filter = max(min(quantile(filter_table, 0.50), 1000), 1000),
                        with_A_sites = TRUE, aligned_position = "center", code = GENETIC_CODE) {
  stopifnot(is(filter_table, "matrix"))
  stopifnot(length(cds) == nrow(filter_table))
  stopifnot(with_A_sites == TRUE)
  stopifnot(filter_cds_mod3 == TRUE)
  message("-- Codon usage analysis")
  cds_filtered <- filter_CDS_by_counts(cds, filter_table,
                       min_counts_cds_filter, filter_cds_mod3,
                       minimum_5p_flank = ifelse(aligned_position == "center", 4,3),
                       mrna = mrna)

  # Get coverage for A sites and P sites
  cds_lengths <- widthPerGroup(cds_filtered, FALSE)
  coverage_list <- coverage_A_and_P(cds_filtered, mrna, reads, cds_lengths, aligned_position = aligned_position)

  # Get sequences and translate
  seqs <- orf_coding_table(cds_filtered, faFile, code)
  seqs[, merged := paste0(AA, ":", codon)]
  # Calc Codon usage
  type <- NULL # avoid bioccheck error
  genes_pos_index <- rep.int(seq_along(cds_filtered), times = cds_lengths)
  codon_Psite <- seq_usage(coverage_list[["P"]], seqs$merged, genes_pos_index, aligned_position = aligned_position)
  codon_Asite <- seq_usage(coverage_list[["A"]], seqs$merged, genes_pos_index,
                seqs.order.table = unique(as.character(codon_Psite$seqs)), aligned_position = aligned_position)
  return(rbindlist(list(codon_Psite[, type := as.factor("P")],
                        codon_Asite[, type := as.factor("A")])))
}

orf_coding_table <- function(grl, faFile, code = GENETIC_CODE, as.factors = TRUE) {
  seqs_nt <- txSeqsFromFa(grl, faFile, is.sorted = TRUE)
  seqs_AA <- as.character(translate(seqs_nt, genetic.code = code))
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
#' This function takes inspiration from the codonDT paper, and among
#' others returns the negative binomial estimates, but in addition many
#' other features.
#'
#' The primary column to use is "mean_txNorm", this is the fair normalized
#' score.
#' @inheritParams outputLibs
#' @inheritParams txSeqsFromFa
#' @param reads either a single library (GRanges, GAlignment, GAlignmentPairs),
#' or a list of libraries returned from \code{outputLibs(df)} with p-sites.
#' If list, the list must have names coresponding to the library names.
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
#' @param min_counts_cds_filter numeric, default:
#'  \code{max(min(quantile(filter_table, 0.50), 100), 100)}.
#'   Minimum number of counts from the 'filter_table'  argument.
#' @param with_A_sites logical, default TRUE. Not used yet, will also return
#' A site scores.
#' @param code a named character vector of size 64. Default: GENETIC_CODE.
#' Change if organism does not use the standard code.
#' @param aligned_position what positions should be taken to calculate per-codon coverage.
#' By default: "center", meaning that positions -1,0,1 will be taken.
#' Alternative: "left", then positions 0,1,2 are taken.
#' @return a data.table of rows per codon / AA. All values are given
#' per library, per site (A or P), sorted by the mean_txNorm_percentage column
#' of the first library in the set, the columns are:
#' \itemize{
#'  \item{variable (character)}{Library name}
#'  \item{seq (character)}{Amino acid:codon}
#'  \item{sum (integer)}{total counts per seq}
#'  \item{sum_txNorm (integer)}{total counts per seq normalized per tx}
#'  \item{var (numeric)}{variance of total counts per seq}
#'  \item{N (integer)}{total number of codons of that type}
#'  \item{mean_txNorm (numeric)}{Default use output, the fair codon usage,
#'  normalized both for gene and genome level for codon and read counts}
#'  \item{...}
#'  \item{alpha (numeric)}{dirichlet alpha MOM estimator (imagine mean and
#'            variance of probability in 1 value, the lower the value,
#'            the higher the variance, mean is decided by the relative
#'            value between samples)}
#'  \item{sum_txNorm (integer)}{total counts per seq normalized per tx}
#'  \item{relative_to_max_score (integer)}{Percentage use of codon}
#'  \item{type (factor(character))}{Either "P" or "A"}
#' }
#'
#' @importFrom data.table frollsum setorderv
#' @importFrom Biostrings GENETIC_CODE translate subseq `subseq<-`
#' @family codon
#' @references https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7196831/
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
#' ## For single library
#' res <- codon_usage_exp(df, fimport(filepath(df[1,], "pshifted")),
#'                  min_counts_cds_filter = 10)
#' # mean_txNorm is adviced scoring column
#' # codon_usage_plot(res, res$mean_txNorm)
#' # Default for plot function is the percentage scaled version of mean_txNorm
#' # codon_usage_plot(res) # This gives check error
#' ## For multiple libs
#' res2 <- codon_usage_exp(df, outputLibs(df, type = "pshifted", output.mode = "list"),
#'                  min_counts_cds_filter = 10)
#' # codon_usage_plot(res2)
codon_usage_exp <- function(df, reads, cds = loadRegion(df, "cds", filterTranscripts(df)),
                            mrna = loadRegion(df, "mrna", names(cds)),
                            filter_cds_mod3 = TRUE, filter_table = assay(countTable(df, type = "summarized")[names(cds)]),
                            faFile = df@fafile,
                            min_counts_cds_filter = max(min(quantile(filter_table, 0.50), 1000), 1000),
                            with_A_sites = TRUE, code = GENETIC_CODE, aligned_position = "center") {
  codon_usage(reads = reads, cds = cds,
              mrna = mrna,
              filter_cds_mod3 = filter_cds_mod3, filter_table = filter_table,
              faFile = faFile, min_counts_cds_filter = min_counts_cds_filter,
              with_A_sites = with_A_sites, aligned_position = aligned_position)
}

#' Plot codon_usage
#'
#' @param res a data.table of output from a codon_usage function
#' @param score_column numeric, default: res$relative_to_max_score.
#'  Which parameter to use as score column.
#' @param ylab character vector, names for libraries to show on Y axis
#' @param legend.position character, default "none", do not display legend.
#' @param limit numeric, 2 values for plot color limits. Default:
#' c(0, max(score_column))
#' @param midpoint numeric, default: limit/2. midpoint of color limit.
#' @param monospace_font logical, default TRUE. Use monospace font, this
#' does not work on systems (require specific font packages), set to
#' FALSE if it crashes for you.
#' @return a ggplot object
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
#' ## For multiple libs
#' res2 <- codon_usage_exp(df, outputLibs(df, type = "pshifted", output.mode = "list"),
#'                  min_counts_cds_filter = 10)
#' # codon_usage_plot(res2, monospace_font = TRUE) # This gives check error
#' codon_usage_plot(res2, monospace_font = FALSE) # monospace font looks better
codon_usage_plot <- function(res, score_column = res$relative_to_max_score,
                             ylab = "Ribo-seq library",
                             legend.position = "none",
                             limit = c(0, max(score_column)),
                             midpoint = limit/2, monospace_font = TRUE) {
  if (is.null(score_column)) stop("score_column can not be NULL!")
  type <- NULL # avoid bioccheck error
  font <- element_text()
  if (monospace_font)  font <- element_text(family = "monospace")
  plot <- ggplot(res, aes(type, seqs, fill = score_column)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "orange", mid = "white",
                         midpoint = midpoint, limit = limit, space = "Lab",
                         name="Rel. Cov.") +
    theme_classic()+ # minimal theme
    theme(legend.position=legend.position,
          axis.text.y = font) +
    ylab(ylab)

  if (length(unique(res$variable)) > 1)
    plot <- plot + facet_wrap(res$variable, ncol = 4)
  return(plot)
}



# library(ORFik); library(data.table); library(ggplot2)
# df <- read.experiment("human_all_merged_l50")
# reads <- fimport(filepath(df, "cov"))
# res <- codon_usage_exp(df, reads)
#
