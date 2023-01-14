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
    codon_sums[, variable := "lib1"]
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

codon_usage_exp <- function(df, reads, cds = loadRegion(df, "cds", filterTranscripts(df)),
                            mrna = loadRegion(df, "mrna", names(cds)),
                            filter_cds_mod3 = TRUE, filter_table = assay(countTable(df, type = "summarized")[names(cds)]),
                            fa = df@fafile, min_counts_cds_filter = 1e4,
                            with_A_sites = TRUE) {
  stopifnot(is(filter_table, "matrix"))
  stopifnot(length(cds) == nrow(filter_table))
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
  # Create merged A and P site ranges
  cds_with_A_site <- windowPerGroup(startSites(cds_filtered, T, T, T), mrna,
                downstream = widthPerGroup(cds_filtered, F) - 1, upstream = 3)
  # Get coverage
  if (is.list(reads)) {
    dt_samp <- lapply(reads, function(lib) {
      coveragePerTiling(cds_with_A_site, reads = lib,
                                   is.sorted = TRUE, as.data.table = T)[,1]
    })
    dt_samp <- setDT(dt_samp)
  } else {
    dt_samp <- coveragePerTiling(cds_with_A_site, reads = reads,
                                 is.sorted = TRUE, as.data.table = T)[,1]
  }
  colnames(dt_samp) <- "score"
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
  seqs <- translate(txSeqsFromFa(cds_filtered, df))
  seqs <- unlist(strsplit(as.character(unlist(seqs, use.names = FALSE)),
                          split = ""))
  # Calc Codon usage
  genes_pos_index <- rep.int(seq_along(cds_filtered), times = cds_lengths)
  ifelse(filter_cds_mod3, "fast", "exact")
  codon <- seq_usage(dt_all_real, seqs, genes_pos_index)
  codon_Asite <- seq_usage(dt_all_real_A, seqs, genes_pos_index,
                                seqs.order.table = as.character(codon$seqs))
  return(rbindlist(list(codon[, type := as.factor("P")],
                        codon_Asite[, type := as.factor("A")])))
}

# library(ORFik); library(data.table); library(ggplot2)
# df <- read.experiment("human_all_merged_l50")
# reads <- fimport(filepath(df, "cov"))
# res <- codon_usage_exp(df, reads)
#
# ggplot(res, aes(seqs, type, fill = relative_to_max_score)) +
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "orange", mid = "white",
#                        midpoint = 0, limit = c(0,100), space = "Lab",
#                        name="Rel. Cov.") +
#   theme_classic()+ # minimal theme
#   theme(legend.position="none") +
#   ylab("Ribo-seq library");

