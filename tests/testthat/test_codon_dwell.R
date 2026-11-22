context("codon dwell")
library(ORFik)

# Simulated genome: random CDSs on one chromosome (both strands, some with two exons), P-site counts drawn from
# E[y_gi] = mu_g * exp(sum_k beta_k[codon at offset k]) with known betas.
sense <- names(GENETIC_CODE)[GENETIC_CODE != "*"]
stops <- names(GENETIC_CODE)[GENETIC_CODE == "*"]

sim_genome <- function(n_genes = 300, n_codons = 200, seed = 1, extra = TRUE) {
  set.seed(seed)
  cds_seqs <- vapply(seq_len(n_genes), function(i)
    paste0("ATG", paste(sample(sense, n_codons - 2, TRUE), collapse = ""),
           sample(stops, 1)), character(1))
  if (extra) { # edge cases, all must be dropped by codon_dwell_design
    cds_seqs <- c(cds_seqs,
                  short = paste0("ATG", paste(sample(sense, 20, TRUE), collapse = ""), "TAA"),
                  not_mod3 = paste0("ATG", paste(sample(sense, 150, TRUE), collapse = ""), "TAAG"),
                  has_N = paste0("ATG", paste(sample(sense, 100, TRUE), collapse = ""), "NNN",
                                 paste(sample(sense, 100, TRUE), collapse = ""), "TAA"),
                  internal_stop = paste0("ATG", paste(sample(sense, 100, TRUE), collapse = ""), "TGA",
                                         paste(sample(sense, 100, TRUE), collapse = ""), "TAA"))
  }
  nms <- c(paste0("gene", seq_len(n_genes)), names(cds_seqs)[-seq_len(n_genes)])
  if (extra) nms[(n_genes + 1):length(nms)] <- c("short", "not_mod3", "has_N", "internal_stop")
  minus <- seq_along(cds_seqs) %% 2 == 0
  two_exons <- seq_along(cds_seqs) %% 3 == 0
  spacer <- 50; intron <- 40
  chr <- character(); pos <- 1; grl <- list()
  for (i in seq_along(cds_seqs)) {
    s <- cds_seqs[i]
    if (minus[i]) s <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
    L <- nchar(s)
    chr <- c(chr, paste(sample(c("A", "C", "G", "T"), spacer, TRUE), collapse = ""))
    start <- pos + spacer
    if (two_exons[i]) {
      cut <- 3 * 17 + 1
      chr <- c(chr, substr(s, 1, cut), paste(sample(c("A", "C", "G", "T"), intron, TRUE), collapse = ""),
               substr(s, cut + 1, L))
      ex <- IRanges(c(start, start + cut + intron), c(start + cut - 1, start + intron + L - 1))
      pos <- start + intron + L
    } else {
      chr <- c(chr, s)
      ex <- IRanges(start, width = L)
      pos <- start + L
    }
    st <- if (minus[i]) "-" else "+"
    gr <- GRanges("chrS", ex, st)
    if (minus[i]) gr <- rev(gr)
    grl[[nms[i]]] <- gr
  }
  genome <- Biostrings::DNAStringSet(c(chrS = paste0(paste(chr, collapse = ""),
                                                     paste(rep("A", spacer), collapse = ""))))
  fa <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(genome, fa)
  Rsamtools::indexFa(fa)
  list(cds = GRangesList(grl), fa = fa, seqs = cds_seqs, n_genes = n_genes)
}

# Draw P-site counts with known betas; returns GRanges with score (one range per codon with reads, frame 0)
sim_reads <- function(g, offsets, beta, mean_per_codon = 8, seed = 2, trim = c(20, 10)) {
  set.seed(seed)
  out <- list()
  for (i in seq_len(g$n_genes)) {
    s <- g$seqs[i]; n <- nchar(s) / 3
    cod <- substring(s, seq(1, 3 * n - 2, 3), seq(3, 3 * n, 3))
    p <- (trim[1] + 1):(n - trim[2]) # P-site codons simulated (as in the model)
    eta <- numeric(length(p))
    for (k in seq_along(offsets)) eta <- eta + beta[[k]][match(cod[p + 1 + offsets[k]], sense)]
    mu_g <- mean_per_codon * exp(rnorm(1, 0, 0.7))
    y <- rpois(length(p), mu_g * exp(eta) / mean(exp(eta)))
    keep <- y > 0
    out[[i]] <- data.table(gene = i, tx = 3 * (p[keep] - 1) + 1, y = y[keep])
  }
  dt <- rbindlist(out)
  gr <- pmapFromTranscriptF(IRanges(dt$tx, width = 1), g$cds[dt$gene], removeEmpty = TRUE)
  gr <- unlist(gr, use.names = FALSE)
  stopifnot(length(gr) == nrow(dt))
  mcols(gr) <- NULL
  gr$score <- dt$y
  sort(gr)
}

g <- sim_genome()
offsets <- -2:0
set.seed(10)
true_beta <- lapply(offsets, function(k) { b <- rnorm(61, 0, c(0.15, 0.25, 0.4)[k + 3]); b - mean(b) })
reads <- sim_reads(g, offsets, true_beta)
design <- codon_dwell_design(g$cds, g$fa, offsets = offsets, trim_start = 20, trim_end = 10)

test_that("codon_dwell_design filters CDSs and builds codes", {
  expect_is(design, "codon_dwell_design")
  expect_equal(length(design$cds), g$n_genes)
  expect_equal(unname(design$dropped), c(1L, 1L, 1L, 1L))
  expect_false(any(c("short", "not_mod3", "has_N", "internal_stop") %in% names(design$cds)))
  expect_equal(nrow(design$codes), g$n_genes * (200 - 30))
  expect_equal(ncol(design$codes), length(offsets))
  # codes match the sequence: gene 2 (minus strand), first kept P-site codon (21), A site = codon 22
  s <- g$seqs[2]
  r <- design$row_start[2] + 1
  expect_equal(design$sense[as.integer(design$codes[r, "0"]) + 1], unname(substr(s, 3 * 21 + 1, 3 * 21 + 3)))
  expect_equal(design$sense[as.integer(design$codes[r, "-2"]) + 1], unname(substr(s, 3 * 19 + 1, 3 * 19 + 3)))
  expect_error(codon_dwell_design(g$cds, g$fa, offsets = -10:5, trim_start = 5), "trim_start")
  expect_error(codon_dwell_design(g$cds, g$fa, offsets = -2:6, trim_end = 5), "trim_end")
})

test_that("P-sites are mapped to the right codons on both strands", {
  cnt <- ORFik:::codon_dwell_counts(reads, design)
  expect_equal(cnt$cds_psites, sum(reads$score))
  expect_equal(sum(cnt$y), sum(reads$score)) # all simulated reads are inside the trimmed window
  # P-sites in frame 1 and 2 of a codon count for that codon
  r2 <- reads; start(r2) <- end(r2) <- ifelse(as.logical(strand(r2) == "+"), start(r2) + 1, start(r2) - 1)
  in_exon <- overlapsAny(r2, design$cds)
  cnt2 <- ORFik:::codon_dwell_counts(r2[in_exon], design)
  expect_equal(cnt2$cds_psites, sum(r2$score[in_exon]))
  expect_gt(cor(cnt$y, cnt2$y), 0.99)
})

test_that("fast P-site mapping (position index) equals findOverlaps + pmapToTranscriptF", {
  slow <- ORFik:::codon_dwell_counts_overlaps(reads, design)
  expect_identical(ORFik:::codon_dwell_counts(reads, design), slow)
  # all frames, reads outside CDSs, unknown chromosome, unstranded reads, GAlignments-like weights absent
  set.seed(7)
  r3 <- GRanges("chrS", IRanges(sample(max(end(unlist(g$cds))) + 100, 20000, TRUE), width = 1),
                sample(c("+", "-"), 20000, TRUE), score = rpois(20000, 3) + 1L)
  expect_identical(ORFik:::codon_dwell_counts(r3, design),
                   ORFik:::codon_dwell_counts_overlaps(r3, design))
  # extra chromosomes in the reads (seqlevels matched by name)
  r3b <- suppressWarnings(c(GRanges("chrOther", IRanges(1:5, width = 1), "+", score = 1L), r3))
  expect_identical(ORFik:::codon_dwell_counts(r3b, design), ORFik:::codon_dwell_counts(r3, design))
  # unstranded reads count on both strands (as findOverlaps)
  ru <- GRanges("chrS", IRanges(start(reads)[1:50], width = 1), "*", score = 2L)
  both <- c(ru, ru); strand(both) <- rep(c("+", "-"), each = 50)
  expect_identical(ORFik:::codon_dwell_counts(ru, design), ORFik:::codon_dwell_counts(both, design))
  r4 <- r3; mcols(r4) <- NULL
  expect_identical(ORFik:::codon_dwell_counts(r4, design),
                   ORFik:::codon_dwell_counts_overlaps(r4, design))
  # overlapping CDSs of different genes (same strand): a read counts once per CDS
  cds_ov <- c(g$cds[1:20], GRangesList(ov1 = g$cds[[1]], ov3 = shift(g$cds[[3]], 30)))
  d_ov <- codon_dwell_design(cds_ov, g$fa, offsets = offsets, trim_start = 20, trim_end = 10)
  expect_gt(length(d_ov$index$layer_end), 1)
  expect_identical(ORFik:::codon_dwell_counts(r3, d_ov), ORFik:::codon_dwell_counts_overlaps(r3, d_ov))
  # designs saved without index / nt are completed
  d_old <- design; d_old$index <- NULL; d_old$nt <- NULL; d_old$nt_start <- NULL
  d_new <- ORFik:::codon_dwell_prepare(d_old, need_nt = TRUE)
  expect_identical(d_new$nt, design$nt); expect_identical(d_new$nt_start, design$nt_start)
  expect_identical(ORFik:::codon_dwell_counts(reads, d_old), slow)
  expect_equal(length(design$nt), 3 * sum(design$ncodons))
  expect_equal(as.integer(design$nt[1:3]), match(strsplit(substr(g$seqs[1], 1, 3), "")[[1]], c("A", "C", "G", "T")) - 1L)
})

test_that("codon_dwell glm recovers simulated dwell betas", {
  res <- codon_dwell(reads = reads, design = design, min_reads_gene = 100, verbose = FALSE)
  expect_true(all(c("library", "offset", "site", "codon", "aa", "beta", "se") %in% colnames(res)))
  expect_equal(nrow(res), 61 * length(offsets))
  s <- attr(res, "summary")
  expect_true(s$converged)
  expect_equal(s$status, "ok")
  truth <- unlist(true_beta)
  expect_gt(cor(res$beta, truth), 0.97)
  z <- (res$beta - truth) / res$se
  expect_lt(max(abs(z)), 5)
  expect_gt(stats::sd(z), 0.6); expect_lt(stats::sd(z), 1.6) # sandwich SEs roughly calibrated
  expect_equal(res[offset == 0, sum(beta)], 0, tolerance = 1e-8) # sum-to-zero
  expect_true(all(res[offset == 0]$site == "A"))
  expect_equal(res[codon == "AAA", unique(aa)], "K")
  expect_gt(s$deviance_explained, 0)
})

test_that("held-out deviance and design reuse work", {
  res <- codon_dwell(reads = list(a = reads, b = reads), design = design,
                     heldout = TRUE, verbose = FALSE)
  s <- attr(res, "summary")
  expect_equal(s$library, c("a", "b"))
  expect_gt(s$heldout_deviance_explained[1], 0)
  expect_equal(res[library == "a"]$beta, res[library == "b"]$beta)
  # building the design inside codon_dwell gives the same result
  res2 <- codon_dwell(reads = reads, cds = g$cds, faFile = g$fa, offsets = offsets, verbose = FALSE)
  expect_equal(res2$beta, res[library == "a"]$beta, tolerance = 1e-8)
})

test_that("profiled fit equals a Poisson GLM with gene dummies", {
  gs <- sim_genome(n_genes = 25, n_codons = 120, seed = 3, extra = FALSE)
  set.seed(4)
  b <- list(rnorm(61, 0, 0.4)); b[[1]] <- b[[1]] - mean(b[[1]])
  rs <- sim_reads(gs, 0L, b, mean_per_codon = 15, seed = 5, trim = c(3, 3))
  d <- codon_dwell_design(gs$cds, gs$fa, offsets = 0L, trim_start = 3, trim_end = 3, min_cds_codons = 50)
  res <- codon_dwell(reads = rs, design = d, min_reads_gene = 1, verbose = FALSE)
  cnt <- ORFik:::codon_dwell_counts(rs, d)
  dat <- data.table(y = cnt$y, gene = factor(d$gene),
                    codon = factor(d$sense[as.integer(d$codes[, 1]) + 1], levels = d$sense))
  fit <- glm(y ~ 0 + gene + codon, family = poisson(), data = dat)
  cf <- coef(fit)[paste0("codon", d$sense[-1])]
  ref <- c(0, cf); names(ref) <- d$sense
  ok <- !is.na(res$beta)
  expect_gt(sum(ok), 55)
  ref <- ref[ok]; ref <- ref - mean(ref)
  expect_equal(unname(res$beta[ok]), unname(ref), tolerance = 1e-5)
})

test_that("ratio method (codon_usage, mean_txNorm) agrees with glm at P and A site", {
  glm_res <- codon_dwell(reads = reads, design = design, verbose = FALSE)
  rat <- codon_dwell(reads = reads, design = design, method = "ratio", verbose = FALSE)
  expect_equal(attr(rat, "summary")$method, "ratio")
  expect_true(all(is.na(rat[offset == -2]$beta)))
  ok <- !is.na(rat$beta)
  expect_equal(sum(ok), 2 * 61)
  expect_gt(cor(glm_res$beta[ok], rat$beta[ok]), 0.9)
  expect_gt(cor(rat$beta[ok], unlist(true_beta)[ok]), 0.9)
  expect_true(all(rat$se[ok] > 0))
  # same numbers as codon_usage() directly
  cnt <- ORFik:::codon_dwell_counts(reads, design)
  keep <- cnt$gene_total >= 100
  cu <- suppressMessages(codon_usage(reads, design$cds[keep],
    extendTrailers(extendLeaders(design$cds[keep], 6L), 6L), g$fa,
    filter_table = matrix(cnt$gene_total[keep], ncol = 1, dimnames = list(names(design$cds)[keep], "l")),
    min_counts_cds_filter = 0))
  cu <- cu[type == "A" & !grepl("^[#*]:", seqs)]
  b <- log(cu$mean_txNorm); b <- b - mean(b)
  expect_equal(rat[offset == 0][match(sub("^.*:", "", cu$seqs), codon)]$beta, b)
})

test_that("libraries without reads or without deep genes give NA rows", {
  empty <- GRanges()
  few <- reads[round(seq(1, length(reads), length.out = 30))]
  res <- codon_dwell(reads = list(none = empty, few = few), design = design, verbose = FALSE)
  s <- attr(res, "summary")
  expect_equal(s$n_genes, c(0L, 0L))
  expect_true(all(grepl("no gene", s$status)))
  expect_true(all(is.na(res$beta)))
  # reads on other chromosomes only
  other <- GRanges("chrOther", IRanges(1:10, width = 1), "+")
  res <- codon_dwell(reads = other, design = design, verbose = FALSE)
  expect_equal(attr(res, "summary")$cds_psites, 0)
})

test_that("codon_dwell works on the ORFik template experiment", {
  df <- ORFik.template.experiment()
  df <- df[df$libtype == "RFP", ][1:2, ]
  res <- codon_dwell(df, offsets = -2:0, trim_start = 5, trim_end = 5,
                     min_cds_codons = 30, min_reads_gene = 10, verbose = FALSE)
  s <- attr(res, "summary")
  expect_equal(s$library, bamVarName(df))
  expect_true(all(s$status == "ok"))
  expect_true(all(s$converged))
  expect_equal(s$n_genes, c(6L, 6L))
  expect_true(all(res[!is.na(beta), abs(sum(beta)) < 1e-8, by = .(library, offset)]$V1))
  rat <- codon_dwell(df, offsets = -2:0, trim_start = 5, trim_end = 5, min_cds_codons = 30,
                     min_reads_gene = 10, method = "ratio", verbose = FALSE)
  expect_gt(cor(res$beta, rat$beta, use = "complete.obs"), 0.2) # 6 short genes: only a sanity check
  # plots
  expect_is(codon_dwell_plot(res), "ggplot")
  expect_is(codon_dwell_plot(res, type = "site", library = s$library), "ggplot")
  expect_error(codon_dwell_plot(res, library = "missing"), "No rows")
})

# Footprint end bias: reads in categories (length 28-30, P-site frame 0-2, offset 12) with capture bias from the 5' and
# 3' end nucleotides; dwell only at E, P, A. Without correction the end bias leaks into offsets -5 and +3..+5.
sim_reads_end_bias <- function(g, offsets, beta, b5, b3, lens = 28:30, pi_f = c(0.7, 0.2, 0.1), off5 = -12,
                               mean_per_codon = 12, seed = 12, trim = c(20, 10)) {
  set.seed(seed)
  acgt <- c("A", "C", "G", "T")
  pi_c <- as.vector(outer(pi_f, c(0.3, 0.45, 0.25))) # frame fastest, then length
  cat_f <- rep(0:2, length(lens)); cat_L <- rep(lens, each = 3)
  out <- list()
  for (i in seq_len(g$n_genes)) {
    s <- g$seqs[i]; n <- nchar(s) / 3
    nts <- match(strsplit(s, "")[[1]], acgt)
    cod <- substring(s, seq(1, 3 * n - 2, 3), seq(3, 3 * n, 3))
    p <- (trim[1] + 1):(n - trim[2])
    eta <- numeric(length(p))
    for (k in seq_along(offsets)) eta <- eta + beta[[k]][match(cod[p + 1 + offsets[k]], sense)]
    mu_g <- mean_per_codon * exp(rnorm(1, 0, 0.5)) * exp(eta) / mean(exp(eta))
    tx <- outer(3 * (p - 1) + 1, cat_f, "+")
    t5 <- tx + off5; t3 <- sweep(t5, 2, cat_L - 1, "+")
    li <- matrix(match(cat_L, lens), length(p), length(cat_L), byrow = TRUE)
    bias <- b5[cbind(nts[t5], as.vector(li))] + b3[cbind(nts[t3], as.vector(li))]
    lam <- mu_g * matrix(pi_c, length(p), length(pi_c), byrow = TRUE) * exp(matrix(bias, length(p)))
    y <- matrix(rpois(length(lam), lam), length(p))
    nz <- which(y > 0)
    out[[i]] <- data.table(gene = i, tx = tx[nz], y = y[nz], size = cat_L[col(y)[nz]])
  }
  dt <- rbindlist(out)
  gr <- unlist(pmapFromTranscriptF(IRanges(dt$tx, width = 1), g$cds[dt$gene], removeEmpty = TRUE),
               use.names = FALSE)
  stopifnot(length(gr) == nrow(dt))
  mcols(gr) <- NULL
  gr$size <- dt$size; gr$score <- dt$y
  gr
}

test_that("end_bias recovers simulated footprint end bias and removes it from end offsets", {
  ge <- sim_genome(n_genes = 400, n_codons = 220, seed = 21, extra = FALSE)
  offs <- c(-5L, -2:0, 3:5)
  set.seed(22)
  tb <- lapply(offs, function(k) {
    if (k %in% -2:0) { b <- rnorm(61, 0, 0.3); b - mean(b) } else rep(0, 61) })
  b5 <- matrix(rnorm(12, 0, 0.4), 4); b5 <- sweep(b5, 2, colMeans(b5))
  b3 <- matrix(rnorm(12, 0, 0.4), 4); b3 <- sweep(b3, 2, colMeans(b3))
  re <- sim_reads_end_bias(ge, offs, tb, b5, b3)
  de <- codon_dwell_design(ge$cds, ge$fa, offsets = offs, trim_start = 20, trim_end = 10)
  sh <- data.frame(fraction = 28:30, offsets_start = -12L)
  raw <- codon_dwell(reads = re, design = de, verbose = FALSE)
  eb <- codon_dwell(reads = re, design = de, end_bias = 1L, shifts = sh, verbose = FALSE)
  s <- attr(eb, "summary")
  expect_equal(s$status, "ok")
  expect_equal(s$end_bias_categories, 9L)
  expect_equal(s$end_bias_fraction, 1)
  cf <- attr(eb, "end_bias")
  expect_equal(nrow(cf), 2 * 2 * 3 * 4) # ends x (flank, terminal) x lengths x nt
  truth <- c(as.vector(b5), as.vector(b3))
  expect_gt(cor(cf[position == 1]$beta, truth), 0.98)
  expect_lt(max(abs(cf[position == 1]$beta - truth)), 0.08)
  expect_lt(max(abs(cf[position == 0]$beta)), 0.08) # no flank bias simulated
  expect_true(all(cf$se > 0))
  # end offsets: spurious signal without correction, gone with it; E/P/A dwell kept
  ends <- raw$offset %in% c(-5, 3:5)
  expect_gt(stats::sd(raw$beta[ends]), 0.1)
  expect_lt(stats::sd(eb$beta[ends]), 0.06)
  epa <- eb$offset %in% -2:0
  expect_gt(cor(eb$beta[epa], unlist(tb[2:4])), 0.97)
  # default (TRUE): 3 nucleotides per end + flank; extra terms near 0, same codon estimates
  eb3 <- codon_dwell(reads = re, design = de, end_bias = TRUE, shifts = sh, verbose = FALSE)
  cf3 <- attr(eb3, "end_bias")
  expect_equal(sort(unique(cf3$position)), 0:3)
  expect_gt(cor(cf3[position == 1]$beta, truth), 0.98)
  expect_lt(max(abs(cf3[position != 1]$beta)), 0.1)
  expect_gt(cor(eb3$beta, eb$beta), 0.99)
  expect_error(codon_dwell(reads = re, design = de, end_bias = 7, shifts = sh, verbose = FALSE), "0..6")
  # same genes, same uncorrected totals
  expect_equal(s$n_genes, attr(raw, "summary")$n_genes)
  expect_equal(s$cds_psites, attr(raw, "summary")$cds_psites)
  # no read lengths -> error status, shifts by name
  re2 <- re; re2$size <- NULL
  r2 <- codon_dwell(reads = re2, design = de, end_bias = TRUE, shifts = sh, verbose = FALSE)
  expect_match(attr(r2, "summary")$status, "size")
  expect_true(all(is.na(r2$beta)))
  expect_identical(ORFik:::codon_dwell_shift_for(list(`/x/ofst/SRR1.ofst` = sh), "SRR1"), sh)
  expect_identical(ORFik:::codon_dwell_shift_for(list(a = sh), "b", "/p/a_pshifted.ofst"), sh)
  expect_error(ORFik:::codon_dwell_shift_for(list(a = sh), "b"), "No shifting table")
  expect_error(codon_dwell(reads = re, design = de, end_bias = TRUE, verbose = FALSE), "shifts")
  expect_error(codon_dwell(reads = re, design = de, end_bias = TRUE, shifts = sh, method = "ratio",
                           verbose = FALSE), "glm")
})
