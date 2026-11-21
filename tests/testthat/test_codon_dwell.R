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
