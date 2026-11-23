# data.table / ggplot2 columns used in this file (see globals.R)
utils::globalVariables(c("seqs", "type", "mean_txNorm", "var_txNorm", "lb", "N",
  "rel", "g", "d", "m", "s", "site", "offset", "aa",
  "codon", "beta", "se", "label", "off", "fill", "status", "seconds", "pos",
  "seqnames", "strand", "w", "st", "n_positions", "iterations", "converged",
  "deviance_explained", "heldout_deviance_explained", "frac", "ck",
  "end_bias_categories", "end_bias_fraction"))

#' Codon dwell time design (sequence features of coding sequences)
#'
#' Builds the sequence part of the codon dwell time model of
#' \code{\link{codon_dwell}}: for every codon position kept in the model
#' (the P-site codon of a footprint), the codon found at each offset
#' around the A site. The design depends only on the coding sequences, not on
#' any library, so build it once per organism / annotation and reuse it for
#' all libraries (pass it to \code{codon_dwell(design = )}, or save it with
#' \code{saveRDS}).
#'
#' Coding sequences are kept when their length is a multiple of 3, they have
#' at least \code{min_cds_codons} codons, contain only A, C, G and T, and have
#' no in-frame stop codon before the last codon. Positions kept per CDS are
#' the P-site codons \code{(trim_start + 1):(ncodons - trim_end)}; the
#' codon at offset k is the codon at \code{P-site codon + 1 + k}
#' (offset 0 = A site, -1 = P site, -2 = E site). The trimming must keep all
#' offsets inside the CDS and away from its last codon:
#' \code{trim_start >= -min(offsets) - 1} and
#' \code{trim_end >= max(offsets) + 2}.
#' @param cds a GRangesList of coding sequences, one per gene is advised
#' (e.g. \code{loadRegion(df, "cds", filterTranscripts(df, 0, 1, 0))}). Must
#' have names.
#' @param faFile a FaFile, path to fasta or ORFik experiment, see
#' \code{\link{findFa}}.
#' @param offsets integer vector, codon offsets relative to the A site
#' (0 = A, -1 = P, -2 = E, negative = towards the start codon / nascent
#' peptide, positive = downstream mRNA). Default \code{-10:5}.
#' @param trim_start integer, default 20. Number of codons removed from the
#' start of each CDS (initiation ramp).
#' @param trim_end integer, default 10. Number of codons removed from the end
#' of each CDS (termination region).
#' @param min_cds_codons integer, default 100. Minimum CDS length in codons.
#' @param code a named character vector of size 64, the genetic code.
#' Default: \code{GENETIC_CODE}. Sense codons are all codons not coding "*".
#' @return a list of class "codon_dwell_design" with elements: codes (raw
#' matrix, positions x offsets, value = index of the sense codon minus one,
#' so 0 is the reference codon, the first sense codon), gene (gene index per
#' position), row_start (first position per gene, 0-based), ncodons (codons per
#' CDS), cds (the kept CDS GRangesList), faFile (fasta path), offsets, sense (sense codons), aa
#' (amino acid per sense codon), trim_start, trim_end, dropped (number of
#' CDSs removed per reason), nt (raw, the CDS nucleotides as codes 0..3 for
#' A, C, G, T, concatenated over CDSs), nt_start (0-based start of each CDS in
#' nt) and index (genome position -> CDS lookup used to map P-sites, see
#' \code{codon_dwell_index}). Designs saved by older versions without nt or
#' index still work: the missing parts are rebuilt when needed.
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' cds <- loadRegion(df, "cds", filterTranscripts(df, 0, 1, 0))
#' design <- codon_dwell_design(cds, findFa(df), offsets = -2:0,
#'                              trim_start = 5, trim_end = 5, min_cds_codons = 30)
#' dim(design$codes)
codon_dwell_design <- function(cds, faFile, offsets = -10:5, trim_start = 20L,
                               trim_end = 10L, min_cds_codons = 100L,
                               code = GENETIC_CODE) {
  stopifnot(is(cds, "GRangesList"), !is.null(names(cds)))
  offsets <- as.integer(offsets)
  stopifnot(length(offsets) > 0, !anyDuplicated(offsets), !anyNA(offsets))
  trim_start <- as.integer(trim_start); trim_end <- as.integer(trim_end)
  if (trim_start < -min(offsets) - 1L)
    stop("trim_start must be >= -min(offsets) - 1 (", -min(offsets) - 1L,
         "), so that all offsets stay inside the CDS")
  if (trim_end < max(offsets) + 2L)
    stop("trim_end must be >= max(offsets) + 2 (", max(offsets) + 2L,
         "), so that offsets never reach the stop codon")
  stopifnot(length(code) == 64, !is.null(names(code)))
  sense <- names(code)[code != "*"]
  stops <- names(code)[code == "*"]
  if (length(sense) > 255) stop("Too many sense codons")

  dropped <- c(not_mod3 = 0L, too_short = 0L, non_ACGT = 0L, internal_stop = 0L)
  cds <- cds[lengths(cds) > 0]
  len <- widthPerGroup(cds, FALSE)
  is_mod3 <- len %% 3L == 0L
  dropped["not_mod3"] <- sum(!is_mod3)
  ncod <- len %/% 3L
  long_enough <- ncod >= max(min_cds_codons, trim_start + trim_end + 1L)
  dropped["too_short"] <- sum(is_mod3 & !long_enough)
  cds <- cds[is_mod3 & long_enough]
  if (length(cds) == 0) stop("No CDS left after length filters")
  cds <- sortPerGroup(cds)
  seqs <- txSeqsFromFa(cds, faFile, is.sorted = TRUE)
  # Codon index (1..64, into names(code)) per codon, via nucleotide codes
  nt <- as.integer(unlist(seqs, use.names = FALSE))
  nt_map <- rep(NA_integer_, max(nt, 16L))
  nt_map[as.integer(Biostrings::DNAString("ACGT"))] <- 0:3
  ni <- nt_map[nt]
  ncod <- widthPerGroup(cds, FALSE) %/% 3L
  n_all <- length(ni) %/% 3L
  idx3 <- 3L * seq_len(n_all)
  acgt <- c("A", "C", "G", "T")
  all_codons <- paste0(rep(acgt, each = 16), rep(rep(acgt, each = 4), 4), rep(acgt, 16))
  cid <- match(all_codons, names(code))[16L * ni[idx3 - 2L] + 4L * ni[idx3 - 1L] + ni[idx3] + 1L]
  gene_of_codon <- rep.int(seq_along(cds), ncod)
  bad_nt <- unique(gene_of_codon[is.na(cid)])
  is_last <- sequence(ncod) == rep.int(ncod, ncod)
  bad_stop <- unique(gene_of_codon[!is_last & !is.na(cid) & names(code)[cid] %in% stops])
  bad_stop <- setdiff(bad_stop, bad_nt)
  dropped["non_ACGT"] <- length(bad_nt); dropped["internal_stop"] <- length(bad_stop)
  keep <- !(seq_along(cds) %in% c(bad_nt, bad_stop))
  if (!any(keep)) stop("No CDS left after sequence filters")
  cid <- cid[keep[gene_of_codon]]
  nt_kept <- as.raw(ni[rep.int(keep, 3L * ncod)])
  cds <- cds[keep]; ncod <- ncod[keep]
  sidx <- match(names(code)[cid], sense) - 1L # 0-based sense codon index, NA = stop codon

  first <- c(0L, cumsum(ncod)[-length(ncod)])
  nrow_g <- ncod - trim_start - trim_end
  gene <- rep.int(seq_along(cds), nrow_g)
  rowcod <- sequence(nrow_g, from = trim_start + 1L)
  codes <- matrix(as.raw(0), nrow = length(gene), ncol = length(offsets))
  for (k in seq_along(offsets)) {
    v <- sidx[first[gene] + rowcod + 1L + offsets[k]]
    if (anyNA(v)) stop("Internal error: offset outside sense codons")
    codes[, k] <- as.raw(v)
  }
  colnames(codes) <- offsets
  structure(list(codes = codes, gene = gene,
                 nt = nt_kept, nt_start = c(0, cumsum(3 * as.numeric(ncod))[-length(ncod)]),
                 index = codon_dwell_index(cds),
                 row_start = c(0L, cumsum(nrow_g)[-length(nrow_g)]),
                 ncodons = ncod, cds = cds, offsets = offsets, sense = sense,
                 aa = unname(code[sense]), trim_start = trim_start,
                 faFile = as.character(BiocGenerics::path(findFa(faFile))),
                 trim_end = trim_end, dropped = dropped),
            class = "codon_dwell_design")
}

#' @export
print.codon_dwell_design <- function(x, ...) {
  cat("codon_dwell_design: ", length(x$cds), " CDSs, ", nrow(x$codes),
      " positions, offsets ", paste(range(x$offsets), collapse = ".."),
      ", trim ", x$trim_start, "/", x$trim_end, " codons\n", sep = "")
  invisible(x)
}

# Genome position -> CDS lookup of a design. The CDS exons are sorted by a numeric key
# (chromosome-strand id * 2^32 + position) and split into layers of non-overlapping exons (overlapping CDSs of
# different genes go to later layers), so P-sites are matched with findInterval per layer instead of findOverlaps,
# and the transcript coordinate comes from the exon's offset in its CDS (no pmapToTranscriptF).
codon_dwell_index <- function(cds) {
  ex <- unlist(cds, use.names = FALSE)
  g <- rep.int(seq_along(cds), lengths(cds))
  w <- as.numeric(width(ex))
  off <- cumsum(w) - w
  tx0 <- as.integer(off - off[!duplicated(g)][g]) # nt of the CDS before this exon (CDS exons sorted 5' -> 3')
  minus <- as.logical(strand(ex) == "-")
  cs <- (as.integer(seqnames(ex)) - 1) * 2 + minus
  ks <- cs * 2^32 + start(ex); ke <- cs * 2^32 + end(ex)
  o <- order(ks, ke)
  ks <- ks[o]; ke <- ke[o]
  layer <- integer(length(o)); rem <- seq_along(o); l <- 0L
  while (length(rem) > 0) {
    l <- l + 1L
    e <- ke[rem]
    pick <- ks[rem] > c(-Inf, cummax(e)[-length(e)])
    layer[rem[pick]] <- l
    rem <- rem[!pick]
  }
  o2 <- order(layer, ks)
  o <- o[o2]
  list(seqlevels = seqlevels(cds), key_start = ks[o2], key_end = ke[o2],
       start = start(ex)[o], end = end(ex)[o], minus = minus[o], gene = g[o], tx0 = tx0[o],
       layer_end = cumsum(tabulate(layer[o2], l)))
}

# Adds the parts of a design that older saved designs lack (index, nt).
codon_dwell_prepare <- function(design, need_nt = FALSE) {
  if (is.null(design$index)) design$index <- codon_dwell_index(design$cds)
  if (need_nt && is.null(design$nt)) {
    seqs <- txSeqsFromFa(design$cds, design$faFile, is.sorted = TRUE)
    nt <- as.integer(unlist(seqs, use.names = FALSE))
    nt_map <- rep(NA_integer_, max(nt, 16L))
    nt_map[as.integer(Biostrings::DNAString("ACGT"))] <- 0:3
    design$nt <- as.raw(nt_map[nt])
    design$nt_start <- c(0, cumsum(3 * as.numeric(design$ncodons))[-length(design$ncodons)])
  }
  design
}

#' Map P-sites to the CDSs of a codon dwell design
#'
#' @param reads P-site reads (GRanges, GAlignments, GAlignmentPairs or a path
#' readable by \code{fimport}). The 5' end of each read is used as P-site.
#' Weights are taken from a "score" metadata column if present (collapsed
#' reads, e.g. from ofst files).
#' @param design a codon_dwell_design
#' @return a list with gene (CDS index), tx (1-based P-site position in the
#' CDS), w (weight) and idx (index of the read), one entry per read and CDS it
#' falls in (a read in two overlapping CDSs counts for both, as with
#' findOverlaps).
#' @keywords internal
codon_dwell_map <- function(reads, design) {
  if (is.character(reads)) reads <- fimport(reads)
  w <- mcols(reads)$score
  if (is.null(w)) w <- rep(1, length(reads))
  if (!is(reads, "GRanges")) reads <- granges(reads)
  ix <- design$index
  if (is.null(ix)) ix <- codon_dwell_index(design$cds)
  empty <- list(gene = integer(), tx = integer(), w = numeric(), idx = integer())
  if (length(reads) == 0) return(empty)
  m <- codon_dwell_map_cpp(as.integer(seqnames(reads)), as.integer(strand(reads)),
                           start(reads), end(reads), match(seqlevels(reads), ix$seqlevels),
                           ix$key_start, ix$key_end, as.integer(ix$layer_end), ix$start, ix$end,
                           ix$minus, ix$gene, ix$tx0)
  list(gene = m$gene, tx = m$tx, w = as.numeric(w)[m$idx], idx = m$idx)
}

#' Count P-sites per codon position of a codon dwell design
#'
#' @inheritParams codon_dwell_map
#' @param map output of \code{codon_dwell_map} (default: computed from reads)
#' @param w weights per entry of map (default map$w), e.g. end bias corrected
#' weights for y.
#' @return a list with y (counts per design position), gene_total
#' (P-sites per CDS, all codons and frames, with the weights map$w) and
#' cds_psites (sum of gene_total).
#' @keywords internal
codon_dwell_counts <- function(reads, design, map = codon_dwell_map(reads, design),
                               w = map$w) {
  G <- length(design$cds)
  y <- numeric(nrow(design$codes)); gene_total <- numeric(G)
  if (length(map$gene) > 0) {
    cc <- codon_dwell_counts_cpp(map$gene, map$tx, as.numeric(map$w), as.numeric(w),
                                 as.integer(design$ncodons), as.numeric(design$row_start),
                                 design$trim_start, design$trim_end, nrow(design$codes))
    y <- cc$y; gene_total <- cc$gene_total
  }
  list(y = y, gene_total = gene_total, cds_psites = sum(gene_total))
}

# Reference implementation of codon_dwell_counts (findOverlaps + pmapToTranscriptF, the path of ORFik 1.x
# before the position index); kept for the tests.
codon_dwell_counts_overlaps <- function(reads, design) {
  if (is.character(reads)) reads <- fimport(reads)
  w <- mcols(reads)$score
  if (is.null(w)) w <- rep(1, length(reads))
  if (!is(reads, "GRanges")) reads <- granges(reads)
  G <- length(design$cds)
  y <- numeric(nrow(design$codes)); gene_total <- numeric(G)
  if (length(reads) > 0) {
    sn <- seqnames(reads); st <- strand(reads)
    minus <- as.logical(st == "-")
    dt <- data.table(s = as.integer(sn), pos = ifelse(minus, end(reads), start(reads)),
                     st = as.integer(st), w = as.numeric(w))
    dt <- dt[, .(w = sum(w)), by = .(s, pos, st)]
    p <- GRanges(factor(levels(sn), levels = levels(sn))[dt$s], IRanges(dt$pos, width = 1L),
                 factor(levels(st), levels = levels(st))[dt$st])
    suppressWarnings(seqlevels(p) <- unique(c(seqlevels(p), seqlevels(design$cds))))
    hits <- findOverlaps(p, design$cds)
    if (length(hits) > 0) {
      tx <- pmapToTranscriptF(p[queryHits(hits)], design$cds[subjectHits(hits)])
      g <- subjectHits(hits)
      cod <- (start(tx) - 1L) %/% 3L + 1L
      ww <- dt$w[queryHits(hits)]
      gene_total <- as.numeric(tabulate_weighted(g, ww, G))
      inside <- cod > design$trim_start & cod <= design$ncodons[g] - design$trim_end
      row <- design$row_start[g[inside]] + cod[inside] - design$trim_start
      y <- as.numeric(tabulate_weighted(row, ww[inside], nrow(design$codes)))
    }
  }
  list(y = y, gene_total = gene_total, cds_psites = sum(gene_total))
}

tabulate_weighted <- function(bin, w, n) {
  out <- numeric(n)
  if (length(bin) == 0) return(out)
  s <- rowsum(w, bin, reorder = FALSE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# Profiled-gene Poisson evaluation (see src/codon_dwell.cpp)
codon_dwell_eval <- function(design, sel, y, beta, hessian = TRUE, scores = FALSE,
                             threads = 1L) {
  codon_dwell_eval_cpp(design$codes, sel$rows, as.double(y), sel$starts,
                       as.double(beta), length(design$sense) - 1L, hessian, scores,
                       as.integer(threads))
}

codon_dwell_check_design <- function(design) {
  stopifnot(is(design, "codon_dwell_design"))
  # codes <= number of sense codons - 1 is checked in C++ (codon_dwell_feature_counts_cpp) before every fit
  stopifnot(length(design$gene) == nrow(design$codes),
            ncol(design$codes) == length(design$offsets))
  invisible(TRUE)
}

# Row selection of genes (logical per design row, grouped by gene)
codon_dwell_rows <- function(design, select) {
  idx <- which(select)
  g <- design$gene[idx]
  br <- if (length(g) > 1) which(g[-1L] != g[-length(g)]) else integer()
  list(rows = idx - 1L, starts = as.integer(c(0L, br, length(g))),
       n_genes = length(br) + as.integer(length(g) > 0))
}

codon_dwell_profiled_hessian <- function(e, active) {
  H <- e$hessian_raw[active, active, drop = FALSE] -
    tcrossprod(e$xp)[active, active, drop = FALSE] # no copy of the P x genes matrix
  diag(H) <- diag(H) + 1e-8
  H
}

codon_dwell_chol <- function(H) {
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R)) {
    diag(H) <- diag(H) + 1e-6 * max(1, mean(diag(H)))
    R <- chol(H)
  }
  R
}

# Fisher scoring on beta only, chord Newton once the step is small.
codon_dwell_fisher <- function(design, sel, y, active, beta = NULL, se = TRUE,
                               iters = 50L, tol = 1e-8, chord = 0.02,
                               max_step = 3, threads = 1L) {
  P <- length(design$offsets) * (length(design$sense) - 1L)
  if (is.null(beta)) beta <- numeric(P)
  Rc <- NULL; step <- Inf; it <- 0L
  for (it in seq_len(iters)) {
    e <- codon_dwell_eval(design, sel, y, beta, hessian = is.null(Rc), threads = threads)
    if (is.null(Rc)) Rc <- codon_dwell_chol(codon_dwell_profiled_hessian(e, active))
    step <- backsolve(Rc, forwardsolve(t(Rc), e$gradient[active]))
    if (max(abs(step)) > max_step) step <- step * (max_step / max(abs(step)))
    beta[active] <- beta[active] + step
    if (max(abs(step)) < tol) break
    if (max(abs(step)) >= chord) Rc <- NULL
  }
  out <- list(beta = beta, iterations = it, converged = max(abs(step)) < tol)
  e <- codon_dwell_eval(design, sel, y, beta, hessian = se, scores = se, threads = threads)
  out$deviance_explained <- 1 - e$deviance / e$null_deviance
  if (se) {
    Hi <- chol2inv(codon_dwell_chol(codon_dwell_profiled_hessian(e, active)))
    out$vcov <- Hi %*% tcrossprod(e$score)[active, active, drop = FALSE] %*% Hi
  }
  out
}

# Per offset: which codons are observed (positions with reads), the anchor codon
# (fixed at 0: the reference codon if observed), the estimated codons and start
# values for the fit (log observed / expected counts under uniform density within
# gene, relative to the anchor; one marginal pass that puts the first Newton step
# close to the optimum).
codon_dwell_features <- function(design, select, y, threads = 1L) {
  S <- length(design$sense); K <- length(design$offsets)
  rows <- which(select)
  yy <- as.double(y[select])
  g <- design$gene[rows]
  m_g <- tabulate(g, length(design$cds))
  n_g <- codon_dwell_tabulate_cpp(g, yy, length(design$cds))
  mu0 <- (n_g / pmax(m_g, 1))[g]
  fc <- codon_dwell_feature_counts_cpp(design$codes, rows - 1L, yy, S - 1L, mu0, threads)
  observed <- fc$n > 0 & fc$y > 0
  anchor <- vapply(seq_len(K), function(k) {
    if (observed[1, k]) 1L else if (any(observed[, k])) which.max(fc$n[, k] * observed[, k]) else 1L
  }, integer(1))
  estimated <- observed; estimated[cbind(anchor, seq_len(K))] <- FALSE
  # parameter c (codon index 2..S) of offset k is at (k - 1) * (S - 1) + c - 1
  active <- which(estimated[-1, , drop = FALSE])
  lr <- log(pmax(fc$y, 1e-300) / pmax(fc$y2, 1e-300))
  lr <- sweep(lr, 2, lr[cbind(anchor, seq_len(K))])
  start <- numeric((S - 1L) * K)
  start[active] <- pmax(pmin(as.vector(lr[-1, , drop = FALSE])[active], 3), -3)
  list(observed = observed, anchor = anchor, active = active, start = start)
}

# Sum-to-zero (over observed codons) betas and their SEs per offset
codon_dwell_sum_to_zero <- function(design, feat, beta, vcov = NULL) {
  S <- length(design$sense); K <- length(design$offsets)
  b_all <- se_all <- rep(NA_real_, S * K)
  for (k in seq_len(K)) {
    obs <- which(feat$observed[, k])
    if (length(obs) == 0) next
    pidx <- (k - 1L) * (S - 1L) + obs - 1L # parameter index; invalid for codon 1
    b <- ifelse(obs == 1L, 0, beta[pmax(pidx, 1L)])
    b[obs == feat$anchor[k]] <- 0
    m <- length(obs)
    b_all[(k - 1L) * S + obs] <- b - mean(b)
    if (!is.null(vcov)) {
      ai <- match(pidx, feat$active)            # position in active set (NA for anchor / ref)
      E <- matrix(0, m, length(feat$active))
      ok <- which(!is.na(ai) & obs != 1L & obs != feat$anchor[k])
      E[cbind(ok, ai[ok])] <- 1
      Tm <- (diag(m) - matrix(1 / m, m, m)) %*% E
      se_all[(k - 1L) * S + obs] <- sqrt(pmax(rowSums((Tm %*% vcov) * Tm), 0))
    }
  }
  list(beta = b_all, se = se_all)
}

codon_dwell_glm <- function(design, counts, keep_gene, heldout = FALSE, threads = 1L) {
  select <- keep_gene[design$gene]
  sel <- codon_dwell_rows(design, select)
  y <- counts$y[select]
  feat <- codon_dwell_features(design, select, counts$y, threads = threads)
  if (length(feat$active) == 0) stop("No codon features observed")
  fit <- codon_dwell_fisher(design, sel, y, feat$active, beta = feat$start, threads = threads)
  stz <- codon_dwell_sum_to_zero(design, feat, fit$beta, fit$vcov)
  summ <- list(iterations = fit$iterations, converged = fit$converged,
               deviance_explained = fit$deviance_explained,
               heldout_deviance_explained = NA_real_)
  if (heldout) {
    kept <- which(keep_gene)
    even_gene <- logical(length(keep_gene))
    even_gene[kept[seq_along(kept) %% 2L == 0L]] <- TRUE
    tr <- select & even_gene[design$gene]; te <- select & !even_gene[design$gene]
    if (any(tr) && any(te)) {
      ftr <- codon_dwell_fisher(design, codon_dwell_rows(design, tr), counts$y[tr],
                                feat$active, beta = fit$beta, se = FALSE, threads = threads)
      e <- codon_dwell_eval(design, codon_dwell_rows(design, te), counts$y[te],
                            ftr$beta, hessian = FALSE, threads = threads)
      summ$heldout_deviance_explained <- 1 - e$deviance / e$null_deviance
    }
  }
  list(beta = stz$beta, se = stz$se, summary = summ, n_positions = length(y))
}

# Ratio estimate: codon_usage() (as codon_usage_exp(), mean_txNorm) on the
# design CDSs and genes, P and A site; beta = log(mean_txNorm), centered over
# internal sense codons.
codon_dwell_ratio <- function(design, counts, keep_gene, reads, mrna, faFile) {
  S <- length(design$sense)
  cds <- design$cds[keep_gene]
  ft <- matrix(counts$gene_total[keep_gene], ncol = 1,
               dimnames = list(names(cds), "library"))
  mrna <- mrna[names(cds)]
  cu <- suppressMessages(codon_usage(reads, cds, mrna, faFile = faFile,
                                     filter_table = ft, min_counts_cds_filter = 0))
  cu <- cu[!grepl("^[#*]:", seqs)]
  cu[, codon := sub("^.*:", "", seqs)]
  cu[, offset := ifelse(type == "A", 0L, -1L)]
  cu <- cu[mean_txNorm > 0 & codon %in% design$sense]
  cu[, lb := log(mean_txNorm)]
  cu[, `:=`(beta = lb - mean(lb), se = sqrt(var_txNorm / N) / mean_txNorm), by = offset]
  k <- match(cu$offset, design$offsets)
  b_all <- se_all <- rep(NA_real_, S * length(design$offsets))
  ok <- !is.na(k)
  idx <- (k[ok] - 1L) * S + match(cu$codon[ok], design$sense)
  b_all[idx] <- cu$beta[ok]; se_all[idx] <- cu$se[ok]
  list(beta = b_all, se = se_all, n_positions = sum(widthPerGroup(cds, FALSE)) / 3,
       summary = list(iterations = NA_integer_, converged = NA,
                      deviance_explained = NA_real_,
                      heldout_deviance_explained = NA_real_))
}

# Centered (mean over the levels) values and SEs of one block of parameters: idx = parameter index per level (NA =
# fixed at 0), active = estimated parameters (columns of vcov).
codon_dwell_center <- function(beta, idx, active, vcov = NULL) {
  b <- ifelse(is.na(idx), 0, beta[pmax(idx, 1L)])
  m <- length(b)
  se <- rep(NA_real_, m)
  if (!is.null(vcov)) {
    ai <- match(idx, active)
    E <- matrix(0, m, length(active))
    ok <- which(!is.na(ai))
    E[cbind(ok, ai[ok])] <- 1
    Tm <- (diag(m) - matrix(1 / m, m, m)) %*% E
    se <- sqrt(pmax(rowSums((Tm %*% vcov) * Tm), 0))
  }
  list(beta = b - mean(b), se = se)
}

# Footprint end bias (choros-style technical terms, first stage of codon_dwell(end_bias = TRUE)).
#
# Reads at codon position i are split into categories c = (read length L, frame f of the P-site in the codon). With
# P-site offset o_L (negative, from the shifting table), the read's 5' end is at CDS nt 3 (i - 1) + 1 + f + o_L and its
# 3' end L - 1 nt further. Model: E[y_ic] = mu_g * d_i * pi_c * exp(b5_L[nt5_ic] + b3_L[nt3_ic]). Conditional on the
# reads at position i, the categories are multinomial with logits a_c + b5 + b3: the dwell d_i and gene flux mu_g cancel,
# so the end bias is estimated without any codon model, with the profiled kernel of codon_dwell (positions as groups,
# category intercepts + 5' and 3' end nucleotide per read length as features; nucleotide A is the reference per length).
# With end_nt > 1, the first end_nt nucleotides from each end get their own terms (position 1 = the terminal nt);
# with flank = TRUE also the nucleotide just outside each end (position 0: the nuclease cut site).
# Second stage: each read gets weight exp(-(b5 + b3)), rescaled so every category keeps its total; the codon
# fit then runs on the corrected counts (E[sum_c corrected y_ic] = mu_g * d_i * const).
codon_dwell_end_bias <- function(map, design, len, shift, end_nt = 1L, flank = TRUE, min_frac = 0.005,
                                 max_positions = 100000L, threads = 1L) {
  end_nt <- as.integer(end_nt)
  stopifnot(length(end_nt) == 1, end_nt >= 1, end_nt <= 6)
  g <- map$gene; tx <- map$tx; w <- map$w
  L <- as.integer(len[map$idx])
  cod <- (tx - 1L) %/% 3L + 1L; f <- (tx - 1L) %% 3L
  inside <- cod > design$trim_start & cod <= design$ncodons[g] - design$trim_end
  o5 <- as.integer(shift$offsets_start)[match(L, as.integer(shift$fraction))]
  cand <- inside & !is.na(o5)
  if (!any(cand)) stop("end_bias: no P-site with a read length in the shifting table")
  ckey <- L * 3L + f
  tab <- data.table(ck = ckey[cand], w = w[cand])[, .(w = sum(w)), by = ck]
  tab[, frac := w / sum(w)]
  tab <- tab[frac >= min_frac][order(-w)]
  len_order <- tab[, .(w = sum(w)), by = .(L = ck %/% 3L)][order(-w)]$L
  lens <- sort(head(len_order, 15L))
  tab <- head(tab[(ck %/% 3L) %in% lens], 61L)
  if (nrow(tab) < 2) stop("end_bias: fewer than 2 read length / frame categories")
  tab <- tab[order(ck)]
  C <- nrow(tab)
  cat_L <- tab$ck %/% 3L; cat_f <- tab$ck %% 3L
  cat_o <- as.integer(shift$offsets_start)[match(cat_L, as.integer(shift$fraction))]
  cat_l <- match(cat_L, lens) - 1L
  nL <- length(lens)
  ci <- match(ckey, tab$ck)
  use <- cand & !is.na(ci)
  # End positions of the used reads (1-based in the CDS), in bounds of the CDS
  tx5 <- tx + o5; tx3 <- tx5 + L - 1L
  J <- c(if (flank) -1L, seq_len(end_nt) - 1L) # nt positions from each end (-1 = flank, outside the read)
  nJ <- length(J)
  use <- use & tx5 - max(0L, -min(J)) >= 1L & tx3 + max(0L, -min(J)) <= 3L * design$ncodons[g]
  nt <- design$nt; nts <- design$nt_start
  base_r <- nts[g[use]]
  # nucleotide codes (0..3) per used read: columns 5' end, 5' end + 1, ..., 3' end, 3' end - 1, ...
  nr <- cbind(sapply(J, function(j) as.integer(nt[base_r + tx5[use] + j])),
              sapply(J, function(j) as.integer(nt[base_r + tx3[use] - j])))
  nr <- matrix(nr, ncol = 2L * nJ)
  row <- design$row_start[g[use]] + cod[use] - design$trim_start
  # ---- stage 1: fit on (a stride subsample of) positions with modelled reads
  upos <- sort(unique(row))
  if (length(upos) > max_positions)
    upos <- upos[unique(round(seq(1, length(upos), length.out = max_positions)))]
  pg <- design$gene[upos]
  pcod <- upos - design$row_start[pg] + design$trim_start
  # per position x category end positions
  t5 <- outer(3L * (pcod - 1L) + 1L, cat_f + cat_o, "+")
  t3 <- sweep(t5, 2, cat_L - 1L, "+")
  okp <- rowSums(t5 + min(J) < 1L | t3 - min(J) > 3L * design$ncodons[pg]) == 0
  upos <- upos[okp]; pg <- pg[okp]; t5 <- t5[okp, , drop = FALSE]; t3 <- t3[okp, , drop = FALSE]
  np <- length(upos)
  if (np < 100) stop("end_bias: fewer than 100 positions for the end bias fit")
  base <- nts[pg]
  pi <- match(row, upos)
  yb <- codon_dwell_tabulate_cpp(ifelse(is.na(pi), NA_integer_, (pi - 1L) * C + ci[use]),
                                 w[use], np * C)
  nf <- max(C - 1L, 4L * nL - 1L, 1L)
  lcat <- matrix(cat_l, np, C, byrow = TRUE)
  ntcol <- function(tp) as.raw(as.vector(t(4L * lcat + matrix(as.integer(nt[base + tp]), np))))
  codes <- do.call(cbind, c(list(as.raw(rep(seq_len(C) - 1L, np))),
                            lapply(J, function(j) ntcol(t5 + j)), lapply(J, function(j) ntcol(t3 - j))))
  K <- ncol(codes)
  fd <- list(codes = codes, offsets = seq_len(K), sense = character(nf + 1L))
  sel <- list(rows = seq_len(np * C) - 1L, starts = as.integer(seq(0, np * C, by = C)), n_genes = np)
  fc <- codon_dwell_feature_counts_cpp(codes, sel$rows, yb, nf)
  is_nt_ref <- c(TRUE, (seq_len(nf) %% 4L) == 0L) # codes 4 * l + 0 (A) are references
  est <- fc$n > 0 & fc$y > 0
  est[1, ] <- FALSE; est[is_nt_ref, -1] <- FALSE
  active <- which(est[-1, , drop = FALSE])
  fit <- codon_dwell_fisher(fd, sel, yb, active, threads = threads)
  # ---- coefficients per end, nucleotide position, read length and nucleotide (centered over the 4 nucleotides)
  acgt <- c("A", "C", "G", "T")
  coef <- rbindlist(lapply(2:K, function(k) rbindlist(lapply(seq_len(nL) - 1L, function(l) {
    code <- 4L * l + 0:3
    idx <- ifelse(code %% 4L == 0L, NA_integer_, (k - 1L) * nf + code)
    idx[!(idx %in% active)] <- NA_integer_
    cc <- codon_dwell_center(fit$beta, idx, active, fit$vcov)
    data.table(end = if (k - 1L <= nJ) "5'" else "3'", position = J[(k - 2L) %% nJ + 1L] + 1L,
               read_length = lens[l + 1L], nt = acgt, beta = cc$beta, se = cc$se)
  }))))
  # ---- stage 2: per-read weights, category totals kept
  lr <- cat_l[ci[use]] + 1L
  B <- numeric(length(lr))
  for (k in seq_len(2L * nJ)) {
    bk <- matrix(coef$beta[(k - 1L) * 4L * nL + seq_len(4L * nL)], 4) # nt x length
    B <- B + bk[cbind(nr[, k] + 1L, lr)]
  }
  wb <- w[use] * exp(-B)
  kc <- codon_dwell_tabulate_cpp(ci[use], w[use], C) / codon_dwell_tabulate_cpp(ci[use], wb, C)
  w_new <- w
  w_new[use] <- wb * kc[ci[use]]
  categories <- data.table(read_length = cat_L, frame = cat_f, offset = cat_o, fraction = tab$frac)
  list(w = w_new, coef = coef, categories = categories, n_positions = np,
       modelled_fraction = sum(w[use]) / sum(w[inside]), iterations = fit$iterations,
       converged = fit$converged)
}

# Shifting table (data.frame with fraction, offsets_start) for one library: shifts is a data.frame (used for all
# libraries) or a list named by library or by file (matched by file base name without extension and "_pshifted").
codon_dwell_shift_for <- function(shifts, library, path = NULL) {
  if (is.data.frame(shifts)) return(shifts)
  if (!is.list(shifts) || is.null(names(shifts))) stop("shifts must be a data.frame or a named list of them")
  stem <- function(x) sub("_pshifted$", "", sub("\\.[^.]*$", "", basename(x)))
  hit <- match(library, names(shifts))
  if (is.na(hit)) hit <- match(stem(library), stem(names(shifts)))
  if (is.na(hit) && !is.null(path)) hit <- match(stem(path), stem(names(shifts)))
  if (is.na(hit)) stop("No shifting table for library ", library)
  shifts[[hit]]
}

# Nucleotides per read end modelled by codon_dwell(end_bias = TRUE)
END_BIAS_NT <- 3L

codon_dwell_one <- function(reads, design, library, method, min_reads_gene,
                            heldout, threads, mrna, faFile, verbose,
                            end_bias = 0L, shifts = NULL) {
  t0 <- Sys.time()
  if (is.function(reads)) reads <- reads()
  path <- if (is.character(reads)) reads else attr(reads, "path")
  if (is.character(reads)) reads <- fimport(reads)
  map <- codon_dwell_map(reads, design)
  counts <- codon_dwell_counts(design = design, map = map)
  eb <- NULL
  if (end_bias > 0) {
    eb <- tryCatch({
      len <- mcols(reads)$size
      if (is.null(len)) stop("end_bias needs read lengths: a 'size' column (pshifted ofst / collapsed reads)")
      codon_dwell_end_bias(map, design, len, codon_dwell_shift_for(shifts, library, path),
                           end_nt = end_bias, threads = threads)
    }, error = function(e) e)
    if (!inherits(eb, "error")) counts$y <- codon_dwell_counts(design = design, map = map, w = eb$w)$y
  }
  rm(map)
  if (method == "glm") rm(reads)
  keep_gene <- counts$gene_total >= min_reads_gene
  S <- length(design$sense); K <- length(design$offsets)
  coef <- data.table(library = library, offset = rep(design$offsets, each = S),
                     site = "", codon = rep(design$sense, K),
                     aa = rep(design$aa, K), beta = NA_real_, se = NA_real_)
  coef[offset == 0L, site := "A"]; coef[offset == -1L, site := "P"]
  coef[offset == -2L, site := "E"]
  summ <- data.table(library = library, method = method,
                     n_genes = sum(keep_gene), n_positions = 0L,
                     cds_psites = counts$cds_psites,
                     psites_used = sum(counts$y[keep_gene[design$gene]]),
                     iterations = NA_integer_, converged = NA,
                     deviance_explained = NA_real_,
                     heldout_deviance_explained = NA_real_, status = "ok")
  if (end_bias > 0) summ[, `:=`(end_bias_categories = NA_integer_, end_bias_fraction = NA_real_)]
  fit <- tryCatch({
    if (inherits(eb, "error")) stop(conditionMessage(eb))
    if (sum(keep_gene) == 0)
      stop("no gene with >= ", min_reads_gene, " CDS P-sites")
    if (method == "glm") {
      codon_dwell_glm(design, counts, keep_gene, heldout = heldout, threads = threads)
    } else codon_dwell_ratio(design, counts, keep_gene, reads, mrna, faFile)
  }, error = function(e) e)
  eb_coef <- NULL
  if (inherits(fit, "error")) {
    summ[, status := conditionMessage(fit)]
    if (verbose) message("- ", library, ": ", conditionMessage(fit))
  } else {
    coef[, `:=`(beta = fit$beta, se = fit$se)]
    summ[, `:=`(n_positions = fit$n_positions,
                iterations = fit$summary$iterations,
                converged = fit$summary$converged,
                deviance_explained = fit$summary$deviance_explained,
                heldout_deviance_explained = fit$summary$heldout_deviance_explained)]
    if (end_bias > 0) {
      summ[, `:=`(end_bias_categories = nrow(eb$categories),
                  end_bias_fraction = eb$modelled_fraction)]
      eb_coef <- cbind(library = library, eb$coef)
    }
  }
  summ[, seconds := round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)]
  if (verbose && summ$status == "ok")
    message("- ", library, ": ", summ$n_genes, " genes, ", summ$n_positions,
            " positions, ", round(summ$seconds, 1), " s")
  list(coef = coef, summary = summ, end_bias = eb_coef)
}

#' Codon dwell times from Ribo-seq P-sites
#'
#' Estimates the relative time ribosomes spend on each codon (log relative
#' dwell time) at codon offsets around the A site, per library.
#'
#' Model (method "glm"): P-site counts per codon position i of gene g follow
#' \code{E[y_gi] = mu_g * exp(sum_k beta_k[codon at offset k])}, a log-linear
#' model with one-hot codon features at each offset. The gene effects
#' \code{mu_g} are profiled out exactly (the likelihood conditional on the
#' gene totals is multinomial within each gene), so no gene parameters are
#' fitted: Fisher scoring iterates only on beta. This gives the same
#' estimates as a Poisson GLM with one dummy variable per gene. Point
#' estimates are quasi-Poisson (consistent under overdispersion), standard
#' errors are gene-clustered sandwich (robust) errors. The design matrix is
#' never built: codon features are stored as codes per position
#' (\code{\link{codon_dwell_design}}), and the information matrix is
#' accumulated in C++ (OpenMP threads, see argument threads). Betas are reported sum-to-zero over the observed codons
#' per offset (beta = log dwell relative to the average codon at that
#' offset); SEs are for these sum-to-zero contrasts. Codons with no position
#' or no read at an offset get NA (and are pooled with the reference codon in
#' the fit; this only happens with very small data).
#'
#' Method "ratio": \code{\link{codon_usage}} (the engine of
#' \code{\link{codon_usage_exp}}) on the same CDSs and genes, with
#' beta = log(mean_txNorm) centered over the internal sense codons, for the P
#' site (offset -1) and A site (offset 0); other offsets are NA. As in
#' codon_usage_exp, all CDS positions are used (no trimming) and start / stop
#' codons are excluded from the table. It is fast and robust, but confounded
#' when codons co-occur (neighbouring codons, amino acids), which the glm
#' separates. SEs are delta-method approximations over genes
#' (\code{sqrt(var_txNorm / N) / mean_txNorm}).
#'
#' Footprint end bias (\code{end_bias}, choros-style technical terms):
#' nuclease cut preferences, ligation and reverse transcription make the
#' nucleotides at and next to the 5' and 3' ends of footprints change how
#' often a footprint is captured. For a read length with a fixed P-site offset
#' the ends sit at fixed codon offsets (with offset 12, a 28 nt read starts at
#' the first nucleotide of the codon at offset -5, its 5' cut site is the last
#' nucleotide of offset -6, and it ends at the first nucleotide of offset +4),
#' so without correction end bias shows up as codon "dwell" at the footprint
#' ends (offsets -6/-5 and +3..+5). With end_bias, P-sites are split into
#' categories (read length, frame of the P-site in its codon; categories with
#' at least 0.5\% of the P-sites, at most 15 read lengths) and the read ends
#' are located with the P-site offsets (\code{shifts}). A first stage fits,
#' within each codon position, the multinomial over categories with category
#' intercepts and, per read length, terms for the nucleotide just outside
#' each end (position 0, the cut site) and the first n nucleotides from each
#' end (positions 1..n, 1 = terminal nucleotide; n = end_bias, 3 for TRUE).
#' Dwell time and gene flux cancel in this conditional model, so the end bias
#' is estimated without a codon model (on up to 100000 positions with reads,
#' a stride subsample). Each read then gets the weight
#' \code{exp(-sum of its end terms)}, rescaled so every category keeps its
#' total, and the codon model is fitted on these corrected counts. Genes are
#' selected on the uncorrected counts. Reads of other categories keep weight
#' 1. The SEs of the codon betas do not include the uncertainty of the end bias
#' weights. On eIF5A / DOHH KO and WT HCT116 libraries (PRJNA545625) the
#' spread of codon betas at offsets -6, -5, +4 and +5 roughly halves, interior
#' offsets (-4..+2) are unchanged and the Pro effects at E and P are kept;
#' part of the -5 signal remains, so interpret footprint-end offsets with care.
#'
#' Prerequisites: reads must be P-sites with correct offsets per read length
#' (see \code{\link{shiftFootprintsByExperiment}}); estimates from libraries
#' with cycloheximide pretreatment or with initiation inhibitors
#' (harringtonine, lactimidomycin) are not comparable to untreated
#' libraries. Libraries with fewer than ~5 million CDS P-sites give noisy
#' estimates. Overlapping CDSs count a read once per CDS.
#' @param df an ORFik \code{\link{experiment}} or NULL. If given, it provides
#' the defaults for reads (the "pshifted" libraries, loaded one at a time),
#' cds and faFile.
#' @param reads NULL (load from df), or P-site reads of one library (GRanges,
#' GAlignments, path) or a named list of them. Weights from a "score" column
#' (collapsed reads) are used.
#' @param cds a GRangesList of CDSs, default (when df is given):
#' \code{loadRegion(df, "cds", filterTranscripts(df, 0, 1, 0,
#' longestPerGene = longestPerGene))}.
#' @param faFile a FaFile or path to the genome fasta, default \code{findFa(df)}
#' (or the fasta the design was built from).
#' @param mrna a GRangesList of transcripts containing the CDSs, only used by
#' method "ratio" (codon_usage needs a few nt of flank). Default:
#' \code{loadRegion(df, "mrna", names(design$cds))} when df is given, else the
#' CDSs extended by 6 nt on both sides.
#' @param design a "codon_dwell_design" from
#' \code{\link{codon_dwell_design}}, default NULL (built from cds, faFile,
#' offsets, trim_start, trim_end and min_cds_codons). Pass a cached design to
#' skip the sequence step when running many libraries; its offsets and
#' trimming are then used.
#' @inheritParams codon_dwell_design
#' @param min_reads_gene numeric, default 100. Genes with fewer P-sites in the
#' CDS (all codons) are not used.
#' @param longestPerGene logical, default TRUE. Use only the longest isoform
#' per gene for the default cds.
#' @param method "glm" (default, profiled-gene quasi-Poisson GLM) or "ratio".
#' @param heldout logical, default FALSE. For "glm", also fit on even genes
#' (rank among the used genes) and report the deviance explained on odd genes.
#' @param type character, default "pshifted". Library type loaded from df.
#' @param end_bias logical or integer, default FALSE (0). Correct for
#' footprint end nucleotide bias (see details): an integer 1..6 is the number
#' of nucleotides modelled from each read end (plus the flanking cut-site
#' nucleotide), TRUE = 3. Needs read lengths (a "size" column, as in
#' pshifted ofst files) and the P-site offsets (shifts).
#' @param shifts the P-site offsets per read length, used when end_bias is
#' TRUE: a data.frame with columns fraction (read length) and offsets_start
#' (negative offset of the P-site from the 5' end), used for all libraries,
#' or a list of them named by library or by file (as from
#' \code{\link{shifts_load}}, matched by file name without extension and
#' "_pshifted"). Default NULL: \code{shifts_load(df)}.
#' @param threads integer, default 4. OpenMP threads used by the glm fit
#' (libraries are processed one at a time, so memory stays at one library plus
#' the design). The dense information matrix products also use the BLAS
#' threads of R.
#' @param verbose logical, default TRUE.
#' @return a data.table with one row per library, offset and sense codon:
#' library, offset, site ("E", "P", "A" for offsets -2, -1, 0), codon, aa,
#' beta (log relative dwell time, sum-to-zero per library and offset) and se.
#' Attribute "summary" is a data.table per library with method, n_genes,
#' n_positions, cds_psites (P-sites in the design CDSs), psites_used,
#' iterations, converged, deviance_explained (in-sample, within-gene profile
#' vs uniform), heldout_deviance_explained, status and seconds. Attribute
#' "design" holds the design dimensions.
#' @references Ribo-DT: https://doi.org/10.1093/nar/gkaa203,
#' choros: https://doi.org/10.1101/2023.02.16.528794 (PMID 36865295)
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RFP", ][1, ]
#' # Tiny example genome: short trimming, E/P/A only, low gene filter
#' res <- codon_dwell(df, offsets = -2:0, trim_start = 5, trim_end = 5,
#'                    min_cds_codons = 30, min_reads_gene = 10)
#' res[site == "A"][order(-beta)]
#' attr(res, "summary")
#' ratio <- codon_dwell(df, offsets = -2:0, trim_start = 5, trim_end = 5,
#'                      min_cds_codons = 30, min_reads_gene = 10, method = "ratio")
#' cor(res$beta, ratio$beta, use = "complete.obs") # A and P site only
#' \dontrun{
#' # Real data: build the design once, reuse it for all libraries
#' design <- codon_dwell_design(loadRegion(df, "cds", filterTranscripts(df, 0, 1, 0)),
#'                              findFa(df))
#' res <- codon_dwell(df, design = design, heldout = TRUE)
#' codon_dwell_plot(res)
#' # Footprint end bias correction (P-site offsets from shifts_load(df))
#' res_eb <- codon_dwell(df, design = design, end_bias = TRUE)
#' attr(res_eb, "end_bias")
#' }
codon_dwell <- function(df = NULL, reads = NULL, cds = NULL, faFile = NULL,
                        design = NULL, mrna = NULL, offsets = -10:5, min_reads_gene = 100,
                        trim_start = 20L, trim_end = 10L, min_cds_codons = 100L,
                        longestPerGene = TRUE, method = c("glm", "ratio"),
                        heldout = FALSE, type = "pshifted", code = GENETIC_CODE,
                        end_bias = FALSE, shifts = NULL,
                        threads = 4L, verbose = TRUE) {
  method <- match.arg(method)
  stopifnot(min_reads_gene >= 1, length(end_bias) == 1, !is.na(end_bias))
  end_bias <- if (is.logical(end_bias)) END_BIAS_NT * end_bias else as.integer(end_bias)
  if (end_bias < 0 || end_bias > 6) stop("end_bias must be TRUE/FALSE or 0..6 nucleotides per end")
  if (end_bias > 0 && method != "glm") stop("end_bias is only available for method 'glm'")
  if (is.null(design)) {
    if (is.null(cds)) {
      if (is.null(df)) stop("Give df, or cds and faFile, or design")
      cds <- loadRegion(df, "cds", filterTranscripts(df, 0, 1, 0,
                                                     longestPerGene = longestPerGene))
    }
    if (is.null(faFile)) {
      if (is.null(df)) stop("Give faFile when df is NULL")
      faFile <- findFa(df)
    }
    design <- codon_dwell_design(cds, faFile, offsets = offsets,
                                 trim_start = trim_start, trim_end = trim_end,
                                 min_cds_codons = min_cds_codons, code = code)
  }
  codon_dwell_check_design(design)
  design <- codon_dwell_prepare(design, need_nt = end_bias > 0)
  if (method == "ratio") {
    if (is.null(faFile)) faFile <- if (!is.null(df)) findFa(df) else design$faFile
    if (is.null(mrna)) {
      mrna <- if (!is.null(df)) {
        loadRegion(df, "mrna", names(design$cds))
      } else extendTrailers(extendLeaders(design$cds, 6L), 6L)
    }
    if (!all(names(design$cds) %in% names(mrna))) stop("mrna must contain all CDS names")
  }
  if (is.null(reads)) {
    if (is.null(df)) stop("Give reads or df")
    libs <- bamVarName(df)
    reads <- lapply(seq_len(nrow(df)), function(i) {
      path <- filepath(df[i, ], type)
      function() fimport(path)
    })
    names(reads) <- libs
    if (end_bias > 0 && is.null(shifts)) {
      sh <- shifts_load(df)
      stem <- function(x) sub("_pshifted$", "", sub("\\.[^.]*$", "", basename(x)))
      hit <- match(stem(filepath(df, type)), stem(names(sh)))
      if (anyNA(hit)) stop("No shifting table (shifts_load(df)) for libraries: ",
                           paste(libs[is.na(hit)], collapse = ", "))
      shifts <- setNames(sh[hit], libs)
    }
  } else if (!is.list(reads) || is(reads, "GRangesList")) {
    if (is(reads, "GRangesList")) reads <- as.list(reads) else reads <- list(reads)
    if (is.null(names(reads))) names(reads) <- paste0("library", seq_along(reads))
  } else if (is.null(names(reads))) names(reads) <- paste0("library", seq_along(reads))
  if (end_bias > 0 && is.null(shifts)) {
    if (is.null(df)) stop("end_bias needs shifts (P-site offsets per read length) or df")
    shifts <- shifts_load(df)
  }
  if (verbose) message("-- Codon dwell (", method, "): ", length(design$cds),
                       " CDSs, ", nrow(design$codes), " positions, ",
                       length(reads), " libraries")
  res <- mapply(codon_dwell_one, reads = reads, library = names(reads),
                MoreArgs = list(design = design, method = method,
                                min_reads_gene = min_reads_gene,
                                heldout = heldout, threads = threads,
                                mrna = mrna, faFile = faFile, verbose = verbose,
                                end_bias = end_bias, shifts = shifts),
                SIMPLIFY = FALSE)
  coef <- rbindlist(lapply(res, `[[`, "coef"))
  attr(coef, "summary") <- rbindlist(lapply(res, `[[`, "summary"), fill = TRUE)
  if (end_bias > 0) attr(coef, "end_bias") <- rbindlist(lapply(res, `[[`, "end_bias"))
  attr(coef, "design") <- list(n_cds = length(design$cds),
                               n_positions = nrow(design$codes),
                               offsets = design$offsets,
                               trim_start = design$trim_start,
                               trim_end = design$trim_end)
  coef[]
}

#' Plot codon dwell times
#'
#' Heatmap of beta (log relative dwell time) per codon and offset, or bars of
#' beta per codon at one offset with 95\% confidence intervals.
#' @param res output of \code{\link{codon_dwell}}
#' @param type "heatmap" (offset x codon) or "site" (bars at one offset).
#' @param library character, the libraries to plot (facets), default: the
#' first library in res.
#' @param offset integer, the offset for type "site", default 0 (A site).
#' @param limit numeric, default NULL. Color scale limit for the heatmap
#' (symmetric), default the 99th percentile of |beta|.
#' @return a ggplot object
#' @family codon
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RFP", ][1, ]
#' res <- codon_dwell(df, offsets = -2:0, trim_start = 5, trim_end = 5,
#'                    min_cds_codons = 30, min_reads_gene = 10)
#' codon_dwell_plot(res)
#' codon_dwell_plot(res, type = "site", offset = 0)
codon_dwell_plot <- function(res, type = c("heatmap", "site"),
                             library = unique(res$library)[1], offset = 0L,
                             limit = NULL) {
  type <- match.arg(type)
  lib_sel <- library
  dt <- copy(res[res$library %in% lib_sel, ])
  if (nrow(dt) == 0) stop("No rows for library: ", paste(lib_sel, collapse = ", "))
  dt[, label := paste0(aa, ":", codon)]
  if (type == "heatmap") {
    lev <- unique(dt[order(aa, codon)]$label)
    dt[, label := factor(label, levels = rev(lev))]
    dt[, off := factor(ifelse(site != "", paste0(offset, " (", site, ")"),
                              as.character(offset)),
                       levels = unique(ifelse(site != "", paste0(offset, " (", site, ")"),
                                              as.character(offset))[order(offset)]))]
    if (is.null(limit)) limit <- stats::quantile(abs(dt$beta), 0.99, na.rm = TRUE)
    dt[, fill := pmax(pmin(beta, limit), -limit)]
    p <- ggplot(dt, aes(x = off, y = label, fill = fill)) + geom_tile() +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-limit, limit),
                           name = "log dwell") +
      labs(x = "Offset from A site (codons)", y = "AA:codon") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 6),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    off_sel <- offset
    dt <- dt[dt$offset == off_sel & !is.na(beta)]
    if (nrow(dt) == 0) stop("No estimates at offset ", off_sel)
    lev <- unique(dt[, .(m = mean(beta)), by = label][order(m)]$label)
    dt[, label := factor(label, levels = lev)]
    p <- ggplot(dt, aes(x = label, y = beta, fill = beta > 0)) + geom_col() +
      geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
                    width = 0.3, na.rm = TRUE) +
      scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                        guide = "none") +
      coord_flip() +
      labs(x = "AA:codon", y = paste0("log relative dwell, offset ", off_sel)) +
      theme_minimal() + theme(axis.text.y = element_text(size = 6))
  }
  if (length(unique(dt$library)) > 1) p <- p + facet_wrap(~ library)
  p
}
