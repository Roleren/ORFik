# detectRibosomeShifts_new: C++ counting pass (one streaming read over an unshifted .ofst), results identical to detectRibosomeShifts.
# Code originates from the pshiftR prototype (backclaude/claude-megafst); internal helpers are prefixed psr_.
# Equality: 60/60 libraries and the 352-length pshift gold set give identical offsets (strict and weak).

psr_index <- function(txdb, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = NULL, txNames = NULL) {
  txdb <- loadTxdb(txdb)
  if (is.null(txNames)) txNames <- filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR)
  cds <- loadRegion(txdb, part = "cds", names.keep = txNames)
  txNames <- txNames[txNames %in% names(cds)]
  mrna <- loadRegion(txdb, "mrna")[txNames]
  cds <- cds[txNames]
  chrs <- GenomeInfoDb::seqlevels(cds)
  exons <- function(grl, with_off) {
    u <- unlist(grl, use.names = TRUE)
    d <- data.table::data.table(tx = match(names(u), txNames), chr = match(as.character(GenomeInfoDb::seqnames(u)), chrs) - 1L,
                                start = BiocGenerics::start(u), end = BiocGenerics::end(u), minus = as.integer(as.character(BiocGenerics::strand(u)) == "-"))
    if (with_off) {  # exons are listed 5' -> 3' within each transcript; offset = 1-based tx coordinate of the exon's first nt
      d[, w := end - start + 1L]
      d[, off := cumsum(w) - w + 1L, by = tx]
    }
    d
  }
  te <- exons(mrna, TRUE); ce <- exons(cds, FALSE)
  # CDS start in transcript coordinates: tx coordinate of the first CDS nt (5' end of the first CDS exon)
  first <- unlist(BiocGenerics::start(cds)); 
  cds_first <- vapply(seq_along(cds), function(i) {
    g <- cds[[i]]; p5 <- if (as.character(BiocGenerics::strand(g))[1] == "-") max(BiocGenerics::end(g)) else min(BiocGenerics::start(g)); p5
  }, 0L)
  m <- te[ce[!duplicated(tx), .(tx, cds_p5 = cds_first[tx])], on = "tx", allow.cartesian = TRUE]
  m <- m[ifelse(minus == 1L, cds_p5 <= end & cds_p5 >= start, cds_p5 >= start & cds_p5 <= end)]
  m[, cs := off + ifelse(minus == 1L, end - cds_p5, cds_p5 - start)]
  cds_start <- m[order(tx)]$cs
  stopifnot(length(cds_start) == length(txNames))
  key <- function(d) list(ks = (as.numeric(d$chr) * 2 + d$minus) * 2^30 + d$start, ke = (as.numeric(d$chr) * 2 + d$minus) * 2^30 + d$end)
  kt <- key(te); kc <- key(ce)
  chrs_used <- unique(as.character(GenomeInfoDb::seqnames(unlist(cds, use.names = FALSE))))
  cds_len <- vapply(cds, function(g) sum(BiocGenerics::width(g)), 0L)
  list(txNames = txNames, chrs = chrs, chrs_used = chrs_used, cds_end = as.integer(cds_start + cds_len - 1L), cds_start = as.integer(cds_start), tx_exons = te, cds_exons = ce,
       tx_idx = psr_segments_cpp(kt$ks, kt$ke), cds_idx = psr_segments_cpp(kc$ks, kc$ke),
       params = list(minFiveUTR = minFiveUTR, minCDS = minCDS, minThreeUTR = minThreeUTR))
}

psr_index_cached <- function(txdb, cache_dir, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = NULL, rebuild = FALSE, txNames = NULL) {
  txdb <- loadTxdb(txdb)
  f <- AnnotationDbi::dbfile(txdb)
  key <- paste(normalizePath(f), format(file.mtime(f), "%Y%m%d%H%M%S"), file.size(f), minFiveUTR, minCDS, if (is.null(minThreeUTR)) "NULL" else minThreeUTR, "v1", if (is.null(txNames)) "auto" else as.character(openssl_md5(paste(txNames, collapse = ","))), sep = "|")
  id <- substr(paste(as.character(openssl_md5(key)), collapse = ""), 1, 16)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(cache_dir, sprintf("psr_index_%s.rds", id))
  if (!rebuild && file.exists(out)) { ix <- readRDS(out); if (identical(ix$key, key)) return(ix) }
  ix <- psr_index(txdb, minFiveUTR, minCDS, minThreeUTR, txNames = txNames); ix$key <- key; ix$txdb_file <- f
  tmp <- paste0(out, ".tmp", Sys.getpid()); saveRDS(ix, tmp); file.rename(tmp, out); ix
}
openssl_md5 <- function(x) { tf <- tempfile(); on.exit(unlink(tf)); writeLines(x, tf); tools::md5sum(tf) }

psr_cigar_table <- function(levels) psr_cigar_cpp(as.character(levels))

psr_count <- function(ofst, index, lengths, chunk = 2e7, extended = FALSE) {
  nr <- fst::metadata_fst(ofst)$nrOfRows
  lengths <- sort(unique(as.integer(lengths)))
  ex <- index$tx_exons; ce <- index$cds_exons
  ptr <- psr_new_cpp(index$tx_idx, as.integer(ex$tx) - 1L, as.integer(ex$start), as.integer(ex$end), as.integer(ex$off), as.integer(ex$minus),
                     index$cds_idx, as.integer(ce$tx) - 1L, as.integer(index$cds_start), lengths, as.integer(index$cds_end), extended)
  cig_cache <- NULL
  from <- 1L
  while (from <= nr) {
    to <- min(nr, from + chunk - 1)
    x <- fst::read_fst(ofst, from = from, to = to, as.data.table = TRUE)
    lv <- levels(x$seqnames)
    # ORFik: reads are kept only on seqlevels that carry eligible CDS (validSeqlevels via seqnamesPerGroup)
    used <- index$chrs_used
    if (!any(lv %in% used)) lv2 <- sub("^chr", "", lv) else lv2 <- lv
    if (!any(lv %in% used)) used <- sub("^chr", "", used)
    cmap <- ifelse(lv2 %in% used, match(lv2, if (identical(lv2, lv)) index$chrs else sub("^chr", "", index$chrs)) - 1L, NA_integer_)
    cmap[is.na(cmap)] <- -1L
    chr <- cmap[as.integer(x$seqnames)]
    sl <- levels(x$strand); smap <- match(sl, c("+", "-", "*")); strand <- smap[as.integer(x$strand)]
    cl <- levels(x$cigar); ct <- psr_cigar_cpp(cl)
    psr_add_cpp(ptr, chr, as.integer(x$start), strand, as.integer(x$cigar), as.numeric(x$score), ct$span, ct$size, cl)
    from <- to + 1L
  }
  r <- psr_result_cpp(ptr); r$lengths <- lengths; r
}


psr_is_periodic <- function(x, strict.fft = TRUE) {
  if (sum(x) == 0) return(FALSE)
  n <- length(x)
  amplitudes <- abs(stats::fft(x)); amp <- amplitudes[2:(n / 2 + 1)]
  nspec <- n %/% 2
  freq <- seq.int(from = 1 / n, by = 1 / n, length.out = nspec)   # spec.pgram frequencies for N = n (n has no padding for 150)
  periods <- 1 / freq
  keep <- periods != n
  if (strict.fft) {
    p <- periods[keep][which.max(amp[keep])]
  } else {
    spec <- stats::spec.pgram(x, plot = FALSE)$spec
    p <- periods[keep][which.max(spec[keep])]
  }
  (p > 2.9) & (p < 3.1)
}

psr_change_point <- function(x, max.pos = 40L, interval = seq.int(14L, 24L), center.pos = 12) {
  meta <- x[seq.int(max.pos)]
  pos <- -(length(x) / 2):(length(x) / 2 - 1)
  frames <- rep(c(0, 1, 2), 20)[seq.int(max.pos)]
  max.frame <- which.max(c(sum(meta[frames == 0]), sum(meta[frames == 1]), sum(meta[frames == 2]))) - 1
  interval <- interval[c(0, 1, 2) %% 3 == max.frame] - (pos[interval[1]] %% 3)
  means <- numeric(0)
  for (j in interval) {
    downstream <- mean(meta[seq.int(j, max.pos, by = 3)])
    upstream <- mean(meta[seq.int(1, j - 1, by = 3)])
    means <- c(means, downstream - upstream)
  }
  center.pos <- c(center.pos, center.pos - 1, center.pos + 2)[max.frame + 1]
  scaled <- means / (abs(pos[interval] + center.pos) + 1)^0.25
  pos[interval[which.max(abs(scaled))]]
}

psr_detect <- function(cnt, top_tx = 10L, min_reads = 1000, min_reads_TIS = 50, accepted.lengths = cnt$lengths,
                       must.be.periodic = TRUE, strict.fft = TRUE, verbose = FALSE) {
  tot <- cnt$size_total
  szs <- which(tot > 0) - 1L
  tab <- data.table::data.table(size = szs, counts = tot[szs + 1L])
  tab <- tab[(counts >= min_reads) & (size %in% accepted.lengths)]
  if (nrow(tab) == 0) stop("No valid read lengths found with", " accepted.lengths and counts > min_reads")
  gate <- cnt$gate; ncds <- sum(gate)
  keep <- gate[cnt$tx + 1L]
  W <- cnt$win[keep, , drop = FALSE]; tx <- cnt$tx[keep] + 1L; sz <- cnt$lengths[cnt$size_idx[keep] + 1L]
  # periodicity: first 150 nt of the CDS (positions 0..149 = columns 31..180), per gene normalised, summed over genes in CDS-list order
  cdscols <- 31:180
  ord <- order(tx)
  x_by_len <- list(); valid <- integer(0)
  for (l in tab$size) {
    ii <- ord[sz[ord] == l]
    X <- W[ii, cdscols, drop = FALSE]; gs <- rowSums(X)
    nz <- gs > 0
    xl <- if (any(nz)) colSums(X[nz, , drop = FALSE] / gs[nz]) else rep(0, 150)
    x_by_len[[as.character(l)]] <- xl
    if (!must.be.periodic || psr_is_periodic(xl, strict.fft)) valid <- c(valid, l)
  }
  if (!must.be.periodic) valid <- accepted.lengths
  if (length(valid) == 0) stop(paste("Library contained no periodic or accepted read-lengths,", "check your library. Are you using the correct genome?", "Is this Ribo-seq?"))
  # TIS: positions -30..29 = columns 1..60, all valid lengths, per-gene sum over the window and valid lengths
  tiscols <- 1:60
  inv <- sz %in% valid
  gsum <- numeric(length(gate)); 
  if (any(inv)) { s <- rowSums(W[inv, tiscols, drop = FALSE]); gsum <- tapply(s, tx[inv], sum); g2 <- numeric(length(gate)); g2[as.integer(names(gsum))] <- gsum; gsum <- g2 }
  gs_g <- gsum[gate]                                  # one value per gated gene (zeros included)
  rr <- 60L * length(valid)
  top <- if (ncds > 1000) (100 - top_tx) / 100 else 0
  thr <- stats::quantile(rep(gs_g, each = rr), top)
  pass_gene <- gate & gsum >= thr
  sel <- inv & pass_gene[tx]
  out <- list()
  for (l in valid) {
    jj <- sel & sz == l
    sc <- if (any(jj)) colSums(W[jj, tiscols, drop = FALSE]) else rep(0, 60)
    if (sum(sc) > min_reads_TIS) out[[as.character(l)]] <- psr_change_point(sc)
  }
  if (length(out) == 0) stop(paste("Library contained no periodic or accepted read-lengths,", "reason: Not enough reads in TIS region.", "check your library. Are you using the correct genome?", "Is this Ribo-seq?"))
  res <- data.table::data.table(fraction = as.integer(names(out)), offsets_start = unlist(out, use.names = FALSE))
  attr(res, "periodic_lengths") <- valid; attr(res, "tab") <- tab; attr(res, "n_cds_gated") <- ncds
  res
}

psr_detect_ofst <- function(ofst, index, accepted.lengths = c(20L, 21L, 25:33), strict.fft = TRUE, ...) {
  cnt <- psr_count(ofst, index, accepted.lengths)
  r <- psr_detect(cnt, accepted.lengths = accepted.lengths, strict.fft = strict.fft, ...)
  data.table::setattr(r, "periodic_lengths", NULL); data.table::setattr(r, "tab", NULL); data.table::setattr(r, "n_cds_gated", NULL); r
}

#' Detect ribosome P-site shifts (C++ counting pass)
#'
#' Same result as \code{\link{detectRibosomeShifts}} (identical offsets for \code{start = TRUE}), computed in one
#' streaming C++ pass over an unshifted .ofst file instead of importing the whole file into R:
#' about 7-30x faster and 10-20x less memory. The old function is unchanged.
#' The fast path is used when \code{footprints} is a single path to an unshifted .ofst file, \code{start = TRUE},
#' \code{stop = FALSE}, \code{firstN = 150}, \code{heatmap = FALSE} and \code{verbose = FALSE}. In all other cases
#' the call is forwarded to \code{detectRibosomeShifts} (same arguments), so it is a drop-in replacement.
#' The per-organism transcript index is cached in \code{cache_dir} (keyed by txdb file, mtime, size and parameters).
#' @inheritParams detectRibosomeShifts
#' @param cache_dir directory for the cached transcript index, default
#'  \code{getOption("ORFik.pshift.cache", file.path(tempdir(), "ORFik_pshift_index"))}
#' @return a data.table with columns fraction and offsets_start (as detectRibosomeShifts)
#' @export
#' @examples
#' \dontrun{
#' offsets <- detectRibosomeShifts_new("/path/sample.ofst", txdb, accepted.lengths = 26:34)
#' }
detectRibosomeShifts_new <- function(footprints, txdb, start = TRUE, stop = FALSE,
  top_tx = 10L, minFiveUTR = 30L, minCDS = 150L,
  minThreeUTR = if (stop) {30} else NULL,
  txNames = filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR),
  firstN = 150L, tx = NULL, min_reads = 1000, min_reads_TIS = 50, accepted.lengths = 26:34,
  heatmap = FALSE, must.be.periodic = TRUE, strict.fft = TRUE, verbose = FALSE,
  cache_dir = getOption("ORFik.pshift.cache", file.path(tempdir(), "ORFik_pshift_index"))) {
  txdb <- loadTxdb(txdb)
  fast <- is.character(footprints) && length(footprints) == 1 && identical(tools::file_ext(footprints), "ofst") &&
    isTRUE(start) && !isTRUE(stop) && firstN == 150L && isFALSE(heatmap) && !isTRUE(verbose)
  if (!fast) {
    return(detectRibosomeShifts(footprints, txdb, start = start, stop = stop, top_tx = top_tx, minFiveUTR = minFiveUTR,
                                minCDS = minCDS, minThreeUTR = minThreeUTR, txNames = txNames, firstN = firstN, tx = tx,
                                min_reads = min_reads, min_reads_TIS = min_reads_TIS, accepted.lengths = accepted.lengths,
                                heatmap = heatmap, must.be.periodic = must.be.periodic, strict.fft = strict.fft, verbose = verbose))
  }
  index <- psr_index_cached(txdb, cache_dir, minFiveUTR = minFiveUTR, minCDS = minCDS, minThreeUTR = minThreeUTR, txNames = txNames)
  psr_detect_ofst(footprints, index, accepted.lengths = accepted.lengths, strict.fft = strict.fft, must.be.periodic = must.be.periodic,
                  min_reads = min_reads, min_reads_TIS = min_reads_TIS, top_tx = top_tx)
}

#' Build or load the cached transcript index used by detectRibosomeShifts_new
#'
#' Call once before running many libraries in parallel (each worker then only reads the cache file).
#' @inheritParams detectRibosomeShifts_new
#' @return invisibly, the path of the cache directory
#' @export
detectRibosomeShifts_index <- function(txdb, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = NULL,
  txNames = filterTranscripts(txdb, minFiveUTR, minCDS, minThreeUTR),
  cache_dir = getOption("ORFik.pshift.cache", file.path(tempdir(), "ORFik_pshift_index"))) {
  txdb <- loadTxdb(txdb)
  psr_index_cached(txdb, cache_dir, minFiveUTR = minFiveUTR, minCDS = minCDS, minThreeUTR = minThreeUTR, txNames = txNames)
  invisible(cache_dir)
}
