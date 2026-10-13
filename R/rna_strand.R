#' Normalize single-end RNA alignments to transcript strand
#'
#' Coordinates, CIGAR and read metadata are preserved. Paired interpretation
#' remains controlled by strandMode; this transformation is for single-end
#' alignments, including a single mate selected before alignment.
#' @param x A GAlignments object.
#' @param read_strand One of alignment, opposite, or unstranded. Unstranded
#' retains alignment strands and requires unstranded downstream coverage.
#' @return The alignment object with RNA orientation metadata.
#' @export
normalize_rna_strand <- function(x, read_strand = c("alignment", "opposite", "unstranded")) {
  read_strand <- match.arg(read_strand)
  if (!is(x, "GAlignments"))
    stop("Single-end RNA normalization requires GAlignments; use strandMode for pairs")
  if (read_strand == "opposite") {
    s <- strand(x)
    # BAMs sorted by genomic position contain long equal-strand runs. Flip
    # the run values, preserving lengths, instead of expanding every read.
    strand(x) <- Rle(factor(c("-", "+", "*")[as.integer(runValue(s))],
                           levels = c("+", "-", "*")), runLength(s))
  }
  attr(x, "rna_strand") <- read_strand
  x
}

# Small sidecars keep the established OFST schema compatible and make the
# orientation survive serialization. They also inform direct BAM imports.
rna_strand_sidecar <- function(path) {
  root <- dirname(path)
  kind <- if (grepl("\\.ofst$", path)) "ofst_all" else "bam"
  if (kind == "ofst_all" && basename(root) == "ofst") {
    root <- dirname(root)
    if (basename(root) == "unique_mappers") {
      root <- dirname(root); kind <- "ofst_unique"
    }
  }
  folder <- file.path(root, ".rna_strand")
  if (kind != "bam") folder <- file.path(folder, kind)
  file.path(folder, paste0(tools::file_path_sans_ext(basename(path)), "_rna_strand.rds"))
}

rna_strand_record <- function(path) {
  if (!is.character(path) || length(path) != 1L) return(NULL)
  sidecar <- rna_strand_sidecar(path)
  if (!file.exists(sidecar)) return(NULL)
  record <- readRDS(sidecar)
  if (!is.null(record$product_signature)) {
    st <- file.info(path)
    if (!identical(c(size = st$size, mtime = as.numeric(st$mtime)), record$product_signature))
      stop("RNA strand provenance is stale for ", path, "; regenerate its strand decision/product")
  }
  record
}

rna_strand_resolve <- function(path, read_strand = "auto") {
  read_strand <- match.arg(read_strand, c("auto", "alignment", "opposite", "unstranded"))
  if (read_strand != "auto") return(read_strand)
  record <- rna_strand_record(path)
  if (is.null(record)) return("alignment")
  if (identical(record$status, "unknown"))
    stop("RNA orientation is unresolved for ", path)
  if (isTRUE(record$normalized)) return(if (record$read_strand == "unstranded") "unstranded" else "alignment")
  record$read_strand
}
