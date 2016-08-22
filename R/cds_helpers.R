#' Creates GRanges object of Transcription Start Sites from CDS.
#'
#' It will return A postions of ATG transcription start sites.
#' @param cds GRanges List object of your CDSs.
#' @return A GRanges object of TSS sites.
#' @export
#' @import GenomicRanges
#' @examples
#' cds <- GRangesList('gene_plus_strand' = GRanges(Rle(c('1'), c(4)),
#'                                                 IRanges(c(925942, 930155, 939040, 939275),
#'                                                         width=c(72, 182, 90, 17)),
#'                                                 Rle(strand(c('+', '+', '+', '+')))),
#'                    'gene_minus_strand' = GRanges(Rle(c('1'), c(4)),
#'                                                  IRanges(c(245929870, 245927931,
#'                                                            245863799, 245858562),
#'                                                          width=c(32, 103, 88, 109)),
#'                                                  Rle(strand(c('-', '-', '-', '-')))))
#' extract_START_sites_from_CDSs(cds)
#'
extract_START_sites_from_CDSs <- function(cds) {
    cdsStarts <- start(cds)
    cdsEnds <- end(cds)
    isPlusStrand <- unlist(runValue(strand(cds)) == "+")

    cdsFixedStarts <- mapply(function(starts, ends, isPlusStrand) {
        if (isPlusStrand) {
            starts[1]
        } else {
            ends[1]
        }
    }, cdsStarts, cdsEnds, isPlusStrand)

    cdsStarts <- GRanges(seqnames = Rle(as.vector(unlist(runValue(seqnames(cds))))), ranges = IRanges(cdsFixedStarts, width = 1),
        strand = strand(as.vector(unlist(runValue(strand(cds))))))
    return(cdsStarts)
}


#' Creates GRanges object of STOP base pairs from CDS.
#'
#' @param cds GRanges List object of your CDSs.
#' @return A GRanges object of last STOP nucletodie.
#' @export
#' @import GenomicRanges
#' @examples
#' #will return G postions as in TAG stop codon
#' cds <- GRangesList('gene_plus_strand' = GRanges(Rle(c('1'), c(4)),
#'                                                 IRanges(c(925942, 930155, 939040, 939275),
#'                                                         width=c(72, 182, 90, 17)),
#'                                                 Rle(strand(c('+', '+', '+', '+')))),
#'                    'gene_minus_strand' = GRanges(Rle(c('1'), c(4)),
#'                                                  IRanges(c(245929870, 245927931,
#'                                                            245863799, 245858562),
#'                                                          width=c(32, 103, 88, 109)),
#'                                                  Rle(strand(c('-', '-', '-', '-')))))
#' extract_STOP_sites_from_CDSs(cds)
#'
extract_STOP_sites_from_CDSs <- function(cds) {
    cdsStarts <- start(cds)
    cdsEnds <- end(cds)
    isPlusStrand <- unlist(runValue(strand(cds)) == "+")

    cdsFixedEnds <- mapply(function(starts, ends, isPlusStrand) {
        if (isPlusStrand) {
            ends[length(ends)]
        } else {
            starts[length(starts)]
        }
    }, cdsStarts, cdsEnds, isPlusStrand)

    cdsEnds <- GRanges(seqnames = Rle(as.vector(unlist(runValue(seqnames(cds))))), ranges = IRanges(cdsFixedEnds, width = 1),
        strand = strand(as.vector(unlist(runValue(strand(cds))))))
    return(cdsEnds)
}
