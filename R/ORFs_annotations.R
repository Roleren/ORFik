#' Create small artificial orfs from cds
#'
#' Usefull to see if short ORFs prediction is dependent on length.\cr
#' Split cds first in two, a start part and stop part.
#' Then say how large the two parts can be and merge them together.
#' It will sample a value in range give.\cr
#' Parts will be forced to not overlap and can not extend outside
#' original cds
#'
#' If artificial cds length is not divisible by 2, like 3 codons,
#' the second codon will always be from the start region etc.\cr
#' Also If there are many very short original cds, the distribution
#' will be skewed towards more smaller artificial cds.
#' @param cds a GRangesList of orfs, must have width
#'  \%\% 3 == 0 and length >= 6
#' @param start5 integer, default: 1 (start of orf)
#' @param end5 integer, default: 4 (max 4 codons from start codon)
#' @param start3 integer, default -4 (max 4 codons from stop codon)
#' @param end3 integer, default: 0 (end of orf)
#' @param bin.if.few logical, default TRUE, instead of per codon,
#' do per 2, 3, 4 codons if you have few samples compared to lengths wanted,
#' If you have 4 cds' and you want 7 different lengths, which is the standard,
#' it will give you possible nt length: 6-12-18-24 instead of original
#' 6-9-12-15-18-21-24.\cr
#' If you have more than 30x cds than lengths wanted this is skipped.
#' (for default arguments this is: 7*30 = 210 cds)
#' @return GRangesList of new ORFs (sorted: + strand increasing start,
#'  - strand decreasing start)
artificial.orfs <- function(cds, start5 = 1, end5 = 4, start3 = -4, end3 = 0,
                            bin.if.few = TRUE) {
  #start5 = 1; end5 = 4; start3 = -4; end3 = 0;bin.if.few = TRUE; groupings <- ORFik:::groupings
  validGRL(class(cds), type = "cds")
  widths <- widthPerGroup(cds)
  names <- names(cds)

  if (start5 > end5) stop("start5 > end5 argument")
  if (start3 > end3) stop("start3 > end3 argument")
  if (!all(widths %% 3 == 0)) stop("not all cds has width moduls 3 = 0")
  if (!all(widths >= 6)) stop("not all cds has width >= 6")
  if (length(cds) < 4) stop("No meaningful test for so few cds!")
  # Make random cuts, distribute them equally
  max.size <- (end5 * 3) - (start3 * 3)
  combinations <- floor((max.size - 1) / 3)
  if (combinations < 1) stop("No possible combination with input")
  bin <- 3
  if (bin.if.few & (length(cds) < (30 * combinations))) {
    while (((max.size / bin) * 30 > length(cds)) & (max.size / bin > 2)) {
      bin <- bin + 3
    }
  }
  possible <- seq.int(6, max.size, by = bin)
  possible <- rep(possible, length(cds) / length(possible))
  miss <- length(cds) - length(possible)
  if (miss > 0) possible <- c(possible, sample(possible, size = miss))
  possible <- sample(possible)

  rand5 <- possible / 2
  dif <- rand5 %% 1
  rand5 <- rand5 + (dif*3)
  rand3 <- (possible / 2) - (dif*3)
  if (!all(possible == (rand5 + rand3))) stop("Bug, report on ORFik issues page!")
  rand3 <- widths - rand3 + 1

  possible_start <- data.table(max = widths - 3, random = rand5)
  possible_start[, min := pmin(max, random)]
  possible_start[, pick := pmin(max, min)]
  start_part <- IRanges(start5, possible_start$pick)

  possible_end <- data.table(min = possible_start$pick + 1, random = rand3)
  possible_end[, max := pmax(min, random)]
  possible_end[, pick := pmax(max, min)]
  end_part <- IRanges(possible_end$pick, widths)

  start_grl <- pmapFromTranscriptF(start_part, cds, removeEmpty = TRUE)
  end_grl <- pmapFromTranscriptF(end_part, cds, removeEmpty = TRUE)
  merged_gr <- c(unlistGrl(start_grl), unlistGrl(end_grl))
  merged_gr <- split(merged_gr, c(groupings(start_grl), groupings(end_grl)))
  names(merged_gr) <- names(cds)
  merged_gr <- reduce(merged_gr)
  new_widths <- widthPerGroup(merged_gr, FALSE)
  if (!all(new_widths > 0) | !all(new_widths %% 3 == 0) |
      length(cds) != length(merged_gr) | !all(names(merged_gr) == names))
    stop("Error during creation of artifical ORFs!")
  message("Grouping statistics:")
  dt <- data.table(widths = new_widths)
  print(summary(dt[, .(group.size = .N), by = widths]))
  message("If minimum group.size is < 30, it can be hard to use for statistical purpose!")
  return(sortPerGroup(merged_gr))
}


#' Overlaps GRanges object with provided annotations.
#'
#' @param rel_orf - GRanges object of your ORF.
#' @param tran - GRanges object of annotation (transcript or cds) that
#' overlapped in some way rel_orf.
#' @param isoform_names - A vector of strings that will be used instead of
#' these defaults:
#' 'perfect_match' - start and stop matches the tran object strand wise
#' 'elong_START_match' - rel_orf is extension from the STOP side of the tran
#' 'trunc_START_match' - rel_orf is truncation from the STOP side of the tran
#' 'elong_STOP_match' - rel_orf is extension from the START side of the tran
#' 'trunc_STOP_match' - rel_orf is truncation from the START side of the tran
#' 'overlap_inside' - rel_orf is inside tran object
#' 'overlap_both' - rel_orf contains tran object inside
#' 'overlap_upstream' - rel_orf is overlaping upstream part of the tran
#' 'overlap_downstream' - rel_orf is overlaping downstream part of the tran
#' 'upstream' - rel_orf is upstream towards the tran
#' 'downstream' - rel_orf is downstream towards the tran
#' 'none' - when none of the above options is true
#' @return A string object of defined isoform towards transcript.
#'
defineIsoform <- function(
  rel_orf, tran, isoform_names = c(
    "perfect_match", "elong_START_match", "trunc_START_match",
    "elong_STOP_match", "trunc_STOP_match", "overlap_inside",
    "overlap_both", "overlap_upstream", "overlap_downstream",
    "upstream", "downstram", "none")) {
    stran <- as.vector(strand(rel_orf))[1]
    if (stran == "+") {

        orf_start <- start(rel_orf)[1]
        orf_end <- end(rel_orf)[length(rel_orf)]
        tran_start <- start(tran)[1]
        tran_end <- end(tran)[length(tran)]

        if (orf_start == tran_start & orf_end == tran_end) {
            return(isoform_names[1])
        }
        if (orf_start == tran_start & orf_end > tran_end) {
            return(isoform_names[2])
        }
        if (orf_start == tran_start & orf_end < tran_end) {
            return(isoform_names[3])
        }
        if (orf_start < tran_start & orf_end == tran_end) {
            return(isoform_names[4])
        }
        if (orf_start > tran_start & orf_end == tran_end) {
            return(isoform_names[5])
        }
        if (orf_start > tran_start & orf_end < tran_end) {
            return(isoform_names[6])
        }
        if (orf_start < tran_start & orf_end > tran_end) {
            return(isoform_names[7])
        }
        if (orf_start < tran_start & orf_end < tran_end) {
            return(isoform_names[8])
        }
        if (orf_start > tran_start & orf_end > tran_end) {
            return(isoform_names[9])
        }
        if (orf_start < tran_start & orf_end < tran_start) {
            return(isoform_names[10])
        }
        if (orf_start > tran_end & orf_end > tran_end) {
            return(isoform_names[11])
        }
    } else {

        orf_start <- end(rel_orf)[length(rel_orf)]
        orf_end <- start(rel_orf)[1]
        tran_start <- end(tran)[length(tran)]
        tran_end <- start(tran)[1]

        if (orf_start == tran_start & orf_end == tran_end) {
            return(isoform_names[1])
        }
        if (orf_start == tran_start & orf_end < tran_end) {
            return(isoform_names[2])
        }
        if (orf_start == tran_start & orf_end > tran_end) {
            return(isoform_names[3])
        }
        if (orf_start > tran_start & orf_end == tran_end) {
            return(isoform_names[4])
        }
        if (orf_start < tran_start & orf_end == tran_end) {
            return(isoform_names[5])
        }
        if (orf_start < tran_start & orf_end > tran_end) {
            return(isoform_names[6])
        }
        if (orf_start > tran_start & orf_end < tran_end) {
            return(isoform_names[7])
        }
        if (orf_start > tran_start & orf_end > tran_end) {
            return(isoform_names[8])
        }
        if (orf_start < tran_start & orf_end < tran_end) {
            return(isoform_names[9])
        }
        if (orf_start > tran_start & orf_end > tran_start) {
            return(isoform_names[10])
        }
        if (orf_start < tran_end & orf_end < tran_end) {
            return(isoform_names[11])
        }
    }
    return(isoform_names[12])
}


#' Overlaps GRanges object with provided annotations.
#'
#' It will return same list of GRanges, but with metdata columns:
#' trainscript_id - id of transcripts that overlap with each ORF
#' gene_id - id of gene that this transcript belongs to
#' isoform - for coding protein alignment in relation to cds on coresponding
#' transcript,
#' for non-coding transcripts alignment in relation to the transcript.
#' @param ORFs - GRanges or GRangesList object of your ORFs.
#' @param con - Path to gtf file with annotations.
#' @return A GRanges object of your ORFs with metadata columns 'gene',
#' 'transcript', isoform' and 'biotype'.
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#'
assignAnnotations <- function(ORFs, con) {

    message("Loading annotations from gtf file")
    txdb <- makeTxDbFromGFF(con, format = "gtf")
    gtf_annot <- rtracklayer::import(con, format = "gtf")
    gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
    transcript_df <- data.frame(gtf_annot$transcript_id, gtf_annot$gene_id,
                                gtf_annot$transcript_biotype)
    names(transcript_df) <- c("transcript_id", "gene_id", "transcript_biotype")

    # remove non unique transcript/biotype combinations
    transcript_df <- transcript_df[!duplicated(transcript_df[
      c("transcript_id", "transcript_biotype")]), ]

    # find out overlaping transcripts
    transcripts <- exonsBy(txdb, by = "tx", use.names = TRUE)
    transcripts_hits <- findOverlaps(ORFs, transcripts, type = "any")
    t_names <- names(transcripts)
    t_names <- t_names[subjectHits(transcripts_hits)]
    # ORFs that have no overlap with transcripts
    notranscript <- which(!(seq_along(ORFs) %in%
                              unique(queryHits(transcripts_hits))))

    # which transcript corespond to our table of annotations
    tran_matches <- match(t_names, transcript_df$transcript_id)

    newORFs <- ORFs[queryHits(transcripts_hits)]
    newORFs$transcript_id <- t_names
    newORFs$gene_id <- transcript_df$gene_id[tran_matches]
    newORFs$transcript_biotype <-
      transcript_df$transcript_biotype[tran_matches]

    notranscriptORFs <- ORFs[notranscript]
    notranscriptORFs$transcript_id <- "none"
    notranscriptORFs$gene_id <- "none"
    notranscriptORFs$transcript_biotype <- "none"
    notranscriptORFs$isoform <- "none"

    newIsoforms_p <- c()
    newIsoforms_o <- c()
    cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
    p_coding <- newORFs[newORFs$transcript_biotype == "protein_coding"]
    other_coding <- newORFs[newORFs$transcript_biotype != "protein_coding"]
    message("Preparing annotations for ORFs overlapping coding regions...")
    for (i in seq_along(p_coding)) {
        newIsoforms_p[i] <- defineIsoform(p_coding[i],
                                          cds[[p_coding[i]$transcript_id]])
    }
    message("Preparing annotations for ORFs overlapping non-coding regions...")
    for (i in seq_along(other_coding)) {
        newIsoforms_o[i] <- defineIsoform(
          other_coding[i], transcripts[[other_coding[i]$transcript_id]])
    }
    p_coding$isoform <- newIsoforms_p
    other_coding$isoform <- newIsoforms_o

    finallyORFs <- c(p_coding, other_coding, notranscriptORFs)
    message("Sorting results...")
    finallyORFs <- sort(finallyORFs)
}
