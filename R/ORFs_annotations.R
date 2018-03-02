#' Overlaps GRanges object with provided annotations.
#'
#' @param rel_orf - GRanges object of your ORF.
#' @param tran - GRanges object of annotation (transcript or cds) that overlapped in
#' some way rel_orf.
#' @param isoform_names - A vector of strings that will be used instead of these defaults:
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
#' @import GenomicRanges
#' @export
#' @examples
#' #define_isoform()
define_isoform <- function(rel_orf, tran, isoform_names = c("perfect_match",
                                                            "elong_START_match",
                                                            "trunc_START_match",
                                                            "elong_STOP_match",
                                                            "trunc_STOP_match",
                                                            "overlap_inside",
                                                            "overlap_both",
                                                            "overlap_upstream",
                                                            "overlap_downstream",
                                                            "upstream",
                                                            "downstram",
                                                            "none")) {
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
#' isoform - for coding protein alignment in relation to cds on coresponding transcript,
#' for non-coding transcripts alignment in relation to the transcript.
#' @param ORFs - GRanges or GRangesList object of your ORFs.
#' @param con - Path to gtf file with annotations.
#' @return A GRanges object of your ORFs with metadata columns 'gene', 'transcript',
#' 'isoform' and 'biotype'.
#' @export
#' @import GenomicRanges
#' @import rtracklayer
#' @import GenomicFeatures
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @examples
#' #assign_annotations()

assign_annotations <- function(ORFs, con) {

    message("Loading annotations from gtf file")
    txdb <- makeTxDbFromGFF(con, format = "gtf")
    gtf_annot <- import(con, format = "gtf")
    gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
    transcript_df <- data.frame(gtf_annot$transcript_id, gtf_annot$gene_id, gtf_annot$transcript_biotype)
    names(transcript_df) <- c("transcript_id", "gene_id", "transcript_biotype")

    # remove non unique transcript/biotype combinations
    transcript_df <- transcript_df[!duplicated(transcript_df[c("transcript_id", "transcript_biotype")]), ]

    # find out overlaping transcripts
    transcripts <- exonsBy(txdb, by = "tx", use.names = T)
    transcripts_hits <- findOverlaps(ORFs, transcripts, type = "any")
    t_names <- names(transcripts)
    t_names <- t_names[subjectHits(transcripts_hits)]
    # ORFs that have no overlap with transcripts
    notranscript <- which(!(1:length(ORFs) %in% unique(queryHits(transcripts_hits))))

    # which transcript corespond to our table of annotations
    tran_matches <- match(t_names, transcript_df$transcript_id)

    newORFs <- ORFs[queryHits(transcripts_hits)]
    newORFs$transcript_id <- t_names
    newORFs$gene_id <- transcript_df$gene_id[tran_matches]
    newORFs$transcript_biotype <- transcript_df$transcript_biotype[tran_matches]

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
    for (i in 1:length(p_coding)) {
        newIsoforms_p[i] <- define_isoform(p_coding[i], cds[[p_coding[i]$transcript_id]])
    }
    message("Preparing annotations for ORFs overlapping non-coding regions...")
    for (i in 1:length(other_coding)) {
        newIsoforms_o[i] <- define_isoform(other_coding[i], transcripts[[other_coding[i]$transcript_id]])
    }
    p_coding$isoform <- newIsoforms_p
    other_coding$isoform <- newIsoforms_o

    finallyORFs <- c(p_coding, other_coding, notranscriptORFs)
    message("Sorting results...")
    finallyORFs <- sort(finallyORFs)
}
