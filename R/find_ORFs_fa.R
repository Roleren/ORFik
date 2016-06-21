#' Creates GRanges with Open Reading Frames from fasta files.
#'
#' Each fasta header is treated separately, and name of the sequence will be used as seqname in
#' returned GRanges object. Frame of the Open Reading Frame is also returned in
#' metadata column 'frame'.
#' @param file - Path to fasta file.
#' @param startCodon Default is c('ATG').
#' @param stopCodon Default is c('TAA', 'TAG', 'TGA').
#' @param longestORF bolean. Default TRUE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping, small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least START + 8*3 [bp] + STOP.
#' @return A List of GRanges of ORFs mapped to fasta file. Each ORF includes START and STOP codon,
#' frame, seqname, and strand.
#' @export
#' @import S4Vectors
#' @import seqinr
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' #find_ORFs_fa()

find_ORFs_fa <- function(file, startCodon = c("ATG"), stopCodon = c("TAA", "TAG", "TGA"), longestORF = T, minimumLength = 0) {
    
    message("Loading fasta file.")
    refseq <- read.fasta(file)
    message("Preparing reverse complementary sequence.")
    refseq_r <- lapply(refseq, function(x) rev(comp(x)))
    message("Defining frames and stop codons.")
    f0 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 0, numcode = 1)
    f1 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 1, numcode = 1)
    f2 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 2, numcode = 1)
    f0r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 0, numcode = 1)
    f1r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 1, numcode = 1)
    f2r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 2, numcode = 1)
    
    stopCodons <- function(x) {
        return(which(x == "*"))
    }
    
    f0 <- lapply(f0, stopCodons)
    f0 <- lapply(f0, function(x) x * 3)
    f1 <- lapply(f1, stopCodons)
    f1 <- lapply(f1, function(x) x * 3 + 1)
    f2 <- lapply(f2, stopCodons)
    f2 <- lapply(f2, function(x) x * 3 + 2)
    f0r <- lapply(f0r, stopCodons)
    f0r <- lapply(f0r, function(x) x * 3)
    f1r <- lapply(f1r, stopCodons)
    f1r <- lapply(f1r, function(x) x * 3 + 1)
    f2r <- lapply(f2r, stopCodons)
    f2r <- lapply(f2r, function(x) x * 3 + 2)
    
    seq_names <- names(refseq)
    finalRanges <- GRanges()
    
    for (seqname in 1:length(seq_names)) {
        message(paste0("Finding ORFs for "), seq_names[seqname])
        # plus strand f0
        f0ORFs <- c()
        for (chunk in 1:length(f0[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f0[[seqname]][chunk - 1] + 1)
            endo <- f0[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f0ORFs <- c(tempORF[end(tempORF) == endo], f0ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f0ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f0ORFs))), f0ORFs, Rle(strand(rep("+", length(f0ORFs)))))
        f0ORFs$Frame <- rep("0", length(f0ORFs))
        finalRanges <- c(f0ORFs, finalRanges)
        
        # f1
        f1ORFs <- c()
        for (chunk in 1:length(f1[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f1[[seqname]][chunk - 1] + 1)
            endo <- f1[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f1ORFs <- c(tempORF[end(tempORF) == endo], f1ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f1ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f1ORFs))), f1ORFs, Rle(strand(rep("+", length(f1ORFs)))))
        f1ORFs$Frame <- rep("1", length(f1ORFs))
        finalRanges <- c(f1ORFs, finalRanges)
        
        # f2
        f2ORFs <- c()
        for (chunk in 1:length(f2[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f2[[seqname]][chunk - 1] + 1)
            endo <- f2[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f2ORFs <- c(tempORF[end(tempORF) == endo], f2ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f2ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f2ORFs))), f2ORFs, Rle(strand(rep("+", length(f2ORFs)))))
        f2ORFs$Frame <- rep("2", length(f2ORFs))
        finalRanges <- c(f2ORFs, finalRanges)
        
        # minus strand
        full_string <- length(refseq_r[[seqname]])
        # f0
        f0ORFs <- c()
        for (chunk in 1:length(f0r[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f0r[[seqname]][chunk - 1] + 1)
            endo <- f0r[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f0ORFs <- c(tempORF[end(tempORF) == endo], f0ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f0ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f0ORFs))), IRanges(full_string - end(f0ORFs) + 1, width = width(f0ORFs)), 
            Rle(strand(rep("-", length(f0ORFs)))))
        f0ORFs$Frame <- rep("0", length(f0ORFs))
        finalRanges <- c(f0ORFs, finalRanges)
        
        # f1
        f1ORFs <- c()
        for (chunk in 1:length(f1r[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f1r[[seqname]][chunk - 1] + 1)
            endo <- f1r[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f1ORFs <- c(tempORF[end(tempORF) == endo], f1ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f1ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f1ORFs))), IRanges(full_string - end(f1ORFs) + 1, width = width(f1ORFs)), 
            Rle(strand(rep("-", length(f1ORFs)))))
        f1ORFs$Frame <- rep("1", length(f1ORFs))
        finalRanges <- c(f1ORFs, finalRanges)
        
        # f2
        f2ORFs <- c()
        for (chunk in 1:length(f2r[[seqname]])) {
            starto <- ifelse(chunk == 1, 1, f2r[[seqname]][chunk - 1] + 1)
            endo <- f2r[[seqname]][chunk]
            seq_to_search <- toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
            tempORF <- find_in_frame_ORFs(seq_to_search, startCodon = startCodon, longestORF = longestORF, minimumLength = minimumLength)
            tempORF <- shift(tempORF, starto - 1)
            f2ORFs <- c(tempORF[end(tempORF) == endo], f2ORFs)  #stop codon must be same as in ORF so we know frame is correct
        }
        f2ORFs <- GRanges(Rle(rep(seq_names[seqname], length(f2ORFs))), IRanges(full_string - end(f2ORFs) + 1, width = width(f2ORFs)), 
            Rle(strand(rep("-", length(f2ORFs)))))
        f2ORFs$Frame <- rep("2", length(f2ORFs))
        finalRanges <- c(f2ORFs, finalRanges)
    }
    
    finalRanges <- sort(finalRanges)
}
