#' Returns start definition according to
#' \url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1}
#' ncbi genetic code number for translation.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of START sies separatd with "|".
#' @export
#'
start_definition <- function(transl_table) {
  STARTdef <- c("ATG|TTG|CTG", #1 The Standard Code
                "ATT|ATC|ATA|ATG|GTG", #2 The Vertebrate Mitochondrial Code
                "ATA|ATG", #3 The Yeast Mitochondrial Code
                "TTA|TTG|CTG|ATT|ATC|ATA|ATG|GTG", #4 Mold/Protozoan/Coelenterate Mitochondrial/Mycoplasma/Spiroplasma
                "TTG|ATT|ATC|ATA|ATG|GTG", #5 The Invertebrate Mitochondrial Code
                "ATG", #6 The Ciliate, Dasycladacean and Hexamita Nuclear Code
                "ATG", #7 ??? No info
                "ATG", #8 ??? No info
                "ATG|GTG", #9 The Echinoderm and Flatworm Mitochondrial Code
                "ATG", #10 The Euplotid Nuclear Code
                "TTG|CTG|ATT|ATC|ATA|ATG|GTG", #11 The Bacterial, Archaeal and Plant Plastid Code
                "CTG|ATG", #12 The Alternative Yeast Nuclear Code
                "TTG|ATA|ATG|GTG", #13 The Ascidian Mitochondrial Code
                "ATG", #14 The Alternative Flatworm Mitochondrial Code
                "ATG", #15 ??? No info
                "ATG", #16 Chlorophycean Mitochondrial Code
                "ATG", #17 ??? No info
                "ATG", #18 ??? No info
                "ATG", #19 ??? No info
                "ATG", #20 ??? No info
                "ATG|GTG", #21 Trematode Mitochondrial Code
                "ATG", #22 Scenedesmus obliquus Mitochondrial Code
                "ATT|ATG|GTG", #23 Thraustochytrium Mitochondrial Code
                "TTG|CTG|ATG|GTG", #24 Pterobranchia Mitochondrial Code
                "TTG|ATG|GTG", #25 Candidate Division SR1 and Gracilibacteria Code
                "CTG|ATG") #26 Pachysolen tannophilus Nuclear Code
  return(STARTdef[transl_table])
}


#' Returns stop definition according to
#' \url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1}
#' ncbi genetic code number for translation.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of STOP sies separatd with "|".
#' @export
#'
stop_definition <- function(transl_table) {
  STOPdef <- c("TAA|TAG|TGA", #1 The Standard Code
               "TAA|TAG|AGA|AGG", #2 The Vertebrate Mitochondrial Code
               "TAA|TAG", #3 The Yeast Mitochondrial Code
               "TAA|TAG", #4 Mold/Protozoan/Coelenterate Mitochondrial/Mycoplasma/Spiroplasma
               "TAA|TAG", #5 The Invertebrate Mitochondrial Code
               "TGA", #6 The Ciliate, Dasycladacean and Hexamita Nuclear Code
               "TGA", #7 ??? No info
               "TGA", #8 ??? No info
               "TAA|TAG", #9 The Echinoderm and Flatworm Mitochondrial Code
               "TAA|TAG", #10 The Euplotid Nuclear Code
               "TAA|TAG|TGA", #11 The Bacterial, Archaeal and Plant Plastid Code
               "TAA|TAG|TGA", #12 The Alternative Yeast Nuclear Code
               "TAA|TAG", #13 The Ascidian Mitochondrial Code
               "TAG", #14 The Alternative Flatworm Mitochondrial Code
               "TGA", #15 ??? No info
               "TGA", #16 Chlorophycean Mitochondrial Code
               "TGA", #17 ??? No info
               "TGA", #18 ??? No info
               "TGA", #19 ??? No info
               "TGA", #20 ??? No info
               "TAA|TAG", #21 Trematode Mitochondrial Code
               "TCA|TAA|TGA", #22 Scenedesmus obliquus Mitochondrial Code
               "TTA|TAA|TAG|TGA", #23 Thraustochytrium Mitochondrial Code
               "TAA|TAG", #24 Pterobranchia Mitochondrial Code
               "TAA|TAG", #25 Candidate Division SR1 and Gracilibacteria Code
               "TAA|TAG") #26 Pachysolen tannophilus Nuclear Code
  return(STOPdef[transl_table])
}

#' Creates IRanges with Open Reading Frames from input sequences.
#'
#' @param seq string. Input sequence.
#' @param seqname string. Name of the sequence, meybe name of the transcript or name of the chromosome.
#' Final GRanges will use this name as seqname.
#' @param startCodon string. Default is NULL. This is gonna overwrite START codons used in numcode. eg. "ATG|GTG".
#' @param strand string. One of "*" which is default, "+" or "-"
#' @param longestORF bolean. Default TRUE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping, small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least START + 8*3 [bp] + STOP.
#' @param numcode numeric.
#' The \url{http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1}
#' ncbi genetic code number for translation.
#' @return A List of GRanges of ORFs mapped to fasta file. Each ORF includes START and STOP codon,
#' frame, seqname, and strand.
#' @export
#' @import S4Vectors
#' @import seqinr
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' #find_ORFs()
#'
find_ORFs <- function(seq,
                      seqname = "unknown",
                      startCodon = NULL,
                      strand = "*",
                      numcode = 1,
                      longestORF = T,
                      minimumLength = 0) {

    startCodon <- ifelse(is.null(startCodon), start_definition(numcode), startCodon)
    stopCodon <- stop_definition(numcode)

    refseq <- s2c(seq)
    finalRanges <- GRanges()

    if (strand != "-") { # strand is + or *
      f0 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 0, numcode = numcode)
      f1 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 1, numcode = numcode)
      f2 <- getTrans(refseq, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 2, numcode = numcode)

      f0 <- which(f0 == "*") * 3
      f1 <- which(f1 == "*") * 3 + 1
      f2 <- which(f2 == "*") * 3 + 2

      #f0
      f0ORFs <- c()
      for (chunk in 1:length(f0)) {
        starto <- ifelse(chunk == 1, 1, f0[chunk - 1] + 1)
        endo <- f0[chunk]
        seq_to_search <- toupper(getSequence(refseq[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f0ORFs <- c(tempORF[end(tempORF) == endo], f0ORFs)
      }
      f0ORFs <- GRanges(Rle(rep(seqname, length(f0ORFs))), f0ORFs, Rle(strand(rep("+", length(f0ORFs)))))
      f0ORFs$Frame <- rep("0", length(f0ORFs))
      finalRanges <- c(f0ORFs, finalRanges)

      # f1
      f1ORFs <- c()
      for (chunk in 1:length(f1)) {
        starto <- ifelse(chunk == 1, 1, f1[chunk - 1] + 1)
        endo <- f1[chunk]
        seq_to_search <- toupper(getSequence(refseq[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f1ORFs <- c(tempORF[end(tempORF) == endo], f1ORFs)
      }
      f1ORFs <- GRanges(Rle(rep(seqname, length(f1ORFs))), f1ORFs, Rle(strand(rep("+", length(f1ORFs)))))
      f1ORFs$Frame <- rep("1", length(f1ORFs))
      finalRanges <- c(f1ORFs, finalRanges)

      # f2
      f2ORFs <- c()
      for (chunk in 1:length(f2)) {
        starto <- ifelse(chunk == 1, 1, f2[chunk - 1] + 1)
        endo <- f2[chunk]
        seq_to_search <- toupper(getSequence(refseq[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f2ORFs <- c(tempORF[end(tempORF) == endo], f2ORFs)
      }
      f2ORFs <- GRanges(Rle(rep(seqname, length(f2ORFs))), f2ORFs, Rle(strand(rep("+", length(f2ORFs)))))
      f2ORFs$Frame <- rep("2", length(f2ORFs))
      finalRanges <- c(f2ORFs, finalRanges)
    }

    if (strand != "+") {
      refseq_r <- rev(comp(refseq))
      f0r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 0, numcode = numcode)
      f1r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 1, numcode = numcode)
      f2r <- getTrans(refseq_r, sens = "F", NAstring = "N", ambiguous = FALSE, frame = 2, numcode = numcode)

      f0r <- which(f0r == "*") * 3
      f1r <- which(f1r == "*") * 3 + 1
      f2r <- which(f2r == "*") * 3 + 2

      full_string <- length(refseq_r)
      # f0
      f0ORFs <- c()
      for (chunk in 1:length(f0r)) {
        starto <- ifelse(chunk == 1, 1, f0r[chunk - 1] + 1)
        endo <- f0r[chunk]
        seq_to_search <- toupper(getSequence(refseq_r[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f0ORFs <- c(tempORF[end(tempORF) == endo], f0ORFs)
      }
      f0ORFs <- GRanges(Rle(rep(seqname, length(f0ORFs))),
                        IRanges(full_string - end(f0ORFs) + 1, width = width(f0ORFs)),
                        Rle(strand(rep("-", length(f0ORFs)))))
      f0ORFs$Frame <- rep("0", length(f0ORFs))
      finalRanges <- c(f0ORFs, finalRanges)

      # f1
      f1ORFs <- c()
      for (chunk in 1:length(f1r)) {
        starto <- ifelse(chunk == 1, 1, f1r[chunk - 1] + 1)
        endo <- f1r[chunk]
        seq_to_search <- toupper(getSequence(refseq_r[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f1ORFs <- c(tempORF[end(tempORF) == endo], f1ORFs)
      }
      f1ORFs <- GRanges(Rle(rep(seqname, length(f1ORFs))),
                        IRanges(full_string - end(f1ORFs) + 1, width = width(f1ORFs)),
                        Rle(strand(rep("-", length(f1ORFs)))))
      f1ORFs$Frame <- rep("1", length(f1ORFs))
      finalRanges <- c(f1ORFs, finalRanges)

      # f2
      f2ORFs <- c()
      for (chunk in 1:length(f2r)) {
        starto <- ifelse(chunk == 1, 1, f2r[chunk - 1] + 1)
        endo <- f2r[chunk]
        seq_to_search <- toupper(getSequence(refseq_r[starto:endo], as.string = T))
        tempORF <- find_in_frame_ORFs(seq_to_search,
                                      startCodon = startCodon,
                                      longestORF = longestORF,
                                      minimumLength = minimumLength)
        tempORF <- shift(tempORF, starto - 1)
        #stop codon must be same as in ORF so we know frame is correct
        f2ORFs <- c(tempORF[end(tempORF) == endo], f2ORFs)
      }
      f2ORFs <- GRanges(Rle(rep(seqname, length(f2ORFs))),
                        IRanges(full_string - end(f2ORFs) + 1, width = width(f2ORFs)),
                        Rle(strand(rep("-", length(f2ORFs)))))
      f2ORFs$Frame <- rep("2", length(f2ORFs))
      finalRanges <- c(f2ORFs, finalRanges)
    }

    finalRanges <- sort(finalRanges)
}
