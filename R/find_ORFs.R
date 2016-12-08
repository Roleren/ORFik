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


#' Creates list of IRanges with Open Reading Frames.
#'
#' @param fastaSeq DNA sequence to search for Open Reading Frames.
#' @param startCodon string. Default is "ATG".
#' @param stopCodon string. Default is "TAA|TAG|TGA".
#' @param longestORF bolean. Default FALSE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least START + 8*3 [bp] + STOP.
#' @return A List of IRanges objects of ORFs.
#' @export
#' @import IRanges
#' @examples
#' #find_in_frame_ORFs()
#'
find_in_frame_ORFs <- function(fastaSeq,
                               startCodon = "ATG",
                               stopCodon = "TAA|TAG|TGA",
                               longestORF = F,
                               minimumLength = 0) {

  codpos <- paste0(if (longestORF) paste0("(?<!(", startCodon, "))"), "(?=(", startCodon, "))(?=(?:([ATGCN]{3}))*?(", stopCodon, "))")
  a <- gregexpr(codpos, fastaSeq, perl = TRUE)[[1]]
  gr <- IRanges(start = as.vector(a), end = attr(a,"capture.start")[, if (longestORF) 4 else 3] + 2)

  gr <- gr[width(gr) >=  6 + minimumLength*3]
  return(gr[order(start(gr), end(gr))])
}


#' Creates GRanges with Open Reading Frames from fasta files.
#'
#' Each fasta header is treated separately, and name of the sequence will be used as seqname in
#' returned GRanges object. Frame of the Open Reading Frame is also returned in
#' metadata column 'frame'.
#' @param file - Path to fasta file.
#' @param startCodon Default is "ATG".
#' @param stopCodon Default is "TAA|TAG|TGA".
#' @param longestORF bolean. Default TRUE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping small ones.
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
#'
find_ORFs_fa <-
  function(file,
           startCodon = c("ATG"),
           stopCodon = c("TAA|TAG|TGA"),
           longestORF = T,
           minimumLength = 0) {
    message("Loading fasta file.")
    refseq <- read.fasta(file)
    message("Preparing reverse complementary sequence.")
    refseq_r <- lapply(refseq, function(x)
      rev(comp(x)))
    message("Defining frames and stop codons.")
    f0 <-
      getTrans(
        refseq,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 0,
        numcode = 1
      )
    f1 <-
      getTrans(
        refseq,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 1,
        numcode = 1
      )
    f2 <-
      getTrans(
        refseq,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 2,
        numcode = 1
      )
    f0r <-
      getTrans(
        refseq_r,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 0,
        numcode = 1
      )
    f1r <-
      getTrans(
        refseq_r,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 1,
        numcode = 1
      )
    f2r <-
      getTrans(
        refseq_r,
        sens = "F",
        NAstring = "N",
        ambiguous = FALSE,
        frame = 2,
        numcode = 1
      )

    stopCodons <- function(x) {
      return(which(x == "*"))
    }

    f0 <- lapply(f0, stopCodons)
    f0 <- lapply(f0, function(x)
      x * 3)
    f1 <- lapply(f1, stopCodons)
    f1 <- lapply(f1, function(x)
      x * 3 + 1)
    f2 <- lapply(f2, stopCodons)
    f2 <- lapply(f2, function(x)
      x * 3 + 2)
    f0r <- lapply(f0r, stopCodons)
    f0r <- lapply(f0r, function(x)
      x * 3)
    f1r <- lapply(f1r, stopCodons)
    f1r <- lapply(f1r, function(x)
      x * 3 + 1)
    f2r <- lapply(f2r, stopCodons)
    f2r <- lapply(f2r, function(x)
      x * 3 + 2)

    seq_names <- names(refseq)
    finalRanges <- GRanges()

    for (seqname in 1:length(seq_names)) {
      message(paste0("Finding ORFs for "), seq_names[seqname])
      # plus strand f0
      f0ORFs <- c()
      for (chunk in 1:length(f0[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f0[[seqname]][chunk - 1] + 1)
        endo <- f0[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <- shift(tempORF, starto - 1)
        f0ORFs <-
          c(tempORF[end(tempORF) == endo], f0ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f0ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f0ORFs))), f0ORFs, Rle(strand(rep(
          "+", length(f0ORFs)
        ))))
      f0ORFs$Frame <- rep("0", length(f0ORFs))
      finalRanges <- c(f0ORFs, finalRanges)

      # f1
      f1ORFs <- c()
      for (chunk in 1:length(f1[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f1[[seqname]][chunk - 1] + 1)
        endo <- f1[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <- shift(tempORF, starto - 1)
        f1ORFs <-
          c(tempORF[end(tempORF) == endo], f1ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f1ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f1ORFs))), f1ORFs, Rle(strand(rep(
          "+", length(f1ORFs)
        ))))
      f1ORFs$Frame <- rep("1", length(f1ORFs))
      finalRanges <- c(f1ORFs, finalRanges)

      # f2
      f2ORFs <- c()
      for (chunk in 1:length(f2[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f2[[seqname]][chunk - 1] + 1)
        endo <- f2[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <- shift(tempORF, starto - 1)
        f2ORFs <-
          c(tempORF[end(tempORF) == endo], f2ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f2ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f2ORFs))), f2ORFs, Rle(strand(rep(
          "+", length(f2ORFs)
        ))))
      f2ORFs$Frame <- rep("2", length(f2ORFs))
      finalRanges <- c(f2ORFs, finalRanges)

      # minus strand
      full_string <- length(refseq_r[[seqname]])
      # f0
      f0ORFs <- c()
      for (chunk in 1:length(f0r[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f0r[[seqname]][chunk - 1] + 1)
        endo <- f0r[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <- shift(tempORF, starto - 1)
        f0ORFs <-
          c(tempORF[end(tempORF) == endo], f0ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f0ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f0ORFs))),
                IRanges(full_string - end(f0ORFs) + 1, width = width(f0ORFs)),
                Rle(strand(rep(
                  "-", length(f0ORFs)
                ))))
      f0ORFs$Frame <- rep("0", length(f0ORFs))
      finalRanges <- c(f0ORFs, finalRanges)

      # f1
      f1ORFs <- c()
      for (chunk in 1:length(f1r[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f1r[[seqname]][chunk - 1] + 1)
        endo <- f1r[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <- shift(tempORF, starto - 1)
        f1ORFs <-
          c(tempORF[end(tempORF) == endo], f1ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f1ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f1ORFs))),
                IRanges(full_string - end(f1ORFs) + 1, width = width(f1ORFs)),
                Rle(strand(rep(
                  "-", length(f1ORFs)
                ))))
      f1ORFs$Frame <- rep("1", length(f1ORFs))
      finalRanges <- c(f1ORFs, finalRanges)

      # f2
      f2ORFs <- c()
      for (chunk in 1:length(f2r[[seqname]])) {
        starto <- ifelse(chunk == 1, 1, f2r[[seqname]][chunk - 1] + 1)
        endo <- f2r[[seqname]][chunk]
        seq_to_search <-
          toupper(getSequence(refseq_r[[seqname]][starto:endo], as.string = T))
        tempORF <-
          find_in_frame_ORFs(
            seq_to_search,
            startCodon = startCodon,
            longestORF = longestORF,
            minimumLength = minimumLength
          )
        tempORF <-
          shift(tempORF, starto - 1)
        f2ORFs <-
          c(tempORF[end(tempORF) == endo], f2ORFs)  #stop codon must be same as in ORF so we know frame is correct
      }
      f2ORFs <-
        GRanges(Rle(rep(seq_names[seqname], length(f2ORFs))),
                IRanges(full_string - end(f2ORFs) + 1, width = width(f2ORFs)),
                Rle(strand(rep(
                  "-", length(f2ORFs)
                ))))
      f2ORFs$Frame <-
        rep("2", length(f2ORFs))
      finalRanges <- c(f2ORFs, finalRanges)
    }

    finalRanges <- sort(finalRanges)
  }
