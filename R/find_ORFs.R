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


#' Replaces string at pos with N.
#'
#' @param string string.
#' @param pos numeric vector.
#' @return string with swapped N at pos
#'
subchar <- function(string, pos) {
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr needed for this function to work. Please install it.", call. = FALSE)
  }

  for (i in pos) {
    stringr::str_sub(string, i, i) <- "N"
    # string <- sub(paste("^(.{", i - 1, "}).", sep = ""), "\\1N", string, perl = TRUE)
  }
  return(string)
}


#' Creates list of IRanges with Open Reading Frames.
#'
#' @param fastaSeq DNA sequence to search for Open Reading Frames.
#' @param startCodon string. Default is "ATG".
#' @param stopCodon string. Default is "TAA|TAG|TGA".
#' @param longestORF bolean. Default TRUE. Defines whether pick longest ORF only.
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
                               longestORF = T,
                               minimumLength = 0) {
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("stringr needed for this function to work. Please install it.", call. = FALSE)
  }

  codpos <- paste0("(?:", startCodon, ")(?:[ATGCN]{3}(?<!", stopCodon, ")){", minimumLength, ",}(?:", stopCodon, ")")

  frame <- stringr::str_locate_all(fastaSeq, codpos)[[1]]
  gr <- IRanges(start = frame[, 1], end = frame[, 2])

  while (dim(frame)[1] != 0 && !longestORF) {
    starts <- c(sapply(frame, as.integer))
    fastaSeq <- subchar(fastaSeq, starts)
    frame <- stringr::str_locate_all(fastaSeq, codpos)[[1]]

    frGR <- IRanges(start = frame[, 1], end = frame[, 2])
    gr <- c(gr, frGR)
  }

  return(gr[order(start(gr), end(gr))])
}
