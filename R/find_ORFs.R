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


#' Creates a GRangeslist of Open Reading Frames mapped to genomic coordinates
#'
#' Input is a Grangeslist of regions to search, together with a DNAStringSet with
#' fastaSequence in same order as the grl.
#' @param grl \code{\link[GenomicRanges]{GRangesList}} of sequences
#'  to search for orfs, in Genomic coordinates
#' @param fastaSeqs DNA sequences to search for Open Reading Frames,
#'  must be DNAStringSet.
#' @param startCodon string. Default is "ATG".
#' @param stopCodon string. Default is "TAA|TAG|TGA".
#' @param longestORF bolean. Default FALSE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least
#'  START + 8*3 [bp] + STOP.
#' @return A GRangesList of ORFs.
#' @export
find_in_frame_ORFs <- function(grl, fastaSeqs, startCodon =  "ATG",
                               stopCodon = "TAA|TAG|TGA", longestORF = F,
                               minimumLength = 0 ){

  if (class(grl) != "GRangesList") stop("Invalid type of grl, must be GRangesList.")
  if (length(fastaSeqs) == 0) stop("Fasta sequences had length 0")
  if (length(fastaSeqs) != length(grl)) stop("Fasta seqs and grl must have same length!")
  result <- ORFs_as_List(fastaSeqs = as.matrix(as.character( fastaSeqs)),
                                 startCodon = startCodon,stopCodon = stopCodon,
                                 longestORF = longestORF,
                                 minimumLength = minimumLength)

  return(map_to_GRanges(grl, result))
}
