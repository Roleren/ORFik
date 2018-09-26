#' Returns start definitions
#'
#' According to:
#' <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
#' index.cgi?chapter=tgencodes#SG1>
#' ncbi genetic code number for translation. This version is a cleaned up
#' version, unknown indices removed.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of START sies separatd with "|".
#' @family findORFs
#' @export
#' @examples
#' startDefinition
#' startDefinition(1)
#'
startDefinition <- function(transl_table) {
  STARTdef <- c("ATG|TTG|CTG", #1 The Standard Code
                "ATT|ATC|ATA|ATG|GTG", #2 The Vertebrate Mitochondrial Code
                "ATA|ATG", #3 The Yeast Mitochondrial Code
                "TTA|TTG|CTG|ATT|ATC|ATA|ATG|GTG",
                #4 Mold/Protozoan/Coelenterate Mitochondrial/Mycoplasma
                "TTG|ATT|ATC|ATA|ATG|GTG",
                #5 The Invertebrate Mitochondrial Code
                "ATG", #6 The Ciliate, Dasycladacean and Hexamita Nuclear Code
                "ATG|GTG", #9 The Echinoderm and Flatworm Mitochondrial Code
                "ATG", #10 The Euplotid Nuclear Code
                "TTG|CTG|ATT|ATC|ATA|ATG|GTG",
                #11 The Bacterial, Archaeal and Plant
                "CTG|ATG", #12 The Alternative Yeast Nuclear Code
                "TTG|ATA|ATG|GTG", #13 The Ascidian Mitochondrial Code
                "ATG", #14 The Alternative Flatworm Mitochondrial Code
                "ATG", #16 Chlorophycean Mitochondrial Code
                "ATG|GTG", #21 Trematode Mitochondrial Code
                "ATG", #22 Scenedesmus obliquus Mitochondrial Code
                "ATT|ATG|GTG", #23 Thraustochytrium Mitochondrial Code
                "TTG|CTG|ATG|GTG", #24 Pterobranchia Mitochondrial Code
                "TTG|ATG|GTG", #25 Candidate Division SR1 and G.bacteria Code
                "CTG|ATG") #26 Pachysolen tannophilus Nuclear Code
  return(STARTdef[transl_table])
}


#' Returns stop definitions
#'
#' According to:
#' <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
#' index.cgi?chapter=tgencodes#SG1>
#' ncbi genetic code number for translation. This version is a cleaned up
#' version, unknown indices removed.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of STOP sies separatd with "|".
#' @family findORFs
#' @export
#' @examples
#' stopDefinition
#' stopDefinition(1)
#'
stopDefinition <- function(transl_table) {
  STOPdef <- c("TAA|TAG|TGA", #1 The Standard Code
               "TAA|TAG|AGA|AGG", #2 The Vertebrate Mitochondrial Code
               "TAA|TAG", #3 The Yeast Mitochondrial Code
               "TAA|TAG", #4 Mold/Protozoan/Coelenterate Mitochondrial/Mplasma/
               "TAA|TAG", #5 The Invertebrate Mitochondrial Code
               "TGA", #6 The Ciliate, Dasycladacean and Hexamita Nuclear Code
               "TAA|TAG", #9 The Echinoderm and Flatworm Mitochondrial Code
               "TAA|TAG", #10 The Euplotid Nuclear Code
               "TAA|TAG|TGA", #11 The Bacterial, Archaeal and Plant Code
               "TAA|TAG|TGA", #12 The Alternative Yeast Nuclear Code
               "TAA|TAG", #13 The Ascidian Mitochondrial Code
               "TAG", #14 The Alternative Flatworm Mitochondrial Code
               "TGA", #16 Chlorophycean Mitochondrial Code
               "TAA|TAG", #21 Trematode Mitochondrial Code
               "TCA|TAA|TGA", #22 Scenedesmus obliquus Mitochondrial Code
               "TTA|TAA|TAG|TGA", #23 Thraustochytrium Mitochondrial Code
               "TAA|TAG", #24 Pterobranchia Mitochondrial Code
               "TAA|TAG", #25 Candidate Division SR1 and Gracilibacteria Code
               "TAA|TAG") #26 Pachysolen tannophilus Nuclear Code
  return(STOPdef[transl_table])
}


#' Find Open Reading Frames.
#'
#' Find all Open Reading Frames (ORFs) on the input sequences
#' in ONLY 5'- 3' direction (+), but within all three possible reading frames.
#' For each sequence of the input vector \code{\link{IRanges}} with START and
#' STOP positions (inclusive) will be returned as
#' \code{\link{IRangesList}}. Returned coordinates are relative to the
#' input sequences.
#'
#' If you want antisence strand too, do:
#' \code{
#' #positive strands
#' pos <- findORFs(seqs)
#' #negative strands (DNAStringSet only if character)
#' neg <- findORFs(reverseComplement(DNAStringSet(seqs)))
#  #merge together
#' relist(c(GRanges(pos, strand = "+"), GRanges(neg, strand = "-")),
#'        skeleton = merge(pos, neg))
#' }
#'
#' @param seqs (DNAStringSet or character) DNA sequences to search for Open
#' Reading Frames. Can be both uppercase or lowercase.
#' @param startCodon (character) Possible START codons to search for. Check
#' \code{\link{startDefinition}} for helper function.
#' @param stopCodon (character) Possible STOP codons to search for. Check
#'  \code{\link{stopDefinition}} for helper function.
#' @param longestORF (logical) Default FALSE. When TRUE will only report the
#' longest ORFs per seqlevel (1 in + and in in - direction),
#' all smaller ORFs will be ignored.
#' When FALSE will report all possible ORFs in all three reading frames.
#' If you want longest ORF per unique (seqname, strand, stopcodon) combination,
#' use \code{\link{longestORFs}}
#' @param minimumLength (integer) Default is 0. Minimum length of ORF, without
#' counting 3bp for START and STOP codons. For example minimumLength = 8 will
#' result in size of ORFs to be at least START + 8*3 (bp) + STOP = 30 bases.
#' Use this param to restrict search.
#' @return (IRangesList) of ORFs locations incuding START and STOP codons
#' grouped by input seqeunces.
#' @export
#' @family findORFs
#' @seealso [findMapORFs()], [findORFsFasta()],
#' [startDefinition()], [stopDefinition()]
#' @examples
#' findORFs("ATGTAA")
#' findORFs("ATGTTAA") # not in frame anymore
#'
#' findORFs("ATGATGTAA") # two ORFs
#' findORFs("ATGATGTAA", longestORF = TRUE) # only longest of two above
#'
#' findORFs(c("ATGTAA", "ATGATGTAA"))
#'
findORFs <- function(
  seqs, startCodon =  startDefinition(1), stopCodon = stopDefinition(1),
  longestORF = FALSE, minimumLength = 0 ){

  if (is.null(seqs) || length(seqs) == 0)
    stop("Fasta sequences had length 0 or is NULL")
  if (is.character(seqs) & substr(seqs[1], 1,1) %in%
      c("a", "t", "c", "g", "n")) {
    startCodon <- tolower(startCodon)
    stopCodon <- tolower(stopCodon)
  }

  result <- orfs_as_List(fastaSeqs = as.character(seqs, use.names = FALSE),
                         startCodon = startCodon, stopCodon = stopCodon,
                         longestORF = longestORF,
                         minimumLength = minimumLength)
  return(split(IRanges(result$orf[[1]], result$orf[[2]]), result$index))
}


#' Find ORFs and immediately map them to their genomic positions.
#'
#' Finds ORFs on the sequences of interest, but returns relative positions to
#' the positions of `grl` argument. For example, `grl` can be exons
#' of known transcripts (with genomic coordinates), and `seq` sequences of
#' those transcripts, in that case, [findMapORFs()] will return
#' genomic coordinates of ORFs found on transcript sequences.
#'
#' This function assumes that `seq` is in widths relative to `grl`,
#' and that their orders match. 1st seq is 1st grl object, etc.
#'
#' @param grl (\code{\link{GRangesList}}) of sequences
#'  to search for ORFs, probably in genomic coordinates
#' @inheritParams findORFs
#' @param groupByTx logical (T), should output GRangesList be grouped by
#' transcripts (T) or by ORFs (F)?
#' @return A GRangesList of ORFs.
#' @export
#' @family findORFs
#' @seealso [findORFs()], [findORFsFasta()],
#' [startDefinition()], [stopDefinition()]
#' @examples
#' # This sequence has ORFs at 1-9 and 4-9
#' seqs <- c("ATGATGTAA") # the dna sequence
#' findORFs(seqs)
#' # lets assume that this sequence comes from two exons as follows
#' gr <- GRanges(seqnames = rep("1", 2), # chromosome 1
#'               ranges = IRanges(start = c(21, 10), end = c(23, 15)),
#'               strand = rep("-", 2), names = rep("tx1", 2))
#' grl <- GRangesList(tx1 = gr)
#' findMapORFs(grl, seqs) # ORFs are properly mapped to its genomic coordinates
#'
#' grl <- c(grl, grl)
#' names(grl) <- c("tx1", "tx2")
#' findMapORFs(grl, c(seqs, seqs))
#'
findMapORFs <- function(
  grl, seqs, startCodon =  startDefinition(1), stopCodon = stopDefinition(1),
  longestORF = FALSE, minimumLength = 0, groupByTx = TRUE){
  validGRL(class(grl))
  if (is.null(seqs) || length(seqs) == 0)
    stop("Fasta sequences had length 0 or is NULL")
  if (length(seqs) != length(grl))
    stop("Fasta seqs and grl must have same length!")

  result <- orfs_as_List(fastaSeqs = as.character(seqs, use.names = FALSE),
                         startCodon = startCodon, stopCodon = stopCodon,
                         longestORF = longestORF,
                         minimumLength = minimumLength)
  return(mapToGRanges(grl, result, groupByTx))
}


#' Finds Open Reading Frames in fasta files.
#'
#' Searches through each fasta header and reports all ORFs found for BOTH
#' sense (+) and antisense strand (-) in all frames. Name of the header will
#' be used as seqnames of reported ORFs.
#' Each fasta header is treated separately, and name of the sequence will
#' be used as seqname in returned GRanges object. This supports circluar
#' genomes.
#'
#' Remember if you have a fasta file of transcripts (transcript coordinates),
#' delete all negative stranded ORFs afterwards by:
#' orfs <- orfs[strandBool(orfs)] # negative strand orfs make no sense then
#' Seqnames are created from header by format: >name info, so name must be
#' first after "biggern than" and space between name and info.
#' @param filePath (character) Path to the fasta file. Can be both uppercase or
#' lowercase.
#' @inheritParams findORFs
#' @param is.circular (logical) Whether the genome in filePath is circular.
#' Prokaryotic genomes are usually circular. Be carefull if you want to
#' extract sequences, remember that seqlengths must be set, else it does not
#' know what last base in sequence is before loop ends!
#' @return (GRanges) object of ORFs mapped from fasta file. Positions are
#' relative to the fasta file.
#' @export
#' @family findORFs
#' @seealso [findORFs()], [findMapORFs()],
#' [startDefinition()], [stopDefinition()]
#' @examples
#' # location of the example fasta file
#' example_genome <- system.file("extdata", "genome.fasta", package = "ORFik")
#' findORFsFasta(example_genome)
#'
findORFsFasta <- function(
  filePath, startCodon =  startDefinition(1), stopCodon = stopDefinition(1),
  longestORF = FALSE, minimumLength = 0, is.circular = FALSE) {

  if (!is(filePath, "character"))
    stop("'filepath' must be of type character.")
  filePath <- path.expand(filePath)
  if(!file.exists(filePath)) stop("'file' does not exist, check working dir!")
  gr <- findORFs_fasta(filePath, startCodon, stopCodon, longestORF,
                       minimumLength, is.circular)
  if (is.circular) {
    isCircular(gr) <- rep(TRUE, length(seqlevels(gr)))
  }
  return(gr)
}
