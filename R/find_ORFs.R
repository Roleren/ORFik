#' Returns start codon definitions
#'
#' According to:
#' <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
#' index.cgi?chapter=tgencodes#SG1>
#' ncbi genetic code number for translation. This version is a cleaned up
#' version, unknown indices removed.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of START sites separatd with "|".
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


#' Returns stop codon definitions
#'
#' According to:
#' <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
#' index.cgi?chapter=tgencodes#SG1>
#' ncbi genetic code number for translation. This version is a cleaned up
#' version, unknown indices removed.
#'
#' @param transl_table numeric.  NCBI genetic code number for translation.
#' @return A string of STOP sites separatd with "|".
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
#' Find all Open Reading Frames (ORFs) on the simple input sequences
#' in ONLY 5'- 3' direction (+), but within all three possible reading frames.
#' Do not use findORFs for mapping to full chromosomes,
#' then use \code{\link{findMapORFs}}!
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
#' @param seqs (DNAStringSet or character vector) - DNA/RNA sequences to search
#' for Open Reading Frames. Can be both uppercase or lowercase. Easiest call to
#' get seqs if you want only regions from a fasta/fasta index pair is:
#' seqs = ORFik:::txSeqsFromFa(grl, faFile), where grl is a GRanges/List of
#' search regions and faFile is a \code{\link{FaFile}}.
#' @param startCodon (character vector) Possible START codons to search for.
#' Check \code{\link{startDefinition}} for helper function.
#' @param stopCodon (character vector) Possible STOP codons to search for.
#' Check \code{\link{stopDefinition}} for helper function.
#' @param longestORF (logical) Default TRUE. Keep only the longest ORF per
#' unique stopcodon: (seqname, strand, stopcodon) combination, Note: Not longest
#' per transcript! You can also use function
#' \code{\link{longestORFs}} after creation of ORFs for same result.
#' @param minimumLength (integer) Default is 0. Which is START + STOP = 6 bp.
#' Minimum length of ORF, without counting 3bps for START and STOP codons.
#' For example minimumLength = 8 will result in size of ORFs to be at least
#' START + 8*3 (bp) + STOP = 30 bases. Use this param to restrict search.
#' @return (IRangesList) of ORFs locations by START and STOP sites
#' grouped by input sequences. In a list of sequences, only the indices of
#' the sequences that had ORFs will be returned, e.g. 3 sequences where only
#' 1 and 3 has ORFs, will return size 2 IRangesList with names c("1", "3").
#' If there are a total of 0 ORFs, an empty IRangesList will be returned.
#' @export
#' @family findORFs
#' @examples
#' ## Simple examples
#' findORFs("ATGTAA")
#' findORFs("ATGTTAA") # not in frame anymore
#'
#' findORFs("ATGATGTAA") # only longest of two above
#' findORFs("ATGATGTAA", longestORF = FALSE) # two ORFs
#'
#' findORFs(c("ATGTAA", "ATGATGTAA")) # 1 ORF per transcript
#'
#' ## Get DNA sequences from ORFs
#' seq <- DNAStringSet(c("ATGTAA", "AAA", "ATGATGTAA"))
#' names(seq) <- c("tx1", "tx2", "tx3")
#' orfs <- findORFs(seq, longestORF = FALSE)
#'
#' #  you can get sequences like this:
#' gr <- unlist(orfs, use.names = TRUE)
#' gr <- GRanges(seqnames = names(seq)[as.integer(names(gr))],
#'  ranges(gr), strand = "+")
#' # Give them some proper names:
#' names(gr) <- paste0("ORF_", seq.int(length(gr)), "_", seqnames(gr))
#' orf_seqs <- getSeq(seq, gr)
#' orf_seqs
#' # Save as .fasta (orf_seqs must be of type DNAStringSet)
#' # writeXStringSet(orf_seqs, "orfs.fasta")
#' ## Reading from file and find ORFs
#' #findORFs(readDNAStringSet("path/to/transcripts.fasta"))
findORFs <- function(seqs, startCodon =  startDefinition(1),
                     stopCodon = stopDefinition(1), longestORF = TRUE,
                     minimumLength = 0) {

  if (is.null(seqs) || length(seqs) == 0)
    stop("Fasta sequences had length 0 or is NULL")
  if (is.character(seqs) & substr(seqs[1], 1, 1) %in%
      c("a", "t", "c", "g", "n")) {
    startCodon <- tolower(startCodon)
    stopCodon <- tolower(stopCodon)
  }

  result <- orfs_as_List(fastaSeqs = as.character(seqs, use.names = FALSE),
                         startCodon = startCodon, stopCodon = stopCodon,
                         minimumLength = minimumLength)

  if (longestORF) return(longestORFs(split(IRanges(result$orf[[1]],
                                                  result$orf[[2]]),
                                          result$index)))
  return(split(IRanges(result$orf[[1]], result$orf[[2]]), result$index))
}


#' Find ORFs and immediately map them to their genomic positions.
#'
#' This function can map spliced ORFs.
#' It finds ORFs on the sequences of interest, but returns relative positions to
#' the positions of `grl` argument. For example, `grl` can be exons
#' of known transcripts (with genomic coordinates), and `seq` sequences of
#' those transcripts, in that case, this function will return
#' genomic coordinates of ORFs found on transcript sequences.
#'
#' This function assumes that `seq` is in widths relative to `grl`,
#' and that their orders match. 1st seq is 1st grl object, etc.
#'
#' See vignette for real life example.
#' @param grl (\code{\link{GRangesList}}) of sequences
#'  to search for ORFs, probably in genomic coordinates
#' @inheritParams findORFs
#' @param groupByTx logical (default: FALSE), should output GRangesList be
#' grouped by exons per ORF (TRUE) or by orfs per transcript (FALSE)?
#' @return A GRangesList of ORFs.
#' @export
#' @family findORFs
#' @examples
#' # First show simple example using findORFs
#' # This sequence has ORFs at 1-9 and 4-9
#' seqs <- DNAStringSet("ATGATGTAA") # the dna transcript sequence
#' findORFs(seqs)
#' # lets assume that this sequence comes from two exons as follows
#' # Then we need to use findMapORFs instead of findORFs,
#' #  for splicing information
#' gr <- GRanges(seqnames = "1", # chromosome 1
#'               ranges = IRanges(start = c(21, 10), end = c(23, 15)),
#'               strand = "-", #
#'               names = "tx1") #From transcript 1 on chr 1
#' grl <- GRangesList(tx1 = gr) # 1 transcript with 2 exons
#' findMapORFs(grl, seqs) # ORFs are properly mapped to its genomic coordinates
#'
#' grl <- c(grl, grl)
#' names(grl) <- c("tx1", "tx2")
#' findMapORFs(grl, c(seqs, seqs))
#' # More advanced example and how to save sequences found in vignette
#'
findMapORFs <- function(grl, seqs, startCodon = startDefinition(1),
                        stopCodon = stopDefinition(1), longestORF = TRUE,
                        minimumLength = 0, groupByTx = FALSE){
  validGRL(class(grl))
  if (is.null(seqs) || length(seqs) == 0)
    stop("Fasta sequences had length 0 or is NULL")
  if (length(seqs) != length(grl))
    stop("Fasta seqs and grl must have same length!")

  result <- orfs_as_List(fastaSeqs = as.character(seqs, use.names = FALSE),
                         startCodon = startCodon, stopCodon = stopCodon,
                         minimumLength = minimumLength)

  result <- split(IRanges(result$orf[[1]], result$orf[[2]]), result$index)
  if (longestORF) result <- longestORFs(result)

  return(mapToGRanges(grl, result, groupByTx))
}


#' Finds Open Reading Frames in fasta files.
#'
#' Should be used for procaryote genomes or transcript sequences as fasta.
#' Makes no sence for eukaryote whole genomes, since those contains splicing
#' (use findMapORFs for spliced ranges).
#' Searches through each fasta header and reports all ORFs found for BOTH
#' sense (+) and antisense strand (-) in all frames. Name of the header will
#' be used as seqnames of reported ORFs.
#' Each fasta header is treated separately, and name of the sequence will
#' be used as seqname in returned GRanges object. This supports circular
#' genomes.
#'
#' Remember if you have a fasta file of transcripts (transcript coordinates),
#' delete all negative stranded ORFs afterwards by:
#' orfs <- orfs[strandBool(orfs)] # negative strand orfs make no sense then.
#' Seqnames are created from header by format: >name info, so name must be
#' first after "biggern than" and space between name and info.
#' Also make sure your fasta file is valid (no hidden spaces etc),
#'  as this might break the coordinate system!
#' @param filePath (character) Path to the fasta file. Can be both uppercase or
#' lowercase. Or a already loaded R object of either types:
#' "BSgenome" or "DNAStringSet" with named sequences
#' @inheritParams findORFs
#' @param is.circular (logical) Whether the genome in filePath is circular.
#' Prokaryotic genomes are usually circular. Be carefull if you want to
#' extract sequences, remember that seqlengths must be set, else it does not
#' know what last base in sequence is before loop ends!
#' @return (GRanges) object of ORFs mapped from fasta file. Positions are
#' relative to the fasta file.
#' @export
#' @importFrom BSgenome getSeq
#' @family findORFs
#' @examples
#' # location of the example fasta file
#' example_genome <- system.file("extdata", "genome.fasta", package = "ORFik")
#' orfs <- findORFsFasta(example_genome)
#' # To store ORF sequences (you need indexed genome .fai file):
#' fa <- FaFile(example_genome)
#' names(orfs) <- paste0("ORF_", seq.int(length(orfs)), "_", seqnames(orfs))
#' orf_seqs <- getSeq(fa, orfs)
#' # You sequences (fa), needs to have isCircular(fa) == TRUE for it to work
#' # on circular wrapping ranges!
#'
#' # writeXStringSet(DNAStringSet(orf_seqs), "orfs.fasta")
findORFsFasta <- function(filePath, startCodon =  startDefinition(1),
                          stopCodon = stopDefinition(1), longestORF = TRUE,
                          minimumLength = 0, is.circular = FALSE) {

  if (is(filePath, "character")) {
    filePath <- path.expand(filePath)
    if (!file.exists(filePath)) stop("'file' does not exist, check working dir!",
                                     "If you wanted to pass character sequences",
                                     "from R, convert to DNAStringSet!")
    gr <- findORFs_fasta(as.character(readDNAStringSet(filePath), use.names = TRUE),
                         startCodon, stopCodon, minimumLength,
                         is.circular)
  } else if (class(filePath) %in% c("DNAStringSet")) {
    if (length(filePath) == 0) return(GRanges())
    if (is.null(names(filePath))) stop("Sequences does not have names, name them first!")

    gr <- findORFs_fasta(as.character(filePath, use.names = TRUE),
                         startCodon, stopCodon, minimumLength,
                         is.circular)
  } else if (class(filePath) %in% c("BSgenome")) {
    message("Finding ORFs for all chromosomes in BSgenome")
    gr <- findORFs_fasta(as.character(BSgenome::getSeq(filePath), use.names = TRUE),
                         startCodon, stopCodon, minimumLength,
                         is.circular)
  }else stop("filePath must either be file path or DNAString/DNAStringSet!")

  if (longestORF) {
    gr <- longestORFs(gr)
  }
  if (is.circular) {
    isCircular(gr) <- rep(TRUE, length(seqlevels(gr)))
  }
  return(gr)
}

#' Find upstream ORFs from transcript annotation
#'
#' Procedure:
#' 1. Create a new search space starting with the 5' UTRs.
#' 2. Redefine TSS with CAGE if wanted.
#' 3. Add the whole of CDS to search space to allow uORFs going into cds.
#' 4. find ORFs on that search space.
#' 5. Filter out wrongly found uORFs, if CDS is included. The CDS,
#'  alternative CDS, uORFs starting within the CDS etc.
#'
#' From default a filtering process is done to remove "fake" uORFs, but only if
#' cds is included, since uORFs that stop on the stop codon on the CDS is not
#' a uORF, but an alternative cds by definition, etc.
#' @inheritParams findMapORFs
#' @param fa a \code{\link{FaFile}}. With fasta sequences corresponding to
#' fiveUTR annotation. Usually loaded from the genome of an organism with
#' fa = ORFik:::findFa("path/to/fasta/genome")
#' @inheritParams uORFSearchSpace
#' @return A GRangesList of uORFs, 1 granges list element per uORF.
#' @export
#' @family findORFs
#' @examples
#' # Load annotation
#' txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                          package = "GenomicFeatures")
#' \dontrun{
#'  txdb <- loadTxdb(txdbFile)
#'  fiveUTRs <- loadRegion(txdb, "leaders")
#'  cds <- loadRegion(txdb, "cds")
#'  if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'    # Normally you would not use a BSgenome, but some custom fasta-
#'    # annotation you  have for your species
#'    findUORFs(fiveUTRs, BSgenome.Hsapiens.UCSC.hg19::Hsapiens, "ATG",
#'              cds = cds)
#'  }
#' }
findUORFs <- function(fiveUTRs, fa, startCodon = startDefinition(1),
                      stopCodon = stopDefinition(1), longestORF = TRUE,
                      minimumLength = 0, cds = NULL,
                      cage = NULL, extension = 1000, filterValue = 1,
                      restrictUpstreamToTx = FALSE, removeUnused = FALSE) {

  uorfSpace <- uORFSearchSpace(fiveUTRs, cage, extension,
                               filterValue, restrictUpstreamToTx, removeUnused,
                               cds)
  seqs <- txSeqsFromFa(uorfSpace, fa)
  uorfs <- findMapORFs(uorfSpace, seqs, startCodon, stopCodon, longestORF,
                       minimumLength, groupByTx = FALSE)
  if(!is.null(cds)) uorfs <- filterUORFs(uorfs, cds)
  return(uorfs)
}
