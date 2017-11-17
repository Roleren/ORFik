
#' Creates GRangesList from the results of get_all_ORFs_as_GRangesList and
#'  a list of group indeces
#'
#' @param grl GRangesList. A GRangesList of the original sequences that gave the orfs
#' @param result List. A List of the results of finding uorfs
#' List syntax is: result[1] contain grouping indeces,
#' result[2] countains two columns of start and stops
#' @return A GRangesList of ORFs.
#' @export
#' @import GenomicRanges
#' @import GenomicFeatures
map_to_GRangesList <- function(grl, result) {
  # Create GRanges object from result tx ranges
  gr = GRanges(seqnames = as.character(names(grl[result$index])),
                   ranges = IRanges(start = unlist(result$orf[1]),
                                    end = unlist(result$orf[2])),
                   strand =as.character(phead(strand(grl[result$index]),1L )))
  names(gr) = names(grl[result$index])
  
  #map from transcript, remove duplicates, remove hit columns
  #syntax for mapping:-> Seqnames(gr) == names(GRL),
  genomicCoordinates = mapFromTranscripts(x =  gr, transcripts =  grl)
  genomicCoordinates = genomicCoordinates[names(gr[genomicCoordinates$xHits]) == names(genomicCoordinates)]
  rm(gr)
  # old version -> genomicCoordinates = genomicCoordinates[!duplicated(genomicCoordinates$xHits)]
  
  genomicCoordinates$xHits = NULL
  genomicCoordinates$transcriptsHits = NULL
  names(genomicCoordinates) = names(grl[result$index])

  #map from transcript, remove duplicates, remove hit columns
  newGRL = split(genomicCoordinates,result$index )
  names(newGRL) = unique(names(genomicCoordinates))

  #Split by exons and create new exon names
  unlNEW = unlist(newGRL, use.names = F)
  unlGRL = unlist(grl[names(newGRL)])
  ol = findOverlaps(query = unlNEW, subject = unlGRL)
  c = ol[names(unlNEW[ol@from]) == names(unlGRL[ol@to])]
  #resize n to c size
  N = unlNEW[c@from]
  ff = unlGRL[c@to]

  # add end of first to first, start of second to second copy
  dups = duplicated(c@from)
  start(N[dups]) = start(ff[dups])
  dupsR = duplicated(rev(c@from))
  rN = rev(N)
  end(rN[dupsR]) = end(rev(ff)[dupsR])
  N = rev(rN)

  #Relist last time !!!CHECK IF +/- is correct
  #Get grouping t by names
  l = Rle(names(N))
  t = unlist(lapply(1:length(l@lengths),function(x){ rep(x,l@lengths[x])}))
  froms = c@from

  # Fast pre initialized for loop for uorf exon names
  Inds =  rep(1, length(N))
  for(x in 2:length(N)){
    if(t[x] != t[x-1]){
      Inds[x] = 1
    }else{
      if(froms[x] != froms[x-1]){
        Inds[x] = Inds[x-1]+1
      }else{
        Inds[x] = Inds[x-1]
      }
    }
  }
  # create orf names and return as grl
  N$names = paste(names(N),"_",Inds, sep = "")
  #choose split by transcript, or orf exons ?
  newGRL = split(N,t)
  names(newGRL) = unique(names(N))

  return(newGRL)
}


#' Creates GRangesList of Open Reading Frames mapped to genomic coordinates
#' Input is a Grangeslist of regions to search, together with a DNAStringSet with
#' fastaSequence in same order as the grl.
#' @param grl GRangesList of sequences to search for orfs, in Genomic coordinates
#' @param fastaSeqs DNA sequences to search for Open Reading Frames, must be DNAStringSet.
#' @param startCodon string. Default is "ATG".
#' @param stopCodon string. Default is "TAA|TAG|TGA".
#' @param longestORF bolean. Default FALSE. Defines whether pick longest ORF only.
#' When FALSE will report all open reaidng frames, even overlapping small ones.
#' @param minimumLength numeric. Default is 0.
#' For example minimumLength = 8 will result in size of ORFs to be at least
#'  START + 8*3 [bp] + STOP.
#' @return A GRangesList of ORFs.
#' @export
get_all_ORFs_as_GRangesList <- function(grl, fastaSeqs,
				 startCodon =  "ATG",stopCodon = "TAA|TAG|TGA",
				  longestORF = F,minimumLength = 0 ){

  result = get_all_ORFs_as_List( fastaSeqs = as.matrix(as.character(fastaSeqs)),
                                 startCodon = startCodon,stopCodon = stopCodon,
                                            longestORF = longestORF,
                                            minimumLength = minimumLength)
  newGRL = map_to_GRangesList(grl,result)
  return(newGRL)
}
