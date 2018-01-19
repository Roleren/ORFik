#' make GRangesList from GRanges, grouped by names or another column(other)
#' ig. if GRanges should be grouped by gene, give gene column as other
#' @param gr a GRanges object
#' @param other a vector of names to group, no 2 groups can have same name
#' @return a GRangesList named after names(Granges) if other is NULL, else
#' names are from unique(other)
#' @export
groupGRangesBy <- function(gr,other = NULL){
  if (class(gr) != "GRanges") stop("gr must be GRanges Object")
  if(is.null(other)){ # if not using other
    if (is.null(names(gr))) stop("gr object have no names")
    l <- S4Vectors::Rle(names(gr))
  } else { # else use other
    if (length(gr) != length(other))
      stop(" in GroupGRangesByOther: lengths of gr and other does not match")
    l <- S4Vectors::Rle(other)
  }
  grouping <- unlist(lapply(1:nrun(l), function(x){ rep(x, runLength(l)[x])}))
  grl <- split(gr, grouping)
  if(is.null(other)){
    names(grl) <- unique(names(gr))
  } else{
    names(grl) <- unique(other)
  }

  return(grl)
}

#' get list of widths per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
widthPerGroup = function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(sum(width(grl)))
  } else {
    return(as.integer(sum(width(grl))))
  }
}

#' get list of seqnames per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
seqnamesPerGroup <- function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(seqnames(phead(grl, 1L)))
  } else {
    return(as.character(seqnames(phead(grl, 1L))))
  }
}

#' A faster reimplementation of GenomicRanges::sort() for GRangesList
#' @param grl a GRangesList
#' @param decreasing should the first in each group have max(start(group))
#'   ->T or min-> default(F) ?
#' @importFrom data.table as.data.table
gSort <- function(grl, decreasing = F){
  DT <- as.data.table(grl)
  if (decreasing){
    asgrl <- makeGRangesListFromDataFrame(
      DT[order(group, -start)], split.field = "group",
      names.field = "group_name", keep.extra.columns = T)
  } else {
    asgrl <- makeGRangesListFromDataFrame(
      DT[order(group, start)], split.field = "group",
      names.field = "group_name", keep.extra.columns = T)
  }
  names(asgrl) <- names(grl)
  return(asgrl)
}
#' sorts a GRangesList object, uses a faster sort than GenomicRanges::sort(),
#' which works poorly for > 10k groups
#' @param grl a GRangesList
#' @param ignore.strand a boolean, should minus strands be sorted from highest
#'  to lowest(T)
#' @export
sortPerGroup <- function(grl, ignore.strand = F){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (!ignore.strand){
    indexesPos <- which(strandPerGroup(grl, F) == "+")

    grl[indexesPos] <- gSort(grl[indexesPos])
    grl[!indexesPos] <- gSort(grl[!indexesPos], decreasing = T)

    return(grl)
  } else {
    return(gSort(grl))
  }
}

#' get list of strands per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
strandPerGroup <- function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(strand(phead(grl, 1L))) # S4Vectors::phead() wrapper
  }else{
    return(as.character(strand(phead(grl, 1L))))
  }
}

#' get first exon per GRangesList group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList
firstExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(phead(grl, 1L)) # S4Vectors::phead() wrapper
}

#' get first exon per GRangesList group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList,
#' @return a GRangesList of last exons per group
lastExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(ptail(grl,1L)) # S4Vectors::ptail() wrapper
}

#' get first start per granges group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
firstStartPerGroup <- function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(start(firstExonPerGroup(grl)))
  }else{
    return(as.integer(start(firstExonPerGroup(grl))))
  }
}

#' get first end per granges group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
firstEndPerGroup <- function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(end(firstExonPerGroup(grl)))
  }else{
    return(as.integer(end(firstExonPerGroup(grl))))
  }
}

#' get last end per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
lastExonEndPerGroup = function(grl,keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(end(lastExonPerGroup(grl)))
  }else{
    return(as.integer(end(lastExonPerGroup(grl))))
  }
}

#' get last end per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
lastExonStartPerGroup = function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(start(lastExonPerGroup(grl)))
  }else{
    return(as.integer(start(lastExonPerGroup(grl))))
  }
}

#' get list of the number of exons per group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
numExonsPerGroup <- function(grl, keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")

  if (!is.null(names(unlist(grl, use.names = F)))){
    exonsPerGroup <- runLength(S4Vectors::Rle(names(unlist(grl,
      use.names = F))))
  } else if (!is.null(names(unlist(grl, use.names = T)))){
    exonsPerGroup <- runLength(S4Vectors::Rle(names(unlist(grl,
      use.names = T))))
  } else stop("no names to group exons")
  if(keep.names){
    if(is.null(names(grl))) stop("grl must have names if keep.names == T")
    names(exonsPerGroup) = names(grl)
  }
  return(exonsPerGroup)
}

#' seqnames cleanup
#' For many datasets, the fa file and the gtf file have different naming
#' This functions tries to fix the naming to the GRanges standard
#' chrX instead of X chr1 instead of 1 etc..
#' @param grl a GRangesList
fixSeqnames <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  temp <- unlist(grl)
  seqnamesTransformed <- as.character(seqnames(temp))
  indexes <- which(nchar(seqnamesTransformed) < 6)
  temp <- temp[indexes]
  seqlevels(temp) <- sub(replacement = "chrY", pattern = "Y", seqlevels(temp))
  seqlevels(temp) <- sub(replacement = "chrX", pattern = "X", seqlevels(temp))
  seqlevels(temp) <- as.character(unique(seqnames(unlist(temp))))
  return(groupGRangesBy(temp))
}

#' make a meta column with exon ranks
#' @param grl a GRangesList
#' @param byTranscript if ORfs are by transcript, check duplicates
makeExonRanks <- function(grl, byTranscript = F){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")

  if (byTranscript){
    oldNames <- names(grl)
    names(grl) <- 1:length(grl)
  }
  l <- Rle(names(unlist(grl)))
  t <- unlist(lapply(1:nrun(l), function(x) {
    rep(x, runLength(l)[x])
  }))

  Inds <- rep(1, length(t))
  if (!byTranscript){
    for (x in 2:length(t)) {
      if (t[x] == t[x - 1]) {
        Inds[x] <- Inds[x - 1] + 1
      }
    }
  } else {
    for (x in 2:length(t)) {
      if (t[x] != t[x - 1]) {
        if (oldNames[t[x]] == oldNames[t[x] - 1]){
          Inds[x] <- Inds[x - 1] + 1
        }
      } else{
        Inds[x] <- Inds[x - 1]
      }
    }
  }
  return(Inds)
}

#' Make ORF names per orf, grl must be grouped by transcript
#' If a list of orfs are grouped by transcripts, but does not have
#' ORF names, then create them and return the new GRangesList
#' @param grl a GRangesList
#' @return a GRangesList now with ORF names
makeORFNames <- function(grl){
  ranks <- makeExonRanks(grl, byTranscript = T)

  asGR <- unlist(grl, use.names = F)
  if (is.null(names(asGR))) asGR <- unlist(grl, use.names = T)
    asGR$names <- paste0(names(asGR), "_", ranks)
  return(groupGRangesBy(asGR))
}

#' Tile a GRangeslist by 1
#' This is not supported originally by GenomicRanges
#' @param grl a GRangesList object
#' @return a GRangesList grouped by original group, tiled to 1
#' @export
tile1 <- function(grl){
  ORFs <- unlist(grl, use.names = F)
  if(is.null(names(ORFs))) stop("unlisted grl have not names,\n
                                try: names(unlist(grl, use.names = F))")
  if(sum(duplicated(names(ORFs)))){
    if(!is.null(ORFs$names)){
      names(ORFs) <- ORFs$names
    } else{
        nametest <- unlist(grl, use.names = T)
        dups <- sum(duplicated(names(nametest)))
        if(dups != sum(duplicated(names(ORFs)))){
          stop("duplicated ORF names,\n
                need a column called 'names' that are unique,\n
                or change names of groups")
        }
      ORFs$names <- names(ORFs)
    }
  }

  tilex <- tile(ORFs, width =  1L)

  names(tilex) <- ORFs$names
  unl <- unlist(tilex, use.names = T)
  tilex <- groupGRangesBy(unl)
  return(sortPerGroup(tilex))
}

#' map genomic to transcript coordinates by reference
#' @param grl a GRangesList of ranges within the reference,
#'  grl must have column called names that gives grouping for result
#' @param reference a GrangesList of ranges
#'  that include and are bigger or equal to grl
#'  ig. cds is grl and gene can be reference
asTX <- function(grl, reference){
  return(pmapToTranscripts(grl, reference[OrfToTxNames(grl)]))
}

#' get transcript sequence from a GrangesList and a faFile
#' @param grl a GRangesList object
#' @param faFile FaFile used to find the transcripts,
#' @param is.sorted a speedup, if you know the ranges are sorted
txSeqsFromFa <- function(grl, faFile, is.sorted = F){
  if(class(faFile) != "FaFile") stop("only FaFile is valid input")
  if(!isSorted) grl <- sortPerGroup(grl)
  return(extractTranscriptSeqs(faFile, transcripts = grl))
}

#' Creates window around GRanged object.
#'
#' It creates window of window_size around input ranges eg.
#' for GRanges starting at 100-100 and window_size of 3 will give
#' 97-103
#' @param GRanges_obj GRanges object of your CDSs start or stop postions.
#' @param window_size Numeric. Default 30. What size of the window to consider.
#' @return A GRanges object of resized by window_size sites.
#' @export
#' @import GenomicRanges
#' @examples
#' window_resize(GRanges(Rle(c('1'), c(4)),
#'                       IRanges(c(100, 200, 200, 100), width=c(1, 1, 1, 1)),
#'                       Rle(strand(c('+', '+', '-', '-')))),
#'              window_size = 50)
#'
window_resize <- function(GRanges_obj, window_size = 30) {
  GRanges_obj <- promoters(GRanges_obj, upstream = window_size, downstream = window_size + 1)
  return(GRanges_obj)
}

#' source bioconductor
#' helperfunction for quick update of bioconductor packages
#' @param packages either NULL if only source and no update/install
#'  or "all" if you want to update all your bioconductor packages
#'  or c(package1, package2, ...)
#'  for specific packages as a character vector
sourceBioC <- function(packages = NULL){
  source("https://bioconductor.org/biocLite.R")
  if(!is.null(packages)){
    if(packages == "all"){
      biocLite()
    } else{
        biocLite(packages)
    }
  }
}
