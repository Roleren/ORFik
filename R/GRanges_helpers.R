#' make GRangesList from GRanges, grouped by names or another column(other)
#' ig. if GRanges should be grouped by gene, give gene column as other
#' @param gr a GRanges object
#' @param other a vector of names to group by
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
  names(grl) <- unique(names(gr))
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
      DT[order(group, -start)],split.field = "group",
      names.field = "group_name", keep.extra.columns = T)
  } else {
    asgrl <- makeGRangesListFromDataFrame(
      DT[order(group, start)],split.field = "group",
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
    indexesMin <- which(strandPerGroup(grl, F) == "-")

    grl[indexesPos] <- gSort(grl[indexesPos])
    grl[indexesMin] <- gSort(grl[indexesMin], decreasing = T)

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

#' get first exon per granges group
#' @param grl a GRangesList
firstExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(phead(grl, 1L)) # S4Vectors::phead() wrapper
}

#' get last exon per granges group
#' @param grl a GRangesList
lastExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(ptail(grl,1L)) # S4Vectors::ptail() wrapper
}

#' get first start per granges group
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

#' get last end per granges group
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
lastExonEndPerGroup = function(grl,keep.names = T){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  if (keep.names){
    return(end(ptail(grl, 1L))) # S4Vectors::ptail() wrapper
  }else{
    return(as.integer(end(ptail(grl, 1L))))
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
  } else if (!is.null(names(unlist(grl)))){
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
#' @param grl a GRangesList
makeORFNames <- function(grl){

  ranks <- makeExonRanks(grl, byTranscript = T)

  asGR <- unlist(grl, use.names = F)
  if (is.null(names(asGR))) asGR <- unlist(grl, use.names = T)
    asGR$names <- paste0(names(asGR), "_", ranks)
  return(groupGRangesBy(asGR))
}
