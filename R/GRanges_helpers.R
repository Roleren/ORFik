
#' make GRangesList from GRanges, grouped by names or another column(other)
#' @description ig. if GRanges should be grouped by gene,
#'  give gene column as other
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
    indecesPos <- strandBool(grl)

    grl[indecesPos] <- gSort(grl[indecesPos])
    grl[!indecesPos] <- gSort(grl[!indecesPos], decreasing = T)

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

#' helper function to get a logical list of True/False, if GRangesList group have
#' + strand = T, if - strand = F
#' Also checks for * strands, so a good check for bugs
#' @param grl a GRangesList or GRanges object
strandBool <- function(grl){
  if(class(grl) == "GRanges"){
    posIndeces <- as.character(strand(grl)) == "+"
  } else {
    posIndeces <- strandPerGroup(grl, F) == "+"
  }

  sums <- sum(posIndeces) + sum(!posIndeces)
  if(is.na(sums)){
    stop("could not get strands from grl object,\n
          most likely NULL object was passed.")
  }
  if(sums != length(grl)){
    stop("grl contains * strands, set them to either + or -")
  }
  return(posIndeces)
}

#' get first exon per GRangesList group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList
firstExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(phead(grl, 1L)) # S4Vectors::phead() wrapper
}

#' get last exon per GRangesList group
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a GRangesList,
#' @return a GRangesList of last exons per group
lastExonPerGroup <- function(grl){
  if (class(grl) != "GRangesList") stop("grl must be GRangesList Object")
  return(ptail(grl, 1L)) # S4Vectors::ptail() wrapper
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

#' get last start per granges group
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
#'  @export
asTX <- function(grl, reference){
  orfNames <- OrfToTxNames(grl)
  if(sum(orfNames %in% names(reference)) != length(orfNames)){
    stop("not all references are present, so can not map to transcripts.")
  }
  return(pmapToTranscripts(grl, reference[orfNames]))
}

#' get transcript sequence from a GrangesList and a faFile
#' @param grl a GRangesList object
#' @param faFile FaFile used to find the transcripts,
#' @param is.sorted a speedup, if you know the ranges are sorted
txSeqsFromFa <- function(grl, faFile, is.sorted = F){
  if(class(faFile) != "FaFile") stop("only FaFile is valid input")
  if(!is.sorted) grl <- sortPerGroup(grl)
  return(extractTranscriptSeqs(faFile, transcripts = grl))
}

#' Reassign the start positions of the first exons per group in grl
#' @description make sure your grl is sorted, since start of "-" strand
#'  objects should be the
#'  max end in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a GRangesList object
#' @param newStarts an integer vector of same length as grl, with new start values
assignFirstExonsStartSite <- function(grl, newStarts){
  if(length(grl) != length(newStarts)) stop("length of grl and newStarts \n
                                            are not equal!")
  posIndeces <- strandBool(grl)

  dt <- as.data.table(grl)
  dt[!duplicated(dt$group),]$start[posIndeces] <- newStarts[posIndeces]
  dt[!duplicated(dt$group),]$end[!posIndeces] <- newStarts[!posIndeces]

  ngrl <- GenomicRanges::makeGRangesListFromDataFrame(dt,
    split.field = "group", names.field = "group_name", keep.extra.columns = T)
  names(ngrl) <- names(grl)

  return(ngrl)
}
#' Reassign the stop positions of the last exons per group
#' @description make sure your grl is sorted,
#'  since stop of "-" strand objects should be the
#'  min start in group, use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a GRangesList object
#' @param newStops an integer vector of same length as grl,
#'  with new start values
assignLastExonsStopSite <- function(grl, newStops){
  if(length(grl) != length(newStops)) stop("length of grl and newStops \n
                                           are not equal!")
  posIndeces <- strandBool(grl)

  dt <- as.data.table(grl)
  idx = dt[, .I[.N], by=group]
  dt[idx$V1]$end[posIndeces] <- newStops[posIndeces]
  dt[idx$V1]$start[!posIndeces] <- newStops[!posIndeces]
  ngrl <- GenomicRanges::makeGRangesListFromDataFrame(dt,
    split.field = "group", names.field = "group_name", keep.extra.columns = T)
  names(ngrl) <- names(grl)

  return(ngrl)
}
#' get rest of objects downstream
#' @description per group get the part downstream of position
#'  defined in downstreamOf
#'  downstreamOf(tx, ORFik:::ORFStopSites(cds, asGR = F))
#'  will return the 3' utrs per transcript as GRangesList,
#'  usually used for interesting
#'  parts of the transcripts, like upstream open reading frames(uorf).
#'  downstreamOf +/- 1 is start/end site
#'  of transformed tx's, depending on strand
#' @param tx a GRangesList, usually of Transcripts to be changed
#' @param downstreamOf a vector of integers, for each group in tx, where
#' is the new start point of first valid exon.
downstreamOfPerGroup <- function(tx, downstreamOf){
  posIndeces <- strandBool(tx)
  posEnds <- end(tx[posIndeces])
  negEnds <- start(tx[!posIndeces])
  posDown <- downstreamOf[posIndeces]
  negDown <- downstreamOf[!posIndeces]
  pos <- posEnds > posDown
  neg <- negEnds < negDown
  posTx <- tx[posIndeces][pos]
  negTx <- tx[!posIndeces][neg]
  downTx <- tx
  downTx[posIndeces] <- posTx
  downTx[!posIndeces] <- negTx
  #check if anyone hits boundary, set those to boundary
  if(anyNA(strandPerGroup(downTx, F))){
    boundaryHits <- which(is.na(strandPerGroup(downTx, F)))
    downTx[boundaryHits] <- firstExonPerGroup(tx[boundaryHits])
    ir <- IRanges(start = downstreamOf[boundaryHits],
                    end = downstreamOf[boundaryHits])
    irl <- split(ir, 1:length(ir))
    names(irl) <- names(tx[boundaryHits])
    ranges(downTx[boundaryHits]) <- irl
  }
  return(assignFirstExonsStartSite(downTx, downstreamOf))
}
#' get rest of objects upstream
#' @description per group get the part upstream of position
#'  defined in upstreamOf
#'  upstream(tx, ORFik:::ORFStopSites(cds, asGR = F))
#'  will return the 5' utrs per transcript, usually used for interesting
#'  parts of the transcripts, like upstream open reading frames(uorf).
#'  downstreamOf +/- 1 is start/end site
#'  of transformed tx's, depending on strand
#' @param tx a GRangesList, usually of Transcripts to be changed
#' @param upstreamOf a vector of integers, for each group in tx, where
#'  is the new start point of first valid exon.
upstreamOfPerGroup <- function(tx, upstreamOf){
  posIndeces <- strandBool(tx)
  posStarts <- start(tx[posIndeces])
  negStarts <- end(tx[!posIndeces])
  posGrlStarts <- upstreamOf[posIndeces]
  negGrlStarts <- upstreamOf[!posIndeces]
  pos <- posStarts < posGrlStarts
  neg <- negStarts > negGrlStarts
  # need to fix pos/neg with possible cage extensions
  outside <- which(sum(pos) == 0)
  pos[outside] = T
  posTx <- tx[posIndeces]
  posTx[outside] <- firstExonPerGroup(posTx[outside])
  outside <- which(sum(neg) == 0)
  neg[outside] = T
  negTx <- tx[!posIndeces]
  negTx[outside] <- firstExonPerGroup(negTx[outside])

  posTx <- posTx[pos]
  negTx <- negTx[neg]
  tx[posIndeces] <- posTx
  tx[!posIndeces] <- negTx

  return(assignLastExonsStopSite(tx, upstreamOf))
}

#' Extend the leaders tss.
#'
#' Either if you have the 5' utr sequences or the
#'  transcript sequences to extend,
#'  any way the 5' part will be extended.
#'  Remember to sort it, if not sorted beforehand.
#'  use ORFik:::sortPerGroup(grl) to get sorted grl.
#' @param grl a GRangesList of 5' utrs or transcripts..
#' @param extension an integer, how much to extend the leaders.
#'  Or a GRangesList where start / stops by strand are the positions
#'  to use as new starts.
#' @param cds If you want to extend 5' leaders downstream, to catch
#'  upstream ORFs going into cds, include it. It will add first
#'  cds exon to grl matched by names.
#'  Do not add for transcripts, as they are already included.
extendLeaders <- function(grl, extension = 1000, cds = NULL){
  if(class(extension) == "numeric" && length(extension) == 1){
    posIndeces <- strandBool(grl)
    promo <- promoters(unlist(firstExonPerGroup(grl), use.names = F),
      upstream = extension)
    newStarts <- rep(NA, length(grl))
    newStarts[posIndeces] <- as.integer(start(promo[posIndeces]))
    newStarts[!posIndeces] <- as.integer(end(promo[!posIndeces]))
  } else if(is.grl(class(grl))){
    starts <- ORFStartSites(extension)
    changedGRL <-downstreamOfPerGroup(grl[names(extension)], starts)
    return(changedGRL)
  } else {
    stop("extension must either be an integer, or a GRangesList")
  }

  extendedLeaders <- assignFirstExonsStartSite(grl, newStarts)
  if(is.null(cds)) return (extendedLeaders)
  return(addFirstCdsOnLeaderEnds(
    extendedLeaders, cds))
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
  GRanges_obj <- promoters(GRanges_obj, upstream = window_size,
                           downstream = window_size + 1)
  return(GRanges_obj)
}

#' Subset GRanges to get desired frame. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’ @export
#' @param x A tiled to size of 1 GRanges object
#' @param frame A numeric indicating which frame to extract
#' @return GRanges object reduced to only first frame
#' @examples
#' #subset_to_frame(x, 1)
#'
subset_to_frame <- function(x, frame){
  if(as.vector(strand(x) == "+")[1]){
    x[seq(frame, length(x), 3)]
  }else{
    x[seq(length(x) + 1 - frame, 1, -3)]
  }
}

#' Subset GRanges to get stop codons. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#' @param x A tiled to size of 1 GRanges object
#' @return GRanges object reduced to only stop codon
#' @export
#' @examples
#' #subset_to_stop(x)
#'
subset_to_stop <- function(x){
  if(as.vector(strand(x))[1] == "+"){
    x[c(length(x) - 3, length(x) - 4, length(x) - 5)]
  } else {
    x[c(4, 5, 6)]
  }
}

#' Helper function to check for GRangesList
#' @param class the class you want to check if is GRL
#' @return a boolean
is.grl <- function(class){
  if(class == "GRangesList"){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Helper Function to check valid GRangesList input
#' @param class as character vector the given class of
#'  supposed GRangesList object
#' @param type a character vector, is it gtf, cds, 5', 3', for messages.
#' @param checNULL should NULL classes be checked and return indeces of these?
validGRL <- function(class, type, checkNULL = FALSE){
  if(length(class) != length(type)) stop("not equal length of classes\n
                                         and types, see validGRL")
  if(checkNULL){
    indeces <-"NULL" == class

    class <- class[!indeces]
    if(length(class) == 0){
      return(rep(T, length(type)))
    }
    type <- type[!indeces]

  }
  for(classI in 1:length(class)){
    if(!is.grl(class[classI])){
      messageI <- paste(type[classI], "must be given
          and be type GRangesList")
      stop(messageI)
    }
  }
  if(checkNULL){
    return(indeces)
  }
}

#' source bioconductor
#' @description helperfunction for quick update of bioconductor packages
#' @param packages either NULL if only source and no update/install
#'  or "all" if you want to update all your bioconductor packages
#'  or c(package1, package2, ...)
#'  for specific packages as a character vector
sourceBioc <- function(packages = NULL){
  source("https://bioconductor.org/biocLite.R")
  if(!is.null(packages)){
    if(packages == "all"){
      biocLite()
    } else {
        biocLite(packages)
    }
  }
}
