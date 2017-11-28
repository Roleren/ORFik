#' Converts different type of files to Granges
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An imported bed-file, to convert to GRanges
#' @param bed6 If bed6, no meta column is added
bedToGR <- function(x, bed6 = TRUE){

  if (!bed6){
    gr <- GRanges(x[, 1], IRanges(x[, 2] - 1, x[, 3]))
    return(gr)
  }
  starts <- ifelse(x[, 6] == "+", x[, 2]-1, x[, 2])
  ends <- ifelse(x[, 6] == "-", x[, 3], x[, 3] - 1)
  gr <- GRanges(x[, 1], IRanges(starts, ends))
  strand(gr) <- x[, 6]
  score(gr) <- x[, 5]
  if (ncol(x) > 6) mcols(gr) <- x[, 7 : ncol(x)]
  return(gr)
}


#' Get Cage-Data From a file-path
#' @param filePath The location of the cage-file
#' @importFrom data.table fread
cageFromFile <- function(filePath){

  if (.Platform$OS.type == "unix") {
    if (file.exists(filePath)) {
      if (gsub(pattern = ".*\\.", "", filePath) == "gzip"){
      rawCage <- bedToGR(as.data.frame(fread(paste("gunzip -c",
        filePath), sep = "\t")))
      } else if (gsub(pattern = ".*\\.", "", filePath) == "bed"){
        rawCage <- bedToGR(as.data.frame(fread(filePath, sep = "\t")))
      } else {stop("Only bed and gzip formats are supported for filePath")}
    } else {stop("Filepath specified does not name existing file.") }
  } else {
    stop("Only unix operating-systems currently support filePath,
         use cage as GRanges argument instead.") }

  message("Loaded cage-file successfully")
  return(rawCage)
}

#' Filter peak of cage-data by value
#' @param rawCage The raw cage-data
#' @param filterValue The number of counts(score) to filter on for a tss to pass as hit
filterCage <- function(rawCage, filterValue = 1){

  if (is.null(rawCage$score)) stop("Found no score column in the cageData-file, bed standard is column 5")
  filteredCage <- rawCage[rawCage$score > filterValue,] #filter on score

  return(filteredCage)
}

#' Check that seqnames of fiveUTRs and cage uses same standard, i.g chr1 vs 1.
#' @param filteredrawCageData Cage-data to check seqnames in
#' @param fiveUTRs The 5' leader sequences as GRangesList
matchSeqnames <- function(filteredCage, fiveUTRs){
  fiveSeqlevels <- seqlevels(unlist(fiveUTRs, use.names = FALSE))
  cageSeqlevels <- seqlevels(filteredCage)
  if (length(grep(pattern = "chr", fiveSeqlevels)) > 0 &&
     length(grep(pattern = "chr", cageSeqlevels)) == 0){
    warning("seqnames use different chromosome naming conventions,
            trying to fix them")
    regexNormalChr <- '(^[a-zA-Z])*([0-9]+)' # <- chr1, chr2, not chrX, chrY etc.
    normalChr <- paste0("chr", grep(regexNormalChr,
                                    cageSeqlevels, value = TRUE))
    normalChrInd <- grep(regexNormalChr, cageSeqlevels)

    for(i in normalChrInd){
      seqlevels(filteredCage)[i] <- sub(regexNormalChr, normalChr[i], cageSeqlevels[i])
    }
    if (length(grep("chrY", fiveSeqlevels)) == 0)
      seqlevels(filteredCage) <- sub("chrY", "Y", seqlevels(filteredCage))
    if (length(grep("chrX", fiveSeqlevels)) == 0)
      seqlevels(filteredCage) <- sub("chrX", "X", seqlevels(filteredCage))
    if (length(grep("chrM", fiveSeqlevels)) == 0)
      seqlevels(filteredCage) <- sub("chrM", "MT", seqlevels(filteredCage))
  }
  return(filteredCage)
}

#' When reassigning Transcript start site, often you want to add downstream too.
#' This is a simple way to do that
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cds If you want to extend 5' leaders downstream, to catch uorfs going into cds, include it.
addFirstCdsOnLeaderEnds <- function(fiveUTRs, cds){

  if (length(cds) == 0){
    warning("cds is empty, returning without using it.")
    return(fiveUTRs)
  }
  if (is.null(names(cds))){
    warning("cds have no names, returning without using it.")
    return(fiveUTRs)
  }
  matchingNames <- names(fiveUTRs) %in% names(cds)
  if (sum(matchingNames) - length(names(fiveUTRs)) != 0){ # <- check valid names
    warning("not all cds names matches fiveUTRs names, returning without using cds.")
    return(fiveUTRs)
  }

  cdsForUTRs <- cds[names(fiveUTRs)] # get only the ones we need

  firstExons <- phead(cdsForUTRs, 1L) #select first in every, they must be sorted!
  gr <- unlist(firstExons, use.names = FALSE)
  mcols(gr) <- NULL # <- remove all mcols
  gr$exon_id <- NA # <- assign the ones we need
  gr$exon_name <- NA
  gr$exon_rank <- NA

  grl <- relist(gr, firstExons)
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)

  return(reduce(fiveUTRsWithCdsExons))
}

#' Extend first exon of each transcript with length specified
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param extension The number of basses upstream to add on transcripts
extendsTSSexons <- function(fiveUTRs, extension = 1000){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  if (is.null(fiveAsgr$exon_rank))
    stop("fiveUTRs need column called exon_rank, see ?makeTranscriptDbFromGFF")
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  posIDs <- firstExons[strand(firstExons) == "+"]
  minIDs <- firstExons[strand(firstExons) == "-"]
  promo <- promoters(firstExons, upstream = extension)

  start(firstExons[names(posIDs)]) <- start(promo[names(posIDs)])
  end(firstExons[names(minIDs)]) <- end(promo[names(minIDs)])

  return(firstExons)
}

#' Find max peak for each transcript,
#' returns as data.table, without names, but with index
#' @param cageOverlaps The cageOverlaps between cage and extended 5' leaders
#' @param filteredrawCageData The filtered raw cage-data used to reassign 5' leaders
#' @importFrom data.table as.data.table
findMaxPeaks = function(cageOverlaps, filteredrawCageData){

  dt <- as.data.table(filteredrawCageData)
  dt <- dt[from(cageOverlaps)]
  dt$to <- to(cageOverlaps)

  maxPeaks <- dt[, max(score), by = to]
  names(maxPeaks) <- c("to", "score")
  maxPeaks <-  merge(maxPeaks, dt)

  return(maxPeaks[!duplicated(maxPeaks$to)])
}



#' Finds max peaks per trancsript from reads in the cagefile
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cageData The location of the cage-file
#' @param extension The number of basses upstream to add on transcripts
findNewTSS = function(fiveUTRs, cageData, extension){

  shiftedfiveUTRs <- extendsTSSexons(fiveUTRs, extension)
  cageOverlaps <- findOverlaps(query = cageData, subject = shiftedfiveUTRs)
  maxPeakPosition <- findMaxPeaks(cageOverlaps, cageData)
  return(maxPeakPosition)
}

#' After all transcript start sites have been updated from cage, grangeslist back together
#' @param firstExons The first exon of every transcript from 5' leaders
#' @param fiveUTRs The 5' leader sequences as GRangesList
makeGrlAndFilter = function(firstExons, fiveUTRs){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  fiveAsgr[fiveAsgr$exon_rank == 1] <- firstExons
  return(relist(fiveAsgr, fiveUTRs))
}

#' add cage max peaks as new transcript start sites for each 5' leader
#' (*) strands are not supported, since direction must be known.
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param maxPeakPosition The max peak for each 5' leader found by cage
addNewTSSOnLeaders = function(fiveUTRs, maxPeakPosition){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  maxPeakPosition$names <- names(firstExons[maxPeakPosition$to])
  posIDs <- maxPeakPosition$to[maxPeakPosition$strand == "+"]
  minIDs <- maxPeakPosition$to[maxPeakPosition$strand == "-"]

  firstExons[posIDs] <- resize(
    firstExons[posIDs], width = end(firstExons[posIDs]) -
      maxPeakPosition$start[maxPeakPosition$strand == "+"] + 1, fix = "end")
  firstExons[minIDs] <- resize(
    firstExons[minIDs], width = maxPeakPosition$end[maxPeakPosition$strand == "-"] -
      start(firstExons[minIDs]) + 1, fix = "end")
  # Might need an chromosome boundary here? current test show no need.

  return( firstExons )
}


#' Reassign all Transcript start sites(tss's) in 5' leaders, that are found by given cage-file
#' A max peak is defined as new tss if it is within boundary of 5'leader range, specified by extension
#' A max peak must also be greater that the cage peak cutoff specified in filterValue
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cage Either a  filePath for cage-file, or already loaded R-object as GRanges
#' @param extension The maximum number of basses upstream the cage-peak can be from original tss
#' @param filterValue The number of counts(score) to filter on for a tss to pass as hit
#' @param cds If you want to extend 5' leaders downstream, to catch upstream ORFs going into cds, include it.
#' @export
reassignTSSbyCage = function(fiveUTRs, cage, extension = 1000, filterValue = 1, cds = NULL){

  ###Read in cage files
  if (class(fiveUTRs) != "GRangesList") stop("fiveUTRs must be type GRangesList!")
  if (!is.null(cds) & class(cds) != "GRangesList") stop("cds must be type GRangesList!")
  if (class(extension) != "numeric") stop("extension must be numeric!")
  if (class(filterValue) != "numeric") stop("filterValue must be numeric!")
  if (is.null(cage)) stop("Cage can not be NULL")

  if (class(cage) == "character"){ # <- get cage file
    filteredCage <- filterCage(cageFromFile(cage), filterValue) # get the cage data
  } else if (class(cage) == "GRanges"){
    filteredCage <- filterCage(cage, filterValue)
  } else {stop("Cage-file must be either a valid character filepath or GRanges object")}
  # check that seqnames match
  filteredCage <- matchSeqnames(filteredCage, fiveUTRs)

  maxPeakPosition <- findNewTSS(fiveUTRs, filteredCage, extension)
  fiveUTRs <- makeGrlAndFilter(addNewTSSOnLeaders(fiveUTRs, maxPeakPosition), fiveUTRs)
  if(!is.null(cds)) fiveUTRs <- addFirstCdsOnLeaderEnds(fiveUTRs, cds)

  return(fiveUTRs)
}
