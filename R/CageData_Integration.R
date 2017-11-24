#' Converts different type of files to Granges
#' Only Accepts bed files for now
#' @param x A file containing column-based counts
#' @param bed6 If bed6, no meta column is added
toGR=function(x,bed6=TRUE){
  require(GenomicRanges)
  if(!bed6){
    gr<- GRanges(x[,1],IRanges(x[,2]-1,x[,3]))
    return(gr)
  }
  starts<- ifelse(x[,6]=="+",x[,2]-1,x[,2])
  ends<- ifelse(x[,6]=="-",x[,3],x[,3]-1)
  gr<- GRanges(x[,1],IRanges(starts,ends))
  strand(gr)<-x[,6]
  score(gr)<-x[,5]
  if(ncol(x)>6)  mcols(gr)=x[,7:ncol(x)]
  return(gr)
}


#' Get CageData as granges
#' @param dataName The location of the cage-file
#' @param filterValue The number of counts(score) to filter on for a tss to pass as hit
#' @import data.table
getFilteredCageData = function(dataName, filterValue = 1){
  rawCageData <- toGR(as.data.frame(fread(paste("gunzip -c",dataName),sep = "\t")))
  message("Loaded cage-file successfully")
  filteredrawCageData <- rawCageData[rawCageData$score > filterValue,] #filter on score

  #chromosome nameing, need to find a smart way for this, it wont work for all
  seqlevels(filteredrawCageData) <- sub("chrY","Y",seqlevels(filteredrawCageData))
  seqlevels(filteredrawCageData) <- sub("chrX","X",seqlevels(filteredrawCageData))
  seqlevels(filteredrawCageData) <- sub("chrM","MT",seqlevels(filteredrawCageData))
  return(filteredrawCageData)
}

#' When reassigning Transcript start site, often you want to add downstream too.
#' This is a simple way to do that
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cds If you want to extend 5' leaders downstream, to catch uorfs going into cds, include it.
#'
addFirstCdsOnLeaderEnds = function(fiveUTRs, cds){

  if(length(cds) == 0){
    message("cds is empty, returning without using it.")
    return(fiveUTRs)
  }
  if(is.null(names(cds))){
    message("cds have no names, returning without using it.")
    return(fiveUTRs)
  }

  if(tryCatch({length(cds[names(fiveUTRs)]) == 0},error = function(err) return(T))){ # <- check valid names
    message("cds have no matching names, returning without using it.")
    return(fiveUTRs)
  }

  cdsForUTRs <- cds[names(fiveUTRs)] # get only the ones we need

  firstExons <- phead(cdsForUTRs, 1L) #select first in every, they must be sorted!
  gr <- unlist(firstExons,use.names = F)
  gr$cds_id <- NULL; gr$cds_name <- NULL; gr$exon_rank <- NULL
  gr$exon_id <- NA;  gr$exon_name <- NA;  gr$exon_rank <- NA
  grl = relist(gr,firstExons)
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl) # ask gunnar why length != fiveUTRs

  return( reduce(fiveUTRsWithCdsExons) )
}

#' Extend first exon of each transcript with length specified
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param extension The number of basses upstream to add on transcripts
extendsTSSExons = function(fiveUTRs, extension = 1000){
  fiveAsgr <- unlist(fiveUTRs)
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  posIDs <- firstExons[strand(firstExons) == "+"]
  minIDs <- firstExons[strand(firstExons)  == "-"]
  promo <- promoters(firstExons, upstream=extension)

  start(firstExons[names(posIDs)]) <- start(promo[names(posIDs)])
  end(firstExons[names(minIDs)]) <- end(promo[names(minIDs)])

  return(firstExons)
}

#' Find max peak for each transcript,
#' returns as data.table, without names, but with index
#' @param cageOverlaps The cageOverlaps between cage and extended 5' leaders
#' @param filteredrawCageData The filtered raw cage-data used to reassign 5' leaders
#' @import data.table
findMaxPeaks = function(cageOverlaps, filteredrawCageData){
  dt <- as.data.table(filteredrawCageData)
  dt <- dt[from(cageOverlaps)]
  dt$to <- to(cageOverlaps)

  test <- dt[,max(score), by = to]
  names(test) <- c("to","score")
  res <-  merge(  test, dt)

  return(res[!duplicated(res$to)]) #check this line!!!
}



#' Finds max peaks per trancsript from reads in the cagefile
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param filePath The location of the cage-file
#' @param extension The number of basses upstream to add on transcripts
#' @param filterValue The number of counts(score) to filter on for a tss to pass as hit
#' @param cageAsGR Alternative for filePath, if you have cage-data already loaded in R as GRanges
findNewTSS = function(fiveUTRs, filePath = NULL, extension, filterValue, cageAsGR = NULL){
  if(!is.null(filePath)){
    filteredrawCageData <- getFilteredCageData(filePath, filterValue) # get the cage data
  }else if(!is.null(cageAsGR)){
    filteredrawCageData = cageAsGR
  }else{ stop("Either filePath or cageAsGR must be specified!")}

  shiftedfiveUTRs <- extendsTSSExons(fiveUTRs, extension)
  cageOverlaps <- findOverlaps(query = filteredrawCageData,subject = shiftedfiveUTRs)
  maxPeakPosition <- findMaxPeaks(cageOverlaps,filteredrawCageData)
  return(maxPeakPosition)
}

#' After all transcript start sites have been updated from cage, grangeslist back together
#' @param firstExons The first exon of every transcript from 5' leaders
#' @param fiveUTRs The 5' leader sequences as GRangesList
makeGrlAndFilter = function(firstExons, fiveUTRs){
  fiveAsgr <- unlist(fiveUTRs)
  fiveAsgr[fiveAsgr$exon_rank == 1] <- firstExons
  return(relist(fiveAsgr,fiveUTRs))
}

#' add cage max peaks as new transcript start sites for each 5' leader
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param maxPeakPosition The max peak for each 5' leader found by cage
addNewTssOnLeaders = function(fiveUTRs, maxPeakPosition){
  fiveAsgr <- unlist(fiveUTRs)
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  maxPeakPosition$names <- names(firstExons[maxPeakPosition$to])
  posIDs <- maxPeakPosition$to[maxPeakPosition$strand == "+"]
  minIDs <- maxPeakPosition$to[maxPeakPosition$strand == "-"]

  firstExons[posIDs] <- resize(firstExons[posIDs], width = end(firstExons[posIDs]) - maxPeakPosition$start[maxPeakPosition$strand == "+"] + 1, fix = "end")
  firstExons[minIDs] <- resize(firstExons[minIDs], width = maxPeakPosition$end[maxPeakPosition$strand == "-"] -  start(firstExons[minIDs]) + 1, fix = "end")
  # Might need an chromosome boundary here? current test show no need.

  return( firstExons )
}


#' Reassign all Transcript start sites(tss's) in 5' leaders, that are found by given cage-file
#' A max peak is defined as new tss if it is within boundary of 5'leader range, specified by extension
#' A max peak must also be greater that the cage peak cutoff specified in filterValue
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param filePath The location of the cage-file
#' @param extension The maximum number of basses upstream the cage-peak can be from original tss
#' @param filterValue The number of counts(score) to filter on for a tss to pass as hit
#' @param cds If you want to extend 5' leaders downstream, to catch upstream ORFs going into cds, include it.
#' @param cageAsGR Alternative for filePath, if you have cage-data already loaded in R as GRanges
#' @export
reassignTSSbyCage = function(fiveUTRs,filePath = NULL, extension = 1000, filterValue = 1, cds = NULL, cageAsGR = NULL){
  ###Read in cage files
  if(class(fiveUTRs) != "GRangesList") stop("fiveUTRs must be type GRangesList!")
  if(!is.null(cds) & class(cds) != "GRangesList") stop("cds must be type GRangesList!")
  if(class(extension) != "numeric") stop("extension must be numeric!")
  if(class(filterValue) != "numeric") stop("filterValue must be numeric!")
  if(!is.null(cageAsGR) & class(cageAsGR) != "GRanges") stop("cageAsGR must be GRanges!")
  if(is.null(filePath) && is.null(cageAsGR)) stop("Either filePath or cageAsGR must be specified!")
  if(!is.null(filePath) & class(filePath) != "character"){ stop("filePath must be character!")
  }else if(is.null(cageAsGR)){
    message(cat("Using cage file:",filePath))
  }


  maxPeakPosition <- findNewTSS(fiveUTRs, filePath, extension, filterValue, cageAsGR)
  fiveUTRs <- makeGrlAndFilter(addNewTssOnLeaders(fiveUTRs, maxPeakPosition), fiveUTRs)
  if(!is.null(cds)) fiveUTRs <- addFirstCdsOnLeaderEnds(fiveUTRs, cds)

  return(fiveUTRs)
}
