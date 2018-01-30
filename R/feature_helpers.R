#' get length of leaders ordered after oldTxNames
#'
#' Normally only a helper function for ORFik
#' @param fiveUTRs a GRangesList object of leaders
#' @param oldTxNames a character vector of names to group fiveUTRs by.
findCageUTRFivelen <- function(fiveUTRs, oldTxNames){
  newfiveprimeLen <- widthPerGroup(fiveUTRs)
  return(newfiveprimeLen[match(oldTxNames,names(newfiveprimeLen))])
}

#' Get transcript lengths
#'
#' A helper function for easy length retrieval
#' @param Gtf a TxDb object of a gtf file
#'  NB! only add this if you dont use tx argument
#' @param changedFiveUTRs a GRangesList object of leaders.
#'  NB! only add this if you used cage data or other things to change the
#'  leaders, therefor we need it to update transcript lengths.
txLen <- function(Gtf = NULL, changedFiveUTRs = NULL){
  tx_len_temp <- transcriptLengths(Gtf)[,c("tx_name","tx_len")]
  tx_len <- tx_len_temp[,"tx_len"]

  if(!is.null(changedFiveUTRs)){
    if(!is.null(Gtf)){
      new5Length <- findCageUTRFivelen(changedFiveUTRs, tx_len_temp$tx_name)
      tx_len_temp <- transcriptLengths(Gtf, T, T, T)
      tx_len_temp$utr5_len <- new5Length
      tx_len <- tx_len_temp$utr5_len +
      tx_len_temp$cds_len + tx_len_temp$utr3_len
    }
  }
  names(tx_len) <- tx_len_temp$tx_name
  return(tx_len)
}

#' Helper Function to check valid RNA input
#' @param class, the given class of RNA object
checkRNA <- function(class){
  if(is.null(class)){
    message("No RNA added, skipping feature te and fpkm of RNA,\n
            also RibosomeReleaseScore will also be not normalized best way possible.")
  } else {
    if(class != "GAlignments" & class != "GRanges"){
      stop("RNA must be either GAlignments or GRanges")
    }
  }
}

#' Helper Function to check valid RFP input
#' @param class, the given class of RFP object
checkRFP <- function(class){
  if(class != "GAlignments" & class != "GRanges"){
    stop("RFP must be either GAlignments or GRanges")
  }
}

#' Helper function to check valid combinations of extension and cageFiveUTRs
#' @param extension a numeric/integer to reassign 5' utrs.
#' @param cageFiveUTRs a GRangesList, if you used cage-data to extend 5' utrs,
validExtension <- function(extension, cageFiveUTRs){
  if(is.null(extension)){
    stop("please specify extension, to avoid bugs\n
                              ,if you did not use cage, set it to 0,\n
                              standard cage extension is 1000")
  } else if(!is.numeric(extension) && !is.integer(extension)){
      stop("extension must be numeric or integer")
  }
  if(extension != 0 && class(cageFiveUTRs) != "GRangesList"){
    stop("if extension is not 0, then cageFiveUTRs must be defined")
  }
}
