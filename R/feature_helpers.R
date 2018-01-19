#' get length of leaders ordered after oldTxNames
#' Normally only a helper function for ORFik
#' @param fiveUTRs a GRangesList object of leaders
#' @param oldTxNames a character vector of names to group fiveUTRs by.
findCageUTRFivelen <- function(fiveUTRs, oldTxNames){
  newfiveprimeLen <- widthPerGroup(fiveUTRs)
  return(newfiveprimeLen[match(oldTxNames,names(newfiveprimeLen))])
}

#' get transcript lengths
#' a helper function for easy length retrieval
#' @param Gtf a TxDb object of a gtf file
#' @param changedFiveUTRs a GRangesList object of leaders.
#'  NB! only add this if you used cage data or other things to change the
#'  leaders, therefor we need it to update transcript lengths.
TxLen <- function(Gtf, changedFiveUTRs = NULL){
  tx_len_temp <- transcriptLengths(Gtf)[,c("tx_name","tx_len")]
  tx_len <- tx_len_temp[,"tx_len"]

  if(!is.null(fiveUTRs)){
    new5Length <- findCageUTRFivelen(changedFiveUTRs, tx_len_temp$tx_name)
    tx_len_temp <- transcriptLengths(Gtf, T, T, T)
    tx_len_temp$utr5_len <- new5Length
    tx_len <- tx_len_temp$utr5_len +
      tx_len_temp$cds_len + tx_len_temp$utr3_len
  }
  names(tx_len) <- tx_len_temp$tx_name
  return(tx_len)
}
