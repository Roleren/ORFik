

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

#' Subset GRanges to get coverage. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’
#' @param cov A coverage object from coverage()
#' @param y GRanges object for which coverage should be extracted
#' @return numeric vector of coverage of input GRanges object
#' @export
#' @examples
#' #subset_coverage(x)
#'
subset_coverage <- function(cov, y) {
  cov1 <- cov[[as.vector(seqnames(y)[1])]]
  return(as.vector(cov1[ranges(y)]))
}

#' Calucalte Entropy value of input reads overlapping per group in grl
#' For example reads per orf.
#' The entropy interval per group is a real number in the interval (0:1)
#' Where 0 is no variance in reads over group.
#' @export
#' @param grl a GRangesList that the reads can map to
#' @param reads a GAlignment object, ig. from ribo seq, or rna seq
#' @return A numeric vector containing entropy per ORF
#'
entropy <- function(grl, reads) {
  # Get count list of overlaps
  tileBy1 <- tile1(grl)
  countsTile <- countOverlaps(unlist(tileBy1, use.names = F), reads)
  names <- names(countsTile)
  names(countsTile) <- NULL
  countList <- split(countsTile, names)
  names(countList) <- NULL
  countList <- IRanges::RleList(countList)

  # generate the entropy variables
  sums <- sum(countList)
  if(sum(sums) == 0){ # no variance in countList, 0 entropy
    return(rep(0, length(tileBy1)))
  }
  N <- unlist(sums)
  # get indeces where entropy is not 0
  validIndeces <- N > 0
  N <- N[validIndeces] # <- sums not 0
  countList <- countList[validIndeces]

  lengths <- lengths(countList)
  reg_len <- rep(3, length(lengths)) # tuplets per orf
  L <- unlist(lengths)
  # which ORFs should not have tuplet length of 3
  bool <- N <= L
  floors <- floor(L / N)
  reg_len[bool] <- floors[bool]
  reg_counts <- as.integer(L / reg_len) # number of tuplets per orf


  # Need to reassign variables to equal length,
  # to be able to vectorize
  # h: The sequences we make the tuplets per orf,
  # if h[1] is: c(1,2,3,4,5,6) and reg_len[1] is: c(3,3)
  # you get int_seqs: ->  1: c(1,2,3 , 4,5,6) <- 2 triplets
  h <- lapply(reg_counts,function(x){
    0:(x-1)
  })

  runLengths <- lengths(IRanges::RleList(h))
  # get sequence from valid indeces
  indeces <- 1:length(runLengths)
  reg_len <- lapply(indeces, function(x){
    rep(reg_len[x], runLengths[x])
  })
  unlh <- unlist(h, use.names = F)
  unlreg_len <- unlist(reg_len, use.names = F)

  # make sequences for reads, start -> stop
  # gives a triplet reading, 1:3, 3,6
  acums <- L
  for(i in 2:length(L)){
    acums[i] <-acums[i-1] + acums[i]
  }
  unlacums <- unlist(lapply(1:(length(runLengths)-1), function(x){
    rep(acums[x], runLengths[x+1])
  }), use.names = F)
  unlacums <- unlist(c(rep(1,runLengths[1]),unlacums))

  which_reads_start <- (unlacums + unlh * unlreg_len)
  which_reads_end <- (unlh * unlreg_len + (unlreg_len + unlacums -1))
  # the actual triplets ->
  int_seqs <- lapply(1:length(which_reads_start), function(x){
    which_reads_start[x]:which_reads_end[x]
  })
  # group int_seqs
  int_grouping <- unlist(lapply(1:length(reg_counts), function(x){
    rep(x, reg_counts[x])
  }), use.names = F)
  int_seqs <- split(int_seqs, int_grouping)
  names(int_seqs) <- NULL
  N <- unlist(lapply(indeces, function(x){
    rep(N[x], reg_counts[x])
  }), use.names = F)

  intcountList <- IntegerList(countList)
  unlintcount <- unlist(unlist(intcountList, use.names = F),
                        use.names = F)


  # get the assigned tuplets per orf, usually triplets
  triplets <- lapply(1:length(int_seqs), function(x){
    unlintcount[int_seqs[[x]]]
  })
  tripletSums <- unlist(lapply(trip, function(x){
    sum(x)
  }), use.names = F)

  # entropy function, interval 0:1 real number
  Xi <-  tripletSums / N
  validXi <- !is.nan(Xi * log2(Xi))
  Hx <- rep(0, length(validXi))
  Hx[validXi] <- Xi[validXi] * log2(Xi[validXi])

  MHx <-rep(0, length(validXi))

  reg_counts<- lapply(indeces, function(x){
    rep(reg_counts[x], runLengths[x])
  })
  unlrg_counts <- unlist(reg_counts, use.names = F)
  MHx <- 1/unlrg_counts * log2(1 / unlrg_counts)

  # sum the mhx to groups

  grouping <- unlist(lapply(indeces, function(x){
    rep(x, runLengths[x])
  }), use.names = F)

  Hx <- sum(NumericList(split(Hx, grouping)))
  MHx <- sum(NumericList(split(MHx, grouping)))

  entropy <- rep(0.0, length(tileBy1))

  # non 0 entropy values set to HX / MHX
  entropy[validIndeces] <- Hx / MHx
  entropy[is.na(entropy)] <- 0
  return(entropy)
}


#' Calculate Coverage of the input reads.
#' What is percent of the vector covered by the values higher than zero?
#'
#' @export
#' @param countsOver A numerc vector
#' @return numeric value in between 0 and 1
#'
calculateCoverage <- function(countsOver) {
  return(sum(countsOver >= 1)/length(countsOver))
}

#' Make a score for each ORFs start region
#' The closer the sequence is to the kozak sequence
#' The higher the score, based on simplification of PWM
#' score system: 4 upstream, 5 downstream of start
#' CACCATGGC, 1+3+1+2, skip ATG, +2+1 = 10
#' CGCCATGGC, 1+!2+1+2, skip ATG, +2+1 = 9
#' Inspired by experimental bit values for each position
#' @param grl a GRangesList grouped by ORF
#' @param faFile a faFile from the fasta file
#' @param species which species to use, currently only support human
kozacSequenceScore <- function(grl, faFile  = NULL, species = "human"){
  firstExons <- firstExonPerGroup(grl)
  kozakLocation <- promoters(firstExons, upstream = 4, downstream = 5)

  sequences <- as.character(txSeqsFromFa(grl, faFile))
  names(sequences) <- NULL
  scores <- rep(0, length(sequences))
  if(species == "human"){
    # split strings and relist as letters of 9 rows
    subSplit <- strsplit(sequences, split = "")
    mat <- matrix(unlist(subSplit), nrow = 9) #this will not when ATG is on start of chr
    pos1 <- as.character(mat[1,])
    pos2 <- as.character(mat[2,])
    pos3 <- as.character(mat[2,])
    pos4 <- as.character(mat[2,])
    pos8 <- as.character(mat[8,])
    pos9 <- as.character(mat[9,])


    match <- grepl(x = pos1, pattern = "C|G|A")
    scores[match] <- scores[match] + 1

    match <- pos2 == "A"
    scores[match] <- scores[match] + 5
    match <- pos2 == "G"
    scores[match] <- scores[match] + 2
    match <- pos2 == "C"
    scores[match] <- scores[match] + 1

    match <- grepl(x = pos3, pattern = "C|A")
    scores[match] <- scores[match] + 1

    match <- pos4 == "C"
    scores[match] <- scores[match] + 2
    match <- pos4 == "G"
    scores[match] <- scores[match] + 1
    match <- pos4 == "A"
    scores[match] <- scores[match] + 1

    match <- pos8 == "G"
    scores[match] <- scores[match] + 2
    match <- pos8 == "A"
    scores[match] <- scores[match] + 1
    match <- pos8 == "T"
    scores[match] <- scores[match] + 1

    match <- grepl(x = pos9, pattern = "C|A")
    scores[match] <- scores[match] + 1


    return(scores)
  }
  else stop("other species are not supported")
}
