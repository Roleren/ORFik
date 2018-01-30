
#' Creates normalizations of the counts, normally,
#' @description used in Translations efficiency calculations
#' @param counts a list of integer counts per object
#' @param lengthSize a list of integer lengths per object
#' @param librarySize a numeric of size 1, the size of the library
#' @family features
#' @export
#' @return a numeric vector
fpkm <- function(counts, lengthSize, librarySize){
  return((as.numeric(counts) * (10^9)) /
           (as.numeric(lengthSize) * as.numeric(librarySize)))
}

#' Get fpkm of ribo-seq data
#' @param grl a GRangesList object with usually ORFs, but can also be
#'  either leaders, cds', 3' utrs or  ORFs are a special case,
#'  see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @return a numeric vector with the fpkm values
fpkmRFP <- function(grl, RFP, pseudoCount = 0){

  libraryRPF <- length(RFP)

  overlapRFP <- countOverlaps(grl, RFP)

  grl_len <- widthPerGroup(grl, F)

  #normalize by grl lengths
  return(fpkm(overlapRFP, grl_len, libraryRPF) + pseudoCount)
}

#' Get fpkm of ribo-seq data
#' @param grl a GRangesList object with usually ORFs, but can also be
#'  either leaders, cds', 3' utrs or  ORFs are a special case,
#'  see argument tx_len
#' @param RNA rna seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param tx_len the transcript lengths of the transcripts
#'  a named (tx names) vector of integers.
#'  If you have the transcripts as GRangesList,
#'  call ORFik:::widthPerGroup(tx, T)
#'  If not:
#'  Normally you call argument as: tx_len = txLen(Gtf)
#'  fractionLength(grl, tx_len = txLen(Gtf))
#'  See ?ORFik:::txLen
#'  NB!!! if you used cage data, then
#'  the tss for the the leaders have changed,
#'  therefor the tx lengths have changed,
#'  if so, call:
#'  te(grl, RNA, RFP, tx_len = TxLen(Gtf, fiveUTRs))
#'  where fiveUTRs are the changed cageLeaders.
#' @return a numeric vector with the fpkm values
fpkmRNA <- function(grl, RNA, tx_len, pseudoCount = 0){

  libraryRna <- length(RNA)

  overlapRNA <- countOverlaps(grl, RNA)

  if(is.null(names(tx_len))) stop("tx_len have no names")
  tx_len <- tx_len[OrfToTxNames(grl)]
  if(sum(is.na(tx_len))) stop("tx_len have na values, check naming")

  #normalize by tx lengths
  return(fpkm(overlapRNA, tx_len, libraryRna) + pseudoCount)
}


#' Subset GRanges to get coverage.
#' @description GRanges object should be beforehand
#'  tiled to size of 1. This subsetting takes account for strand.
#' @param cov A coverage object from coverage()
#' @param y GRanges object for which coverage should be extracted
#' @family features
#' @export
#' @examples
#' #subset_coverage(x)
#' @return numeric vector of coverage of input GRanges object
subset_coverage <- function(cov, y) {
  cov1 <- cov[[as.vector(seqnames(y)[1])]]
  return(as.vector(cov1[ranges(y)]))
}

#' Calucalte Entropy value of input reads overlapping
#' @description For example reads per orf.
#' The entropy interval per group is a real number in the interval (0:1)
#' Where 0 is no variance in reads over group.
#' @export
#' @param grl a GRangesList that the reads can map to
#' @param reads a GAlignment object, ig. from ribo seq, or rna seq,
#'  can also be GRanges or GRangesList
#' @family features
#' @export
#' @return A numeric vector containing entropy per ORF
#'
entropy <- function(grl, reads) {
  # Get count list of overlaps
  tileBy1 <- tile1(grl)
  countsTile <- countOverlaps(unlist(tileBy1, use.names = F), reads)
  names <- names(countsTile)
  names(countsTile) <- NULL
  countList <- split(countsTile, names)

  countList <- IRanges::RleList(countList)
  countList <- countList[names(grl)]
  names(countList) <- NULL
  # generate the entropy variables
  sums <- sum(countList)
  if(sum(as.numeric(sums)) == 0){ # no variance in countList, 0 entropy
    return(rep(0, length(tileBy1)))
  }
  N <- unlist(sums, use.names = F)
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
  unlacums <- unlist(c(rep(1,runLengths[1]), unlacums))

  which_reads_start <- (unlacums + unlh * unlreg_len)
  which_reads_end <- (unlh * unlreg_len + (unlreg_len + unlacums -1))
  # the actual triplets ->
  int_seqs <- lapply(1:length(which_reads_start), function(x){
    which_reads_start[x]:which_reads_end[x]
  })

  N <- unlist(lapply(indeces, function(x){
    rep(N[x], reg_counts[x])
  }), use.names = F)

  intcountList <- IntegerList(countList)
  unlintcount <- unlist(unlist(intcountList, use.names = F),
                        use.names = F)

  # get the assigned tuplets per orf, usually triplets
  triplets <- lapply(int_seqs, function(x){
    unlintcount[x]
  })
  tripletSums <- unlist(lapply(triplets, function(x){
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
#' @description What is percent of the vector covered by the values higher than zero?
#' See article: 10.1002/embj.201488411
#' @param countsOver A numerc vector
#' @family features
#' @export
#' @return numeric value in between 0 and 1
#'
calculateCoverage <- function(countsOver) {
  return(sum(countsOver >= 1)/length(countsOver))
}

#' Fragment Length Organization Similarity Score
#' @description See article: 10.1016/j.celrep.2014.07.045
#' @param grl a GRangesList object with ORFs
#' @param RFP ribozomal footprints, given as Galignment or GRanges object
#' @param cds a GRangesList of coding sequences
#' @param start usually 26, the start of the floss interval
#' @param end usually 34, the end of the floss interval
#' @family features
#' @export
#' @return a vector of FLOSS of length same as grl
floss <- function(grl, RFP, cds, start = 26, end = 34){
  #TODO: add 0 to the ones that dont hit
  if(start > end) stop("start is bigger than end")
  if(class(RFP) == "GRangesList"){
    stop("RFP must be either GAlignment or GRanges type")
  }
  # for orfs
  overlaps <- findOverlaps(grl, RFP)
  if(class(RFP) == "GRanges"){
    rfpWidth <- width(RFP[to(overlaps)])
  } else {
    rfpWidth <- qwidth(RFP[to(overlaps)])
  }
  rfpPassFilter <- (rfpWidth >= start) & (rfpWidth <= end)
  rfpValidMatch <- rfpWidth[rfpPassFilter]
  ORFGrouping <- from(overlaps)[rfpPassFilter]
  if(sum(as.numeric(ORFGrouping)) == 0){
    return(as.numeric(rep(0, length(grl))))
  }
  whichNoHit <- NULL # which ribo-seq did not hit grl
  if(length(unique(ORFGrouping)) != length(grl)){
    whichNoHit <- S4Vectors::setdiff.Vector(1:length(grl),
                                            unique(ORFGrouping))
  }
  orfFractions <- split(rfpValidMatch, ORFGrouping)
  listing<- IRanges::RleList(orfFractions)
  tableFracs <- table(listing)
  colnames(tableFracs) <- NULL
  orfFractions <- lapply(1:length(listing), function(x){
    tableFracs[x,]/sum(tableFracs[x,])
  })

  # for cds
  overlapsCds <- findOverlaps(cds, RFP)
  if(class(RFP) == "GRanges"){
    rfpWidth <- width(RFP[to(overlapsCds)])
  } else {
    rfpWidth <- qwidth(RFP[to(overlapsCds)])
  }

  rfpPassFilterCDS <- ((rfpWidth >= start) & (rfpWidth <= end))
  rfpValidMatchCDS <- rfpWidth[rfpPassFilterCDS]
  cdsFractions <- split(rfpValidMatchCDS, rfpValidMatchCDS)
  totalLength <- length(rfpValidMatchCDS)
  cdsFractions <- sapply(cdsFractions, function(x){
    length(x)/totalLength
  })
  cdsFractions <- as.double(cdsFractions)
  # floss score ->
  score <- sapply(1:length(orfFractions), function(x){
    sum(abs(orfFractions[[x]] - cdsFractions)) * 0.5
  })
  if(!is.null(whichNoHit)){
    tempScores <- as.numeric(rep(NA, length(grl)))
    tempScores[unique(ORFGrouping)] <- score
    tempScores[whichNoHit] <- 0.
    score <- tempScores
  }
  if(length(score) != length(grl) || anyNA(score)) stop("could not find\n
                                        floss-score for all objects, most\n
                                        likely objects are wrongly annotated.")
  return(score)
}

#' Translational efficiency
#' @description Uses Rna-seq and Ribo-seq to get te of grl
#' is defined as (density of RPFs within ORF)/(RNA expression).
#' See article: 10.1126/science.1168978
#' @param grl a GRangesList object with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs is a special case, see argument tx_len
#' @param RNA rna seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param tx_len the transcript lengths of the transcripts
#'  a named (tx names) vector of integers.
#'  If you have the transcripts as GRangesList,
#'  call ORFik:::widthPerGroup(tx, T)
#'  If not:
#'  Normally you call argument as: tx_len = txLen(Gtf)
#'  fractionLength(grl, tx_len = txLen(Gtf))
#'  See ?ORFik:::txLen
#'  NB!!! if you used cage data, then
#'  the tss for the the leaders have changed,
#'  therefor the tx lengths have changed,
#'  if so, call:
#'  te(grl, RNA, RFP, tx_len = TxLen(Gtf, fiveUTRs))
#'  where fiveUTRs are the changed cageLeaders.
#' @param with.fpkm logical F, if true return the fpkm values together with te
#' @importFrom data.table data.table
#' @family features
#' @export
#' @return a numeric vector of fpkm ratios, if with.fpkm is TRUE
#'  ,return a data.table with te and fpkm values
te <- function(grl, RNA, RFP, tx_len, with.fpkm = F){

  libraryRna <- length(RNA)
  libraryRPF <- length(RFP)

  overlapRNA <- countOverlaps(grl, RNA)
  overlapRFP <- countOverlaps(grl, RFP)

  if(is.null(names(tx_len))) stop("tx_len have no names")
  tx_len <- tx_len[OrfToTxNames(grl)]
  if(sum(is.na(tx_len))) stop("tx_len have na values, check naming")
  grl_len <- widthPerGroup(grl, F)

  #normalize by tx lengths
  fpkmRNA <- fpkm(overlapRNA, tx_len, libraryRna)
  #normalize by grl lengths
  fpkmRFP <- fpkm(overlapRFP, grl_len, libraryRPF)
  if(with.fpkm){
    return(data.table(fpkmRFP = fpkmRFP, fpkmRNA = fpkmRNA,
                      te = fpkmRFP / fpkmRNA))
  }
  return(fpkmRFP / fpkmRNA)
}

#' Fraction Length
#' @description is defined as (length of grls)/(length of transcripts)
#' see article: 10.1242/dev.098343
#' @param grl a GRangesList object with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param tx_len the transcript lengths of the transcripts
#'  a named (tx names) vector of integers.
#'  If you have the transcripts as GRangesList,
#'  call ORFik:::widthPerGroup(tx, T)
#'  If not:
#'  Normally you call argument as: tx_len = txLen(Gtf)
#'  fractionLength(grl, tx_len = txLen(Gtf))
#'  See ?ORFik:::txLen
#'  NB!!! if you used cage data, then
#'  the tss for the the leaders have changed,
#'  therefor the tx lengths have changed,
#'  if so, call:
#'  te(grl, RNA, RFP, tx_len = TxLen(Gtf, fiveUTRs))
#'  where fiveUTRs are the changed cageLeaders.
#' @family features
#' @export
#' @return a numeric vector of ratios
fractionLength <- function(grl, tx_len){
  grl_len <- widthPerGroup(grl, F)
  tx_len <- tx_len[OrfToTxNames(grl)]
  names(tx_len) <- NULL
  return(grl_len / tx_len)
}

#' Disengagement score (DS)
#' @description is defined as (RPFs over ORF)/(RPFs downstream to tx end).
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1242/dev.098344
#' @param grl a GRangesList object with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param GtfOrTx if Gtf: a TxDb object of a gtf file,
#'  if tx: a GrangesList of transcripts, called from:
#'  exonsBy(Gtf, by = "tx", use.names = T)
#' @family features
#' @export
#' @return a named vector of numeric values of scores
disengagementScore <- function(grl, RFP, GtfOrTx){
  overlapGrl <- countOverlaps(grl, RFP) + 1

  if(class(GtfOrTx) == "TxDb"){
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = T)
  } else if(class(GtfOrTx) == "GRangesList") {
    tx <- GtfOrTx
  } else {
    stop("GtfOrTx is neithter of type TxDb or GRangesList")
  }

  tx <- tx[OrfToTxNames(grl, F)]

  grlStops <- ORFStopSites(grl,asGR = F)
  downstreamTx <- downstreamOfPerGroup(tx, grlStops)

  overlapDownstream <- countOverlaps(downstreamTx, RFP) + 1
  score <- overlapGrl / overlapDownstream
  names(score) <- NULL
  return(score)
}

#' Ribosome Release Score (RRS)
#' @description is defined as (RPFs over ORF)/(RPFs over 3' utrs).
#' normalized by lengths
#' , if RNA is added as argument, normalize by RNA counts over areas too
#' to justify location of 3' utrs
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1016/j.cell.2013.06.009
#' @param grl a GRangesList object with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param GtfOrThreeUtrs if Gtf: a TxDb object of a gtf file,
#'  if ThreeUtrs: a GrangesList of transcripts, called from:
#'  threeUTRsByTranscript(Gtf, use.names = T)
#' @param RNA rna seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @family features
#' @export
#' @return a named vector of numeric values of scores
RibosomeReleaseScore <- function(grl, RFP, GtfOrThreeUtrs, RNA = NULL){
  overlapGrl <- countOverlaps(grl, RFP) + 1

  if(class(GtfOrThreeUtrs) == "TxDb"){
    threeUtrs <- threeUTRsByTranscript(GtfOrThreeUtrs, by = "tx",
                                        use.names = T)
  } else if(class(GtfOrThreeUtrs) == "GRangesList") {
    threeUtrs <- GtfOrThreeUtrs
  } else {
    stop("GtfOrThreeUtrs is neithter of type TxDb or GRangesList")
  }
  threeUtrs <- threeUtrs[OrfToTxNames(grl, F)]


  overlapThreeUtrs <- countOverlaps(threeUtrs, RFP) + 1
  rrs <- (overlapGrl / widthPerGroup(grl)) /
    (overlapThreeUtrs /  widthPerGroup(threeUtrs))

  if(!is.null(RNA)){ # normalize by rna ratio
    rnaRatio <- (countOverlaps(grl, RNA) + 1) /
      (countOverlaps(threeUtrs, RNA) + 1)
    rrs <- rrs / rnaRatio
  }
  names(rrs) <- NULL
  return(rrs)
}

#' Get all possible feature in ORFik
#' @description If you want to get all the features easily, run this function.
#' Each feature have a link to an article describing feature,
#' try ?floss
#' @param grl a GRangesList object with usually ORFs, but can also be
#'  either leaders, cds', 3' utrs or  ORFs are a special case,
#'  see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param RNA rna seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param Gtf a TxDb object of a gtf file,
#' @param tx a GrangesList of transcripts,
#'  normally called from: exonsBy(Gtf, by = "tx", use.names = T)
#'  only add this if you are not including Gtf file
#'  You do not need to reassign these to the cage peaks, it will do it for you.
#' @param fiveUTRs fiveUTRs as GRangesList, must be original unchanged fiveUTRs
#' @param cds a GRangesList of coding sequences
#' @param threeUTRs  a GrangesList of transcript 3' utrs,
#'  normally called from: threeUTRsByTranscript(Gtf, use.names = T)
#' @param faFile a FaFile from the fasta file, see ?FaFile
#' @param riboStart usually 26, the start of the floss interval, see ?floss
#' @param riboStop usually 34, the end of the floss interval
#' @param extension a numeric/integer needs to be set! set to 0 if you did not
#'  use cage, if you used cage to change tss' when finding the orfs,
#'  standard cage extension is 1000
#' @param orfFeatures a logical,  is the grl a list of orfs? Must be assigned.
#' @param cageFiveUTRs a GRangesList, if you used cage-data to extend 5' utrs,
#'  include this, also extension must match with the extension used for these.
#' @importFrom data.table data.table
#' @family features
#' @export
#' @return a data.table with scores, each column is one score type, name of
#'  columns are the names of the scores, i.g floss or fpkmRFP.
#' @examples
#'  \dontrun{
#'  #The easiest way to run the method is to include as this:
#'  allFeatures(grl = grl, orfFeatures =  T, RFP = RFP, RNA = RNA, Gtf = Gtf,
#'   faFile = faFile, extension = 0)
#'  #The other arguments will then be found by the Gtf file,
#'  #and extension = 0 means
#'  #you did not use cage to extend 5' utrs
#'  }
allFeatures <- function(grl, RFP, RNA = NULL,  Gtf = NULL, tx = NULL,
                        fiveUTRs = NULL, cds = NULL, threeUTRs = NULL,
                        faFile = NULL, riboStart = 26, riboStop = 34,
                        extension = NULL, orfFeatures = T,
                        cageFiveUTRs = NULL){

  validExtension(extension, cageFiveUTRs)
  validGRL(class(grl), "grl")
  checkRFP(class(RFP))
  checkRNA(class(RNA))
  tx_len <- NULL
  if(is.null(Gtf)){
    validGRL(c(class(fiveUTRs), class(cds), class(threeUTRs)),
             c("fiveUTRs", "cds", "threeUTRs"))
    if(class(tx) != "GRangesList"){ stop("if Gtf is not given,\n
                             tx must be specified as a GRangesList")
    }
  } else {
    if(class(Gtf) != "TxDb") stop("gtf must be TxDb object")

    notIncluded <- validGRL(c(class(fiveUTRs), class(cds),
                              class(threeUTRs), class(tx)), c("fiveUTRs",
                                "cds", "threeUTRs", "tx"), TRUE)
    if(notIncluded[1]){
      fiveUTRs <- fiveUTRsByTranscript(Gtf, use.names = TRUE)
    }
    if(notIncluded[2]){
      cds <- cdsBy(Gtf, by = "tx", use.names = TRUE)
    }
    if(notIncluded[3]){
      threeUTRs <- threeUTRsByTranscript(Gtf, use.names = TRUE)
    }
    if(notIncluded[4]){
      tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
    }
  }
  if(!is.null(cageFiveUTRs)){
    tx <- ORFik:::extendLeaders(tx, extension = cageFiveUTRs)
  }
  tx_len <- ORFik:::widthPerGroup(tx, TRUE)

  floss <- ORFik:::floss(grl, RFP, cds, riboStart, riboStop)
  entropyRFP <- ORFik:::entropy(grl, RFP)
  fractionLengths <- ORFik:::fractionLength(grl, tx_len)
  disengagementScores <- ORFik:::disengagementScore(grl, RFP, tx)
  RRS <- ORFik:::RibosomeReleaseScore(grl, RFP, threeUTRs, RNA)
  scores <- data.table(floss = floss, entropyRFP = entropyRFP,
                       fractionLengths = fractionLengths,
                       disengagementScores = disengagementScores,
                       RRS = RRS)
  if(!is.null(RNA)){ # if rna seq is included
    # TODO: add fpkm functions specific for rfp and rna
    te <- ORFik:::te(grl, RNA, RFP, tx_len, with.fpkm = T)
    scores$te <- te$te
    scores$fpkmRFP <- te$fpkmRFP
    scores$fpkmRNA <- te$fpkmRNA
  }
  if(orfFeatures){ # if features are found for orfs
    scores$ORFScores <- ORFik:::ORFScores(grl, RFP)$ORFscore

    if(class(faFile) == "FaFile"){
      scores$kozak <-  ORFik:::grl = ORFik:::kozakSequenceScore(grl, faFile)
    } else { message("faFile not included, skipping kozak sequence score")}
    distORFCDS <- ORFik:::distOrfToCds(grl, fiveUTRs, cds, extension)
    scores$distORFCDS <- distORFCDS
    scores$inFrameCDS <- ORFik:::inFrameWithCDS(distORFCDS)
    scores$isOverlappingCds <- ORFik:::isOverlappingCds(distORFCDS)
    scores$rankInTx <- ORFik:::OrfRankOrder(grl)
  } else {
    message("orfFeatures set to False, dropping all orf features.")
  }
  return(scores)
}
