
#' Create normalizations of counts
#'
#' Normally use function \code{\link{fpkm}}, if you want unusual normalization
#' , you can use this.
#' Short for: Fragments per kilobase of transcript per million fragments
#' Normally used in Translations efficiency calculations
#' see article: 10.1038/nbt.1621
#' @param counts a list of integer counts per object
#' @param lengthSize a list of integer lengths per object
#' @param librarySize a numeric of size 1, the size of the library
#' @family features
#' @export
#' @return a numeric vector
fpkm_calc <- function(counts, lengthSize, librarySize){
  return((as.numeric(counts) * (10^9)) /
           (as.numeric(lengthSize) * as.numeric(librarySize)))
}

#' Create normalizations of counts
#'
#' Short for: Fragments per kilobase of transcript per million fragments
#'  For ribo seq on orfs, send orfs as grl,
#'  for rna seq always send transcripts as grl
#' see article: 10.1038/nbt.1621
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs are a special case (uORFs, potential new cds' etc).
#' @param reads a GAlignment, GRanges or GRangesList object,
#'  usually of ribo-seq, rna seq, cage seq, etc.
#' @param pseudoCount an integer, 0,
#'  set it to 1 if you want to avoid NA and inf values.
#' @family features
#' @export
#' @return a numeric vector with the fpkm values
fpkm <- function(grl, reads, pseudoCount = 0){
  grl_len <- widthPerGroup(grl, FALSE)

  overlaps <- countOverlaps(grl, reads)

  librarySize <- length(reads)

  return(fpkm_calc(overlaps, grl_len, librarySize) + pseudoCount)
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} that the reads can map to
#' @param reads a GAlignment object, ig. from ribo seq, or rna seq,
#'  can also be GRanges or GRangesList
#' @family features
#' @export
#' @return A numeric vector containing entropy per ORF
#'
entropy <- function(grl, reads) {
  # Get count list of overlaps
  tileBy1 <- tile1(grl)
  unlTile <- unlist(tileBy1, use.names = FALSE)
  if(!is.null(unlTile$names)){ # TODO: check if this is safe enough
    names(unlTile) <- unlTile$names
  }

  countsTile <- countOverlaps(unlTile, reads)
  names <- names(countsTile)
  names(countsTile) <- NULL
  countList <- split(countsTile, names)

  countList <- IRanges::RleList(countList)
  countList <- countList[names(grl)]
  names(countList) <- NULL
  # generate the entropy variables
  sums <- sum(countList)
  if (sum(as.numeric(sums)) == 0) { # no variance in countList, 0 entropy
    return(rep(0, length(tileBy1)))
  }
  N <- unlist(sums, use.names = FALSE)
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
  unlh <- unlist(h, use.names = FALSE)
  unlreg_len <- unlist(reg_len, use.names = FALSE)

  # make sequences for reads, start -> stop
  # gives a triplet reading, 1:3, 3,6
  acums <- L
  for(i in 2:length(L)){
    acums[i] <-acums[i-1] + acums[i]
  }
  unlacums <- unlist(lapply(1:(length(runLengths)-1), function(x){
    rep(acums[x], runLengths[x+1])
  }), use.names = FALSE)
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
  unlintcount <- unlist(unlist(intcountList, use.names = FALSE),
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
  }), use.names = FALSE)

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
#'
#' Defined as: for riboseq reads between size start and stop,
#' sum the fraction of riboseq with width in the interval start:stop
#' for orfs.
#' Divided by the fraction of rriboseq with width in the interval
#' start:stop for coding sequences.
#' A sum of ratios
#' See article: 10.1016/j.celrep.2014.07.045
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object with ORFs
#' @param RFP ribozomal footprints, given as Galignment or GRanges object
#' @param cds a \code{\link[GenomicRanges]{GRangesList}} of coding sequences
#' @param start usually 26, the start of the floss interval
#' @param end usually 34, the end of the floss interval
#' @family features
#' @export
#' @return a vector of FLOSS of length same as grl
floss <- function(grl, RFP, cds, start = 26, end = 34){
  #TODO: add 0 to the ones that dont hit
  if (start > end) stop("start is bigger than end")
  if (class(RFP) == "GRangesList"){
    stop("RFP must be either GAlignment or GRanges type")
  }
  # for orfs
  overlaps <- findOverlaps(grl, RFP)
  rfpWidth <- riboSeqReadWidths(RFP[to(overlaps)])

  rfpPassFilter <- (rfpWidth >= start) & (rfpWidth <= end)
  rfpValidMatch <- rfpWidth[rfpPassFilter]
  ORFGrouping <- from(overlaps)[rfpPassFilter]
  if (sum(as.numeric(ORFGrouping)) == 0){
    return(as.numeric(rep(0, length(grl))))
  }
  whichNoHit <- NULL # which ribo-seq did not hit grl
  if (length(unique(ORFGrouping)) != length(grl)) {
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
  rfpWidth <- riboSeqReadWidths(RFP[to(overlapsCds)])

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
  if (!is.null(whichNoHit)) {
    tempScores <- as.numeric(rep(NA, length(grl)))
    tempScores[unique(ORFGrouping)] <- score
    tempScores[whichNoHit] <- 0.
    score <- tempScores
  }
  if (length(score) != length(grl) || anyNA(score)) stop("could not find\n
                                        floss-score for all objects, most\n
                                        likely objects are wrongly annotated.")
  return(score)
}

#' Translational efficiency
#' @description Uses Rna-seq and Ribo-seq to get te of grl
#' is defined as (density of RPFs within ORF)/(RNA expression).
#' See article: 10.1126/science.1168978
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  can be either transcripts, 5' utrs, cds', 3' utrs or
#'  ORFs are a special case (uORFs, potential new cds' etc).
#' @param RNA rna seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @param tx a GRangesList of the transcripts,
#'  NB!!! if you used cage data, then
#'  the tss for the the leaders have changed,
#'  therefor the tx lengths have changed,
#'  if so, call:
#'  te(grl, RNA, RFP, tx= TxLen(Gtf, fiveUTRs))
#'  where fiveUTRs are the changed cageLeaders.
#' @param with.fpkm logical F, if true return the fpkm values together with te
#' @param pseudoCount an integer, 0,
#'  set it to 1 if you want to avoid NA and inf values.
#' @importFrom data.table data.table
#' @family features
#' @export
#' @return a numeric vector of fpkm ratios, if with.fpkm is TRUE
#'  ,return a data.table with te and fpkm values
te <- function(grl, RNA, RFP, tx, with.fpkm = F, pseudoCount = 0){
  tx <- tx[OrfToTxNames(grl)]
  #normalize by tx lengths
  fpkmRNA <- fpkm(tx, RNA, pseudoCount)
  #normalize by grl lengths
  fpkmRFP <- fpkm(grl, RFP, pseudoCount)
  if (with.fpkm) {
    return(data.table(fpkmRFP = fpkmRFP, fpkmRNA = fpkmRNA,
                      te = fpkmRFP / fpkmRNA))
  }
  return(fpkmRFP / fpkmRNA)
}

#' Fraction Length
#' @description is defined as (length of grls)/(length of transcripts)
#' see article: 10.1242/dev.098343
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param tx_len the transcript lengths of the transcripts
#'  a named (tx names) vector of integers.
#'  If you have the transcripts as GRangesList,
#'  call ORFik:::widthPerGroup(tx, T)
#'  If not:
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
  grl_len <- widthPerGroup(grl, FALSE)
  tx_len <- tx_len[OrfToTxNames(grl)]
  names(tx_len) <- NULL
  return(grl_len / tx_len)
}

#' Disengagement score (DS)
#' @description is defined as (RPFs over ORF)/(RPFs downstream to tx end).
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1242/dev.098344
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually either leaders,
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

  if (class(GtfOrTx) == "TxDb") {
    tx <- exonsBy(GtfOrTx, by = "tx", use.names = T)
  } else if (class(GtfOrTx) == "GRangesList") {
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
#'
#' Is defined as (RPFs over ORF)/(RPFs over 3' utrs),
#' normalized by lengths
#' , if RNA is added as argument, normalize by RNA counts over areas too
#' to justify location of 3' utrs
#' It can be seen as a ribosome stalling feature.
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1016/j.cell.2013.06.009
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually either leaders,
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
#' @return a named vector of numeric values of scores, NA means that
#'  no 3' utr was found for that transcript.
RibosomeReleaseScore <- function(grl, RFP, GtfOrThreeUtrs, RNA = NULL){
  if (class(GtfOrThreeUtrs) == "TxDb") {
    threeUTRs <- threeUTRsByTranscript(GtfOrThreeUtrs, by = "tx",
                                        use.names = T)
  } else if (class(GtfOrThreeUtrs) == "GRangesList") {
    threeUTRs <- GtfOrThreeUtrs
  } else {
    stop("GtfOrThreeUtrs is neithter of type TxDb or GRangesList")
  }
  # check that naming is correct, else change it.
  orfNames <- OrfToTxNames(grl, F)
  validNamesThree <- names(threeUTRs) %in% orfNames
  validNamesGRL <- orfNames %in% names(threeUTRs)
  rrs <- rep(NA,length(grl))
  if (sum(validNamesGRL) != length(grl)) {
    threeUTRs <- threeUTRs[validNamesThree]
    grl <- grl[validNamesGRL]
  }
  overlapGrl <- countOverlaps(grl, RFP) + 1
  threeUTRs <- threeUTRs[orfNames[validNamesGRL]]
  overlapThreeUtrs <- countOverlaps(threeUTRs, RFP) + 1

  rrs[validNamesGRL] <- (overlapGrl / widthPerGroup(grl)) /
    (overlapThreeUtrs / widthPerGroup(threeUTRs))

  if (!is.null(RNA)) { # normalize by rna ratio
    rnaRatio <- (countOverlaps(grl, RNA) + 1) /
      (countOverlaps(threeUTRs, RNA) + 1)
    rrs[validNamesGRL] <- rrs[validNamesGRL] / rnaRatio
  }
  names(rrs) <- NULL
  return(rrs)
}

#' Ribosome Stalling Score (RSS)
#'
#' Is defined as (RPFs over ORF stop sites)/(RPFs over ORFs),
#' normalized by lengths
#' A pseudo-count of one was added to both the ORF and downstream sums.
#' See article: 10.1016/j.cels.2017.08.004
#'
#' For a more accurate analysis see article:  10.1016/j.cels.2017.08.004.
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually either leaders,
#'  cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len
#' @param RFP ribo seq reads as GAlignment, GRanges
#'  or GRangesList object
#' @family features
#' @export
#' @return a named vector of numeric values of scores
RibosomeStallingScore <- function(grl, RFP){
  grl_len <- widthPerGroup(grl, FALSE)
  overlapGrl <- countOverlaps(grl, RFP)
  stopCodons <- ORFStopCodons(grl, TRUE)
  overlapStop <- countOverlaps(stopCodons, RFP)

  rss <- ((overlapStop + 1) / 3) / ((overlapGrl + 1) / grl_len)
  names(rss) <- NULL
  return(rss)
}

#' Get all possible features in ORFik
#'
#' If you want to get all the features easily, run this function.
#' Each feature have a link to an article describing feature,
#' try ?floss
#'
#' Remember to sort ORFs, use: sortPerGroup(grl)
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#'  with usually ORFs, but can also be
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
#' @param includeNonVarying a logical T, if TRUE, include all features not dependent on
#'  Ribo-seq data and RNA-seq data, that is: Kozak, fractionLengths, distORFCDS,
#'  inFrameCDS, isOverlappingCds and rankInTx
#' @importFrom data.table data.table
#' @family features
#' @export
#' @return a data.table with scores, each column is one score type, name of
#'  columns are the names of the scores, i.g \code{\link{floss}}
#'  or \code{\link{fpkm}}
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
                        cageFiveUTRs = NULL, includeNonVarying = T){

  validExtension(extension, cageFiveUTRs)
  validGRL(class(grl), "grl")
  checkRFP(class(RFP))
  checkRNA(class(RNA))
  tx_len <- NULL
  if (is.null(Gtf)) {
    validGRL(c(class(fiveUTRs), class(cds), class(threeUTRs)),
             c("fiveUTRs", "cds", "threeUTRs"))
    if (class(tx) != "GRangesList") { stop("if Gtf is not given,\n
                             tx must be specified as a GRangesList")
    }
  } else {
    if(class(Gtf) != "TxDb") stop("gtf must be TxDb object")

    notIncluded <- validGRL(c(class(fiveUTRs), class(cds),
                              class(threeUTRs), class(tx)), c("fiveUTRs",
                                "cds", "threeUTRs", "tx"), TRUE)
    if (notIncluded[1]) {
      fiveUTRs <- fiveUTRsByTranscript(Gtf, use.names = TRUE)
    }
    if (notIncluded[2]) {
      cds <- cdsBy(Gtf, by = "tx", use.names = TRUE)
    }
    if (notIncluded[3]) {
      threeUTRs <- threeUTRsByTranscript(Gtf, use.names = TRUE)
    }
    if (notIncluded[4]) {
      tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
    }
  }
  if (!is.null(cageFiveUTRs)) {
    tx <- extendLeaders(tx, extension = cageFiveUTRs)
  }
  tx_len <- widthPerGroup(tx, TRUE)

  floss <- floss(grl, RFP, cds, riboStart, riboStop)
  entropyRFP <- entropy(grl, RFP)

  disengagementScores <- disengagementScore(grl, RFP, tx)
  RRS <- RibosomeReleaseScore(grl, RFP, threeUTRs, RNA)
  RSS <- RibosomeStallingScore(grl, RFP)
  scores <- data.table(floss = floss, entropyRFP = entropyRFP,
                       disengagementScores = disengagementScores,
                       RRS = RRS, RSS = RSS)
  if (includeNonVarying) {
    scores$fractionLengths <- fractionLength(grl, tx_len)
  }

  if (!is.null(RNA)) { # if rna seq is included
    # TODO: add fpkm functions specific for rfp and rna
    te <- te(grl, RNA, RFP, tx, with.fpkm = T)
    scores$te <- te$te
    scores$fpkmRFP <- te$fpkmRFP
    scores$fpkmRNA <- te$fpkmRNA
  } else {
    scores$fpkmRFP <- fpkm(grl, RFP)
  }
  if (orfFeatures) { # if features are found for orfs
    scores$ORFScores <- ORFScores(grl, RFP)$ORFscore
    scores$ioScore <- insideOutsideORF(grl, RFP, tx)

    if (includeNonVarying) {
      if (class(faFile) == "FaFile") {
        scores$kozak <- kozakSequenceScore(grl, faFile)
      } else { message("faFile not included, skipping kozak sequence score")}

      distORFCDS <- distOrfToCds(grl, fiveUTRs, cds, extension)
      scores$distORFCDS <- distORFCDS
      scores$inFrameCDS <- inFrameWithCDS(distORFCDS)
      scores$isOverlappingCds <- isOverlappingCds(distORFCDS)
      scores$rankInTx <- OrfRankOrder(grl)
    }
  } else {
    message("orfFeatures set to False, dropping all orf features.")
  }
  return(scores)
}
