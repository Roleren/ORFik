pmapToTranscriptF <-
function (x, transcripts, ignore.strand = FALSE, x.is.sorted = TRUE,
    tx.is.sorted = TRUE, reduce.ranges = TRUE, set.seqlengths = TRUE)
{
    if ((length(x) == 0))
        return(x)
    if (length(x) != length(transcripts)) {
        if (length(x) == 1) {
            x <- rep(x, length(transcripts))
        }
        else if (length(transcripts) == 1) {
            transcripts <- rep(transcripts, length(x))
        }
        else stop("recycling is supported when length(x) == 1 or,",
            "length(transcripts) == 1; otherwise the lengths must match")
    }
    if (is.grl(x) & !x.is.sorted)
        x <- sortPerGroup(x, ignore.strand)
    xOriginal <- x
    oldNames <- names(x)
    xClass <- class(x)
    oldTxNames <- names(transcripts)
    txWidths <- if (is.rl(transcripts)) {
        widthPerGroup(transcripts, FALSE)
    }
    else if (is(transcripts, "IRanges") | is(transcripts, "GRanges")) {
        if (is(transcripts, "IRanges")) {
            transcripts@width
        }
        else {
            transcripts@ranges@width
        }
    }
    else stop("transcripts must either be IRanges, IRangesList, GRanges or GRangesList")
    xPartitionWidths <- if (is.rl(xOriginal)) {
        ends <- xOriginal@partitioning@end
        ends - c(0L, head(ends, -1L))
    }
    else NULL
    if (is.rl(x)) {
        indices <- groupings(x)
        xWidths <- widthPerGroup(x)
        x <- .unlistGrl(x)
    }
    else if (is(x, "IRanges") | is(x, "GRanges")) {
        xWidths <- width(x)
        indices <- seq.int(1, length(x))
    }
    else stop("x must either be IRanges, IRangesList, GRanges or GRangesList")
    names(x) <- NULL
    xStrandOriginal <- if (is.grl(xOriginal)) {
        as.character(xOriginal@unlistData@strand)
    }
    else if (is(xOriginal, "GRanges")) {
        as.character(strand(xOriginal))
    }
    else NULL
    if (!is.logical(ignore.strand))
        stop("ignore.strand must be logical")
    ignore.strand <- as.logical(ignore.strand)
    if (length(indices) && max(indices) > length(transcripts))
        stop("invalid names of IRanges")
    if (length(x) != length(indices))
        stop("length of ranges != indices")
    notEqualSeqnames <- is.gr_or_grl(xOriginal) & is.gr_or_grl(transcripts) &
        !all(seqlevels(xOriginal) %in% seqlevels(transcripts))
    if (notEqualSeqnames)
        stop("subscript contains out-of-bounds indices")
    if (is.grl(transcripts) & !tx.is.sorted)
        transcripts <- sortPerGroup(transcripts, ignore.strand)
    if (!all(xWidths <= txWidths)) {
        stop("Invalid ranges to map, check them.", " One has width bigger than its reference")
    }
    if (is.grl(transcripts) | is(transcripts, "IRangesList")) {
        tx <- .unlistGrl(transcripts)
        txGroupings <- groupings(transcripts)
        ends <- transcripts@partitioning@end
        exonN <- ends - c(0L, head(ends, -1L))
    }
    else {
        tx <- transcripts
        txGroupings <- seq.int(1, length(transcripts))
        exonN <- rep.int(1L, length(transcripts))
    }
    if (!is(tx, "IRanges"))
        tx <- tx@ranges
    names(tx) <- NULL
    txStrand <- if (ignore.strand) {
        rep(TRUE, length(transcripts))
    }
    else strandBool(transcripts)
    xStrand <- if (ignore.strand) {
        rep(TRUE, length(indices))
    }
    else txStrand[indices]
    if (is.grl(xOriginal) | is(xOriginal, "IRangesList")) {
        posWidths <- xPartitionWidths[txStrand]
        negWidths <- xPartitionWidths[!txStrand]
        indicesPos <- rep.int(seq_along(posWidths), posWidths)
        indicesNeg <- rep.int(seq_along(negWidths), negWidths)
    }
    else {
        indicesPos <- seq_along(xOriginal[txStrand])
        indicesNeg <- seq_along(xOriginal[!txStrand])
    }
    exonCumSumPos <- c(0, cumsum(exonN[txStrand]))[indicesPos]
    exonCumSumNeg <- c(0, cumsum(exonN[!txStrand]))[indicesNeg]
    txStrand <- txStrand[txGroupings]
    if (any(txStrand)) {
        pos <- pmapToTranscriptsCPP(start(x)[xStrand], end(x)[xStrand],
            start(tx)[txStrand], end(tx)[txStrand], txGroupings[txStrand],
            "+", exonCumSumPos)
    }
    else {
        pos <- list(ranges = list(vector("integer"), vector("integer")),
            index = vector("integer"))
    }
    if (any(!txStrand)) {
        neg <- pmapToTranscriptsCPP(start(x)[!xStrand], end(x)[!xStrand],
            start(tx)[!txStrand], end(tx)[!txStrand], txGroupings[!txStrand],
            "-", exonCumSumNeg)
    }
    else {
        neg <- list(ranges = list(vector("integer"), vector("integer")),
            index = vector("integer"))
    }
    xStart <- vector("integer", length(indices))
    xStart[xStrand] <- pos$ranges[[1]]
    xStart[!xStrand] <- neg$ranges[[1]]
    xEnd <- vector("integer", length(indices))
    xEnd[xStrand] <- pos$ranges[[2]]
    xEnd[!xStrand] <- neg$ranges[[2]]
    if (is.gr_or_grl(xOriginal) && is.gr_or_grl(transcripts)) {
        query_seqnames <- as.character(if (is.grl(xOriginal)) {
            xOriginal@unlistData@seqnames
        } else xOriginal@seqnames)
        subject_seqnames <- if (is.grl(transcripts)) {
            as.character(transcripts@unlistData@seqnames)[transcripts@partitioning@end]
        }
        else as.character(transcripts@seqnames)
        wrong_sequence <- query_seqnames != subject_seqnames[indices]
        xStart[wrong_sequence] <- 0L
        xEnd[wrong_sequence] <- -1L
    }
    result <- IRanges(xStart, xEnd)
    if (is.gr_or_grl(xClass)) {
        newSeqnames <- if (!is.null(oldTxNames)) {
            oldTxNames
        }
        else seq.int(length(result))
        newStrand <- if (ignore.strand) {
            "*"
        }
        else if (is.grl(transcripts)) {
            strandPerGroup(transcripts, FALSE)[indices]
        }
        else as.character(strand(transcripts))[indices]
        result <- GRanges(seqnames = newSeqnames[indices], ranges = result,
            strand = newStrand)
        unmapped <- start(result) == 0
        if (!ignore.strand)
            unmapped <- unmapped | (strand(result) != xStrandOriginal)
        if (any(unmapped)) {
            strand(result)[unmapped] <- "*"
        }
        if (set.seqlengths) {
            seqlevels.used <- seqlevels(result)
            if (is.null(names(transcripts))) {
                seqlevels.used <- as.integer(seqlevels.used)
            }
            suppressWarnings({
                seqlengths(result) <- if (is.grl(transcripts)) {
                  txWidths[seqlevels.used]
                }
                else txWidths[seqlevels.used]
            })
        }
    }
    if (is.grl(xClass) | is(xOriginal, "IRangesList")) {
        result <- relist(result, xOriginal)
        names(result) <- oldNames
        if (reduce.ranges && any(xPartitionWidths > 1L)) {
            result <- reduce(result, drop.empty.ranges = FALSE)
        }
    }
    else names(result) <- oldNames
    return(result)
}
pmapFromTranscriptF <-
function (x, transcripts, removeEmpty = FALSE)
{
    if (!length(x))
        return(transcripts[integer()])
    if (is.null(names(x))) {
        if (length(x) != length(transcripts))
            stop("invalid names of ranges, when length(x) != length(transcripts)")
    }
    if (is(x, "GRanges"))
        x <- ranges(x)
    if (is(x, "IRangesList")) {
        x = unlist(x, use.names = TRUE)
        indices <- strtoi(names(x))
        names(x) <- NULL
    }
    else if (is(x, "IRanges")) {
        if (is.null(names(x))) {
            indices <- seq.int(1, length(x))
        }
        else {
            indices <- strtoi(names(x))
        }
        names(x) <- NULL
    }
    else stop("x must either be IRanges, IRangesList or GRanges")
    if (!is.logical(removeEmpty))
        stop("removeEmpty must be logical")
    removeEmpty <- as.logical(removeEmpty)
    if (max(indices) > length(transcripts))
        stop("invalid names of IRanges")
    if (length(x) != length(indices))
        stop("length of ranges != indices")
    tx <- ranges(transcripts)
    names <- names(tx)
    names(tx) <- NULL
    tx <- tx[indices]
    if (!all(width(x) <= as.integer(sum(width(tx))))) {
        stop("Invalid ranges to map, check them. One is bigger than its reference")
    }
    groupings <- groupings(tx)
    tx <- unlist(tx, use.names = FALSE)
    xStrand <- strandBool(transcripts)[indices]
    txStrand <- xStrand[groupings]
    if (any(txStrand)) {
        pos <- pmapFromTranscriptsCPP(start(x)[xStrand], end(x)[xStrand],
            start(tx)[txStrand], end(tx)[txStrand], groupings[txStrand],
            "+", removeEmpty)
    }
    else {
        pos <- list(ranges = list(vector("integer"), vector("integer")),
            index = vector("integer"))
    }
    if (any(!txStrand)) {
        neg <- pmapFromTranscriptsCPP(start(x)[!xStrand], end(x)[!xStrand],
            start(tx)[!txStrand], end(tx)[!txStrand], groupings[!txStrand],
            "-", removeEmpty)
    }
    else {
        neg <- list(ranges = list(vector("integer"), vector("integer")),
            index = vector("integer"))
    }
    if (removeEmpty) {
        txStrand <- order(c(pos$index, neg$index)) <= length(pos$index)
    }
    temp <- indices
    nIndices <- vector("integer", length(txStrand))
    nIndices[txStrand] <- pos$index
    nIndices[!txStrand] <- neg$index
    indices <- indices[nIndices]
    xStart <- vector("integer", length(txStrand))
    xStart[txStrand] <- unlist(pos$ranges[1], use.names = FALSE)
    xStart[!txStrand] <- unlist(neg$ranges[1], use.names = FALSE)
    xEnd <- vector("integer", length(txStrand))
    xEnd[txStrand] <- unlist(pos$ranges[2], use.names = FALSE)
    xEnd[!txStrand] <- unlist(neg$ranges[2], use.names = FALSE)
    result <- GRanges(seqnames = seqnamesPerGroup(transcripts,
        FALSE)[indices], ranges = IRanges(xStart, xEnd), strand = strandPerGroup(transcripts,
        FALSE)[indices])
    names(result) <- names[indices]
    result <- split(result, nIndices)
    names(result) <- names(transcripts)[temp]
    seqlevels(result) <- seqlevels(transcripts)
    seqinfo(result) <- seqinfo(transcripts)
    return(result)
}
