parse_cigar <- function(cigar, shift, is_plus_strand) {
    c_signs <- unlist(explodeCigarOps(cigar))
    c_counts <- unlist(explodeCigarOpLengths(cigar))

    i = ifelse(is_plus_strand, 0, length(c_signs) + 1)
    increment = ifelse(is_plus_strand, 1, -1)
    limit = 0
    plusShift = 0
    while (shift >= limit) {
        i = i + increment
        if (c_signs[i] == "M") {
            limit = limit + c_counts[i]
        } else if (c_signs[i] == "N" || c_signs[i] == "D") {
            plusShift = plusShift + c_counts[i]
        } else if (c_signs[i] == "I") {
            plusShift = plusShift - c_counts[i]
        } else {
            warning(paste0("Not supported sign:", c_signs[i]))
        }
    }
    shift = shift + plusShift
    return(shift)
}

#' Shifts footprints
#'
#' Function shifts footprints from bam file and loads them to memory as GRanges.
#' Resizes them to single base in 5' end fashion,
#' treated as p site. Takes account for junctions in cigars. Length of footprint is saved in
#' 'size' parameter of GRanges output. Footprints are also sorted according to their genomic
#' position, ready for saving e.g. as bed file.
#' @param bam_input Path to bam file. Bam file should have .bai file generated.
#' @param selected_lengths Numeric vector of lengths of footprints you select for shifting.
#' @param selected_shifts Numeric vector of shifts for coresponding selected_lengths.
#' eg. c(10, -10) with selected_lengths of c(31, 32) means length of 31 will be shifted left by 10.
#' Footprints of length 32 will be shifted right by 10.
#' @return A GRanges object of shifted footprints.
#' @export
#' @import GenomicAlignments
#' @import GenomeInfoDb
#' @examples
#' #shift_footprints()

shift_footprints <- function(bam_input, selected_lengths, selected_shifts) {

    selected_shifts <- -1 * selected_shifts
    allFootrpintsShifted <- GRanges()

    message("Loading reads into memory...")
    riboReads <- readGAlignments(bam_input)
    message("Reads loaded.")

    # check for input correctness
    if (length(selected_lengths) != length(selected_shifts)) {
        stop("Incorrect input. Not equal number of elements in selected_lengths and selected_shifts!")
    }
    if (sum(selected_lengths > abs(selected_shifts)) != length(selected_shifts)) {
        stop("Incorrect input. selected_shifts cant be bigger than selected_lengths!")
    }

    for (i in 1:length(selected_lengths)) {
        message("Shifting footprints of length ", selected_lengths[i])
        riboReadsW <- riboReads[qwidth(riboReads) == selected_lengths[i]]
        if (length(riboReadsW) == 0) {
            next
        }
        is_cigar <- width(riboReadsW) != qwidth(riboReadsW)
        cigar_strings <- cigar(riboReadsW[is_cigar])
        sizes <- qwidth(riboReadsW)

        riboReadsW <- granges(riboReadsW, use.mcols = TRUE)
        riboReadsW$size <- sizes  #move footprint length to size
        riboReadsW <- resize(riboReadsW, 1)  #resize to 5' only

        cigars <- riboReadsW[is_cigar]
        notCigars <- riboReadsW[!is_cigar]

        # Not Cigars - shift 5' ends, + shift right, - shift left
        if (length(notCigars) != 0) {
            is_plus <- as.vector(strand(notCigars) == "+")
            shiftedNotCigarsPlus <- shift(notCigars[is_plus], selected_shifts[i])
            shiftedNotCigarsMinus <- shift(notCigars[!is_plus], -1 * selected_shifts[i])
            allFootrpintsShifted <- c(allFootrpintsShifted, shiftedNotCigarsPlus, shiftedNotCigarsMinus)
        }
        # Cigars
        if (length(cigars) != 0) {
            is_plus <- as.vector(strand(cigars) == "+")
            shift_by <- rep(selected_shifts[i], length(cigars))
            shift_by <- mapply(parse_cigar, cigar_strings, shift_by, is_plus)
            shift_by[!is_plus] <- -1 * shift_by[!is_plus]
            shifted_cigars <- shift(cigars, shift_by)
            allFootrpintsShifted <- c(allFootrpintsShifted, shifted_cigars)
        }
    }

    message("Sorting shifted footprints...")
    allFootrpintsShifted <- sortSeqlevels(allFootrpintsShifted)
    allFootrpintsShifted <- sort(allFootrpintsShifted)
    message("Process completed successfuly.")
    return(allFootrpintsShifted)
}
