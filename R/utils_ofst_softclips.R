# Canonicalize only in-memory merge inputs, never the source OFST files.
.ofst_remove_softclips <- function(dt) {
  for (col in intersect(c("cigar", "cigar1", "cigar2"), names(dt))) {
    values <- unique(as.character(dt[[col]]))
    if (anyNA(values) || any(!nzchar(values)))
      .ofst_abort("remove_softclips requires non-missing, non-empty ", col, " values.")
    changed <- grepl("S", values, fixed = TRUE)
    if (!any(changed)) next
    clipped <- values[changed]
    ops <- gsub("[0-9]+", "", clipped)
    if (any(!grepl("^([0-9]+[MIDNSHP=X])+$", clipped)) ||
        any(!grepl("^H?S?[MIDNP=X]*S?H?$", ops)))
      .ofst_abort("remove_softclips: invalid CIGAR or internal soft clip in ", col, ".")
    cleaned <- values
    cleaned[changed] <- gsub("[0-9]+S", "", clipped)
    if (any(cleaned == ""))
      .ofst_abort("remove_softclips: a CIGAR contains only soft clips in ", col,
                  "; stripping would leave an empty CIGAR. Remove unmapped/clip-only records separately.")
    index <- match(as.character(dt[[col]]), values)
    data.table::set(dt, j = col, value = cleaned[index])
    query_col <- switch(col, cigar = "qwidth", cigar1 = "qwidth1", cigar2 = "qwidth2")
    if (query_col %in% names(dt))
      data.table::set(dt, j = query_col,
        value = GenomicAlignments::cigarWidthAlongQuerySpace(cleaned)[index])
  }
  dt
}
