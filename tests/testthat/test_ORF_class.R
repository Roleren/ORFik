library(ORFik)
context("ORF class")

test_that("ORF class slots are proper class", {
  grangesObj <- GRanges(Rle(c("chr2", "chr2", "chr2", "chr2"), c(1, 3, 2, 4)),
                        IRanges(1:10, width=10:1))
  trailer <- GRanges(Rle(c("chr2", "chr2", "chr2", "chr2"), c(1, 3, 2, 4)),
                     IRanges(1:10, width=10:1))

  artificialORF <- new("ORF",
                       orfRanges = grangesObj,
                       trailerRanges = trailer,
                       ORFlength = sum(width(grangesObj)),
                       trailerLength = sum(width(trailer)),
                       transcript = "t_name",
                       start = "ATG",
                       family = "CDS")

  expect_is(slot(artificialORF, "orfRanges"), "GRanges")
  expect_is(slot(artificialORF, "trailerRanges"), "GRanges")
  expect_is(slot(artificialORF, "start"), "character")
})
