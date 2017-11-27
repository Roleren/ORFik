library(ORFik)
context("ORF helpers")

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = c(1, 10, 20, 30, 40), end = c(5, 15, 25, 35, 45)),
                            strand = Rle(strand(rep("+", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(10, 20, 30), end = c(15, 25, 35)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))

test_that("define_trailer works as intended for plus strand", {

  #at the start
  trailer <- define_trailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_equal(start(trailer), c(30, 40))
  expect_equal(end(trailer), c(35, 45))

  #middle
  trailer2 <- define_trailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 40)
  expect_equal(end(trailer2), 45)

  #at the end
  trailer3 <- define_trailer(ORFranges3, transcriptRanges)
  expect_is(trailer3, "GRanges")
  expect_equal(length(trailer3), 0)

  #trailer size 3
  trailer4 <- define_trailer(ORFranges2, transcriptRanges, 3)
  expect_equal(start(trailer4), 40)
  expect_equal(end(trailer4), 42)
})


transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = rev(c(1, 10, 20, 30, 40)), end = rev(c(5, 15, 25, 35, 45))),
                            strand = Rle(strand(rep("-", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = rev(c(1, 10, 20)), end = rev(c(5, 15, 25))),
                     strand = Rle(strand(rep("-", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(10, 20, 30)), end = rev(c(15, 25, 35))),
                      strand = Rle(strand(rep("-", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(20, 30, 40)), end = rev(c(25, 35, 45))),
                      strand = Rle(strand(rep("-", 3))))

test_that("define_trailer works as intended for minus strand", {

  #at the end
  trailer <- define_trailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_is(trailer, "GRanges")
  expect_equal(length(trailer), 0)

  #middle
  trailer2 <- define_trailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 1)
  expect_equal(end(trailer2), 5)

  #at the start
  trailer3 <- define_trailer(ORFranges3, transcriptRanges)
  expect_equal(start(trailer3), c(1, 10))
  expect_equal(end(trailer3), c(5, 15))

  #trailer size 3
  trailer4 <- define_trailer(ORFranges2, transcriptRanges, 3)
  expect_equal(start(trailer4), 3)
  expect_equal(end(trailer4), 5)
})

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 4)),
                            ranges = IRanges(start = rev(c(10, 20, 30, 40)), end = rev(c(15, 25, 35, 45))),
                            strand = Rle(strand(rep("-", 4))))



test_that("find_in_frame_ORFs works as intended for plus strand", {

  #longestORF F with different frames
  test_ranges <- orfs_as_IRanges("ATGGGTAATA",
                                  "ATG|TGG|GGG",
                                  "TAA|AAT|ATA",
                                  longestORF = FALSE,
                                  minimumLength = 0)
  expect_is(test_ranges, "IRanges")
  expect_equal(start(test_ranges), c(1, 2, 3))
  expect_equal(end(test_ranges), c(9, 10, 8))

  #longestORF T
  test_ranges <- orfs_as_IRanges("ATGATGTAATAA",
                                  "ATG|TGA|GGG",
                                  "TAA|AAT|ATA",
                                  longestORF = T,
                                  minimumLength = 0)
  expect_is(test_ranges, "IRanges")
  expect_equal(start(test_ranges), c(1, 2))
  expect_equal(end(test_ranges), c(9, 10))

  #longestORF F with minimum size 12 -> 6 + 3*2
  test_ranges <- orfs_as_IRanges("ATGTGGAATATGATGATGATGTAATAA",
                                    "ATG|TGA|GGG",
                                    "TAA|AAT|ATA",
                                    longestORF = F,
                                    minimumLength = 2)
  expect_is(test_ranges, "IRanges")
  expect_equal(start(test_ranges), c(10, 13, 11, 14))
  expect_equal(end(test_ranges), c(24, 24, 25, 25))

  #longestORF T with minimum size 12 -> 6 + 3*2
  test_ranges <- orfs_as_IRanges("ATGTGGAATATGATGATGATGTAATAA",
                                  "ATG|TGA|GGG",
                                  "TAA|AAT|ATA",
                                  longestORF = T,
                                  minimumLength = 2)
  expect_is(test_ranges, "IRanges")
  expect_equal(start(test_ranges), c(10, 11))
  expect_equal(end(test_ranges), c(24, 25))

  #find nothing
  test_ranges <- orfs_as_IRanges("B",
                                  "ATG|TGA|GGG",
                                  "TAA|AAT|ATA",
                                  longestORF = T,
                                  minimumLength = 2)
  expect_is(test_ranges, "IRanges")
  expect_equal(length(test_ranges), 0)

})



# Create data for get_all_ORFs_as_GRangesList test_that#1
seqname <- c("tx1","tx2","tx3","tx4")
seqs <- c("ATGGGTATTTATA","ATGGGTAATA","ATGGG", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
grIn1 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(10, 20), end = c(19, 22)),
                 strand = rep("-", 2), names = rep(seqname[1],2))
grIn2 <- GRanges(seqnames = rep("1", 1),
                 ranges = IRanges(start = c(1010), end = c(1019)),
                 strand = rep("-", 1), names = rep(seqname[2],1))

grIn3 <- GRanges(seqnames = rep("1", 1),
                 ranges = IRanges(start = c(2000), end = c(2004)),
                 strand = rep("-", 1), names = rep(seqname[3],1))

grIn4 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(3000,3030), end = c(3029,3036)),
                 strand = rep("-", 2), names = rep(seqname[4],2))

grl <- GRangesList(grIn1,grIn2,grIn3,grIn4)
names(grl) <- seqname

test_that("get_all_ORFs_as_GRangesList works as intended", {

  #longestORF F with different frames
  test_ranges <-find_in_frame_ORFs(grl,seqs,
                                    "ATG|TGG|GGG",
                                    "TAA|AAT|ATA",
                                    longestORF = F,
                                    minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(as.integer(unlist(start(test_ranges))), c(10, 20, 1011, 1010,1012))
  expect_equal(as.integer(unlist(end(test_ranges))), c(19, 21, 1019, 1018,1017))


})

# Create data for get_all_ORFs_as_GRangesList test_that#2
seqname <- c("tx1","tx1")
seqs <- c("ATGATGTAATAA")
grIn1 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(1, 2), end = c(1, 12)),
                 strand = rep("+", 2), names = rep(seqname[1],2))

grl <- GRangesList(grIn1)
names(grl) <- "tx1"
test_that("map_to_GRanges works as intended for strange exons", {

  #longestORF F with different frames
  test_ranges <- find_in_frame_ORFs(grl,seqs,
                                             "ATG|TGG|GGG",
                                             "TAA|AAT|ATA",
                                             longestORF = F,
                                             minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(as.integer(unlist(start(test_ranges))), c(1, 2, 4))
  expect_equal(as.integer(unlist(end(test_ranges))), c(1, 9, 9))

})


