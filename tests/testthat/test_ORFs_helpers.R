library(ORFik)
context("ORF helpers")

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = c(1, 10, 20, 30, 40),
                              end = c(5, 15, 25, 35, 45)),
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
                            ranges = IRanges(start = rev(c(1, 10, 20, 30, 40)),
                              end = rev(c(5, 15, 25, 35, 45))),
                                strand = Rle(strand(rep("-", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(1, 10, 20)),
                        end = rev(c(5, 15, 25))),
                          strand = Rle(strand(rep("-", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(10, 20, 30)),
                        end = rev(c(15, 25, 35))),
                          strand = Rle(strand(rep("-", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(20, 30, 40)),
                        end = rev(c(25, 35, 45))),
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
                            ranges = IRanges(start = rev(c(10, 20, 30, 40)),
                              end = rev(c(15, 25, 35, 45))),
                                strand = Rle(strand(rep("-", 4))))



test_that("findORFs works as intended for plus strand", {

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
seqname <- c("tx1", "tx2", "tx3", "tx4")
seqs <- c("ATGGGTATTTATA", "ATGGGTAATA",
          "ATGGG", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
grIn1 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(21, 10), end = c(23, 19)),
                 strand = rep("-", 2), names = rep(seqname[1], 2))
grIn2 <- GRanges(seqnames = rep("1", 1),
                 ranges = IRanges(start = c(1010), end = c(1019)),
                 strand = rep("-", 1), names = rep(seqname[2], 1))

grIn3 <- GRanges(seqnames = rep("1", 1),
                 ranges = IRanges(start = c(2000), end = c(2004)),
                 strand = rep("-", 1), names = rep(seqname[3], 1))

grIn4 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(3030, 3000), end = c(3036, 3029)),
                 strand = rep("-", 2), names = rep(seqname[4], 2))

grl <- GRangesList(grIn1, grIn2, grIn3, grIn4)
names(grl) <- seqname

test_that("findORFs works as intended for minus strand", {

  #longestORF F with different frames
  test_ranges <-findORFs(grl,seqs,
                                    "ATG|TGG|GGG",
                                    "TAA|AAT|ATA",
                                    longestORF = F,
                                    minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "-")
  expect_equal(as.integer(unlist(start(test_ranges))), c(10, 21, 1011, 1010, 1012))
  expect_equal(as.integer(unlist(end(test_ranges))), c(19, 22, 1019, 1018, 1017))
  expect_equal(as.integer(unlist(width(test_ranges))), c(10, 2, 9, 9, 6))
  expect_equal(sum(widthPerGroup(test_ranges) %% 3), 0)
})


# Create data for get_all_ORFs_as_GRangesList test_that#2
namesTx <- c("tx1", "tx2")
seqs <- c("ATGATGTAATAA", "ATGTAA")
grIn1 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(1, 3), end = c(1, 13)),
                 strand = rep("+", 2), names = rep(namesTx[1], 2))
grIn2<- GRanges(seqnames = rep("1", 6),
                 ranges = IRanges(start = c(1, 1000, 2000, 3000, 4000, 5000),
                                  end = c(1, 1000, 2000, 3000, 4000, 5000)),
                 strand = rep("+", 6), names = rep(namesTx[2], 6))
grl <- GRangesList(grIn1, grIn2)
names(grl) <- namesTx
test_that("mapToGRanges works as intended for strange exons positive strand", {

  #longestORF F with different frames
  test_ranges <- findORFs(grl,seqs,
                                             "ATG|TGG|GGG",
                                             "TAA|AAT|ATA",
                                             longestORF = F,
                                             minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(1, 3, 5,1, 1000, 2000,
                                                         3000, 4000, 5000))
  expect_equal(as.integer(unlist(end(test_ranges))), c(1, 10, 10,1, 1000, 2000,
                                                       3000, 4000, 5000))
  expect_equal(sum(widthPerGroup(test_ranges) %% 3), 0)
  expect_equal(unlist(grl)$names,c("tx1", "tx1", "tx2", "tx2", "tx2", "tx2",
                                   "tx2", "tx2"))
  expect_equal(unlist(test_ranges)$names,c("tx1_1", "tx1_1", "tx1_2", "tx2_1",
                                           "tx2_1", "tx2_1", "tx2_1", "tx2_1",
                                           "tx2_1"))
})

# Create data for get_all_ORFs_as_GRangesList test_that#3
ranges(grIn1) <- rev(ranges(grIn1))
strand(grIn1) <- rep("-", length(grIn1))
ranges(grIn2) <- rev(ranges(grIn2))
strand(grIn2) <- rep("-", length(grIn2))

grl <- GRangesList(grIn1, grIn2)
names(grl) <- namesTx

test_that("mapToGRanges works as intended for strange exons negative strand", {

  #longestORF F with different frames
  test_ranges <- findORFs(grl,seqs,
                                    "ATG|TGG|GGG",
                                    "TAA|AAT|ATA",
                                    longestORF = F,
                                    minimumLength = 0)

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "-")
  expect_equal(as.integer(unlist(start(test_ranges))), c(5, 5, 5000, 4000, 3000,
                                                         2000, 1000, 1))
  expect_equal(as.integer(unlist(end(test_ranges))), c(13, 10, 5000, 4000, 3000,
                                                       2000, 1000, 1))
  expect_equal(sum(widthPerGroup(test_ranges) %% 3), 0)
  expect_equal(unlist(grl)$names,c("tx1", "tx1", "tx2", "tx2", "tx2",
                                   "tx2", "tx2", "tx2"))
  expect_equal(unlist(test_ranges)$names,c("tx1_1","tx1_2", "tx2_1", "tx2_1",
                                   "tx2_1", "tx2_1", "tx2_1", "tx2_1"))
})

namesTx <- c("tx1", "tx2", "tx3", "tx4")
seqs <- c("ATGATGTAATAA", "ATGTAA", "AAAATGAAATAAA", "AAAATGAAATAA")


grIn3 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(2000, 2008), end = c(2004, 2015)),
                 strand = rep("+", 2), names = rep(namesTx[3], 2))
grIn4 <- GRanges(seqnames = rep("1", 2),
                 ranges = IRanges(start = c(3030, 3000), end = c(3036, 3004)),
                 strand = rep("-", 2), names = rep(namesTx[4], 2))
grl <- GRangesList(grIn1, grIn2, grIn3, grIn4)
names(grl) <- namesTx

test_that("mapToGRanges works as intended for strange exons both strands", {

  #longestORF F with different frames
  test_ranges <- findORFs(grl,seqs,
                                    "ATG|TGG|GGG",
                                    "TAA|AAT|ATA",
                                    longestORF = F,
                                    minimumLength = 0)

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "-")
  expect_equal(as.integer(unlist(start(test_ranges))), c(5, 5, 5000, 4000,
                                                         3000, 2000, 1000, 1, 2003,
                                                         2008, 3030, 3000))
  expect_equal(as.integer(unlist(end(test_ranges))), c(13, 10, 5000, 4000,
                                                       3000, 2000, 1000, 1, 2004,
                                                       2014, 3033, 3004))
  expect_equal(sum(widthPerGroup(test_ranges) %% 3), 0)
})


test_that("GRangesList sorting works as intended", {

  test_ranges <- grl[3:4]

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(2000,
                                                         2008, 3030, 3000))
  expect_equal(as.integer(unlist(end(test_ranges))), c(2004,
                                                       2015, 3036, 3004))

  test_ranges <- sortPerGroup(test_ranges, ignore.strand = T)
  expect_equal(as.integer(unlist(start(test_ranges))), c(2000,
                                                         2008, 3000, 3030))
  expect_equal(as.integer(unlist(end(test_ranges))), c(2004,
                                                       2015, 3004, 3036))
})

test_that("ORFStartCodons works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                        strand = Rle(strand(rep("+", 3))))

  ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                        strand = Rle(strand(rep("+", 3))))
  ORFranges4 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(50, 40, 30), end = c(55, 45, 35)),
                        strand = Rle(strand(rep("-", 3))))
  ORFranges5 <- GRanges(seqnames = Rle(rep("1", 4)),
                        ranges = IRanges(start = c(1000, 1002, 1004, 1006),
                                         end = c(1000, 1002, 1004, 1006)),
                        strand = Rle(strand(rep("+", 4))))
  ORFranges6 <- GRanges(seqnames = Rle(rep("1", 4)),
                        ranges = IRanges(start = c(1002, 1004, 1005, 1006),
                                         end = c(1002, 1004, 1005, 1006)),
                        strand = Rle(strand(rep("+", 4))))
  ORFranges4 <- sort(ORFranges4, decreasing = TRUE)
  names(ORFranges) <- rep("tx1_1" ,3)
  names(ORFranges2) <- rep("tx1_2", 3)
  names(ORFranges3) <- rep("tx1_3", 3)
  names(ORFranges4) <- rep("tx4_1", 3)
  names(ORFranges5) <- rep("tx1_4", 4)
  names(ORFranges6) <- rep("tx1_5", 4)
  grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2,
                     tx1_3 = ORFranges3, tx4_1 = ORFranges4,
                     tx1_4 = ORFranges5, tx1_5 = ORFranges6)


  test_ranges <- ORFStartCodons(grl, TRUE)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(1,
                                                         20, 30, 53, 1000,
                                                         1002, 1004, 1002, 1004))
  expect_equal(as.integer(unlist(end(test_ranges))), c(3,
                                                       22, 32, 55, 1000,
                                                       1002, 1004, 1002, 1005))

})

test_that("ORFStopCodons works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                        strand = Rle(strand(rep("+", 3))))

  ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                        strand = Rle(strand(rep("+", 3))))
  ORFranges4 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                        strand = Rle(strand(rep("-", 3))))
  ORFranges5 <- GRanges(seqnames = Rle(rep("1", 4)),
                        ranges = IRanges(start = c(1000, 1002, 1004, 1006),
                                         end = c(1000, 1002, 1004, 1006)),
                        strand = Rle(strand(rep("+", 4))))
  ORFranges6 <- GRanges(seqnames = Rle(rep("1", 4)),
                        ranges = IRanges(start = c(1002, 1003, 1004, 1006),
                                         end = c(1002, 1003, 1004, 1006)),
                        strand = Rle(strand(rep("+", 4))))
  ORFranges4 <- sort(ORFranges4, decreasing = TRUE)
  names(ORFranges) <- rep("tx1_1" ,3)
  names(ORFranges2) <- rep("tx1_2", 3)
  names(ORFranges3) <- rep("tx1_3", 3)
  names(ORFranges4) <- rep("tx4_1", 3)
  names(ORFranges5) <- rep("tx1_4", 4)
  names(ORFranges6) <- rep("tx1_5", 4)
  grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2,
                     tx1_3 = ORFranges3, tx4_1 = ORFranges4,
                     tx1_4 = ORFranges5, tx1_5 = ORFranges6)


  test_ranges <- ORFStopCodons(grl, TRUE)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(23,
                                                         43, 53, 30, 1002,
                                                         1004, 1006, 1003, 1006))
  expect_equal(as.integer(unlist(end(test_ranges))), c(25,
                                                       45, 55, 32, 1002,
                                                       1004, 1006, 1004, 1006))

  # check with meta columns
  ORFranges$names <- rep("tx1_1" ,3)
  ORFranges2$names <- rep("tx1_2", 3)
  ORFranges3$names <- rep("tx1_3", 3)
  ORFranges4$names <- rep("tx4_1", 3)
  ORFranges5$names <- rep("tx1_4", 4)
  ORFranges6$names <- rep("tx1_5", 4)
  grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2,
                     tx1_3 = ORFranges3, tx4_1 = ORFranges4,
                     tx1_4 = ORFranges5, tx1_5 = ORFranges6)

  test_ranges <- ORFStopCodons(grl, TRUE)
})

test_that("uniqueORFs works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))),
                       names = paste0("tx1_", rep(1, 3)))
  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))),
                       names = paste0("tx2_", rep(1, 3)))
  ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                        strand = Rle(strand(rep("-", 3))),
                        names = paste0("tx3_", rep(1, 3)))

  names(ORFranges) <- rep("tx1_1" ,3)
  names(ORFranges2) <- rep("tx2_1", 3)
  names(ORFranges3) <- rep("tx3_1", 3)
  grl <- GRangesList(tx1_1 = ORFranges, tx2_1 = ORFranges2,
                     tx3_1 = ORFranges3)

  test_ranges <- uniqueORFs(grl)

  expect_is(test_ranges, "GRangesList")
  expect_equal(strandPerGroup(test_ranges,F)[1], "+")
  expect_equal(length(test_ranges), 2)
  expect_equal(names(test_ranges),  c("1", "2"))

})

