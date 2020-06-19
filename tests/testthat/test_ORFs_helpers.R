context("ORF helpers")
library(ORFik)

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = c(1, 10, 20, 30, 40),
                                             end = c(5, 15, 25, 35, 45)),
                            strand = Rle(strand(rep("+", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(10, 20, 30),
                                       end = c(15, 25, 35)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40),
                                       end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))

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


test_that("defineTrailer works as intended for plus strand", {

  #at the start
  trailer <- defineTrailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_equal(start(trailer), c(30, 40))
  expect_equal(end(trailer), c(35, 45))

  #middle
  trailer2 <- defineTrailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 40)
  expect_equal(end(trailer2), 45)

  #at the end
  trailer3 <- defineTrailer(ORFranges3, transcriptRanges)
  expect_is(trailer3, "GRanges")
  expect_equal(length(trailer3), 0)

  #trailer size 3
  trailer4 <- defineTrailer(ORFranges2, transcriptRanges, 3)
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

test_that("defineTrailer works as intended for minus strand", {

  #at the end
  trailer <- defineTrailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_is(trailer, "GRanges")
  expect_equal(length(trailer), 0)

  #middle
  trailer2 <- defineTrailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 1)
  expect_equal(end(trailer2), 5)

  #at the start
  trailer3 <- defineTrailer(ORFranges3, transcriptRanges)
  expect_equal(start(trailer3), c(1, 10))
  expect_equal(end(trailer3), c(5, 15))

  #trailer size 3
  trailer4 <- defineTrailer(ORFranges2, transcriptRanges, 3)
  expect_equal(start(trailer4), 3)
  expect_equal(end(trailer4), 5)
})

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 4)),
                            ranges = IRanges(start = rev(c(10, 20, 30, 40)),
                                             end = rev(c(15, 25, 35, 45))),
                            strand = Rle(strand(rep("-", 4))))

test_that("findORFsFasta works as intended", {
  filePath <- system.file("extdata", "genome.fasta",
                          package = "ORFik")

  test_result <- findORFsFasta(filePath, longestORF = FALSE)
  expect_is(test_result, "GRanges")
  expect_equal(length(test_result), 3990)

  ## allow circular
  test_result <- findORFsFasta(filePath, longestORF = FALSE,
                               is.circular = TRUE)
  expect_is(test_result, "GRanges")
  expect_equal(length(test_result), 3998)

})



test_that("findORFs works as intended for plus strand", {

  #longestORF F with different frames
  test_ranges <- findORFs("ATGGGTAATA", "ATG|TGG|GGG", "TAA|AAT|ATA",
                          longestORF = FALSE, minimumLength = 0)
  expect_is(test_ranges, "IRangesList")
  expect_equal(unlist(start(test_ranges), use.names = FALSE), c(1, 2, 3))
  expect_equal(unlist(end(test_ranges), use.names = FALSE), c(9, 10, 8))

  #longestORF T
  test_ranges <- findORFs("ATGATGTAATAA", "ATG|TGA|GGG", "TAA|AAT|ATA",
                          longestORF = TRUE, minimumLength = 0)
  expect_is(test_ranges, "IRangesList")
  expect_equal(unlist(start(test_ranges), use.names = FALSE), c(1, 2))
  expect_equal(unlist(end(test_ranges), use.names = FALSE), c(9, 10))

  #longestORF F with minimum size 12 -> 6 + 3*2
  test_ranges <- findORFs("ATGTGGAATATGATGATGATGTAATAA", "ATG|TGA|GGG",
                          "TAA|AAT|ATA", longestORF = FALSE, minimumLength = 2)
  expect_is(test_ranges, "IRangesList")
  expect_equal(unlist(start(test_ranges), use.names = FALSE), c(10, 13, 11, 14))
  expect_equal(unlist(end(test_ranges), use.names = FALSE), c(24, 24, 25, 25))

  #longestORF T with minimum size 12 -> 6 + 3*2
  test_ranges <- findORFs("ATGTGGAATATGATGATGATGTAATAA", "ATG|TGA|GGG",
                          "TAA|AAT|ATA", longestORF = TRUE, minimumLength = 2)
  expect_is(test_ranges, "IRangesList")
  expect_equal(unlist(start(test_ranges), use.names = FALSE), c(10, 11))
  expect_equal(unlist(end(test_ranges), use.names = FALSE), c(24, 25))

  #find nothing
  test_ranges <- findORFs("B", "ATG|TGA|GGG", "TAA|AAT|ATA", minimumLength = 2)
  expect_is(test_ranges, "IRangesList")
  expect_equal(length(test_ranges), 0)

})


test_that("findMapORFs works as intended for minus strand", {

  #longestORF F with different frames
  test_ranges <-findMapORFs(grl,seqs,
                            "ATG|TGG|GGG",
                            "TAA|AAT|ATA",
                            longestORF = FALSE,
                            minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "-")
  expect_equal(as.integer(unlist(start(test_ranges))),
               c(21, 10, 1011, 1010, 1012))
  expect_equal(as.integer(unlist(end(test_ranges))),
               c(22, 19, 1019, 1018, 1017))
  expect_equal(as.integer(unlist(width(test_ranges))), c(2, 10, 9, 9, 6))
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
  test_ranges <- findMapORFs(grl,seqs,
                             "ATG|TGG|GGG",
                             "TAA|AAT|ATA",
                             longestORF = FALSE,
                             minimumLength = 0)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges,FALSE)[1], "+")
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
  test_ranges <- findMapORFs(grl,seqs,
                             "ATG|TGG|GGG",
                             "TAA|AAT|ATA",
                             longestORF = FALSE,
                             minimumLength = 0)

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "-")
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
  test_ranges <- findMapORFs(grl,seqs,
                             "ATG|TGG|GGG",
                             "TAA|AAT|ATA",
                             longestORF = FALSE,
                             minimumLength = 0)

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "-")
  expect_equal(as.integer(unlist(start(test_ranges))),
               c(5, 5, 5000, 4000, 3000, 2000, 1000, 1, 2003,
                 2008, 3030, 3000))
  expect_equal(as.integer(unlist(end(test_ranges))),
               c(13, 10, 5000, 4000, 3000, 2000, 1000, 1, 2004,
                 2014, 3033, 3004))
  expect_equal(sum(widthPerGroup(test_ranges) %% 3), 0)
})

test_that("pmapFromTranscriptsF works as intended", {
            xStart = c(1, 5, 10, 1000, 5, 6, 1, 1)
            xEnd = c(6, 8, 12, 2000, 10, 10, 3, 1)
            TS = c(1,5, 1000, 1005, 1008, 2000, 2003, 4000, 5000, 7000, 85, 70,
                   101, 9)
            TE = c(3, 9, 1003, 1006, 1010, 2001, 2020, 4500, 6000, 8000, 89,
                   82, 105, 9)
            indices = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 7)
            strand = c(rep("+", 10), rep("-", 3), "+")
            seqnames = rep("1", length(TS))
            result <- split(IRanges(xStart, xEnd), c(seq.int(1, 5), 5, 6, 7))
            transcripts <- split(GRanges(seqnames, IRanges(TS, TE), strand),
                                 indices)
            test_ranges <- pmapFromTranscriptF(result, transcripts,  TRUE)
            expect_is(test_ranges, "GRangesList")
            expect_equal(start(unlistGrl(test_ranges)),
                         c(1, 5, 1005, 1008, 2010, 5498, 7000, 85, 78, 78,
                           103, 9))
            expect_equal(end(unlistGrl(test_ranges)),
                         c(3, 7, 1006, 1009, 2012, 6000, 7497, 85, 82, 82,
                           105, 9))
})



test_that("GRangesList sorting works as intended", {

  test_ranges <- grl[3:4]

  test_ranges <- sortPerGroup(test_ranges)
  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(2000,
                                                         2008, 3030, 3000))
  expect_equal(as.integer(unlist(end(test_ranges))), c(2004,
                                                       2015, 3036, 3004))

  test_ranges <- sortPerGroup(test_ranges, ignore.strand = TRUE)
  expect_equal(as.integer(unlist(start(test_ranges))), c(2000,
                                                         2008, 3000, 3030))
  expect_equal(as.integer(unlist(end(test_ranges))), c(2004,
                                                       2015, 3004, 3036))
})

test_that("startCodons works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20),
                                        end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(20, 30, 40),
                                         end = c(25, 35, 45)),
                        strand = Rle(strand(rep("+", 3))))

  ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50),
                                         end = c(35, 45, 55)),
                        strand = Rle(strand(rep("+", 3))))
  ORFranges4 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(50, 40, 30),
                                         end = c(55, 45, 35)),
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


  test_ranges <- startCodons(grl, TRUE)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))),
               c(1, 20, 30, 53, 1000, 1002, 1004, 1002, 1004))
  expect_equal(as.integer(unlist(end(test_ranges))),
               c(3, 22, 32, 55, 1000, 1002, 1004, 1002, 1005))

})

test_that("stopCodons works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 10, 20),
                                        end = c(5, 15, 25)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(20, 30, 40),
                                         end = c(25, 35, 45)),
                        strand = Rle(strand(rep("+", 3))))

  ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50),
                                         end = c(35, 45, 55)),
                        strand = Rle(strand(rep("+", 3))))
  ORFranges4 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(30, 40, 50),
                                         end = c(35, 45, 55)),
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


  test_ranges <- stopCodons(grl, TRUE)

  expect_is(test_ranges, "GRangesList")
  expect_is(strand(test_ranges),"CompressedRleList")
  expect_is(seqnames(test_ranges),"CompressedRleList")
  expect_equal(strandPerGroup(test_ranges, FALSE)[1], "+")
  expect_equal(as.integer(unlist(start(test_ranges))), c(23,43, 53, 30, 1002,
                                                         1004, 1006, 1003,
                                                         1006))
  expect_equal(as.integer(unlist(end(test_ranges))), c(25,45, 55, 32, 1002,
                                                       1004, 1006, 1004,
                                                       1006))

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

  test_ranges <- stopCodons(grl, TRUE)
  negStopss <- GRangesList(tx1_1 = GRanges("1", c(7, 5, 3, 1), "-"),
                           tx1_2 = GRanges("1", c(15, 13, 11, 9), "-"))
  expect_equal(stopSites(stopCodons(negStopss, FALSE), is.sorted = TRUE),
               c(1,9))

  negStopss <- GRangesList(tx1_1 = GRanges("1", IRanges(c(9325,8012),
                                                        c(9418, 8013)),
                                           "-"))
  expect_equal(startSites(stopCodons(negStopss, FALSE), is.sorted = TRUE),
               9325)
})

ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(10, 20, 30),
                                       end = c(15, 25, 35)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40),
                                       end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges$names <- rep("tx1_1" ,3)
ORFranges2$names <- rep("tx1_2", 3)
ORFranges3$names <- rep("tx1_3", 3)
orfs <- c(ORFranges,ORFranges2,ORFranges3)
grl <- groupGRangesBy(orfs, orfs$names)

test_that("startRegion works as intended", {
  transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                              ranges = IRanges(start = c(1, 10, 20, 30, 40),
                                               end = c(5, 15, 25, 35, 45)),
                              strand = Rle(strand(rep("+", 5))))
  transcriptRanges <- groupGRangesBy(transcriptRanges,
                                     rep("tx1", length(transcriptRanges)))

  test_ranges <- startRegion(grl, transcriptRanges)
  expect_equal(as.integer(unlist(start(test_ranges))), c(1, 4, 10, 14, 20))
  expect_equal(as.integer(unlist(end(test_ranges))), c(3, 5, 12, 15, 22))
  test_ranges <- startRegion(grl)
  expect_equal(as.integer(unlist(end(test_ranges))), c(3, 12, 22))
})

test_that("uniqueGroups works as intended", {
  grl[3] <- grl[1]
  test_ranges <- uniqueGroups(grl)

  expect_is(test_ranges, "GRangesList")
  expect_equal(strandPerGroup(test_ranges, FALSE), c("+", "+"))
  expect_equal(length(test_ranges), 2)
  expect_equal(names(test_ranges),  c("1", "2"))

})

test_that("uniqueOrder works as intended", {
  gr1 <- GRanges("1", IRanges(1,10), "+")
  gr2 <- GRanges("1", IRanges(20, 30), "+")
  # make a grl with duplicated ORFs (gr1 twice)
  grl <- GRangesList(tx1_1 = gr1, tx2_1 = gr2, tx3_1 = gr1)
  test_result <- uniqueOrder(grl) # remember ordering
  expect_equal(test_result, as.integer(c(1,2,1)))
})

test_that("findUORFs works as intended", {
   # Load annotation
   txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                           package = "GenomicFeatures")
   txdb <- loadTxdb(txdbFile)
   fiveUTRs <- loadRegion(txdb, "leaders")
   cds <- loadRegion(txdb, "cds")
   if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
     # Normally you would not use a BSgenome, but some custome fasta-
     # annotation you  have for your species
     uorfs <- findUORFs(fiveUTRs["uc001bum.2"],
                        BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                        "ATG", cds = cds)
     expect_equal(names(uorfs[1]), "uc001bum.2_5")
     expect_equal(length(uorfs), 1)
   }

})

test_that("artificial.orfs works as intended", {
  cds <- GRangesList(tx1 = GRanges("chr1", IRanges(start = c(100), end = 150),"+"),
                     tx2 = GRanges("chr1", IRanges(200, 205), "+"),
                     tx3 = GRanges("chr1", IRanges(300, 311), "+"),
                     tx4 = GRanges("chr1", IRanges(400, 999), "+"),
                     tx5 = GRanges("chr1", IRanges(500, 511), "-"))
  res <- artificial.orfs(cds)
  expect_equal(100, startSites(res[1]))
})
