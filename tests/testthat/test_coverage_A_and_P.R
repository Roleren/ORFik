context("coverage_A_and_P")
library(ORFik)

test_that("A-site coverage is the P-site coverage of the previous codon, including the stop codon", {
  # one transcript: 10 nt leader, 12 nt CDS (4 codons incl. stop), 10 nt trailer, + strand
  mrna <- GRangesList(tx1 = GRanges("chr1", IRanges(1, 32), "+"))
  cds <- GRangesList(tx1 = GRanges("chr1", IRanges(11, 22), "+"))
  # one read per position, score = transcript position, so every value is unique
  reads <- GRanges("chr1", IRanges(1:32, width = 1), "+", score = 1:32)
  for (ap in c("center", "start")) {
    res <- ORFik:::coverage_A_and_P(cds, mrna, reads, aligned_position = ap)
    P <- res$P$score; A <- res$A$score
    expect_equal(length(P), 12); expect_equal(length(A), 12)
    expect_equal(A, P - 3)            # every position, incl. the stop codon's 3rd nt (was P[12] there)
  }
})
