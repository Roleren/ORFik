coverage_chunk_fixture <- function() {
  x <- GAlignments(seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr1", "chr2"),
                                         levels = c("chr1", "chr2", "empty"))),
                   pos = c(2L, 5L, 8L, 2L, 15L),
                   cigar = c("3M4N2M", "2S3M1I2M", "2M2D3M", "3M4N2M", "2M"),
                   strand = Rle(factor(c("+", "-", "*", "+", "-"), levels = c("+", "-", "*"))))
  seqlengths(x) <- c(40L, 30L, 10L)
  mcols(x)$score <- c(2L, 3L, 5L, 7L, 11L)
  x
}

coverage_chunk_values <- function(x) {
  list(forward = lapply(as.list(f(x)), as.numeric),
       reverse = lapply(as.list(r(x)), as.numeric), seqinfo = seqinfo(x), mode = strandMode(x))
}

test_that("alignment coverage chunks preserve blocks, strands, weights and Seqinfo", {
  x <- coverage_chunk_fixture()
  for (ignore in c(FALSE, TRUE)) for (w in list("AUTO", "score", 1L, 0.5, mcols(x)$score, mcols(x)$score/2)) {
    expected <- covRleFromGR(x, w, ignore)
    for (size in c(1, 2, 4, 10)) {
      actual <- suppressMessages(covRleFromGR(x, w, ignore, chunk.size = size))
      expect_equal(coverage_chunk_values(actual), coverage_chunk_values(expected))
    }
  }
  # Independent single-end oracle, not the implementation under test.
  gr <- unlist(grglist(x), use.names = FALSE)
  weights <- rep(mcols(x)$score, lengths(grglist(x)))
  for (ignore in c(FALSE, TRUE)) {
    got <- suppressMessages(covRleFromGR(x, ignore.strand = ignore, chunk.size = 2))
    keep <- if (ignore) rep(TRUE, length(gr)) else as.character(strand(gr)) %in% c("+", "*")
    expect_equal(lapply(as.list(f(got)), as.numeric),
                 lapply(as.list(coverage(gr[keep], weight = as.numeric(weights[keep]))), as.numeric))
  }
})

test_that("coverage retries the real expansion stage and halves failed chunks", {
  x <- coverage_chunk_fixture()
  expected <- covRleFromGR(x)
  original <- ORFik:::.coverage_cigar_ranges
  calls <- integer()
  testthat::local_mocked_bindings(.coverage_cigar_ranges = function(x, ...) {
    calls <<- c(calls, length(x))
    if (length(x) > 1) stop("long vectors not supported yet: memory.c")
    original(x, ...)
  }, .package = "ORFik")
  got <- suppressMessages(covRleFromGR(x))
  expect_equal(coverage_chunk_values(got), coverage_chunk_values(expected))
  expect_equal(calls, c(5L, 2L, rep(1L, 5)))
})

test_that("coverage propagates unrelated errors and terminates irreducible failures", {
  x <- coverage_chunk_fixture()
  calls <- 0L
  testthat::local_mocked_bindings(.coverage_cigar_ranges = function(...) {
    calls <<- calls + 1L
    stop("invalid CIGAR operation")
  }, .package = "ORFik")
  expect_error(covRleFromGR(x), "invalid CIGAR")
  expect_equal(calls, 1L)
  testthat::local_mocked_bindings(.coverage_cigar_ranges = function(...) stop("negative length vectors are not allowed"),
                                .package = "ORFik")
  expect_error(suppressMessages(covRleFromGR(x)), "cannot expand even one alignment/pair")
})

test_that("coverage chunks handle empty input and count overflow", {
  x <- coverage_chunk_fixture()
  for (ignore in c(FALSE, TRUE))
    expect_equal(coverage_chunk_values(covRleFromGR(x[FALSE], ignore.strand = ignore, chunk.size = 1)),
                 coverage_chunk_values(covRleFromGR(x[FALSE], ignore.strand = ignore)))
  x <- x[c(1,1,1)]
  mcols(x)$score <- rep(2^30L, 3)
  got <- suppressMessages(covRleFromGR(x, chunk.size = 1))
  expect_equal(as.numeric(f(got)[["chr1"]])[2:4], rep(3 * 2^30, 3))
  expect_false(anyNA(runValue(f(got)[["chr1"]])))
})

test_that("coverage chunk and weight inputs are validated before processing", {
  x <- coverage_chunk_fixture()
  for (size in list(0, -1, NA_real_, Inf, 1.5, c(1,2), "2", 2^31))
    expect_error(covRleFromGR(x, chunk.size = size), "chunk.size")
  for (w in list(1:2, numeric(), "missing", NA_character_, matrix(1,5,1), factor(1:5)))
    expect_error(covRleFromGR(x, w, chunk.size = 2), "weight")
  withr::local_options(ORFik.coverage.chunk.size = 1)
  expect_message(covRleFromGR(x), "completed alignments/pairs")
})

test_that("paired-end chunking retains pairs and strand modes", {
  first <- coverage_chunk_fixture()[c(1,1,2)]
  last <- GAlignments(seqnames = seqnames(first), pos = start(first) + 1L,
                      cigar = cigar(first), strand = strand(invertStrand(first)),
                      seqinfo = seqinfo(first))
  for (mode in 0:2) {
    x <- GAlignmentPairs(first, last, strandMode = mode)
    mcols(x)$score <- c(2L, 3L, 5L)
    for (ignore in c(FALSE, TRUE)) {
      expected <- covRleFromGR(x, ignore.strand = ignore)
      got <- suppressMessages(covRleFromGR(x, ignore.strand = ignore, chunk.size = 1))
      expect_equal(coverage_chunk_values(got), coverage_chunk_values(expected))
    }
  }
})

test_that("direct IRanges coverage preserves unusual CIGAR operations", {
  x <- GAlignments(seqnames = Rle(factor(rep("chr1", 9))),
                   pos = rep(2L, 9),
                   cigar = c("2H3S2=1X2I3M2D2N1M", "10S", "1X", "2I", "3M",
                             "2M10N2M10N2M", "1M1P2M", "3M", "3M"),
                   strand = Rle(factor(rep(c("+", "-", "*"), 3), levels=c("+","-","*"))))
  seqlengths(x) <- 40L
  mcols(x)$score <- c(0,0.5,2,-1,1,3,4,2^30,2^30)
  bad <- x
  bad@cigar[3] <- "*"
  expect_error(covRleFromGR(bad), "cigar")
  grl <- grglist(x)
  gr <- unlist(grl, use.names=FALSE)
  w <- rep(mcols(x)$score, lengths(grl))
  for(ignore in c(FALSE,TRUE)) {
    got <- covRleFromGR(x, ignore.strand=ignore)
    for(st in if(ignore) "both" else c("+","-")) {
      keep <- if(ignore) rep(TRUE,length(gr)) else as.character(strand(gr)) %in% c(st,"*")
      expected <- coverage(gr[keep], weight=w[keep])
      actual <- if(st=="-") r(got) else f(got)
      expect_equal(lapply(as.list(actual),as.numeric),lapply(as.list(expected),as.numeric))
    }
    expect_equal(coverage_chunk_values(suppressMessages(covRleFromGR(x,ignore.strand=ignore,chunk.size=2))),
                 coverage_chunk_values(got))
  }
})

test_that("compressed-list cumulative endpoint failures trigger recovery", {
  # Base R's compact integer sequences make this a real >2^31-element
  # structural test without a large allocation. Each IRanges is valid:
  # its largest end is 2*n-1, still below INT_MAX. Do not materialize it.
  n <- as.integer(2^30 - 1)
  ranges <- S4Vectors::new2("IRanges", start=seq_len(n), width=seq_len(n), check=FALSE)
  simple <- IRangesList(ranges,ranges,ranges,compress=FALSE)
  expect_length(simple,3)
  expect_equal(sum(as.double(lengths(simple))),3*as.double(n))
  error <- tryCatch(PartitioningByEnd(simple),error=identity)
  expect_s3_class(error,"error")
  expect_true(ORFik:::.coverage_size_error(error))
  for(msg in c("long vectors not supported yet", "negative length vectors are not allowed",
               "cannot allocate vector of size 20 Gb",
               "IRangesList object 'x' is too big (the cumulated length of its list\n  elements is >= 2^32)"))
    expect_true(ORFik:::.coverage_size_error(simpleError(msg)))
  for(msg in c("invalid CIGAR", "seqlengths contains NA", "invalid weight column"))
    expect_false(ORFik:::.coverage_size_error(simpleError(msg)))
})

test_that("coverage-stage failure retries without adding partial chromosome results", {
  x <- coverage_chunk_fixture()
  expected <- covRleFromGR(x)
  testthat::local_mocked_bindings(coverage = function(x, ...) {
    if (is(x, "IRanges") && length(x) > 2L)
      stop("long vectors not supported yet: coverage")
    IRanges::coverage(x, ...)
  }, .package = "ORFik")
  got <- suppressMessages(covRleFromGR(x))
  expect_equal(coverage_chunk_values(got), coverage_chunk_values(expected))
})

test_that("covRleList and coverage export inherit forced chunking", {
  x <- coverage_chunk_fixture()
  expected <- ORFik:::covRleListFromGR(x, verbose=FALSE)
  withr::local_options(ORFik.coverage.chunk.size=1)
  got <- suppressMessages(ORFik:::covRleListFromGR(x,verbose=FALSE))
  expect_equal(names(got@list),names(expected@list))
  for(i in seq_along(got@list))
    expect_equal(coverage_chunk_values(got@list[[i]]),coverage_chunk_values(expected@list[[i]]))
  path <- tempfile("chunked-coverage-",tmpdir=tempdir())
  suppressMessages(ORFik:::export.cov(x,path,seqinfo=seqinfo(x),format="rds"))
  expect_equal(coverage_chunk_values(readRDS(paste0(path,".covrds"))),
               coverage_chunk_values(suppressMessages(covRleFromGR(x))))
})

test_that("contig indexing preserves non-lexical Seqinfo order and unused contigs", {
  set.seed(923)
  n <- 500L
  chromosomes <- paste0("chr",sample.int(100))
  x <- GAlignments(seqnames=Rle(factor(sample(chromosomes[1:60],n,TRUE),levels=chromosomes)),
                   pos=sample.int(50L,n,TRUE),cigar=rep("2M3N2M",n),
                   strand=Rle(factor(sample(c("+","-","*"),n,TRUE),levels=c("+","-","*"))))
  seqlengths(x) <- rep(100L,100)
  mcols(x)$score <- sample.int(10,n,TRUE)
  expanded <- grglist(x)
  gr <- unlist(expanded,use.names=FALSE)
  w <- as.numeric(rep(mcols(x)$score,lengths(expanded)))
  got <- suppressMessages(covRleFromGR(x,chunk.size=101))
  for(st in c("+","-")) {
    keep <- as.character(strand(gr)) %in% c(st,"*")
    ref <- coverage(gr[keep],weight=w[keep])
    ans <- if(st=="+") f(got) else r(got)
    expect_equal(lapply(as.list(ans),as.numeric),lapply(as.list(ref),as.numeric))
  }
  expect_identical(seqinfo(got),seqinfo(x))
})
