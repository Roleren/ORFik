# Optional integration checks use existing binaries only (never download tools).
star_test_binary <- function(name, fallback) {
  binary <- Sys.which(name)
  if (!nzchar(binary)) binary <- path.expand(fallback)
  skip_if(!file.exists(binary), paste(name, "is not available"))
  binary
}

test_that("real fastp corrects the low-quality overlapping base only when enabled", {
  skip_on_os("windows")
  fastp <- star_test_binary("fastp", "~/bin/fastp")
  root <- withr::local_tempdir()
  read1 <- "ACGTACGATCGTACCGATGCTAGCTAGGCTACGATCGACTGATCGTAGCTACGATCGATGCTAGCATCGATCGTAGCTAGCA"
  correct2 <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(read1)))
  wrong2 <- correct2
  substr(wrong2, 35, 35) <- if (substr(correct2, 35, 35) == "A") "C" else "A"
  qual2 <- strrep("I", nchar(read1)); substr(qual2, 35, 35) <- "!"
  files <- file.path(root, c("read_1.fastq", "read_2.fastq"))
  writeLines(c("@pair/1", read1, "+", strrep("I", nchar(read1))), files[1])
  writeLines(c("@pair/2", wrong2, "+", qual2), files[2])
  for (correction in c(FALSE, TRUE)) {
    out <- file.path(root, as.character(correction))
    STAR.align.single(files[1], files[2], out, index.dir = root,
                      star.path = NULL, fastp = fastp, steps = "tr",
                      adapter.sequence = "disable", base.correction = correction,
                      max.cpus = 1, verbose = FALSE)
    result <- readLines(file.path(out, "trim/trimmed2_read_1.fastq"))
    expect_identical(result[2], if (correction) correct2 else wrong2)
    expect_identical(substr(result[4], 35, 35), if (correction) "I" else "!")
  }
})

test_that("real STAR maps a novel junction only when allowed and retains indexed splicing", {
  skip_on_os("windows")
  star <- star_test_binary("STAR", "~/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR")
  root <- withr::local_tempdir()
  withr::local_seed(81)
  genome <- paste(sample(c("A", "C", "G", "T"), 10000, TRUE), collapse = "")
  for (start in c(1001, 4001)) {
    substr(genome, start, start + 1) <- "GT"
    substr(genome, start + 98, start + 99) <- "AG"
  }
  fasta <- file.path(root, "genome.fa")
  writeLines(c(">chrTest", genome), fasta)
  # The first junction is indexed; the second exists only in the reads.
  gtf <- file.path(root, "annotation.gtf")
  writeLines(c('chrTest\ttest\texon\t921\t1000\t.\t+\t.\tgene_id "g"; transcript_id "t";',
               'chrTest\ttest\texon\t1101\t1180\t.\t+\t.\tgene_id "g"; transcript_id "t";'), gtf)
  index <- file.path(root, "index")
  dir.create(file.path(index, "genomeDir"), recursive = TRUE)
  status <- system2(star, shQuote(c("--runMode", "genomeGenerate",
                   "--genomeDir", file.path(index, "genomeDir"),
                   "--genomeFastaFiles", fasta, "--sjdbGTFfile", gtf,
                   "--sjdbOverhang", "79", "--genomeSAindexNbases", "5",
                   "--genomeChrBinNbits", "10", "--runThreadN", "1",
                   "--outFileNamePrefix", paste0(root, "/index_"))),
                   stdout = file.path(root, "index.stdout"), stderr = file.path(root, "index.stderr"))
  expect_identical(status, 0L)
  reads <- file.path(root, "junctions.fastq")
  seqs <- vapply(c(1001, 4001), function(start)
    paste0(substr(genome, start - 40, start - 1), substr(genome, start + 100, start + 139)), "")
  writeLines(c("@indexed", seqs[1], "+", strrep("I", 80),
               "@novel", seqs[2], "+", strrep("I", 80)), reads)
  for (allow in c(FALSE, TRUE)) {
    out <- file.path(root, as.character(allow))
    STAR.align.single(reads, output.dir = out, index.dir = index,
                      star.path = star, fastp = NULL, steps = "ge", allow.introns = allow,
                      max.cpus = 1, keep.index.in.memory = "noShared", verbose = FALSE)
    sj <- read.table(file.path(out, "aligned/junctions_SJ.out.tab"))
    expect_true(any(sj$V2 == 1001 & sj$V3 == 1100 & sj$V6 == 1))
    expect_identical(any(sj$V2 == 4001 & sj$V3 == 4100 & sj$V6 == 0), allow)
  }
})
