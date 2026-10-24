context("STAR index failure handling")
library(ORFik)

test_that("STAR_MAKE_INDEX.sh fails when the genome index fails, even if later indices succeed", {
  skip_on_os("windows")
  dir <- tempfile("star_index_test_"); dir.create(dir)
  fake <- file.path(dir, "STAR")
  writeLines(c("#!/usr/bin/env bash",
               "dir=''; fasta=''; prev=''",
               "for a in \"$@\"; do [ \"$prev\" = '--genomeDir' ] && dir=$a; [ \"$prev\" = '--genomeFastaFiles' ] && fasta=$a; prev=$a; done",
               "case \"$fasta\" in *genome.fa) exit 137;; esac",
               "mkdir -p \"$dir\"; touch \"$dir/Genome\" \"$dir/SA\" \"$dir/SAindex\" \"$dir/genomeParameters.txt\"; exit 0"), fake)
  Sys.chmod(fake, "755")
  genome <- file.path(dir, "genome.fa"); contam <- file.path(dir, "contam.fa")
  writeLines(c(">chr1", "ACGTACGTACGTACGTACGTACGTACGT"), genome)
  writeLines(c(">c1", "ACGTTTTTGGGGCCCCAAAA"), contam)
  script <- system.file("STAR_Aligner", "STAR_MAKE_INDEX.sh", package = "ORFik")
  status <- system2("bash", c(script, "-o", file.path(dir, "out/"), "-S", fake, "-m", "1", "-R", "1000000000",
                              "-a", "1", "-T", "-", "-f", genome, "-c", contam), stdout = FALSE, stderr = FALSE)
  expect_false(status == 0)
  expect_false(file.exists(file.path(dir, "out", "genomeDir", "genomeParameters.txt")))
})
