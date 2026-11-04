context("Assembly accessions for refseq/genbank downloads")
library(ORFik)

test_that("is_assembly_accession detects NCBI assembly accessions only", {
  expect_true(ORFik:::is_assembly_accession("GCF_000196035.1"))
  expect_true(ORFik:::is_assembly_accession("GCA_000005845"))
  expect_false(ORFik:::is_assembly_accession("Listeria monocytogenes"))
  expect_false(ORFik:::is_assembly_accession("GCF_123"))
  expect_false(ORFik:::is_assembly_accession(c("GCF_000196035.1", "GCA_000005845")))
})

test_that("loadTxdb gives GFF transcripts with missing or duplicated names their feature ID", {
  gff <- tempfile(fileext = ".gff3")
  writeLines(c("##gff-version 3",
    "chr1\tRefSeq\tgene\t101\t400\t.\t+\t.\tID=gene-a1;Name=prs",
    "chr1\tRefSeq\tCDS\t101\t400\t.\t+\t0\tID=cds-1;Parent=gene-a1",
    "chr1\tRefSeq\tgene\t1001\t1300\t.\t+\t.\tID=gene-a2;Name=prs",
    "chr1\tRefSeq\tCDS\t1001\t1300\t.\t+\t0\tID=cds-2;Parent=gene-a2",
    "chr1\tRefSeq\tgene\t2001\t2300\t.\t-\t.\tID=gene-a3;Name=ispD",
    "chr1\tRefSeq\tCDS\t2001\t2300\t.\t-\t0\tID=cds-3;Parent=gene-a3",
    "chr1\tRefSeq\tmRNA\t3001\t3300\t.\t+\t.\tID=rna-x4;Parent=gene-a3",
    "chr1\tRefSeq\texon\t3001\t3300\t.\t+\t.\tID=exon-4;Parent=rna-x4"), gff)
  txdb <- suppressWarnings(loadTxdb(gff))
  tx <- GenomicFeatures::transcripts(txdb)$tx_name
  expect_false(anyNA(tx))
  expect_false(anyDuplicated(tx) > 0)
  expect_true(all(c("gene-a1", "gene-a2", "ispD", "rna-x4") %in% tx))
  expect_equal(length(suppressMessages(loadRegion(txdb, "cds"))), 3)
})

test_that("loadTxdb names GFF transcripts without transcript_id when the file uses transcript_id", {
  gff <- tempfile(fileext = ".gff3")
  writeLines(c("##gff-version 3",
    "chr1\tRefSeq\tgene\t101\t400\t.\t+\t.\tID=gene-A;Name=A",
    "chr1\tRefSeq\tmRNA\t101\t400\t.\t+\t.\tID=rna-XM_1;Parent=gene-A;Name=XM_1;transcript_id=XM_1",
    "chr1\tRefSeq\texon\t101\t400\t.\t+\t.\tID=exon-XM_1-1;Parent=rna-XM_1;transcript_id=XM_1",
    "chr1\tRefSeq\tCDS\t101\t400\t.\t+\t0\tID=cds-1;Parent=rna-XM_1",
    "chr1\ttRNAscan-SE\tgene\t1001\t1073\t.\t+\t.\tID=gene-TRNAC;Name=TRNAC",
    "chr1\ttRNAscan-SE\ttRNA\t1001\t1073\t.\t+\t.\tID=rna-TRNAC;Parent=gene-TRNAC",
    "chr1\ttRNAscan-SE\texon\t1001\t1073\t.\t+\t.\tID=exon-TRNAC-1;Parent=rna-TRNAC",
    "chr1\tcmsearch\tpseudogene\t2001\t2117\t.\t+\t.\tID=gene-LOC1;Name=LOC1;pseudo=true",
    "chr1\tcmsearch\texon\t2001\t2117\t.\t+\t.\tID=id-LOC1;Parent=gene-LOC1"), gff)
  tx <- GenomicFeatures::transcripts(suppressWarnings(suppressMessages(loadTxdb(gff))))$tx_name
  expect_false(anyNA(tx))
  expect_true(all(c("XM_1", "rna-TRNAC", "gene-LOC1") %in% tx))
})
