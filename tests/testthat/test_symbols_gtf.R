context("Gene symbol GTF fallback")
library(ORFik)

test_that("symbols_fill_from_gtf fills missing/id-only symbols, keeps biomaRt ones and column positions", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    'chr1\tt\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "Abc1";',
    'chr1\tt\tgene\t200\t300\t.\t+\t.\tgene_id "G2"; gene_name "Def2";',
    'chr1\tt\tgene\t400\t500\t.\t+\t.\tgene_id "G3"; gene_name "G3";',
    'chr1\tt\tgene\t600\t700\t.\t+\t.\tgene_id "G4";'), gtf)
  s <- data.table::data.table(ensembl_gene_id = c("G1", "G2", "G3", "G4", "G1"),
                              external_gene_name = c("Bio1", NA, "", "G4", "Bio1"),
                              ensembl_tx_name = paste0("T", 1:5), uniprot_id = NA)
  r <- symbols_fill_from_gtf(s, gtf)
  expect_identical(names(r), c(names(s), "symbol_source"))
  expect_identical(r$external_gene_name, c("Bio1", "Def2", "", "G4", "Bio1"))
  expect_identical(r$symbol_source, c("biomart", "gtf", "none", "none", "biomart"))
  expect_identical(s$external_gene_name[2], NA_character_)   # input not modified
})

test_that("symbols_fill_from_gtf reads NCBI GFF3 (gene- prefix, gene attribute, locus tag not a symbol)", {
  gff <- tempfile(fileext = ".gff3")
  writeLines(c("##gff-version 3",
    "c1\tGenbank\tgene\t1\t300\t.\t+\t.\tID=gene-ZMO1_ZMO0001;Name=rpsL;gene=rpsL;locus_tag=ZMO1_ZMO0001",
    "c1\tGenbank\tgene\t400\t700\t.\t+\t.\tID=gene-ZMO1_ZMO0002;Name=ZMO1_ZMO0002;locus_tag=ZMO1_ZMO0002"), gff)
  s <- data.table::data.table(ensembl_gene_id = c("ZMO1_ZMO0001", "ZMO1_ZMO0002"), external_gene_name = NA_character_,
                              ensembl_tx_name = c("t1", "t2"), uniprot_id = NA)
  r <- symbols_fill_from_gtf(s, gff)
  expect_identical(r$external_gene_name, c("rpsL", NA))
  expect_identical(r$symbol_source, c("gtf", "none"))
})
