test_that("genome download errors retain the underlying service failure", {
  testthat::local_mocked_bindings(
    getENSEMBL.Seq = function(...) stop("Ensembl service unavailable"),
    .package = "biomartr"
  )
  expect_error(
    ORFik:::get_genome_fasta(TRUE, tempdir(), "test_species", "test_species",
                            "toplevel", "ensembl", TRUE),
    "Genome download failed for test_species: Ensembl service unavailable"
  )
})

test_that("unsuccessful genome downloads are diagnosed before indexing", {
  for (result in list(FALSE, character(), NA_character_, "", "Not available")) {
    testthat::local_mocked_bindings(
      getENSEMBL.Seq = function(...) result,
      .package = "biomartr"
    )
    expect_error(
      ORFik:::get_genome_fasta(TRUE, tempdir(), "test_species", "test_species",
                              "toplevel", "ensembl", TRUE),
      "No genome file was returned for test_species"
    )
  }
})

test_that("successful downloads and local FASTA inputs are indexed", {
  reference <- system.file("extdata/references/homo_sapiens",
                           "Homo_sapiens_dummy.fasta", package = "ORFik")
  output <- withr::local_tempdir()
  genome <- file.path(output, "genome.fa")
  file.copy(reference, genome)
  testthat::local_mocked_bindings(
    getENSEMBL.Seq = function(...) genome,
    .package = "biomartr"
  )
  expect_identical(
    ORFik:::get_genome_fasta(TRUE, output, "Homo_sapiens", "Homo_sapiens",
                            "toplevel", "ensembl", FALSE), genome
  )
  expect_true(file.exists(paste0(genome, ".fai")))
  expect_identical(
    ORFik:::get_genome_fasta(genome, output, "Homo_sapiens", "Homo_sapiens",
                            "toplevel", "ensembl", TRUE), genome
  )
})
