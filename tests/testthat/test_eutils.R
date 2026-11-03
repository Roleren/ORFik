context("NCBI E-utilities access")
library(ORFik)

test_that("eutils_read_xml reads esearch and sends long id lists as POST", {
  skip_on_cran(); skip_on_bioc(); skip_if_offline("eutils.ncbi.nlm.nih.gov")
  url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1179332&retmode=xml"
  ids <- unlist(xml2::as_list(ORFik:::eutils_read_xml(url))$eSearchResult$IdList, use.names = FALSE)
  expect_gt(length(ids), 0)
  # an efetch URL longer than the GET limit is sent as POST and still returns every record
  long <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=fasta&retmode=xml&id=",
                 paste(c(ids, rep(ids[1], 300)), collapse = ","))
  expect_gt(nchar(long), 1500)
  doc <- ORFik:::eutils_read_xml(long)
  expect_equal(length(xml2::xml_find_all(doc, "//EXPERIMENT_PACKAGE")), length(ids))
})
