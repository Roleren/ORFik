context("NCBI E-utilities access")
library(ORFik)

test_that("eutils_read_xml reads esearch and sends long id lists as POST", {
  skip_on_cran(); skip_on_bioc(); skip_if_offline("eutils.ncbi.nlm.nih.gov")
  url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA521877&retmax=1000&retmode=xml"
  ids <- unlist(xml2::as_list(ORFik:::eutils_read_xml(url))$eSearchResult$IdList, use.names = FALSE)
  expect_gt(length(ids), 0)
  # PRJNA521877 has 512 SRA experiments: the efetch URL is far longer than the GET limit, so it is sent as POST
  long <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=fasta&retmode=xml&id=",
                 paste(ids, collapse = ","))
  expect_gt(nchar(long), 1500)
  doc <- ORFik:::eutils_read_xml(long)
  expect_equal(length(xml2::xml_find_all(doc, "//EXPERIMENT_PACKAGE")), length(ids))
})
