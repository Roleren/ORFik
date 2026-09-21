context("SRA metadata")
library(ORFik)

test_that("sample_info_single rescues run attributes when Statistics node is missing", {
  xml <- xml2::read_xml('<EXPERIMENT_PACKAGE>
    <EXPERIMENT>
      <IDENTIFIERS><PRIMARY_ID>SRX207952</PRIMARY_ID></IDENTIFIERS>
      <STUDY_REF><IDENTIFIERS><EXTERNAL_ID>PRJNA182756</EXTERNAL_ID></IDENTIFIERS></STUDY_REF>
      <DESIGN><LIBRARY_DESCRIPTOR>
        <LIBRARY_STRATEGY>OTHER</LIBRARY_STRATEGY>
        <LIBRARY_SELECTION>other</LIBRARY_SELECTION>
        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
        <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR></DESIGN>
      <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 2000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
    </EXPERIMENT>
    <STUDY><IDENTIFIERS><PRIMARY_ID>SRP017378</PRIMARY_ID><EXTERNAL_ID>PRJNA182756</EXTERNAL_ID></IDENTIFIERS></STUDY>
    <SAMPLE alias="GSM1047591">
      <IDENTIFIERS><PRIMARY_ID>SRS377853</PRIMARY_ID><EXTERNAL_ID>SAMN01821863</EXTERNAL_ID></IDENTIFIERS>
      <SAMPLE_NAME><TAXON_ID>9606</TAXON_ID><SCIENTIFIC_NAME>Homo sapiens</SCIENTIFIC_NAME></SAMPLE_NAME>
      <TITLE>Tr.rp</TITLE>
      <SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>source_name</TAG><VALUE>Immortalized primary fibroblasts</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES>
    </SAMPLE>
    <SUBMISSION center_name="GEO" accession="SRA062117"/>
    <Organization><Contact><Name><Last>Loayza-Puch</Last></Name></Contact></Organization>
    <RUN_SET><RUN accession="SRR627627" total_spots="178192873" total_bases="9087836523" size="6128138982" published="2013-03-25 08:27:41">
      <IDENTIFIERS><PRIMARY_ID>SRR627627</PRIMARY_ID></IDENTIFIERS>
      <EXPERIMENT_REF accession="SRX207952"/>
    </RUN></RUN_SET>
  </EXPERIMENT_PACKAGE>')

  dt <- ORFik:::sample_info_single(xml2::as_list(xml)$EXPERIMENT_PACKAGE)

  expect_equal(dt$Run, "SRR627627")
  expect_equal(dt$spots, 178192873L)
  expect_equal(dt$bases, 9087836523)
  expect_equal(dt$avgLength, 51L)
  expect_equal(dt$size_MB, floor(6128138982 / 1024^2))
  expect_equal(dt$Experiment, "SRX207952")
})

test_that("ebi_paths_to_ascp handles FTP and native ENA Aspera paths", {
  paths <- c(
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR627/SRR627626/SRR627626.fastq.gz",
    "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR627/SRR627626/SRR627626.fastq.gz"
  )

  expect_equal(ORFik:::ebi_paths_to_ascp(paths), c(
    "era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR627/SRR627626/SRR627626.fastq.gz",
    "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR627/SRR627626/SRR627626.fastq.gz"
  ))
})
