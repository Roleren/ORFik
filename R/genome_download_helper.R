#' @inherit getGenomeAndAnnotation
get_genome_fasta <- function(genome, output.dir, organism,
                             assembly_type, db, gunzip) {
  if (genome != FALSE) { # fasta genome of organism
    if (db == "ensembl") {
      message(paste("Starting", assembly_type, "genome retrieval of",
                    organism, "from ensembl: "))
      genome <- biomartr:::getENSEMBL.Seq(organism, type = "dna",
                                          release = NULL,
                                          id.type = assembly_type,
                                          path = output.dir)[1]
      if (genome == FALSE) {
        message("Remember some small genome organisms like yeast,",
                " does not have primary assemblies, ",
                "then change assembly_type to toplevel or use",
                "db = refseq.")
      }
      if (gunzip) # unzip gtf file
        genome <- R.utils::gunzip(genome, overwrite = TRUE)
    } else {
      genome  <- biomartr::getGenome(db = db, organism,
                                     path = output.dir, gunzip = gunzip)
    }

    message("Making .fai index of genome")
    indexFa(genome)
  } else { # check if it already exists
    genome <- grep(pattern = organism,
                   x = list.files(output.dir, full.names = TRUE),
                   value = TRUE)
    genome <- grep(pattern = "\\.fa", x = genome, value = TRUE)
    genome <- grep(pattern = "\\.fai", x = genome, value = TRUE, invert = TRUE)
    if (length(genome) != 1) genome <- FALSE
  }
  return(genome)
}

#' @inherit getGenomeAndAnnotation
get_genome_gtf <- function(GTF, output.dir, organism, assembly_type, gunzip) {
  if (GTF) { # gtf of organism
    gtf <- biomartr:::getENSEMBL.gtf(organism = organism,
                                     type = "dna",
                                     id.type = assembly_type,
                                     path = output.dir)

    if (gunzip) # unzip gtf file
      gtf <- R.utils::gunzip(gtf, overwrite = TRUE)
    message("Making txdb of GTF")
    organismCapital <- paste0(toupper(substr(organism, 1, 1)),
                              substr(organism, 2, nchar(organism)))
    organismCapital <- gsub("_", " ", organismCapital)
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, organism = organismCapital)
    if (genome != FALSE)
      seqlevelsStyle(txdb) <- seqlevelsStyle(FaFile(genome))[1]
    txdb_file <- paste0(gtf, ".db")
    AnnotationDbi::saveDb(txdb, txdb_file)
  } else { # check if it already exists
    gtf <- grep(pattern = organism,
                x = list.files(output.dir, full.names = TRUE),
                value = TRUE)
    gtf <- grep(pattern = "\\.gtf", x = gtf, value = TRUE)
    gtf <- grep(pattern = "\\.db", x = gtf, value = TRUE, invert = TRUE)
    if (length(gtf) != 1) gtf <- FALSE
  }
  return(gtf)
}

#' @inherit getGenomeAndAnnotation
get_noncoding_rna <- function(ncRNA, output.dir, organism, gunzip) {
  if (is.logical(ncRNA)) return(ncRNA)

  if (ncRNA != "") {
    if (ncRNA == "auto") {
      a <- biomartr::getENSEMBLInfo()
      ncRNA <- a[grep(organism, a$name, TRUE),]$common_name
      if (length(ncRNA) == 0) stop("ncRNA was auto,",
                                   "but could not find organism")
      message(paste0("ncRNA auto: organism common name:",
                     ncRNA))
    }
    message("Downloading ncRNA's")
    file <- "http://www.noncode.org/datadownload/NONCODEv5_"
    org <- ncRNA
    extension <- ".fa.gz"
    out <- paste0(output.dir, "/NONCODE_ncRNA_",org, extension)
    download.file(paste0(file, org, extension), destfile = out)
    ncRNA <- out
    if (gunzip) # unzip gtf file
      ncRNA <- R.utils::gunzip(ncRNA, overwrite = TRUE)
  }
  return(ncRNA)
}

#' @inherit getGenomeAndAnnotation
get_phix_genome <- function(phix, output.dir, gunzip) {
  if (phix) {
    message("Downloading phix genome")
    if (Sys.info()[1] == "Linux") { # Faster version for Linux
      phix.url <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Escherichia_virus_phiX174/all_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
      phix <- paste0(output.dir, "/Escherichia_virus_phiX174.fa.gz")
      download.file(phix.url, destfile = phix,
                    method = "wget", extra = "--passive-ftp")
    } else {
      phix <- biomartr::getGenome(db = "refseq", "Escherichia virus phiX174",
                                  path = output.dir, gunzip = FALSE)
    }
    if (gunzip) # unzip phix file
      phix <- R.utils::gunzip(phix, overwrite = TRUE)
  }
  return(phix)
}
