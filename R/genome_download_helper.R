#' @inherit getGenomeAndAnnotation
#' @keywords internal
get_genome_fasta <- function(genome, output.dir, organism,
                             assembly_type, db, gunzip) {
  if (genome != FALSE) { # fasta genome of organism
    if (db == "ensembl") {
      message(paste("Starting", assembly_type, "genome retrieval of",
                    organism, "from ensembl: "))
      genome <- tryCatch(biomartr:::getENSEMBL.Seq(organism, type = "dna",
                                                   release = NULL,
                                                   id.type = assembly_type,
                                                   path = output.dir)[1],
                         error = function(e) {
                           return(e)
                         }
      )
      if (inherits(genome, "error")) {
        if (genome$message[1] == "Given file does not exist") {
        message("Could not find assembly_type: ", assembly_type)
        assembly_types_cand <- c("toplevel", "primary_assembly")
        assembly_type <- assembly_types_cand[!(assembly_types_cand %in% assembly_type)]
        message("Switching to search for assembly_type: ", assembly_type)
        genome <- biomartr:::getENSEMBL.Seq(organism, type = "dna",
                                            release = NULL,
                                            id.type = assembly_type,
                                            path = output.dir)[1]
        } else stop(genome)
      }

      if (is.logical(genome)) {
        if (genome == FALSE) {
          message("Remember some small genome organisms like yeast,",
                  " does not have primary assemblies, ",
                  "then change assembly_type to toplevel and/or use",
                  "db = refseq.")
        }
      }
      if (gunzip) # unzip gtf file
        genome <- R.utils::gunzip(genome, overwrite = TRUE)
    } else {
      genome  <- biomartr::getGenome(db = db, organism,
                                     path = output.dir, gunzip = gunzip)
    }

    if (any(genome == "Not available")) stop("Could not find genome, check spelling!")
    message("Making .fai index of genome")
    indexFa(genome)
    message("Genome download and fasta indexing complete")
  } else { # check if it already exists
    genome <- grep(pattern = organism,
                   x = list.files(output.dir, full.names = TRUE),
                   value = TRUE)
    if (db == "refseq") {
      genome <- grep(pattern = "\\.fna", x = genome, value = TRUE)
    } else genome <- grep(pattern = "\\.fa", x = genome, value = TRUE)
    genome <- grep(pattern = "\\.fai", x = genome, value = TRUE, invert = TRUE)
    if (length(genome) != 1) {
      if (length(genome) > 1) {
        if (db == "refseq") {
          genome <- grep(pattern = "_refseq", x = genome, value = TRUE)
        } else genome <- grep(pattern = "\\.dna", x = genome, value = TRUE)
      }
      if (length(genome) > 1) {
        warning("Found multiple candidates for pre downloaded genome,
                 setting to FALSE!
                 You can update path manually in the returned object")
        genome <- FALSE
      } else if (length(genome) == 0) genome <- FALSE
    }
  }

  return(genome)
}

#' @inherit getGenomeAndAnnotation
#' @param genome character path, default NULL.
#' Path to fasta genome, corresponding to the gtf. must be indexed
#' (.fai file must exist there).
#' If you want to make sure chromosome naming of the GTF matches the genome
#' and correct seqlengths. If value is NULL or FALSE, it will be ignored.
#' @inheritParams makeTxdbFromGenome
#' @keywords internal
get_genome_gtf <- function(GTF, output.dir, organism, assembly_type, db,
                           gunzip, genome, optimize = FALSE) {
  if (GTF != FALSE) { # gtf of organism
    if (db == "ensembl") {
      gtf <- biomartr:::getENSEMBL.gtf(organism = organism,
                                       type = "dna",
                                       id.type = assembly_type,
                                       path = output.dir)
    } else {
      message("Some refseq gffs are malformed, like Arabidopsis thaliana,",
      " and might crash during gff reading step!",
      " Until this is fixed in rtracklayer, check out example",
      " in ?getGenomeAndAnnotation on how to rescue those gffs")
      gtf <-biomartr::getGFF(db = db, organism = organism,
                        path = output.dir, reference = TRUE)
    }

    if (gunzip) # unzip gtf file
      gtf <- R.utils::gunzip(gtf, overwrite = TRUE)
    makeTxdbFromGenome(gtf, genome, organism, optimize)
  } else { # check if it already exists
    gtf <- grep(pattern = organism,
                x = list.files(output.dir, full.names = TRUE),
                value = TRUE)
    if (db == "ensembl") {
      gtf <- grep(pattern = "\\.gtf", x = gtf, value = TRUE)
    } else gtf <- grep(pattern = "\\.gff", x = gtf, value = TRUE)

    gtf <- grep(pattern = "\\.db", x = gtf, value = TRUE, invert = TRUE)
    if (length(gtf) != 1) {
      warning("Found multiple candidates for pre downloaded gtf,
              setting to FALSE!")
      gtf <- FALSE
    }
  }
  return(gtf)
}

#' @inherit getGenomeAndAnnotation
#' @keywords internal
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
#' @keywords internal
get_phix_genome <- function(phix, output.dir, gunzip) {
  if (phix) {
    message("Downloading phix genome")
    if (Sys.info()[1] == "Linux") { # Faster version for Linux

      phix <- paste0(output.dir, "/Escherichia_virus_phiX174.fa.gz")
      tryCatch(download.file(phix.url, destfile = phix,
                    method = "wget", extra = "--passive-ftp"),
               error = function(e, phix.url) {
                 phix.url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
                 download.file(phix.url, destfile = phix,
                               method = "wget")
               })
    } else {
      phix <- biomartr::getGenome(db = "refseq", "Escherichia phage phiX174",
                                  path = output.dir, gunzip = FALSE)
    }
    if (gunzip) # unzip phix file
      phix <- R.utils::gunzip(phix, overwrite = TRUE)
  }
  return(phix)
}

#' Download Silva SSU & LSU sequences
#'
#' Version downloaded is 138.1. NR99_tax (non redundant)
#'
#' If it fails from timeout, set higher timeout: options(timeout = 200)
#' @inheritParams get_phix_genome
#' @return filepath to downloaded file
#' @export
#' @examples
#' output.dir <- tempdir()
#' # get_silva_rRNA(output.dir)
get_silva_rRNA <- function(output.dir) {
  silva <- paste0(output.dir, "/rRNA_SSU&LSU_silva_138.1.fasta.gz")
  if (file.exists(silva)) {
    message("Silva rRNA file already found, will not download")
    return(silva)
  }
  message("Downloading silva rRNA SSU & LSU (version 138.1)")

  silva.ssu.url <- "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
  download.file(silva.ssu.url, destfile = silva)

  silva.lsu.url <- "https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
  download.file(silva.lsu.url, destfile = silva, mode = "a")

  return(silva)
}
