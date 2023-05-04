#' @inherit getGenomeAndAnnotation
#' @keywords internal
get_genome_fasta <- function(genome, output.dir, organism,
                             assembly_type, db, gunzip) {
  if (genome != FALSE) { # If not auto detect
    if (is.logical(genome)) { # Then download
      message("- Download genome (fasta)")
      if (db == "ensembl") {
        message(paste("Starting", assembly_type, "genome retrieval of",
                      organism, "from ensembl: "))
        genome <- tryCatch(getENSEMBL.Seq(organism, type = "dna",
                                          release = NULL,
                                          id.type = assembly_type,
                                          path = output.dir)[1],
                           error = function(e) {
                             return(e)
                           }
        )
        if (inherits(genome, "error") | is.logical(genome)) {
          swap_assembly_type <- is.logical(genome)
          if (!swap_assembly_type)
            swap_assembly_type <- genome$message[1] == "Given file does not exist"
          if (swap_assembly_type) {
          message("Could not find assembly_type: ", assembly_type)
          assembly_types_cand <- c("toplevel", "primary_assembly")
          assembly_type <- assembly_types_cand[!(assembly_types_cand %in% assembly_type)]
          message("Switching to search for assembly_type: ", assembly_type)
          genome <- getENSEMBL.Seq(organism, type = "dna",
                                   release = NULL, id.type = assembly_type,
                                   path = output.dir)[1]
          } else stop(genome)
        }

        if (is.logical(genome)) {
          if (genome == FALSE) {
            message("Remember some small genome organisms like yeast,",
                    " does not have primary assemblies, ",
                    "then change assembly_type to toplevel and/or use:",
                    " db = refseq.")
          }
        } else {
          if (gunzip) # unzip gtf file
            genome <- R.utils::gunzip(genome, overwrite = TRUE)
        }

      } else {
        genome  <- biomartr::getGenome(db = db, organism,
                                       path = output.dir, gunzip = gunzip)
      }
    } else if (is.character(genome)) { # Use hard drive file
      message("- Using pre-existing genome (fasta)")
      stopifnot(file.exists(genome))
      is_compressed <- tools::file_ext(genome) == "gz"
      if (is_compressed) genome <- R.utils::gunzip(genome, overwrite = TRUE)
    }

    if (any(genome == "Not available") | is.logical(genome))
      stop("Could not find genome, check spelling!")
    message("Making .fai index of genome")
    indexFa(genome)
    message("Genome fetch and fasta indexing complete")
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
                           gunzip, genome, optimize = FALSE,
                           uniprot_id = FALSE,
                           gene_symbols = FALSE,
                           pseudo_5UTRS_if_needed = NULL,
                           remove_annotation_outliers = TRUE) {
  if (GTF != FALSE) { # If GTF should be processed
    if (is.logical(GTF)) { # If download
      message("- Download annotation (gtf/gff3)")
      is_compressed <- gunzip
      if (db == "ensembl") {
        gtf <- getENSEMBL.gtf(organism = organism, type = "dna",
                              id.type = assembly_type, path = output.dir)
      } else {
        message("Some refseq gffs are malformed, like Arabidopsis thaliana,",
                " and might crash during gff reading step!",
                " Until this is fixed in rtracklayer, check out example",
                " in ?getGenomeAndAnnotation on how to rescue those gffs")
        gtf <-biomartr::getGFF(db = db, organism = organism,
                               path = output.dir, reference = TRUE,
                               remove_annotation_outliers = remove_annotation_outliers)
      }
    } else if (is.character(GTF)) { # File already exists
      message("- Using pre-existing annotation (gtf/gff3)")
      gtf <- GTF # It already exists
      stopifnot(file.exists(gtf))
      is_compressed <- tools::file_ext(gtf) == "gz"
    }

    if (is_compressed) # unzip gtf file
      gtf <- R.utils::gunzip(gtf, overwrite = TRUE)
    makeTxdbFromGenome(gtf, genome, organism, optimize, gene_symbols,
                       uniprot_id, pseudo_5UTRS_if_needed)
  } else { # Try to auto detect, else set to FALSE
    gtf <- auto_detect_gtf_in_dir(organism, output.dir, db)
  }
  return(gtf)
}

auto_detect_gtf_in_dir <- function(organism, output.dir, db) {
  gtf <- grep(pattern = organism,
              x = list.files(output.dir, full.names = TRUE),
              value = TRUE)
  if (db == "ensembl") {
    gtf <- grep(pattern = "\\.gtf|\\.gff3", x = gtf, value = TRUE)
  } else gtf <- grep(pattern = "\\.gff", x = gtf, value = TRUE)

  gtf <- grep(pattern = "\\.db$", x = gtf, value = TRUE, invert = TRUE)
  if (length(gtf) != 1) {
    warning("Found multiple candidates for pre downloaded gtf,
              setting to FALSE!")
    gtf <- FALSE
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

get_rRNA <- function(rRNA, output.dir) {
  if (!(rRNA %in% c("", FALSE, TRUE, "silva"))) {
    if (!file.exists(rRNA)) stop(paste("local rRNA file is specified, but does not exist:",
                                       rRNA))
  } else if (rRNA == "silva") rRNA <- get_silva_rRNA(output.dir)
  return(rRNA)
}

get_tRNA <- function(tRNA) {
  if (!(tRNA %in% c("", FALSE, TRUE))) {
    if (!file.exists(tRNA)) stop(paste("local tRNA file is specified, but does not exist:",
                                       tRNA))
  }
  return(tRNA)
}

contaminants_download <- function(tRNA, rRNA, phix, ncRNA, output.dir,
                                  organism, gunzip) {
  tRNA <- get_tRNA(tRNA)
  rRNA <- get_rRNA(rRNA, output.dir)
  phix <- get_phix_genome(phix, output.dir, gunzip)
  ncRNA <- get_noncoding_rna(ncRNA, output.dir, organism, gunzip)
  conts <- as.character(c(phix, ncRNA, tRNA, rRNA))
  names(conts) <- c("phix", "ncRNA", "tRNA", "rRNA")
  return(conts)
}

contaminants_processing <- function(conts, gtf, genome,
                                   merge_contaminants, output.dir) {
  # If any defined contaminants
  any_contaminants <- any(!(conts %in% c("", FALSE)))
  if (any_contaminants) {
    message("- Processing contaminants")
    # Find which contaminants to find from gtf:
    non_gtf_contaminants <- conts[!(conts %in% c(TRUE, FALSE, ""))]
    gtf_contaminants <- conts[conts == TRUE]
    total_seqs <- DNAStringSet()
    if (length(gtf_contaminants) > 0) {
      if (is.logical(gtf) | is.logical(genome))
        stop("gtf or genome not specified, so impossible to find gtf contaminants!")
      # Make fasta file of those contaminants
      gtf.imp <- importGtfFromTxdb(gtf)
      txdb <- loadTxdb(paste0(gtf, ".db"))
      tx <- loadRegion(txdb)
      # Loop through the gtf contaminants
      for (t in names(gtf_contaminants)) {
        valids <- gtf.imp[grep(x = gtf.imp$transcript_biotype, pattern = t)]
        if (length(gtf.imp) == 0) {
          warning(paste("Found no transcripts in gtf of type:", t))
          next
        }
        valid_tx <- tx[unique(valids$transcript_id)]
        seqs <- txSeqsFromFa(valid_tx, genome, TRUE, TRUE)
        path.cont <- paste0(output.dir, "/", t, ".fasta")
        writeXStringSet(seqs, path.cont)
        total_seqs <- c(total_seqs, seqs)
        gtf_contaminants[t] <- conts[t] <- path.cont
      }
    }
  }

  if (merge_contaminants & any_contaminants) {
    message("Merging contaminant genomes:")
    for (cont in non_gtf_contaminants)
      total_seqs <- c(total_seqs, readDNAStringSet(cont))

    all_contaminants <- paste(c(names(non_gtf_contaminants),
                                names(gtf_contaminants)), collapse = "_")
    all_cont <- paste0("merged_contaminants_", all_contaminants, ".fasta")
    all_cont <- file.path(output.dir, all_cont)
    writeXStringSet(total_seqs, filepath = all_cont)

    output <- c(gtf, genome, all_cont)
    names(output) <- c("gtf", "genome", "contaminants")
  } else {
    output <- c(gtf, genome, conts)
    names(output) <- c("gtf", "genome", names(conts))
  }
  output <- output[!(output %in% c("", "TRUE","FALSE"))]
  return(output)
}
