#' Download genome (fasta), annotation (GTF) and contaminants
#'
#' This function automatically downloads (if files not already exists)
#' genomes and contaminants specified for genome alignment.
#' Will create a R transcript database (TxDb object) from the annotation. \cr
#' It will also index the genome for you\cr
#' If you misspelled something or crashed, delete wrong files and
#' run again.\cr
#' Do remake = TRUE, to do it all over again.
#'
#' If you want custom genome or gtf from you hard drive, assign it
#' after you run this function, like this:\cr
#' annotation <- getGenomeAndAnnotation(GTF = FALSE, genome = FALSE)\cr
#' annotation["genome"] = "path/to/genome.fasta"\cr
#' annotation["gtf"] = "path/to/gtf.gtf"
#' @param organism scientific name of organism, Homo sapiens,
#' Danio rerio, Mus musculus, etc. See \code{biomartr:::get.ensembl.info()}
#' for full list of supported organisms.
#' @param output.dir directory to save downloaded data
#' @param db database to use for genome and GTF,
#' default adviced: "ensembl" (remember to set assembly_type to "primary_assembly",
#' else it will contain haplotypes, very large file!).
#' Alternatives: "refseq" (primary assembly) and "genbank" (mix)
#' @param GTF logical, default: TRUE, download gtf of organism specified
#' in "organism" argument. If FALSE, check if the downloaded
#' file already exist. If you want to use a custom gtf from you hard drive,
#' set GTF = FALSE,
#' and assign: \cr annotation <- getGenomeAndAnnotation(gtf = FALSE)\cr
#' annotation["gtf"] = "path/to/gtf.gtf".\cr
#' If db is not "ensembl", you will instead get a gff file.
#' @param genome logical, default: TRUE, download genome of organism
#' specified in "organism" argument. If FALSE, check if the downloaded
#' file already exist. If you want to use a custom gtf from you hard drive,
#' set \code{GTF = FALSE},
#' and assign: \cr \code{annotation <- getGenomeAndAnnotation(genome = FALSE)}\cr
#' \code{annotation["genome"] = "path/to/genome.fasta"}.\cr
#' Will download the primary assembly from Ensembl.
#' @param merge_contaminants logical, default TRUE. Will merge
#' the contaminants specified into one fasta file, this considerably
#' saves space and is much quicker to align with STAR than each contaminant
#' on it's own. If no contaminants are specified, this is ignored.
#' @param phix logical, default FALSE, download phiX sequence to filter
#'  out Illumina control reads. ORFik defines Phix as a contaminant genome.
#' Phix is used in Illumina sequencers for sequencing quality control.
#' Genome is: refseq, Escherichia phage phiX174.
#' If sequencing facility created fastq files with the command \code{bcl2fastq},
#' then there should be very few phix reads left in the fastq files recieved.
#' @param ncRNA logical or character, default FALSE (not used, no download),
#' ncRNA is used as a contaminant genome.
#' If TRUE, will try to find ncRNA sequences from the gtf file, usually represented as
#' lncRNA (long noncoding RNA's). Will let you know if no ncRNA sequences were found in
#' gtf.\cr If not found try character input:\cr
#' Alternatives: "auto" or manual assign like "human".
#' If "auto" will try to find ncRNA file on NONCODE from organism,
#' Homo sapiens -> human etc. "auto" will not work for all,
#' then you must specify the name used by
#' NONCODE, go to the link below and find it.
#' If not "auto" / "" it must be a character vector
#' of species common name (not scientific name) Homo sapiens is human,
#' Rattus norwegicus is rat etc, download ncRNA sequence to filter out with.
#' From NONCODE online server, if you cant find
#' common name see: http://www.noncode.org/download.php/
#' @param tRNA logical or character, default FALSE (not used, no download),
#' tRNA is used as a contaminant genome.
#' If TRUE, will try to find tRNA sequences from the gtf file, usually represented as
#' Mt_tRNA (mature tRNA's). Will let you know if no tRNA sequences were found in
#' gtf. If not found try character input:\cr
#' if not "" it must be a character vector  to valid path of mature
#' tRNAs fasta file to remove as contaminants on your disc. Find and download
#' your wanted mtRNA at: http://gtrnadb.ucsc.edu/, or run trna-scan on
#' you genome.
#' @param rRNA logical or character, default FALSE (not used, no download),
#' rRNA is used as a contaminant genome.
#' If TRUE, will try to find rRNA sequences from the gtf file, usually represented as
#' rRNA (ribosomal RNA's). Will let you know if no rRNA sequences were found in
#' gtf. If not found you can try character input:\cr If "silva" will download silva SSU & LSU
#' sequences for all species (250MB file) and use that. If you want a smaller file go to
#' https://www.arb-silva.de/ \cr
#' If not "" or "silva" it must be a character vector to valid path of mature
#' rRNA fasta file to remove as contaminants on your disc.
#' @param gunzip logical, default TRUE, uncompress downloaded files
#' that are zipped when downloaded, should be TRUE!
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @param assembly_type a character string specifying from which assembly type
#' the genome shall be retrieved from (ensembl only, else this argument is ignored):
#' Default is \code{assembly_type = "primary_assembly")}.
#' This will give you no haplotypes (copies of the same chromosome with small variations).
#' As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed.\cr
#' \code{assembly_type = "toplevel")}.
#' This will give you all haplotypes.
#' As an example the toplevel fasta genome in human is over 70 GB uncompressed.
#' @param remove_annotation_outliers logical, default TRUE. Only for refseq.
#'  shall outlier lines be removed from the input annotation_file?
#'  If yes, then the initial annotation_file will be overwritten and
#'  the removed outlier lines will be stored at tempdir for further
#'  exploration. Among others Aridopsis refseq contains malformed lines,
#'  where this is needed
#' @param notify_load_existing logical, default TRUE. If annotation exists
#' (defined as: locally (a file called outputs.rds) exists in outputdir),
#' print a small message notifying the user it is not redownloading. Set to
#' FALSE, if this is not wanted
#' @inheritParams makeTxdbFromGenome
#' @importFrom biomartr getGTF getGenome getENSEMBLInfo
#' @importFrom Rsamtools indexFa
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom AnnotationDbi saveDb
#' @importFrom Biostrings DNAStringSet writeXStringSet readDNAStringSet
#' @return a named character vector of path to genomes and gtf downloaded,
#'  and additional contaminants if used. If merge_contaminants is TRUE, will not
#'  give individual fasta files to contaminants, but only the merged one.
#' @family STAR
#' @export
#' @examples
#'
#' ## Get Saccharomyces cerevisiae genome and gtf (create txdb for R)
#' #getGenomeAndAnnotation("Saccharomyces cerevisiae", tempdir(), assembly_type = "toplevel")
#' ## Download and add pseudo 5' UTRs
#' #getGenomeAndAnnotation("Saccharomyces cerevisiae", tempdir(), assembly_type = "toplevel",
#' #  pseudo_5UTRS_if_needed = 100)
#' ## Get Danio rerio genome and gtf (create txdb for R)
#' #getGenomeAndAnnotation("Danio rerio", tempdir())
#'
#' output.dir <- "/Bio_data/references/zebrafish"
#' ## Get Danio rerio and Phix contamints to deplete during alignment
#' #getGenomeAndAnnotation("Danio rerio", output.dir, phix = TRUE)
#'
#' ## Optimize for ORFik (speed up for large annotations like human or zebrafish)
#' #getGenomeAndAnnotation("Danio rerio", tempdir(), optimize = TRUE)
#'
#' ## How to save malformed refseq gffs:
#' ## First run function and let it crash:
#' #annotation <- getGenomeAndAnnotation(organism = "Arabidopsis thaliana", output.dir = "~/Desktop/test_plant/",
#' #  assembly_type = "primary_assembly", db = "refseq")
#' ## Then apply a fix (example for linux, too long rows):
#' # \code{system("cat ~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq.gff | awk '{ if (length($0) < 32768) print }' > ~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq_trimmed2.gff")}
#' ## Then updated arguments:
#' annotation <- c("~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq_trimmed.gff",
#'  "~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq.fna")
#' names(annotation) <- c("gtf", "genome")
#' # Then make the txdb (for faster R use)
#' # makeTxdbFromGenome(annotation["gtf"], annotation["genome"], organism = "Arabidopsis thaliana")
getGenomeAndAnnotation <- function(organism, output.dir, db = "ensembl",
                                   GTF = TRUE, genome = TRUE,
                                   merge_contaminants = TRUE, phix = FALSE,
                                   ncRNA = FALSE, tRNA = FALSE, rRNA = FALSE,
                                   gunzip = TRUE, remake = FALSE,
                                   assembly_type = "primary_assembly",
                                   optimize = FALSE,
                                   gene_symbols = FALSE,
                                   pseudo_5UTRS_if_needed = NULL,
                                   remove_annotation_outliers = TRUE,
                                   notify_load_existing = TRUE) {
  # Pre checks
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) & !remake) {
    if (notify_load_existing) message("Loading premade Genome files,
                                       do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }
  if (!(assembly_type %in% c("toplevel", "primary_assembly")))
    stop("Please select one the available assembly types: \ntoplevel, primary_assembly")
  conts <- as.character(c(phix, ncRNA, tRNA, rRNA))
  any_contaminants <- any(!(conts %in% c("", FALSE)))
  # Start process
  dir.create(output.dir, recursive = TRUE)
  organism <- gsub(" ", "_", organism)
  if (db == "refseq") {
    organism <- gsub("_", " ", organism)
  }
  ## Go through all contaminants:
  tRNA <- get_tRNA(tRNA)
  rRNA <- get_rRNA(rRNA)
  phix <- get_phix_genome(phix, output.dir, gunzip)
  ncRNA <- get_noncoding_rna(ncRNA, output.dir, organism, gunzip)
  # Get species fasta genome and gtf
  genome <- get_genome_fasta(genome, output.dir, organism,
                             assembly_type, db, gunzip)
  gtf <- get_genome_gtf(GTF, output.dir, organism, assembly_type, db,
                        gunzip, genome, optimize = optimize,
                        gene_symbols = gene_symbols,
                        pseudo_5UTRS_if_needed = pseudo_5UTRS_if_needed,
                        remove_annotation_outliers = remove_annotation_outliers)

  if (any_contaminants) {
    # Find which contaminants to find from gtf:
    conts <- as.character(c(phix, ncRNA, tRNA, rRNA))
    names(conts) <- c("phix", "ncRNA", "tRNA", "rRNA")
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
        gtf_contaminants[t] <- path.cont
      }
    }
  }

  if (merge_contaminants & any_contaminants) {
    message("Merging contaminant genomes:")
    for (cont in non_gtf_contaminants)
      total_seqs <- c(total_seqs, readDNAStringSet(cont))

    all_contaminants <- paste(c(names(non_gtf_contaminants),
                                names(gtf_contaminants)), collapse = "_")

    all_cont <- paste0(output.dir, "/", "merged_contaminants_", all_contaminants, ".fasta")
    writeXStringSet(total_seqs, filepath = all_cont)
    output <- c(gtf, genome, all_cont)
    names(output) <- c("gtf", "genome", "contaminants")
  } else {
    output <- c(gtf, genome, phix, ncRNA, tRNA, rRNA)
    names(output) <- c("gtf", "genome", "phix", "ncRNA", "tRNA", "rRNA")
  }
  output <- output[!(output %in% c("", "TRUE","FALSE"))]
  message("All data downloaded and ready at:")
  message(output.dir)


  saveRDS(object = output, finished.file)
  return(output)
}
