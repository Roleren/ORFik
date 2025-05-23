#' Download genome (fasta), annotation (GTF) and contaminants
#'
#' This function automatically downloads (if files not already exists)
#' genomes and contaminants specified for genome alignment.
#' By default, it will use ensembl reference,
#' upon completion, the function will store
#' a file called \code{file.path(output.dir, "outputs.rds")} with
#' the output paths of your completed genome/annotation downloads.
#' For most non-model nonvertebrate organisms, you need
#' my fork of biomartr for it to work:
#' remotes::install_github("Roleren/biomartr)
#' If you misspelled something or crashed, delete wrong files and
#' run again.\cr
#' Do remake = TRUE, to do it all over again.\cr
#'
#' Some files that are made after download:\cr
#' - A fasta index for the genome\cr
#' - A TxDb to speed up GTF/GFF reading\cr
#' - Seperat of merged contaminant files\cr
#' Files that can be made:\cr
#' - Gene symbols (hgnc, etc)\cr
#' - Uniprot ids (For name of protein structures)\cr
#' If you want custom genome or gtf from you hard drive, assign existing
#' paths like this: \cr
#' annotation <- getGenomeAndAnnotation(GTF = "path/to/gtf.gtf",
#' genome = "path/to/genome.fasta")\cr
#' @inheritParams biomartr::getGenome
#' @param organism scientific name of organism, Homo sapiens,
#' Danio rerio, Mus musculus, etc. See \code{biomartr:::get.ensembl.info()}
#' for full list of supported organisms.
#' @param output.dir directory to save downloaded data
#' @param db database to use for genome and GTF,
#' default adviced: "ensembl" (remember to set assembly_type to "primary_assembly",
#' else it will contain haplotypes, very large file!).
#' Alternatives: "refseq" (reference assemblies) and "genbank" (all assemblies)
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
#' if TRUE or defned path, ncRNA is used as a contaminant reference.
#' If TRUE, will try to find ncRNA sequences from the gtf file, usually represented as
#' lncRNA (long noncoding RNA's). Will let you know if no ncRNA sequences were found in
#' gtf.\cr If not found try character input:\cr
#' Alternatives; "auto":
#' Will try to find ncRNA file on NONCODE from organism,
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
#' rRNA is used as a contaminant reference
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
#' @param assembly character, default is assembly = organism, which means getting
#' the first assembly in list, otherwise the name of the assembly wanted, like
#' "GCA_000005845" will get ecoli substrain k12, which is the most used ones for
#' references. Usually ignore this for non bacterial species.
#' @param refseq_genbank_format = c("gtf", "gff3")[1] Gtf format files are usually
#' more secure from bugs downstream, so we highly advice to use them. GFF3 files
#' can sometimes include information you might not find in the gtf, so sometimes
#' it makes sense to use it.
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
#' @references https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4919035/
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
#' # Drosophila melanogaster (toplevel exists only)
#' #getGenomeAndAnnotation("drosophila melanogaster", output.dir = file.path(config["ref"],
#' # "Drosophila_melanogaster_BDGP6"), assembly_type = "toplevel")
#' ## How to save malformed refseq gffs:
#' ## First run function and let it crash:
#' #annotation <- getGenomeAndAnnotation(organism = "Arabidopsis thaliana",
#' #  output.dir = "~/Desktop/test_plant/",
#' #  assembly_type = "primary_assembly", db = "refseq")
#' ## Then apply a fix (example for linux, too long rows):
#' # fixed_gff <- fix_malformed_gff("~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq.gff")
#' ## Then updated arguments:
#' # annotation <- c(fixed_gff, "~/Desktop/test_plant/Arabidopsis_thaliana_genomic_refseq.fna")
#' # names(annotation) <- c("gtf", "genome")
#' # Then make the txdb (for faster R use)
#' # makeTxdbFromGenome(annotation["gtf"], annotation["genome"], organism = "Arabidopsis thaliana")
getGenomeAndAnnotation <- function(organism, output.dir, db = "ensembl",
                                   GTF = TRUE, genome = TRUE,
                                   merge_contaminants = TRUE, phix = FALSE,
                                   ncRNA = FALSE, tRNA = FALSE, rRNA = FALSE,
                                   gunzip = TRUE, remake = FALSE,
                                   assembly_type = c("primary_assembly", "toplevel"),
                                   optimize = FALSE, gene_symbols = FALSE,
                                   uniprot_id = FALSE,
                                   pseudo_5UTRS_if_needed = NULL,
                                   remove_annotation_outliers = TRUE,
                                   notify_load_existing = TRUE,
                                   assembly = organism,
                                   refseq_genbank_format = c("gtf", "gff3")[1]) {
  # Pre checks
  stopifnot(is(organism, "character"))
  stopifnot(is(output.dir, "character"))
  finished.file <- file.path(output.dir, "outputs.rds")
  if (file.exists(finished.file) & !remake) {
    if (notify_load_existing) message("Loading premade Genome files,",
                                  " do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }
  if (!all(assembly_type %in% c("toplevel", "primary_assembly")))
    stop("Please select one the available assembly types: \ntoplevel, primary_assembly")
  dir.create(output.dir, recursive = TRUE)

  # Start process
  organism <- gsub(" ", "_", organism)
  if (db == "refseq") {
    organism <- gsub("_", " ", organism)
  }
  ## Download all contaminants wanted:
  conts <- contaminants_download(tRNA, rRNA, phix, ncRNA, output.dir, organism,
                                 gunzip)
  # Get species fasta genome and gtf
  genome <- get_genome_fasta(genome, output.dir, organism, assembly,
                             assembly_type, db, gunzip)
  gtf <- get_genome_gtf(GTF, output.dir, organism, assembly,
                        db, gunzip, genome, optimize = optimize,
                        uniprot_id = uniprot_id,
                        gene_symbols = gene_symbols,
                        pseudo_5UTRS_if_needed = pseudo_5UTRS_if_needed,
                        remove_annotation_outliers = remove_annotation_outliers,
                        refseq_genbank_format = refseq_genbank_format)
  output <- contaminants_processing(conts, gtf, genome, merge_contaminants,
                                    output.dir)

  message("All data downloaded and ready at:")
  message(output.dir)
  saveRDS(object = output, finished.file)
  return(output)
}
