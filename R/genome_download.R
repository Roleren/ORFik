#' Download genome (fasta), annotation (GTF) and contaminants
#'
#' This function downloads the genome, annotation, and contaminant references
#' used for alignment if they are not already present. By default it uses the
#' Ensembl reference. When the run completes, it stores a file called
#' \code{file.path(output.dir, "outputs.rds")} with the output paths for the
#' downloaded files.
#'
#' For most non-model non-vertebrate organisms, you need
#' my fork of biomartr for it to work:
#' remotes::install_github("Roleren/biomartr")
#' If a download fails or you used the wrong organism or assembly, delete the
#' incorrect files and run again.\cr
#' Set \code{remake = TRUE} to rebuild everything from scratch.\cr
#'
#' Files created during the download process include:\cr
#' - A fasta index for the genome\cr
#' - A TxDb to speed up GTF/GFF reading\cr
#' - Separate or merged contaminant files\cr
#' Optional files that can also be created:\cr
#' - Gene symbols (hgnc, etc)\cr
#' - Uniprot ids (For name of protein structures)\cr
#' To use references already on disk, pass their paths as \code{genome} and
#' \code{GTF}. Supply \code{organism} and \code{output.dir} as well.
#' A completed run reuses the paths in \code{outputs.rds}; choose a separate
#' output directory or set \code{remake = TRUE} when changing references.
#' @inheritParams biomartr::getGenome
#' @param organism scientific name of organism, Homo sapiens,
#' Danio rerio, Mus musculus, etc. See \code{biomartr:::get.ensembl.info()}
#' for full list of supported organisms.
#' @param output.dir directory to save downloaded data
#' @param db database to use for the genome and annotation. The recommended
#' choice is `"ensembl"`. When using Ensembl, remember to set
#' \code{assembly_type = "primary_assembly"} unless you want haplotypes
#' included, which can make the files much larger. Alternatives are `"refseq"`
#' for reference assemblies and `"genbank"` for all available assemblies.
#' @param GTF logical or character, default TRUE. TRUE downloads the annotation
#' for \code{organism}; FALSE searches \code{output.dir} for an existing file.
#' A character path uses that local GTF/GFF file and creates its TxDb.
#' The argument name is case-sensitive: use \code{GTF}, not \code{gtf}.
#' For RefSeq and GenBank, see \code{refseq_genbank_format}.
#' @param genome logical or character, default TRUE. TRUE downloads the genome;
#' FALSE searches \code{output.dir} for an existing file. A character path uses
#' that local FASTA and creates its index. For Ensembl, \code{assembly_type}
#' selects the sequence set. FASTA and annotation chromosome names must match.
#' @param merge_contaminants logical, default TRUE. Merge the requested
#' contaminant references into one fasta file. This usually saves space and is
#' faster to align with STAR than keeping each contaminant as a separate
#' reference. Ignored if no contaminants are requested.
#' @param phix logical, default FALSE, download phiX sequence to filter
#'  out Illumina control reads. ORFik defines Phix as a contaminant genome.
#' Phix is used in Illumina sequencers for sequencing quality control.
#' Genome is: refseq, Escherichia phage phiX174.
#' If sequencing facility created fastq files with the command \code{bcl2fastq},
#' then there should be very few phix reads left in the fastq files received.
#' @param ncRNA logical or character, default FALSE. If `TRUE` or a path is
#' supplied, ncRNA is used as a contaminant reference. If `TRUE`, ORFik tries
#' to find ncRNA sequences from the GTF file, usually represented as lncRNA
#' entries, and reports if none are found.\cr If that does not work, try a
#' character input.\cr A value of `"auto"` tells ORFik to look for an ncRNA
#' file for the organism on NONCODE, using the organism's common name
#' (for example, Homo sapiens -> human). `"auto"` does not work for every
#' species, so you may need to provide the common name used by NONCODE
#' directly. If the value is not `"auto"` or `""`, it must be a character
#' vector of common species names, not scientific names. See
#' http://www.noncode.org/download.php/ if you need to look up the correct
#' NONCODE name.
#' @param tRNA logical or character, default FALSE. If enabled, tRNA is used
#' as a contaminant reference. If `TRUE`, ORFik tries to extract tRNA
#' sequences from the GTF file, usually represented as `Mt_tRNA` entries, and
#' reports if none are found. If that is not sufficient, provide a character
#' vector with a valid path to a fasta file containing mature tRNAs on disk.
#' You can obtain these from http://gtrnadb.ucsc.edu/ or by running tRNAscan
#' on the genome.
#' @param rRNA logical or character, default FALSE. If enabled, rRNA is used
#' as a contaminant reference. If `TRUE`, ORFik tries to extract rRNA
#' sequences from the GTF file and reports if none are found. If set to
#' `"silva"`, ORFik downloads the Silva SSU and LSU reference files for all
#' species. If you need a smaller reference, download one manually from
#' https://www.arb-silva.de/.\cr
#' If the value is neither `""` nor `"silva"`, it must be a character vector
#' giving a valid path to an rRNA fasta file on disk.
#' @param gunzip logical, default TRUE, uncompress downloaded files
#' that are zipped when downloaded, should be TRUE!
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @param remove_annotation_outliers logical, default TRUE. Only used for
#' RefSeq annotations. If `TRUE`, malformed outlier lines are removed from the
#' input annotation file. The cleaned annotation overwrites the original file,
#' and the removed lines are stored in `tempdir()` for inspection. This is
#' needed for some RefSeq annotations, including Arabidopsis.
#' @param notify_load_existing logical, default TRUE. If a previous download is
#' already present in `output.dir` and recorded in `outputs.rds`, print a short
#' message instead of silently reusing it. Set to `FALSE` to suppress this.
#' @param assembly character, default `assembly = organism`, which selects the
#' first assembly returned for that organism. You can also supply a specific
#' assembly name, for example `"GCA_000005845"` for the commonly used
#' E. coli K-12 reference. This argument is usually only needed for bacterial
#' genomes.
#' @param refseq_genbank_format = c("gtf", "gff3")[1] Annotation format to
#' request from RefSeq or GenBank. GTF is usually the safer choice for
#' downstream compatibility. GFF3 may contain information that is not present
#' in the GTF, so use it when you specifically need that extra annotation.
#' @inheritParams makeTxdbFromGenome
#' @importFrom biomartr getGTF getGenome getENSEMBLInfo
#' @importFrom Rsamtools indexFa
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom AnnotationDbi saveDb
#' @importFrom Biostrings DNAStringSet writeXStringSet readDNAStringSet
#' @return A named character vector with the paths to the downloaded genome,
#' annotation, and any contaminant references. If `merge_contaminants = TRUE`,
#' only the merged contaminant fasta is returned, not the individual files.
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
