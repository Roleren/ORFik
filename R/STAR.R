#' Create STAR genome index
#'
#' Used as reference when aligning data \cr
#' Get genome and gtf by running getGenomeAndFasta()
#'
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.
#' @param arguments a named character vector containing paths wanted to
#' use for index creation. They must be named correctly:
#' names must be a subset of:
#' c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA")
#' @param output.dir directory to save indices, default:
#' paste0(dirname(arguments[1]), "/STAR_index/"), where arguments is the
#' arguments input for this function.
#' @param star.path path to STAR, default: STAR.install(),
#' if you don't have STAR installed at default location, it will install it there,
#' set path to a runnable star if you already have it.
#' @param max.cpus integer, default: min(90, detectCores() - 1),
#'  number of threads to use. Default is minimum of 90 and maximum cores - 1
#' @param script location of STAR index script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @inheritParams base::system
#' @return output.dir, can be used as as input for STAR.align..
#' @family STAR
#' @export
#' @examples
#' ## Manual way, specify all paths yourself.
#' #arguments <- c(path.GTF, path.genome, path.phix, path.rrna, path.trna, path.ncrna)
#' #names(arguments) <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA")
#' #STAR.index(arguments, "output.dir")
#'
#' ## Or use ORFik way:
#' output.dir <- "/Bio_data/references/Human"
#' # arguments <- getGenomeAndAnnotation("Homo sapiens", output.dir)
#' # STAR.index(arguments, output.dir)
STAR.index <- function(arguments, output.dir = paste0(dirname(arguments[1]), "/STAR_index/"),
                       star.path = STAR.install(), max.cpus = min(90, detectCores() - 1),
                       wait = TRUE, remake = FALSE,
                       script = system.file("STAR_Aligner",
                                            "STAR_MAKE_INDEX.sh",
                                            package = "ORFik")) {
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) & !remake) {
    message("Loading premade files information,
            do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }

  if (!file.exists(script))
    stop("STAR index script not found, check path of script!")
  if (is.null(names(arguments)))
    stop("arguments must have names, see ?STAR.index")
  possible <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA")
  if (!all(names(arguments) %in% possible))
    stop("At least one of arguments with invalid name!")

  # match which indices to make
  exts <- c("g", "f", "p", "r", "t", "n")
  names(exts) <- possible
  hits <- paste0("-", exts[names(arguments)], " ", arguments)

  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  out <- paste("-o", output.dir)
  max.cpus <- paste("-m", max.cpus)
  full <- paste(script, out, star.path, max.cpus,
                paste(hits, collapse = " "))
  message("STAR indexing:\n")
  print(full); print("\n")
  if (.Platform$OS.type == "unix") {

    message("Starting indexing at time:")
    print(Sys.time())
    out <- system(command = full, wait = wait)
    out <- ifelse(out == 0 & wait, "Index done", "Index process failed!")
    if (!wait)
      out <- "Wait for index to be complete before you run Alignment!"
    message(out)
  } else stop("STAR is not supported on windows!")
  saveRDS(object = output.dir, finished.file)
  return(output.dir)
}

#' Align all libraries in folder with STAR
#'
#' Does either all files as paired end or single end,
#' so if you have mix, split them in two different folders.\cr
#' #' If STAR halts at .... loading genome, it means the STAR
#' index was aborted early, then you need to run:
#' STAR.remove.crashed.genome(), with the genome that crashed, and rerun.
#'
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.
#' The trimmer used is fastp (the fastest I could find), works on mac and linux.
#' If you want to use your own trimmer set file1/file2 to the location of
#' the trimmed files from your program.\cr
#' A note on trimming from creator of STAR about trimming:
#' "adapter trimming it definitely needed for short RNA sequencing.
#' For long RNA-seq, I would agree with Devon that in most cases adapter trimming
#' is not advantageous, since, by default, STAR performs local (not end-to-end) alignment,
#' i.e. it auto-trims." So trimming can be skipped for longer reads.
#' @param input.dir path to fast files to align, can either be
#' fasta files (.fastq, .fq, .fa etc) or compressed files with .gz.
#' Also either paired end or single end reads.
#' @param index.dir path to STAR index folder. Path returned from ORFik function
#' STAR.index, when you created the index folders.
#' @param fastp path to fastp trimmer, default: install.fastp(), if you
#' have it somewhere else already installed, give the path. Only works for
#' unix (linux or Mac OS), if not on unix, use your favorite trimmer and
#' give the output files from that trimmer as input.dir here.
#' @param paired.end default "no", alternative "yes". Will auto detect
#'  pairs by names. If yes running on a folder:
#'  The folder must then contain an even number of files
#'  and they must be named with the same prefix and sufix of either
#'   _1 and _2, 1 and 2, etc.
#' @param steps a character, default: "tr-ge", trimming then genome alignment\cr
#'  steps of depletion and alignment wanted:
#'  The posible candidates you can use are:
#'  tr: trim reads, ph: phix depletion, rR: rrna depletion,
#'  nc: ncrna depletion, tR: trna depletion, ge: genome alignment,
#'  all: run all steps)\cr
#'  If not "all", a subset of these ("tr-ph-rR-nc-tR-ge")\cr
#'  In bash script it is reformated to this style:
#'  (trimming and genome do: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge")
#'  the step where you align to the genome is usually always included, unless you
#'  are doing pure contaminant analysis.
#'  For Ribo-seq and TCP(RCP-seq) you should do rR (ribosomal RNA depletion),
#'  so when you made the
#'  STAR index you need the rRNA step (usually just download a Silva rRNA database
#'  for SSU&LSU at: https://www.arb-silva.de/)
#' @param adapter.sequence character, default: "auto" (auto detect adapter, is not
#' very reliable for Ribo-seq, so then you must include,
#' else alignment will most likely fail!). Else manual assigned adapter like:
#' "ATCTCGTATGCCGTCTTCTGCTTG" or "AAAAAAAAAAAAA".
#' @param min.length 15, minimum length of reads to pass filter.
#' @param trim.front 0, default trim 0 bases 5'. For Ribo-seq set use 0.
#' Ignored if tr (trim) is not one of the arguments in "steps"
#' @param alignment.type default: "Local": standard local alignment with soft-clipping allowed,
#' "EndToEnd" (global): force end-to-end read alignment, does not soft-clip.
#' @param include.subfolders "n" (no), do recursive search downwards for fast files if "y".
#' @param resume default: NULL, continue from step, lets say steps are "tr-ph-ge":
#'  (trim, phix depletion, genome alignment) and resume is "ph", you will then use
#'  the assumed already trimmed data and start / continue from there starting at phix,
#'  usefull if something crashed. Like if you specified wrong STAR version, but the trimming
#'  step was completed.
#' @param script.folder location of STAR index script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @param script.single location of STAR single file alignment script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @inheritParams STAR.index
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
#' @examples
#' # Use your own paths for annotation or the ORFik way
#'
#' ## use ORFik way:
#' output.dir <- "/Bio_data/references/Human"
#' # arguments <- getGenomeAndAnnotation("Homo sapiens", output.dir)
#' # index <- STAR.index(arguments, output.dir)
#' # STAR.align.folder("data/raw_data/human_rna_seq", "data/processed/human_rna_seq",
#' #                    index, paired.end = "no")
STAR.align.folder <- function(input.dir, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              paired.end = "no",
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 15, trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE,
                              include.subfolders = "n", resume = NULL,
                              script.folder = system.file("STAR_Aligner",
                                                          "RNA_Align_pipeline_folder.sh",
                                                          package = "ORFik"),
                              script.single = system.file("STAR_Aligner",
                                                          "RNA_Align_pipeline.sh",
                                                          package = "ORFik")) {
  if (!file.exists(script.folder))
    stop("STAR folder alignment script not found, check path of script!")
  if (!file.exists(script.single))
    stop("STAR single file alignment script not found, check path of script!")
  cleaning <- system.file("STAR_Aligner", "cleanup_folders.sh",
                         package = "ORFik", mustWork = TRUE)
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))

  full <- paste(script.folder, "-f", input.dir, "-o", output.dir,
                "-p", paired.end,
                "-l", min.length, "-g", index.dir, "-s", steps, resume,
                "-a", adapter.sequence, "-t", trim.front,
                "-A", alignment.type, "-m", max.cpus, "-i", include.subfolders,
                star.path, fastp, "-I",script.single, "-C", cleaning)
  if (.Platform$OS.type == "unix") {
    print(full)
    message("Starting alignment at time:")
    print(Sys.time())
    out <- system(command = full, wait = wait)
    out <- ifelse(out == 0, "Alignment done", "Alignment process failed!")
    message(out)
  } else stop("STAR is not supported on windows!")
  return(output.dir)
}

#' Align single or paired end pair with STAR
#'
#' If you want more than two files use: STAR.align.folder\cr
#' If genome aligner halts at .... loading genome, it means the star
#' index was aborted early, then you need to run:
#' STAR.remove.crashed.genome(), with the genome that crashed, and rerun.
#'
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.\cr
#' The trimmer used is fastp (the fastest I could find), works on mac and linux.
#' If you want to use your own trimmer set file1/file2 to the location of
#' the trimmed files from your program.
#' @inheritParams STAR.align.folder
#' @param file1 library file, if paired must be R1 file
#' @param file2 default NULL, set if paired end to R2 file
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
#' @examples
#' # Use your own paths for annotation or the ORFik way
#'
#' ## use ORFik way:
#' output.dir <- "/Bio_data/references/Human"
#' # arguments <- getGenomeAndAnnotation("Homo sapiens", output.dir)
#' # index <- STAR.index(arguments, output.dir)
#' # STAR.align.single("data/raw_data/human_rna_seq/file1.bam", "data/processed/human_rna_seq",
#' #                    index)
STAR.align.single <- function(file1, file2 = NULL, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 15, trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE, resume = NULL,
                              script.single = system.file("STAR_Aligner",
                                                   "RNA_Align_pipeline.sh",
                                                   package = "ORFik")
) {
  if (!file.exists(script.single))
    stop("STAR single file alignment script not found, check path of script!")
  cleaning <- system.file("STAR_Aligner", "cleanup_folders.sh",
                          package = "ORFik", mustWork = TRUE)

  file2 <- ifelse(is.null(file2), "", paste("-F", file2))
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))

  full <- paste(script.single, "-f", file1, file2, "-o", output.dir,
                "-l", min.length, "-g", index.dir, "-s", steps,
                resume, "-a", adapter.sequence, "-t", trim.front,
                "-A", alignment.type, "-m", max.cpus, star.path, fastp,
                "-C", cleaning)
  if (.Platform$OS.type == "unix") {
    print(full)
    message("Starting alignment at time:")
    print(Sys.time())
    out <- system(command = full, wait = wait)
    out <- ifelse(out == 0, "Alignment done", "Alignment process failed!")
    message(out)
  } else stop("STAR is not supported on windows!")
  return(output.dir)
}

# STAR.merge.tcp <- function(input.dir) {
#   dir.create(paste0(input.dir, "/merged"))
#   system("samtools --help")
# }

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
#' Danio rerio, Mus musculus, etc.
#' @param output.dir directory to save downloaded data
#' @param db database to use for genome and GTF,
#' default adviced: "ensembl" (will contain haplotypes, large file!).
#' Alternatives: "refseq" (primary assembly) and "genbank" (mix)
#' @param GTF logical, default: TRUE, download gtf of organism specified
#' in "organism" argument. If FALSE, check if the downloaded
#' file already exist. If you want to use a custom gtf from you hard drive,
#' set GTF = FALSE,
#' and assign: \cr annotation <- getGenomeAndAnnotation(gtf = FALSE)\cr
#' annotation["gtf"] = "path/to/gtf.gtf".\cr
#' Only db = "ensembl" allowed for GTF.
#' @param genome logical, default: TRUE, download genome of organism
#' specified in "organism" argument. If FALSE, check if the downloaded
#' file already exist. If you want to use a custom gtf from you hard drive,
#' set GTF = FALSE,
#' and assign: \cr annotation <- getGenomeAndAnnotation(genome = FALSE)\cr
#' annotation["genome"] = "path/to/genome.fasta".\cr
#' Will download the primary assembly for ensembl
#' @param merge_contaminants logical, default TRUE. Will merge
#' the contaminants specified into one fasta file, this considerably
#' saves space and is much quicker to align with STAR than each contamint
#' on it's own. If no contaminants are specified, this is ignored.
#' @param phix logical, default FALSE, download phix sequence to filter
#'  out with. Phix is used as a contaminant genome.
#' Only use if illumina sequencing. Phix is used in Illumina sequencers for
#' sequencing quality control. Genome is: refseq, Escherichia virus phiX174
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
#' gtf. If not found try character input:\cr
#' if not "" it must be a character vector to valid path of mature
#' rRNA fasta file to remove as contaminants on your disc. Find and download
#' your wanted rRNA at: https://www.arb-silva.de/
#' @param gunzip logical, default TRUE, uncompress downloaded files
#' that are zipped when downloaded, should be TRUE!
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @param assembly_type a character string specifying from which assembly type
#' the genome shall be retrieved from (ensembl only, else this argument is ignored):
#' Default is
#' \code{assembly_type = "primary_assembly")}.
#' This will give you all no copies of any chromosomes.
#' As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed.\cr
#' \code{assembly_type = "toplevel")}.
#' This will give you all multi-chromosomes (copies of the same chromosome with small variations).
#' As an example the toplevel fasta genome in human is over 70 GB uncompressed.
#' To get primary assembly with 1 chromosome variant per chromosome:
#' @importFrom biomartr getGTF getGenome getENSEMBLInfo
#' @importFrom Rsamtools indexFa
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom AnnotationDbi saveDb
#' @return a named character vector of path to genomes and gtf downloaded,
#'  and additional contaminants if used. If merge_contaminants is TRUE, will not
#'  give individual fasta files to contaminants, but only the merged one.
#' @family STAR
#' @export
#' @examples
#' output.dir <- "/Bio_data/references/zebrafish"
#' #getGenomeAndAnnotation("Danio rerio", output.dir)
#'
#' ## Get Phix contamints to deplete during alignment
#' #getGenomeAndAnnotation("Danio rerio", output.dir, phix = TRUE)
#'
getGenomeAndAnnotation <- function(organism, output.dir, db = "ensembl",
                            GTF = TRUE, genome = TRUE,
                            merge_contaminants = TRUE, phix = FALSE,
                            ncRNA = FALSE, tRNA = FALSE, rRNA = FALSE,
                            gunzip = TRUE, remake = FALSE,
                            assembly_type = "primary_assembly") {
  # Pre checks
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) & !remake) {
    message("Loading premade Genome files,
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
  ## Go through all contaminants:
  if (!(tRNA %in% c("", FALSE, TRUE))) {
    if (!file.exists(tRNA)) stop(paste("tRNA file is given and does not exist:",
                                       tRNA))
  }
  if (!(rRNA %in% c("", FALSE, TRUE))) {
    if (!file.exists(rRNA)) stop(paste("rRNA file is given and does not exist:",
                                       tRNA))
  }
  phix <- get_phix_genome(phix, output.dir, gunzip)
  ncRNA <- get_noncoding_rna(ncRNA, output.dir, organism, gunzip)

  # Get species fasta genome and gtf
  genome <- get_genome_fasta(genome, output.dir, organism,
                             assembly_type, db, gunzip)
  gtf <- get_genome_gtf(GTF, output.dir, organism, assembly_type, gunzip)

  if (any_contaminants) {
    # Find which contaminants to find from gtf:
    conts <- as.character(c(phix, ncRNA, tRNA, rRNA))
    names(conts) <- c("phix", "ncRNA", "tRNA", "rRNA")
    non_gtf_contaminants <- conts[!(conts %in% c(TRUE, FALSE, ""))]
    gtf_contaminants <- conts[conts == TRUE]
    if (length(gtf_contaminants) > 0) {
      if (is.logical(gtf) | is.logical(genome))
        stop("gtf or genome not specified, so impossible to find gtf contaminants!")
      # Make fasta file of those contaminants
      total_seqs <- DNAStringSet()
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

#' Download and prepare STAR
#'
#' Will not run "make", only use precompiled STAR file.\cr
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.
#' @param folder path to folder for download, fille will be named
#' "STAR-version", where version is version wanted.
#' @param version default "2.7.4a"
#' @importFrom utils download.file untar tar
#' @return path to runnable STAR
#' @export
#' @references https://www.ncbi.nlm.nih.gov/pubmed/23104886
#' @family STAR
#' @examples
#' #STAR.install("~/bin", version = "2.7.4a")
STAR.install <- function(folder = "~/bin", version = "2.7.4a") {
  if (.Platform$OS.type != "unix")
    stop("STAR does not work on Windows, try RSubread")
  url <- paste0("https://github.com/alexdobin/STAR/archive/",
                version,".tar.gz")
  folder <- path.expand(folder)
  path <- paste0(folder, "/STAR-", version)
  pathgz <- paste0(path, ".tar.gz")
  bin <- ifelse(Sys.info()[1] == "Linux",
                paste0(path, "/bin/Linux_x86_64/STAR"),
                paste0(path, "/bin/MacOSX_x86_64/STAR"))
  names(bin) <- "STAR"
  if (file.exists(bin)) {
    message(paste("Using STAR at location:",
                  bin))
    return(bin)
  }
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  utils::download.file(url, destfile = pathgz)
  untar(pathgz, exdir = folder)
  return(bin)
}

#' Download and prepare fastp trimmer
#'
#' On Linux, will not run "make", only use precompiled fastp file.\cr
#' On Mac OS it will use precompiled binaries.\cr
#' Does not work yet for Windows!
#' @param folder path to folder for download, file will be named
#' "fastp", this should be most recent version. On mac it will search
#' for a folder called fastp-master inside folder given. Since there
#' is no precompiled version of fastp for Mac OS.
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @return path to runnable fastp
#' @export
#' @references https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/
#' @family STAR
#' @examples
#' #install.fastp()
install.fastp <- function(folder = "~/bin") {
  if (.Platform$OS.type != "unix")
    stop("fastp does not work on Windows, try RSubread")

  is_linux <- Sys.info()[1] == "Linux" # else it is mac
  url <- ifelse(is_linux, # else it is mac
                "http://opengene.org/fastp/fastp",
                "https://github.com/OpenGene/fastp/archive/master.zip")
  path <- ifelse(is_linux, # else it is mac
                 paste0(folder, "/fastp"),
                 paste0(folder, "/fastp-master/fastp"))

  if (file.exists(path)) {
    message(paste("Using fastp at location:",
                  path))
    return(path)
  }
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  if (!is_linux) path <- paste0(folder, "/fastp.zip")
  utils::download.file(url, destfile = path)
  if (!is_linux) { # For mac os
    message("On mac OS, must build fastp, since no precompiled binaries exists")
    message("This will only be done once")
    utils::unzip(path, exdir = folder)
    system(paste0("make -C ", folder, "/fastp-master/"))
    path <- paste0(folder, "/fastp-master/fastp")
  }
  # Update access rights
  system(paste("chmod a+x", path))
  return(path)
}

#' Remove crashed STAR genome
#'
#' This happens if you abort STAR run early, and it halts at: ..... loading genome
#' @param index.path path to index folder of genome
#' @inheritParams STAR.index
#' @return return value from system, 0 if all good.
#' @export
#' @family STAR
#' @examples
#' # STAR.remove.crashed.genome(index.path = "/home/data/human_index/phix/)
STAR.remove.crashed.genome <- function(index.path, star.path = STAR.install()) {
  message("Trying to remove loaded genome:")
  out <- paste(star.path, "--genomeDir", index.path, "--genomeLoad Remove")
  system(out)
}
