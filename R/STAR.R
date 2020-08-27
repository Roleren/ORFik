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
#' @inheritParams base::system
#' @return output.dir, can be used as as input for STAR.align..
#' @family STAR
#' @export
#' @examples
#' # In this argument specify the one you want
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
                       wait = TRUE,
                       script = system.file("STAR_Aligner",
                                            "STAR_MAKE_INDEX.sh",
                                            package = "ORFik")) {
  if (!file.exists(script)) stop("STAR index script not found, check path of script!")
  if (is.null(names(arguments))) stop("arguments must have names, see ?STAR.index")
  possible <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA")
  if (!all(names(arguments) %in% possible)) stop("At least one of arguments with invalid name!")

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
    if (!wait) out <- "Wait for index to be complete before you run Alignment!"
    message(out)
  } else stop("STAR is not supported on windows!")
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
#' @param input.dir path to fast files to align
#' @param index.dir path to STAR index
#' @param fastp path to fastp trimmer, default: install.fastp(), if you
#' have it somewhere else already installed, give the path. If you are not on linux
#' and you want to trim, use your favorite trimmer and give the output files from that
#' trimmer as input.dir here.
#' @param paired.end default "no", alternative "yes". Will auto detect
#'  pairs by names.
#' @param steps a character, default: "tr-ge", trimming --> genome alignment\cr
#'  steps of depletion and alignment wanted:
#'  If not "all", a subset of these ("tr-ph-rR-nc-tR-ge")\cr
#'  In bash script it it reformated to this style:
#'  (trimming and genome do: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge")
#'  tr: trim, ph: phix, rR: rrna, nc: ncrna, tR: trna, ge: genome)\cr
#'  genome: the step where you align to the genome, this is usually always included.
#' @param adapter.sequence "auto"
#' @param min.length 15, minimum length of reads to pass filter.
#' @param trim.front 0, default trim 0 bases 5'. For Ribo-seq set use 0.
#' Ignored if tr (trim) is not one of the arguments in "steps"
#' @param alignment.type default: "Local": standard local alignment with soft-clipping allowed,
#' "EndToEnd" (global): force end-to-end read alignment, does not soft-clip.
#' @param include.subfolders "n" (no), do recursive search downwards for fast files if "y".
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
                              include.subfolders = "n",
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
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))

  full <- paste(script.folder, "-f", input.dir, "-o", output.dir,
                "-p", paired.end,
                "-l", min.length, "-g", index.dir, "-s", steps,
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
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.
#' @inheritParams STAR.align.folder
#' @param file1 library file, if paired must be R1 file
#' @param file2 default NULL, set if paired end to R2 file
#' @param resume default: NULL, continue from step, lets say steps are "tr-ph-ge":
#'  (trim, phix depletion, genome alignment) and resume is "ph", you will use the trimmed
#'  data and continue from there starting at phix, usefull if something crashed.
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
#' Will create a R transcript database (TxDb object) from the annotation. \cr
#' It will also index the genome \cr
#' If you misspelled something or crashed, delete wrong files and
#' run again.\cr
#' Do remake = TRUE, to do it all over again.
#'
#' If you want custom genome or gtf from you hard drive, assign it back
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
#' @param phix logical, default FALSE, download phix sequence to filter
#'  out with.
#' Only use if illumina sequencing. Phix is used in Illumina sequencers for
#' sequencing quality control. Genome is: refseq, Escherichia virus phiX174
#' @param ncRNA character, default "" (no download), a contaminant genome.
#' Alternatives: "auto" or manual assign like "human".
#' If "auto" will try to find ncRNA file from organism, Homo sapiens -> human
#' etc. "auto" will not work for all, then you must specify the name used by
#' NONCODE, go to the link below and find it.
#' If not "auto" / "" it must be a character vector
#' of species common name (not scientific name) Homo sapiens is human,
#' Rattus norwegicus is rat etc, download ncRNA sequence to filter out with.
#' From NONCODE online server, if you cant find
#' common name see: http://www.noncode.org/download.php/
#' @param tRNA chatacter, default "" (not used),
#' if not "" it must be a character vector  to valid path of mature
#' tRNAs fasta file to remove as contaminants on your disc. Find and download
#' your wanted mtRNA at: http://gtrnadb.ucsc.edu/, or run trna-scan on
#' you genome.
#' @param rRNA chatacter, default "" (not used),
#' if not "" it must be a character vector to valid path of mature
#' rRNA fasta file to remove as contaminants on your disc. Find and download
#' your wanted rRNA at: https://www.arb-silva.de/
#' @param gunzip logical, default TRUE, uncompress downloaded files
#' that are zipped when downloaded, should be TRUE!
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @param assembly_type a character string specifying from which assembly type the genome
#' shall be retrieved from (ensembl only, else this argument is ignored):
#' Default is
#' \code{assembly_type = "toplevel")}.
#' This will give you all multi-chromosomes (copies of the same chromosome with small variations).
#' As an example the toplevel fasta genome in human is over 70 GB uncompressed.
#' To get primary assembly with 1 chromosome variant per chromosome:
#' \code{assembly_type = "primary_assembly")}.
#' As an example, the  primary_assembly fasta genome in human is only a few GB uncompressed:
#' @importFrom biomartr getGTF getGenome getENSEMBLInfo
#' @importFrom Rsamtools indexFa
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom AnnotationDbi saveDb
#' @return a character vector of path to genomes and gtf downloaded,
#'  and additional contaminants if used.
#' @family STAR
#' @export
#' @examples
#' output.dir <- "/Bio_data/references/zebrafish"
#' #getGenomeAndAnnotation("Danio rerio", output.dir)
getGenomeAndAnnotation <- function(organism, output.dir, db = "ensembl",
                            GTF = TRUE, genome = TRUE, phix = FALSE,
                            ncRNA = "", tRNA = "", rRNA = "",
                            gunzip = TRUE, remake = FALSE,
                            assembly_type = "primary_assembly") {
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) & !remake) {
    message("Loading premade files information,
            do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }
  if (!(assembly_type %in% c("toplevel", "primary_assembly")))
    stop("Please select one the available assembly types: \ntoplevel, primary_assembly")
  dir.create(output.dir, recursive = TRUE)
  if (tRNA != "") {
    if (!file.exists(tRNA)) stop(paste("tRNA file is given and does not exist:",
                                       tRNA))
  }
  if (rRNA != "") {
    if (!file.exists(rRNA)) stop(paste("rRNA file is given and does not exist:",
                                       tRNA))
  }

  if (ncRNA != "") {
    if (ncRNA == "auto") {
      a <- biomartr::getENSEMBLInfo()
      ncRNA <- a[grep(gsub(" ", "_", organism),
                   a$name, TRUE),]$common_name
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

  if (genome != FALSE) { # fasta genome of organism
    if (db == "ensembl") {
      message(paste("Starting primary assembly genome retrieval of",
                    organism, "from ensembl: "))
      genome <- biomartr:::getENSEMBL.Seq(organism, type = "dna",
                                          release = NULL,
                                          id.type = assembly_type,
                                          path = output.dir)[1]
      if (gunzip) # unzip gtf file
        genome <- R.utils::gunzip(genome, overwrite = TRUE)
    } else {
      genome  <- biomartr::getGenome(db = db, organism,
                                     path = output.dir, gunzip = gunzip)
    }

    message("Making .fai index of genome")
    indexFa(genome)
  } else { # check if it already exists
    genome <- grep(pattern = gsub(" ", "_", organism),
                  x = list.files(output.dir, full.names = TRUE),
                  value = TRUE)
    genome <- grep(pattern = "\\.fa", x = genome, value = TRUE)
    genome <- grep(pattern = "\\.fai", x = genome, value = TRUE, invert = TRUE)
    if (length(genome) != 1) genome <- FALSE
  }
  if (GTF) { # gtf of organism
    # gtf <- biomartr::getGTF(db = db, organism,
    #                         path = output.dir,
    #                         assembly_type = assembly_type)
    gtf <- biomartr:::getENSEMBL.gtf(organism, type = "dna",
                                     id.type = assembly_type,
                                     output.dir)

    if (gunzip) # unzip gtf file
      gtf <- R.utils::gunzip(gtf, overwrite = TRUE)
    message("Making txdb of GTF")
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, organism = organism)
    if (genome != FALSE)
      seqlevelsStyle(txdb) <- seqlevelsStyle(FaFile(genome))[1]
    txdb_file <- paste0(gtf, ".db")
    AnnotationDbi::saveDb(txdb, txdb_file)
  } else { # check if it already exists
    gtf <- grep(pattern = gsub(" ", "_", organism),
                x = list.files(output.dir, full.names = TRUE),
                value = TRUE)
    gtf <- grep(pattern = "\\.gtf", x = gtf, value = TRUE)
    gtf <- grep(pattern = "\\.db", x = gtf, value = TRUE, invert = TRUE)
    if (length(gtf) != 1) gtf <- FALSE
  }
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

  message("All data downloaded and ready at:")
  message(output.dir)
  output <- c(gtf, genome, phix, ncRNA, tRNA, rRNA)
  names(output) <- c("gtf", "genome", "phix", "ncRNA", "tRNA", "rRNA")
  output <- output[!(output %in% c("", "FALSE"))]

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
#' Will not run "make", only use precompiled fastp file.\cr
#' Only works for Linux
#' @param folder path to folder for download, fille will be named
#' "fastp", this should be most recent version
#' @importFrom utils download.file
#' @return path to runnable fastp
#' @export
#' @references https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/
#' @family STAR
#' @examples
#' #install.fastp()
install.fastp <- function(folder = "~/bin") {
  if (.Platform$OS.type != "unix")
    stop("fastp does not work on Windows, try RSubread")
  url <- "http://opengene.org/fastp/fastp"
  path <- paste0(folder, "/fastp")
  if (file.exists(path)) {
    message(paste("Using fastp at location:",
                  path))
    return(path)
  }
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  utils::download.file(url, destfile = path)
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
