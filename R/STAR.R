#' Create STAR genome index
#'
#' Used as reference when aligning data \cr
#' Get genome and gtf by running getGenomeAndFasta()
#'
#' Can only run on unix systems (Linux and Mac)
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
#' @param script location of STAR script
#' @inheritParams base::system
#' @return output.dir, can be used as as input for STAR.align..
#' @family STAR
#' @export
#' @examples
#' # In this argument specify the one you want
#' #arguments <- c(path.GTF, path.genome, path.phix, path.rrna, path.trna, path.ncrna)
#' #names(arguments) <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA")
#' #STAR.index(arguments, "output.dir")
STAR.index <- function(arguments, output.dir = paste0(dirname(arguments[1]), "/STAR_index/"),
                       star.path = STAR.install(), max.cpus = min(90, detectCores() - 1),
                       wait = TRUE,
                       script = "/export/valenfs/projects/Pipelines/STAR_Aligner/STAR_MAKE_INDEX.sh") {
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
#' Can only run on unix systems (Linux and Mac)
#' @param input.dir path to fast files to align
#' @param index.dir path to STAR index
#' @param fastp path to fastp trimmer, default: install.fastp(), if you
#' have it somewhere else already installed, give the path.
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
#' @inheritParams STAR.index
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
STAR.align.folder <- function(input.dir, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              paired.end = "no",
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 15, trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE,
                              include.subfolders = "n",
                              script = "/export/valenfs/projects/Pipelines/STAR_Aligner/RNA_Align_pipeline_folder.sh") {
  if (!file.exists(script)) stop("STAR alignment script not found, check path of script!")
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))

  full <- paste(script, "-f", input.dir, "-o", output.dir, "-p", paired.end,
                "-l", min.length, "-g", index.dir, "-s", steps,
                "-a", adapter.sequence, "-t", trim.front,
                "-A", alignment.type, "-m", max.cpus, "-i", include.subfolders,
                star.path)
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
#' @inheritParams STAR.align.folder
#' @param file1 library file, if paired must be R1 file
#' @param file2 default NULL, set if paired end to R2 file
#' @param resume default: NULL, continue from step, lets say steps are "tr-ph-ge":
#'  (trim, phix depletion, genome alignment) and resume is "ph", you will use the trimmed
#'  data and continue from there starting at phix, usefull if something crashed.
#' @family STAR
#' @export
STAR.align.single <- function(file1, file2 = NULL, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 15, trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE, resume = NULL,
                              script = "/export/valenfs/projects/Pipelines/STAR_Aligner/RNA_Align_pipeline.sh"
) {
  if (!file.exists(script)) stop("STAR alignment script not found, check path of script!")
  file2 <- ifelse(is.null(file2), "", paste("-F", file2))
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))

  full <- paste(script, "-f", file1, file2, "-o", output.dir,
                "-l", min.length, "-g", index.dir, "-s", steps,
                resume, "-a", adapter.sequence, "-t", trim.front,
                "-A", alignment.type, "-m", max.cpus, star.path)
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

#' Download genome (fasta), GTF and contaminants
#'
#' Will create a R transcript database (TxDb object) from the genome. \cr
#' It will also index the genome \cr
#' If you misspelled something or crashed, delete wrong files and run again.\cr
#' Do remake = TRUE, to do it all over again.
#' @param organism scientific name of organism
#' @param output.folder folder to save downloaded data
#' @param db database to use, default adviced: "ensembl"
#' @param GTF logical, default: TRUE, download gtf of organism
#' @param genome logical, default: TRUE, download genome of organism
#' @param fasta path to genome fast file
#' @param phix logical, default FALSE, download phix sequence to filter out with.
#' Only use if illumina sequencing. Phix is used in Illumina sequencers for
#' sequencing quality control. Genome is: refseq, Escherichia virus phiX174
#' @param ncRNA chatacter, default "" (no download), a contaminant genome,
#' if not "" it must be a character vector
#' of species common name (not scientific) Homo sapiens is human, Rattus norwegicus is rat etc,
#' download ncRNA sequence to filter out with. From NONCODE online server, if you cant find
#' common name see: http://www.noncode.org/download.php/
#' @param tRNA chatacter, default "" (not used),
#' if not "" it must be a character vector  to valid path of mature
#' tRNAs fasta file to remove as contaminants. Find and download your wanted mtRNA at:
#' http://gtrnadb.ucsc.edu/, or run trna-scan on you genome.
#' @param rRNA chatacter, default "" (not used),
#' if not "" it must be a character vector to valid path of mature
#' rRNA fasta file to remove as contaminants. Find and download your wanted rRNA at:
#' https://www.arb-silva.de/
#' @param gunzip logical, default TRUE, uncompress downloaded files
#' that are zipped as default
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @importFrom biomartr getGTF
#' @importFrom biomartr getGenome
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @return a character vector of genomeas and gtf downloaded
#' @family STAR
#' @export
getGenomeAndFasta <- function(organism, output.dir, db = "ensembl",
                              GTF = TRUE, genome = TRUE, phix = FALSE,
                              ncRNA = "", tRNA = "", rRNA = "",
                              gunzip = TRUE, remake = FALSE) {
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) | !remake) {
    message("Loading premade files information,
            do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }

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
    genome  <- biomartr::getGenome(db = db, organism, path = output.dir, gunzip = gunzip)
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
    gtf <- biomartr::getGTF(   db = db, organism, path = output.dir)
    if (gunzip) # unzip gtf file
      gtf <- R.utils::gunzip(gtf, overwrite = TRUE)
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, organism = organism)
    if (genome != FALSE)
      seqlevelsStyle(txdb) <- seqlevelsStyle(FaFile(genome))[1]
    txdb_file <- paste0(gtf, ".db")
    saveDb(txdb, txdb_file)
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
    phix <- biomartr::getGenome(db = "refseq", "Escherichia virus phiX174",
                                path = output.dir, gunzip = FALSE)
    if (gunzip) # unzip gtf file
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
#' Only works for Linux and Mac (unix systems)
#' @param folder path to folder for download, fille will be named
#' "STAR-version", where version is version wanted.
#' @param version default "2.7.4a"
#' @importFrom utils download.file
#' @importFrom utils tar
#' @return path to runnable STAR
#' @export
#' @references https://www.ncbi.nlm.nih.gov/pubmed/23104886
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
  os <- Sys.info()
  bin <- ifelse(os[1] == "Linux",
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
#' @examples
#' # STAR.remove.crashed.genome(index.path = "/home/data/human_index/phix/)
STAR.remove.crashed.genome <- function(index.path, star.path = STAR.install()) {
  message("Trying to remove loaded genome:")
  out <- paste(star.path, "--genomeDir", index.path, "--genomeLoad Remove")
  system(out)
}
