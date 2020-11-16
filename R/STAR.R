#' Create STAR genome index
#'
#' Used as reference when aligning data \cr
#' Get genome and gtf by running getGenomeAndFasta()
#'
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.\cr
#' If for some reason the internal STAR index bash script will not work for you,
#' like if you have a very small genome. You can copy the internal index script,
#' edit it and give that as the Index script used for this function.
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
#'  number of threads to use. Default is minimum of 90 and maximum cores - 1. So if you
#'  have 8 cores it will use 7.
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
  possible <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA", "contaminants")
  if (!all(names(arguments) %in% possible))
    stop("At least one of arguments with invalid name!")

  # match which indices to make
  exts <- c("g", "f", "p", "r", "t", "n", "c")
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
    #memory_GB <- as.integer(gsub(" |MemTotal:|kB", replacement = "", a)) / 1e6

    message("Starting indexing at time:")
    print(Sys.time())
    out <- system(command = full, wait = wait)
    out <- ifelse(out == 0 & wait, "Index done", "Index process failed!")
    if (!wait)
      out <- "Wait for index to be complete before you run Alignment!"
    message(out)
    saveRDS(object = output.dir, finished.file)
  } else stop("STAR is not supported on windows!")

  return(output.dir)
}

#' Align all libraries in folder with STAR
#'
#' Does either all files as paired end or single end,
#' so if you have mix, split them in two different folders.\cr
#' If STAR halts at .... loading genome, it means the STAR
#' index was aborted early, then you need to run:
#' STAR.remove.crashed.genome(), with the genome that crashed, and rerun.
#'
#' Can only run on unix systems (Linux and Mac), and requires
#' minimum 30GB memory on genomes like human, rat, zebrafish etc.\cr
#' If for some reason the internal STAR alignment bash script will not work for you,
#' like if you have a very small genome. You can copy the internal alignment script,
#' edit it and give that as the Index script used for this function.\cr
#' The trimmer used is fastp (the fastest I could find), works on mac and linux.
#' If you want to use your own trimmer set file1/file2 to the location of
#' the trimmed files from your program.\cr
#' A note on trimming from creator of STAR about trimming:
#' "adapter trimming it definitely needed for short RNA sequencing.
#' For long RNA-seq, I would agree with Devon that in most cases adapter trimming
#' is not advantageous, since, by default, STAR performs local (not end-to-end) alignment,
#' i.e. it auto-trims." So trimming can be skipped for longer reads.
#' @param input.dir path to fast files to align, the valid input files will be search for from formats:
#' fast files (.fasta, .fastq, .fq, or.fa) with or without compression of .gz.
#' Also either paired end or single end reads. Pairs will automatically be detected from
#' similarity of naming, usualy with a .1 and .2 in the end. If files are renamed, where pairs
#' are not similarily named, this process will fail to find correct pairs.
#' @param index.dir path to STAR index folder. Path returned from ORFik function
#' STAR.index, when you created the index folders.
#' @param fastp path to fastp trimmer, default: install.fastp(), if you
#' have it somewhere else already installed, give the path. Only works for
#' unix (linux or Mac OS), if not on unix, use your favorite trimmer and
#' give the output files from that trimmer as input.dir here.
#' @param paired.end a logical: default FALSE, alternative TRUE. If TRUE, will auto detect
#'  pairs by names. If yes running on a folder:
#'  The folder must then contain an even number of files
#'  and they must be named with the same prefix and sufix of either
#'   _1 and _2, 1 and 2, etc. If SRR numbers are used, it will start on lowest and
#'   match with second lowest etc.
#' @param steps a character, default: "tr-ge", trimming then genome alignment\cr
#'  steps of depletion and alignment wanted:
#'  The posible candidates you can use are:\cr
#' \itemize{
#'  \item{tr : }{trim reads}
#'  \item{co : }{contamination merged depletion}
#'  \item{ph : }{phix depletion}
#'  \item{rR : }{rrna depletion}
#'  \item{nc : }{ncrna depletion}
#'  \item{tR : }{trna depletion}
#'  \item{ge : }{genome alignment}
#'  \item{all: }{run steps: "tr-co-ge" or "tr-ph-rR-nc-tR-ge", depending on if you
#'  have merged contaminants or not}
#' }
#'  If not "all", a subset of these ("tr-co-ph-rR-nc-tR-ge")\cr
#'  If co (merged contaminants) is used, non of the specific contaminants can be specified,
#'  since they should be a subset of co.\cr
#'  The step where you align to the genome is usually always included, unless you
#'  are doing pure contaminant analysis.
#'  For Ribo-seq and TCP(RCP-seq) you should do rR (ribosomal RNA depletion),
#'  so when you made the
#'  STAR index you need the rRNA step, either use rRNA from .gtf or manual download.
#'  (usually just download a Silva rRNA database
#'  for SSU&LSU at: https://www.arb-silva.de/) for your species.
#' @param adapter.sequence character, default: "auto". Auto detect adapter using fastp
#' adapter auto detection, checking first 1.5M reads. (auto detect adapter, is not
#' very reliable for Ribo-seq, so then you must include a manually specified,
#' else alignment will most likely fail!). If already trimmed or trimming not wanted:
#' adapter.sequence = "disable" .You can manually assign adapter like:
#' "ATCTCGTATGCCGTCTTCTGCTTG" or "AAAAAAAAAAAAA". You can also specify one of the three
#' presets:\cr
#' \itemize{
#'  \item{illumina (standard for 100 bp sequencing): }{AGATCGGAAGAGC}
#'  \item{small_RNA (standard for ~50 bp sequencing): }{TGGAATTCTCGG}
#'  \item{nextera: }{CTGTCTCTTATA}
#' }
#' @param min.length 20, minimum length of aligned read without mismatches
#' to pass filter.
#' @param mismatches 3, max non matched bases. Excludes soft-clipping, this only
#' filters reads that have defined mismatches in STAR.
#' Only applies for genome alignment step.
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
#' @param max.multimap numeric, default 10. If a read maps to more locations than specified,
#' will skip the read. Set to 1 to only get unique mapping reads. Only applies for
#' genome alignment step. The depletions are allowing for multimapping.
#' @param script.folder location of STAR index script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @param script.single location of STAR single file alignment script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @param multiQC logical, default TRUE. Do mutliQC comparison of STAR
#' alignment between all the samples. Outputted in aligned/LOGS folder.
#' See ?STAR.multiQC
#' @inheritParams STAR.index
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
#' @examples
#' # First specify directories wanted
#' annotation.dir <- "~/Bio_data/references/Human"
#' fastq.input.dir <- "~/Bio_data/raw_data/Ribo_seq_subtelny/"
#' bam.output.dir <- "~/Bio_data/processed_data/Ribo_seq_subtelny_2014/"
#'
#' ## Download some SRA data and metadata
#' # info <- download.SRA.metadata("DRR041459", fastq.input.dir)
#' # download.SRA(info, fastq.input.dir, rename = FALSE)
#' ## Now align 2 different ways, without and with contaminant depletion
#'
#' ## No contaminant depletion:
#' # annotation <- getGenomeAndAnnotation("Homo sapiens", annotation.dir)
#' # index <- STAR.index(annotation)
#' # STAR.align.folder(fastq.input.dir, bam.output.dir,
#' #                   index, paired.end = FALSE)
#'
#' ## All contaminants merged:
#' # annotation <- getGenomeAndAnnotation(
#' #    organism = "Homo_sapiens",
#' #    phix = TRUE, ncRNA = TRUE, tRNA = TRUE, rRNA = TRUE,
#' #    output.dir = annotation.dir
#' #    )
#' # index <- STAR.index(annotation)
#' # STAR.align.folder(fastq.input.dir, bam.output.dir,
#' #                   index, paired.end = FALSE,
#' #                   steps = "tr-ge")
STAR.align.folder <- function(input.dir, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              paired.end = FALSE,
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 20, mismatches = 3,
                              trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE,
                              include.subfolders = "n", resume = NULL,
                              max.multimap = 10, multiQC = TRUE,
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
  if (is.logical(paired.end)) {
    paired.end <- ifelse(paired.end, "yes", "no")
  } else if(is.character(paired.end)) {
    if (!(paired.end %in% c("yes", "no"))) stop("Argument 'paired.end' must be yes/no")
  } else stop("Argument 'paired.end' must be logical or character yes/no")

  cleaning <- system.file("STAR_Aligner", "cleanup_folders.sh",
                         package = "ORFik", mustWork = TRUE)
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))

  full <- paste(script.folder, "-f", input.dir, "-o", output.dir,
                "-p", paired.end,
                "-l", min.length, "-T", mismatches, "-g", index.dir,
                "-s", steps, resume, "-a", adapter.sequence,
                "-t", trim.front, "-M", max.multimap,
                "-A", alignment.type, "-m", max.cpus, "-i", include.subfolders,
                star.path, fastp, "-I",script.single, "-C", cleaning)
  if (.Platform$OS.type == "unix") {
    print(full)
    message("Starting alignment at time:")
    print(Sys.time())
    out <- system(command = full, wait = wait)
    out <- ifelse(out == 0, "Alignment done", "Alignment process failed!")
    message(out)
    if (multiQC & wait & (out == "Alignment done") &
          dir.exists(paste0(output.dir,"/aligned/"))) {
      STAR.allsteps.multiQC(output.dir, steps = steps)
    }

  } else stop("STAR is not supported on windows!")
  return(output.dir)
}

#' Align single or paired end pair with STAR
#'
#' Given a single NGS fastq/fasta library, or a paired setup of 2 mated
#' libraries. Run alignment and optionally remove contaminants.
#' @inherit STAR.align.folder
#' @inheritParams STAR.align.folder
#' @param file1 library file, if paired must be R1 file. Allowed formats are:
#' (.fasta, .fastq, .fq, or.fa) with or without compression of .gz. This filename usually
#'  contains a suffix of .1
#' @param file2 default NULL, set if paired end to R2 file. Allowed formats are:
#' (.fasta, .fastq, .fq, or.fa) with or without compression of .gz. This filename usually
#'  contains a suffix of .2
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
#' @examples
#'
#' ## Specify output libraries:
#' output.dir <- "/Bio_data/references/Human"
#' bam.dir <- "data/processed/human_rna_seq"
#' # arguments <- getGenomeAndAnnotation("Homo sapiens", output.dir)
#' # index <- STAR.index(arguments, output.dir)
#' # STAR.align.single("data/raw_data/human_rna_seq/file1.bam", bam.dir,
#' #                    index)
STAR.align.single <- function(file1, file2 = NULL, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              steps = "tr-ge", adapter.sequence = "auto",
                              min.length = 20, mismatches = 3, trim.front = 0,
                              alignment.type = "Local", max.cpus = min(90, detectCores() - 1),
                              wait = TRUE, resume = NULL, max.multimap = 10,
                              script.single = system.file("STAR_Aligner",
                                                   "RNA_Align_pipeline.sh",
                                                   package = "ORFik")
) {
  if (!file.exists(script.single))
    stop("STAR single file alignment script not found, check path of script!")
  # TODO, decide to add this in or not.
  # cleaning <- system.file("STAR_Aligner", "cleanup_folders.sh",
  #                         package = "ORFik", mustWork = TRUE)

  file2 <- ifelse(is.null(file2), "", paste("-F", file2))
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))

  full <- paste(script.single, "-f", file1, file2, "-o", output.dir,
                "-l", min.length, "-T", mismatches, "-g", index.dir,
                "-s", steps, resume, "-a", adapter.sequence,
                "-t", trim.front, "-A", alignment.type, "-m", max.cpus,
                "-M", max.multimap,
                star.path, fastp) # "-C", cleaning
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
