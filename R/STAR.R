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
#' It is recommended to run through the RStudio local job tab, to give full info
#' about the run. The system console will not stall, as can happen in happen in
#' normal RStudio console.
#' @param arguments a named character vector containing paths wanted to
#' use for index creation. They must be named correctly:
#' names must be a subset of:
#' c("gtf", "genome", "contaminants", "phix", "rRNA", "tRNA","ncRNA")
#' @param output.dir directory to save indices, default:
#' paste0(dirname(arguments[1]), "/STAR_index/"), where arguments is the
#' arguments input for this function.
#' @param star.path path to STAR, default: STAR.install(),
#' if you don't have STAR installed at default location, it will install it there,
#' set path to a runnable star if you already have it.
#' @param max.cpus integer, default: \code{min(90, BiocParallel:::bpparam()$workers)},
#'  number of threads to use. Default is minimum of 90 and maximum cores - 2. So if you
#'  have 8 cores it will use 6. Note: FASTP will use maximum 16 threads as from testing
#'  I see performance actually degrades using anything higher. From testing I also see
#'  STAR gets no performance gain after ~50 threads. I do suspect this will change
#'  when hard drives gets better in the future.
#' @param max.ram integer, default 30, in Giga Bytes (GB).
#' Maximum amount of RAM allowed for STAR limitGenomeGenerateRAM argument. RULE:
#' idealy 10x genome size, but do not set too close to machine limit. Default fits
#' well for human genome size (3 GB * 10 = 30 GB)
#' @param SAsparse int > 0,  default 1. If you do not have at least 64GB RAM,
#' you might need to set this to 2.
#' suffux array sparsity, i.e.  distance between indices:
#' use bigger numbers to decrease needed RAM at the cost of mapping
#' speed reduction. Only applies to genome, not conaminants.
#' @param tmpDirStar character, default "-". STAR automatic temp folder creation,
#' deleted when done. The directory can not exists, as a safety STAR must make it!.
#' If you are on a NFS file share drive, and you have a non NFS tmp dir,
#' set this to \code{tempfile()} or the manually specified folder to get a
#' considerable speedup!
#' @param script location of STAR index script,
#' default internal ORFik file. You can change it and give your own if you
#' need special alignments.
#' @param remake logical, default: FALSE, if TRUE remake everything specified
#' @param notify_load_existing logical, default TRUE. If annotation exists
#' (defined as: locally (a file called outputs.rds) exists in outputdir),
#' print a small message notifying the user it is not redownloading. Set to
#' FALSE, if this is not wanted
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
                       star.path = STAR.install(),
                       max.cpus = min(90, BiocParallel::bpparam()$workers),
                       max.ram = 30, SAsparse = 1, tmpDirStar = "-",
                       remake = FALSE,
                       script = system.file("STAR_Aligner",
                                            "STAR_MAKE_INDEX.sh",
                                            package = "ORFik"),
                       notify_load_existing = TRUE) {
  finished.file <- paste0(output.dir, "/outputs.rds")
  if (file.exists(finished.file) & !remake) {
    if (notify_load_existing) message("Loading premade index files, ",
                                 "do remake = TRUE if you want to run again")
    return(readRDS(finished.file))
  }
  if (.Platform$OS.type != "unix") stop("STAR is not supported on windows, run through R in WSL!")

  if (!file.exists(script))
    stop("STAR index script not found, check path of script!")
  if (is.null(names(arguments)))
    stop("arguments must have names, see ?STAR.index")
  possible <- c("gtf", "genome", "phix", "rRNA", "tRNA","ncRNA", "contaminants")
  if (!all(names(arguments) %in% possible))
    stop("At least one of arguments with invalid name!")

  sufficient_memory_to_run_this_check(ref_path = output.dir, max.ram)

  # match which indices to make
  exts <- c("g", "f", "p", "r", "t", "n", "c")
  names(exts) <- possible
  hits <- paste0("-", exts[names(arguments)], " ", arguments)

  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  out <- paste("-o", output.dir)
  max.cpus <- paste("-m", max.cpus)
  max.ram <- paste("-R", format(max.ram*1e9, scientific = FALSE))
  SAsparse <- paste("-a", SAsparse)
  tmpDirStar <- paste("-T", tmpDirStar)

  call <- paste(script, out, star.path, max.cpus, max.ram, SAsparse,
                tmpDirStar, paste(hits, collapse = " "))

  message("STAR indexing:\n")
  print(call); print("\n")
  message("Starting indexing at time:")
  print(Sys.time())
  ret <- system(command = call)
  if (ret != 0) stop("STAR INDEX step failed, see error above for more info.")
  message("- Index done")
  saveRDS(object = output.dir, finished.file)


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
#' Can only run on unix systems (Linux, Mac and WSL (Windows Subsystem Linux)),
#' and requires a minimum of 30GB memory on genomes like human, rat, zebrafish etc.\cr
#' If for some reason the internal STAR alignment bash script will not work for you,
#' like if you want more customization of the STAR/fastp arguments.
#' You can copy the internal alignment script,
#' edit it and give that as the script used for this function.\cr
#' The trimmer used is fastp (the fastest I could find), also works on
#' (Linux, Mac and WSL (Windows Subsystem Linux)).
#' If you want to use your own trimmer set file1/file2 to the location of
#' the trimmed files from your program.\cr
#' A note on trimming from creator of STAR about trimming:
#' "adapter trimming it definitely needed for short RNA sequencing.
#' For long RNA-seq, I would agree with Devon that in most cases adapter trimming
#' is not advantageous, since, by default, STAR performs local (not end-to-end) alignment,
#' i.e. it auto-trims." So trimming can be skipped for longer reads.
#' @param input.dir path to fast files to align, the valid input files will be search for from formats:
#' (".fasta", ".fastq", ".fq", or ".fa") with or without compression of .gz.
#' Also either paired end or single end reads. Pairs will automatically be detected from
#' similarity of naming, separated by something as .1 and .2 in the end. If files are renamed, where pairs
#' are not similarily named, this process will fail to find correct pairs!
#' @param index.dir path to STAR index folder. Path returned from ORFik function
#' STAR.index, when you created the index folders.
#' @param fastp path to fastp trimmer, default: install.fastp(), if you
#' have it somewhere else already installed, give the path. Only works for
#' unix (linux or Mac OS), if not on unix, use your favorite trimmer and
#' give the output files from that trimmer as input.dir here.
#' @param paired.end a logical: default FALSE, alternative TRUE. If TRUE, will auto detect
#'  pairs by names. Can not be a combination of both TRUE and FALSE!\cr
#'  If running in folder mode:
#'  The folder must then contain an even number of files
#'  and they must be named with the same prefix and sufix of either
#'   _1 and _2, 1 and 2, etc. If SRR numbers are used, it will start on lowest and
#'   match with second lowest etc.
#' @param steps a character, default: "tr-ge", trimming then genome alignment\cr
#'  steps of depletion and alignment wanted:
#'  The posible candidates you can use are:\cr
#' \itemize{
#'  \item{tr : trim reads}
#'  \item{co : contamination merged depletion}
#'  \item{ph : phix depletion}
#'  \item{rR : rrna depletion}
#'  \item{nc : ncrna depletion}
#'  \item{tR : trna depletion (Mature tRNA, so no intron checks done)}
#'  \item{ge : genome alignment}
#'  \item{all: run steps: "tr-co-ge" or "tr-ph-rR-nc-tR-ge", depending on if you
#'  have merged contaminants or not}
#' }
#'  If not "all", a subset of these ("tr-co-ph-rR-nc-tR-ge")\cr
#'  If co (merged contaminants) is used, non of the specific contaminants can be specified,
#'  since they should be a subset of co.\cr
#'  The step where you align to the genome is usually always included, unless you
#'  are doing pure contaminant analysis or only trimming.
#'  For Ribo-seq and TCP(RCP-seq) you should do rR (ribosomal RNA depletion),
#'  so when you made the
#'  STAR index you need the rRNA step, either use rRNA from .gtf or manual download.
#'  (usually just download a Silva rRNA database
#'  for SSU&LSU at: https://www.arb-silva.de/) for your species.
#' @param adapter.sequence character, default: "auto". Auto detect adapter using fastp
#' adapter auto detection, checking first 1.5M reads. (Auto detection of adapter will
#' not work 100\% of the time (if the library is of low quality), then you must rerun
#' this function with specified adapter from fastp adapter analysis.
#' , using FASTQC or other adapter detection tools, else alignment will most likely fail!).
#' If already trimmed or trimming not wanted:
#' adapter.sequence = "disable" .You can manually assign adapter like:
#' "ATCTCGTATGCCGTCTTCTGCTTG" or "AAAAAAAAAAAAA". You can also specify one of the three
#' presets:\cr
#' \itemize{
#'  \item{illumina (TrueSeq ~75/100 bp sequencing) : AGATCGGAAGAGC}
#'  \item{small_RNA (standard for ~50 bp sequencing): TGGAATTCTCGG}
#'  \item{nextera: CTGTCTCTTATA}
#' }
#' Paired end auto detection uses overlap sequence of pairs, to use the slower
#' more secure paired end adapter detection, specify as: "autoPE".
#' @param quality.filtering logical, default FALSE. Not needed for modern
#' library prep of RNA-seq, Ribo-seq etc (usually < ~ 0.5% of reads are removed).
#' If you are aligning bad quality data, set this to TRUE.\cr
#' These filters will then be applied (default of fastp), filter if:
#' \itemize{
#'  \item{Number of N bases in read : > 5}
#'  \item{Read quality : > 40\% of bases in the read are <Q15}
#' }
#' @param min.length 20, minimum length of aligned read without mismatches
#' to pass filter. Anything under 20 is dangerous, as chance of random hits will
#' become high!
#' @param mismatches 3, max non matched bases. Excludes soft-clipping, this only
#' filters reads that have defined mismatches in STAR.
#' Only applies for genome alignment step.
#' @param trim.front 0, default trim 0 bases on 5' ends.
#' Ignored if tr (trim) is not one of the arguments in "steps".
#' For Ribo-seq use default 0, unless you have 5' end custom barcodes to remove.
#' Alignment to STAR might fail if you have large barcodes, which are not removed!
#' @param trim.tail 0, default trim 0 bases on 3' ends.
#' Ignored if tr (trim) is not one of the arguments in "steps".
#' For Ribo-seq use default 0, unless you have 3' end custom barcodes to remove.
#' Alignment to STAR might fail if you have large barcodes, which are not removed!
#' @param alignment.type default: "Local": standard local alignment with soft-clipping allowed,
#' "EndToEnd" (global): force end-to-end read alignment, does not soft-clip.
#' @param allow.introns logical, default TRUE. Allow large gaps of N in reads
#' during genome alignment, if FALSE:
#' sets --alignIntronMax to 1 (no introns). NOTE: You will still get some spliced reads
#' if you assigned a gtf at the index step.
#' @param include.subfolders "n" (no), do recursive search downwards for fast files if "y".
#' @param resume default: NULL, continue from step, lets say steps are "tr-ph-ge":
#'  (trim, phix depletion, genome alignment) and resume is "ge", you will then use
#'  the assumed already trimmed and phix depleted data and start at genome alignment,
#'  useful if something crashed. Like if you specified wrong STAR version, but the trimming
#'  step was completed. Resume mode can only run 1 step at the time.
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
#' @param keep.contaminants logical, default FALSE. Create and keep
#' contaminant aligning bam files, default is to only keep unaliged fastq reads,
#' which will be further processed in "ge" genome alignment step. Useful if you
#' want to do further processing on contaminants, like specific coverage of
#' specific tRNAs etc.
#' @param keep.contaminants.type logical, default "bam".
#' If aligned files of contaminants are kept, which format to output as,
#' only supports "bam" for now. Fasta / Fastq will be implemented later.
#' @param keep.unaligned.genome logical, default FALSE. Create and keep
#' reads that did not align at the genome alignment step,
#' default is to only keep the aliged bam file. Useful if you
#' want to do further processing on plasmids/custom sequences.
#' @param keep.index.in.memory logical or character, default FALSE (i.e. LoadAndRemove).
#' For STAR.align.single:\cr
#' If TRUE, will keep index in memory, useful if you need to loop over single calls,
#' instead of using STAR.align.folder (remember last run should use FALSE, to remove index).
#' For STAR.align.folder:\cr
#' Only applies to last library, will always keep for all libraries before last.
#' Alternative useful for MAC machines especially is "noShared", for machines
#' that do not support shared memory index, usually gives error: "abort trap 6".
#' @inheritParams STAR.index
#' @return output.dir, can be used as as input in ORFik::create.experiment
#' @family STAR
#' @export
#' @examples
#' # First specify directories wanted (temp directory here)
#' config_file <- tempfile()
#' #config.save(config_file, base.dir = tempdir())
#' #config <- ORFik::config(config_file)
#'
#' ## Yeast RNA-seq samples (small genome)
#' #project <- ORFik::config.exper("chalmers_2012", "Saccharomyces_cerevisiae", "RNA-seq", config)
#' #annotation.dir <- project["ref"]
#' #fastq.input.dir <- project["fastq RNA-seq"]
#' #bam.output.dir <- project["bam RNA-seq"]
#'
#' ## Download some SRA data and metadata (subset to 50k reads)
#' # info <- download.SRA.metadata("SRP012047", outdir = conf["fastq RNA-seq"])
#' # info <- info[1:2,] # Subset to 2 first libraries
#' # download.SRA(info, fastq.input.dir, rename = FALSE, subset = 50000)
#'
#' ## No contaminant depletion:
#' # annotation <- getGenomeAndAnnotation("Saccharomyces cerevisiae", annotation.dir)
#' # index <- STAR.index(annotation)
#' # STAR.align.folder(fastq.input.dir, bam.output.dir,
#' #                   index, paired.end = FALSE) # Trim, then align to genome
#'
#' ## Human Ribo-seq sample (NB! very large genome and libraries!)
#' ## Requires >= 32 GB memory
#' #project <- ORFik::config.exper("subtelny_2014", "Homo_sapiens", "Ribo-seq", config)
#' #annotation.dir <- project["ref"]
#' #fastq.input.dir <- project["fastq Ribo-seq"]
#' #bam.output.dir <- project["bam Ribo-seq"]
#'
#' ## Download some SRA data and metadata (full libraries)
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
                              quality.filtering = FALSE, min.length = 20, mismatches = 3,
                              trim.front = 0, trim.tail = 0, max.multimap = 10,
                              alignment.type = "Local", allow.introns = TRUE,
                              max.cpus = min(90, BiocParallel::bpparam()$workers),
                              include.subfolders = "n", resume = NULL,
                              multiQC = TRUE, keep.contaminants = FALSE,
                              keep.contaminants.type = c("bam", "fastq")[1],
                              keep.unaligned.genome = FALSE,
                              keep.index.in.memory = FALSE,
                              script.folder = system.file("STAR_Aligner",
                                                          "RNA_Align_pipeline_folder.sh",
                                                          package = "ORFik"),
                              script.single = system.file("STAR_Aligner",
                                                          "RNA_Align_pipeline.sh",
                                                          package = "ORFik")) {


  if (is.logical(paired.end)) {
    paired.end <- ifelse(paired.end, "yes", "no")
  } else if(is.character(paired.end)) {
    if (!(paired.end %in% c("yes", "no"))) stop("Argument 'paired.end' must be yes/no")
  } else stop("Argument 'paired.end' must be logical or character yes/no")

  validate_star_input(script.single, index.dir, keep.contaminants.type, alignment.type,
                      allow.introns, keep.index.in.memory, trim.front, trim.tail,
                      steps, include.subfolders, script.folder, mode = "folder")

  cleaning <- system.file("STAR_Aligner", "cleanup_folders.sh",
                         package = "ORFik", mustWork = TRUE)
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))
  quality.filtering <- ifelse(quality.filtering, "-q default", "")
  keep.index.in.memory <- ifelse(is.logical(keep.index.in.memory),
                                 ifelse(keep.index.in.memory, "y", "n"),
                                 keep.index.in.memory)
  keep.contaminants <- ifelse(keep.contaminants, "-K yes", "-K no")
  keep.contaminants.type <- paste("-X", keep.contaminants.type)
  keep.unaligned.genome <- ifelse(keep.unaligned.genome, "-u Fastx", "-u None")

  call <- paste(script.folder, "-f", input.dir, "-o", output.dir,
                "-p", paired.end,
                "-l", min.length, "-T", mismatches, "-g", index.dir,
                "-s", steps, resume, "-a", adapter.sequence,
                "-t", trim.front, "-z", trim.tail,
                "-M", max.multimap, quality.filtering,
                "-A", alignment.type, "-B", allow.introns,"-m", max.cpus,
                "-i", include.subfolders, "-k", keep.index.in.memory,
                keep.contaminants, keep.contaminants.type,
                keep.unaligned.genome, star.path, fastp, "-I",script.single,
                "-C", cleaning)

  return(STAR.align.internal(call, output.dir, multiQC))
}

#' Align single or paired end pair with STAR
#'
#' Given a single NGS fastq/fasta library, or a paired setup of 2 mated
#' libraries. Run either combination of fastq trimming, contamination removal and
#' genome alignment. Works for (Linux, Mac and WSL (Windows Subsystem Linux))
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
#' ## Specify output libraries (using temp config)
#' config_file <- tempfile()
#' #config.save(config_file, base.dir = tempdir())
#' #config <- ORFik::config(config_file)
#' #project <- ORFik::config.exper("yeast_1", "Saccharomyces_cerevisiae", "RNA-seq", config)
#' # Get genome of yeast (quite small)
#' # arguments <- getGenomeAndAnnotation("Saccharomyces cerevisiae", project["ref"])
#' # index <- STAR.index(arguments)
#'
#' ## Make fake reads
#' #genome <- readDNAStringSet(arguments["genome"])
#' #which_chromosomes <- sample(seq_along(genome), 1000, TRUE, prob = width(genome))
#' #nt50_windows <- lapply(which_chromosomes, function(x)
#' # {window <- sample(width(genome[x]) - 51, 1); genome[[x]][seq(window, window+49)]})
#' #nt50_windows <- DNAStringSet(nt50_windows)
#' #names(nt50_windows) <- paste0("read_", seq_along(nt50_windows))
#' #dir.create(project["fastq RNA-seq"], recursive = TRUE)
#' #fake_fasta <- file.path(project["fastq RNA-seq"], "fake-RNA-seq.fasta")
#' #writeXStringSet(nt50_windows, fake_fasta, format = "fasta")
#' ## Align the fake reads and import bam
#' # STAR.align.single(fake_fasta, NULL, project["bam RNA-seq"], index, steps = "ge")
#' #bam_file <- list.files(file.path(project["bam RNA-seq"], "aligned"),
#' #  pattern = "\\.bam$", full.names = TRUE)
#' #fimport(bam_file)
STAR.align.single <- function(file1, file2 = NULL, output.dir, index.dir,
                              star.path = STAR.install(), fastp = install.fastp(),
                              steps = "tr-ge", adapter.sequence = "auto",
                              quality.filtering = FALSE, min.length = 20,
                              mismatches = 3, trim.front = 0, trim.tail = 0,
                              max.multimap = 10, alignment.type = "Local",
                              allow.introns = TRUE,
                              max.cpus = min(90, BiocParallel::bpparam()$workers),
                              resume = NULL, multiQC = FALSE, keep.contaminants = FALSE,
                              keep.unaligned.genome = FALSE,
                              keep.index.in.memory = FALSE,
                              script.single = system.file("STAR_Aligner",
                                                   "RNA_Align_pipeline.sh",
                                                   package = "ORFik")
) {
  validate_star_input(script.single, index.dir, keep.contaminants.type = "bam",
                      alignment.type, allow.introns, keep.index.in.memory,
                      trim.front, trim.tail, steps)

  file2 <- ifelse(is.null(file2), "", paste("-F", file2))
  resume <- ifelse(is.null(resume), "", paste("-r", resume))
  star.path <- ifelse(is.null(star.path), "", paste("-S", star.path))
  fastp <- ifelse(is.null(fastp), "", paste("-P", fastp))
  quality.filtering <- ifelse(quality.filtering, "-q default", "")
  keep.index.in.memory <- ifelse(is.logical(keep.index.in.memory),
                                 ifelse(keep.index.in.memory, "y", "n"),
                                 keep.index.in.memory)
  keep.contaminants <- ifelse(keep.contaminants, "-K yes", "-K no")
  keep.unaligned.genome <- ifelse(keep.unaligned.genome, "-u Fastx", "-u None")
  call <- paste(script.single, "-f", file1, file2, "-o", output.dir,
                "-l", min.length, "-T", mismatches, "-g", index.dir,
                "-s", steps, resume, "-a", adapter.sequence,
                "-t", trim.front, "-z", trim.tail,
                "-A", alignment.type, "-m", max.cpus,
                "-M", max.multimap,
                "-k", keep.index.in.memory, quality.filtering,
                keep.contaminants, keep.unaligned.genome,
                star.path, fastp)

  return(STAR.align.internal(call, output.dir, multiQC))
}

STAR.align.internal <- function(call, output.dir, multiQC = FALSE, steps = "auto") {
  print(paste("Starting time:", Sys.time()))
  print("Full system call:")
  print(call)
  ret <- system(command = call)
  if (ret != 0) stop("STAR alignment step failed, see error above for more info.")
  message("Alignment done")
  if (multiQC && dir.exists(paste0(output.dir,"/aligned/"))) {
    STAR.allsteps.multiQC(output.dir, steps = steps)
  }
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
#'
#' ORFik for now only uses precompiled STAR binaries, so if you already have
#' a STAR version it is adviced to redownload the same version, since
#' STAR genome indices usually does not work between STAR versions.
#' @param folder path to folder for download, fille will be named
#' "STAR-version", where version is version wanted.
#' @param version default "2.7.4a"
#' @importFrom utils download.file untar tar
#' @return path to runnable STAR
#' @export
#' @references https://www.ncbi.nlm.nih.gov/pubmed/23104886
#' @family STAR
#' @examples
#' ## Default folder install:
#' #STAR.install()
#' ## Manual set folder:
#' folder <- "/I/WANT/IT/HERE"
#' #STAR.install(folder, version = "2.7.4a")
#'
STAR.install <- function(folder = "~/bin", version = "2.7.4a") {
  if (.Platform$OS.type != "unix")
    stop("STAR does not work on Windows, install R/ORFik using WSL")

  url <- paste0("https://git", "hub.com/alexdobin/STAR/archive/",
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
#' For windows must be installed through WSL (Windows Subsystem Linux)
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
#' ## With default folder:
#' #install.fastp()
#'
#' ## Or set manual folder:
#' folder <- "~/I/WANT/IT/HERE/"
#' #install.fastp(folder)
install.fastp <- function(folder = "~/bin") {
  if (.Platform$OS.type != "unix")
    stop("On windows OS, install R/ORFik using WSL")

  is_linux <- Sys.info()[1] == "Linux" # else it is mac
  url <- ifelse(is_linux, # else it is mac
                "https://opengene.org/fastp/fastp",
                paste0("https://git", "hub.com/OpenGene/fastp/archive/master.zip"))

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
  message("Downloading fastp, this will be done only once!")
  utils::download.file(url, destfile = path)
  if (!is_linux) { # For mac os
    message("On mac OS, must build fastp, since no precompiled binaries exists")
    message("This will only be done once")
    message("If this fails, please install using conda, and use the path to the conda installation")
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
#' @return return value from system call, 0 if all good.
#' @export
#' @family STAR
#' @examples
#' index.path = "/home/data/human_GRCh38/STAR_INDEX/genomeDir/"
#' # STAR.remove.crashed.genome(index.path = index.path)
#' ## If you have the index argument from STAR.index function:
#' # index.path <- STAR.index()
#' # STAR.remove.crashed.genome(file.path(index.path, "genomeDir"))
#' # STAR.remove.crashed.genome(file.path(index.path, "contaminants_genomeDir"))
STAR.remove.crashed.genome <- function(index.path, star.path = STAR.install()) {
  message("Trying to remove loaded genome:")
  tempdir.used <- file.path(tempdir(), "remove_")
  message(paste("Log files for removal placed in:", tempdir()))
  out <- paste(star.path, "--genomeDir", index.path, "--genomeLoad Remove",
               "--outFileNamePrefix", tempdir.used)
  status <- system(out)
  if (status == 0) {
    message("Genome removed, if still not working",
            " wait 5 minutes or restart system")
  }
  return(status)
}

validate_star_input <- function(script.single, index.dir, keep.contaminants.type, alignment.type,
                                allow.introns, keep.index.in.memory, trim.front, trim.tail, steps,
                                include.subfolders = NULL, script.folder = NULL, mode = "single") {
  if (.Platform$OS.type != "unix") stop("For Windows OS, run through WSL!")

  stopifnot(keep.contaminants.type %in% c("bam", "fastq") & length(keep.contaminants.type) == 1)
  if (keep.contaminants.type == "fastq") stop("Contaminant as fastq not yet implemented,",
                                              " for now use: 'samtools fastq path_to_bam.bam'")
  if (!file.exists(script.single))
    stop("STAR single file alignment script not found, check path of script!")
  if (is.null(steps) || !is.character(steps) || length(steps) > 1)
    stop("steps argument must be length 1 character!")
  if (!dir.exists(index.dir) & steps != "tr")
    stop("STAR index path must be a valid directory called /STAR_index")
  stopifnot(alignment.type %in% c("Local", "EndToEnd"))
  stopifnot(is.logical(allow.introns) & length(allow.introns) == 1)

  stopifnot(is.logical(keep.index.in.memory) | is.character(keep.index.in.memory))
  stopifnot(all(is.numeric(trim.front) & is.numeric(trim.tail)))
  if (mode == "folder") {
    stopifnot(include.subfolders %in% c("y", "n"))
    if (is.null(script.folder) || !file.exists(script.folder))
      stop("STAR folder alignment script not found, check path of script!")
  }
}

