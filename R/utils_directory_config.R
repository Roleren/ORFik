
#' Set directories for experiment
#'
#' Defines a folder for:
#' 1. fastq files (raw_data)\cr
#' 2. bam files (processed data)\cr
#' 3. references (organism annotation and STAR index)\cr
#' 4. Experiment (name of experiment)
#' @param experiment short name of experiment (must be valid as a folder name)
#' @param assembly name of organism and assembly (must be valid as a folder name)
#' @param type name of sequencing type, Ribo-seq, RNA-seq, CAGE.. Can be more than
#' one.
#' @param config a named character vector of length 3,
#'  default: \code{ORFik::config()}
#' @param sub_dir_single character, path. Default: \code{file.path(type, experiment, "")}
#' The subdirectory relative to config defined main locations. If defined location
#' should be used directly without making subdirectories, set to "".
#' @param name_with_type_suffix logical, default TRUE. Make fastq name
#' like 'fastq RNA-seq', setting it to FALSE gives name 'fastq'. Only allowed
#' when length(type) == 1
#' @export
#' @return named character vector of paths for experiment
#' @examples
#' # Where should files go in general?
#' ORFik::config()
#' # Paths for project: "Alexaki_Human" containing Ribo-seq and RNA-seq:
#' #config.exper("Alexaki_Human", "Homo_sapiens_GRCh38_101", c("Ribo-seq", "RNA-seq"))
config.exper <- function(experiment, assembly, type,
                         config = ORFik::config(),
                         sub_dir_single = file.path(type, experiment, ""),
                         name_with_type_suffix = TRUE) {
  stopifnot(is(sub_dir_single, "character"))
  if (length(type) > 1 && all(sub_dir_single == ""))
    stop("sub_dir_single = '' only allowed when length of type is 1")
  stopifnot(is(name_with_type_suffix, "logical"))
  if (length(type) > 1 && !name_with_type_suffix)
    stop("name_with_type_suffix = FALSE only allowed when length of type is 1")
  # Create fastq dir
  dirs <- file.path(config["fastq"], sub_dir_single)
  # Create bam dir
  dirs <- c(dirs, file.path(config["bam"], sub_dir_single))
  # Create reference dir
  dirs <- c(dirs, file.path(config["ref"], assembly, ""))
  # Add experiment name
  dirs <- c(dirs, paste(experiment, type, sep = "_"))
  suffix_name <- if(name_with_type_suffix) {paste0(" ", type)} else rep("", length(type))
  names(dirs) <- c(paste0("fastq", suffix_name), paste0("bam", suffix_name),
                   "ref", paste0("exp", suffix_name))
  return(dirs)
}

#' Read directory config for ORFik experiments
#'
#' Defines a folder for:
#' 1. fastq files (raw data)\cr
#' 2. bam files (processed data)\cr
#' 3. references (organism annotation and STAR index)\cr
#' 4. experiments (Location to store and load all \code{\link{experiment}} .csv files)
#' Update or use another config using \code{config.save()} function.
#' @param file location of config csv, default:
#'  config_file(old_config_location = old_config_location)
#' @inheritParams config_file
#' @return a named character vector of length 3
#' @import BiocFileCache
#' @export
#' @examples
#' ## Make with default config path
#' #config()
config <- function(file = config_file(old_config_location = old_config_location),
                   old_config_location = "~/Bio_data/ORFik_config.csv") {
  make_config <- (!is.null(names(file)) && names(file) == "new") |
    !file.exists(file)

  if (make_config) {
    message("--------------------------------")
    message("Setting up config file for ORFik")
    message("This config defines where to save:
            fastq, bam files, references and experiments")
    message(paste("Saving to: ", file))
    if (!dir.exists(dirname(file))) {
      dir.create(dirname(file))
    }

    if (file.exists(old_config_location)) {
      message("ORFik now uses BiocFileCache for config, moving your old config")
      old_config <- data.table::fread(old_config_location, header = TRUE)
      config.save(file, conf = old_config)
      file.remove(old_config_location)
    } else config.save(file)
  }
  dt <- data.table::fread(file, header = TRUE)
  stopifnot(colnames(dt) == c("type", "directory"))
  old_config_format <- nrow(dt) == 3
  if (old_config_format) {
    exp_dir <- file.path(dirname(dt$directory[1]), "ORFik_experiments/")
    message("Setting new exp dir in config to: ", exp_dir)
    dt <- rbind(dt, data.table(type = "exp", directory = exp_dir))
    config.save(file, conf = dt)
  }

  res <- dt$directory
  names(res) <- dt$type
  return(res)
}

#' Get path for ORFik config in cache
#'
#' @param cache path to bioc cache directory with rname from query argument.
#' Default is: \code{BiocFileCache::getBFCOption("CACHE")}
#' For info, see: [BiocFileCache::BiocFileCache()]
#' @param query default: "ORFik_config". Exact rname of the file in cache.
#' @param ask logical, default interactive().
#' @param old_config_location path, old config location before BiocFileCache
#' implementation.
#' Will copy this to cache directory and delete old version.
#' This is done to follow bioc rules on not writing to user home directory.
#' @return a file path in cache
#' @export
#' @examples
#' config_file()
#' # Another config path
#' config_file(query = "ORFik_config_2")
#'
config_file <- function(cache = BiocFileCache::getBFCOption("CACHE"),
                        query = "ORFik_config", ask = interactive(),
                        old_config_location = "~/Bio_data/ORFik_config.csv") {
  old_location_format <- file.exists(old_config_location)
  cache <- BiocFileCache::BiocFileCache(cache, ask)
  bfcq <- bfcquery(cache, query, exact = TRUE)
  bfcq_nrow <- nrow(bfcq)
  name <- "new"
  if (old_location_format) {
    if (nrow(bfcq) > 0) {
      BiocFileCache::bfcremove(cache, bfcq$rid)
    }
    file <- BiocFileCache::bfcnew(cache, query, ext = ".csv")
  } else {
    if (nrow(bfcq) > 0) {
      name <- "existing"
      file <- file.path(BiocFileCache::bfccache(cache), bfcq$rpath)[1]
    } else {
      file <- BiocFileCache::bfcnew(cache, query, ext = ".csv")
    }
  }
  names(file) <- name
  return(file)
}

#' Save/update directory config for ORFik experiments
#'
#' Defines a folder for fastq files (raw_data), bam files (processed data) and
#' references (organism annotation and STAR index)
#' @inheritParams config
#' @param fastq.dir directory where ORFik puts fastq file directories,
#'  default: file.path(base.dir, "raw_data"), which is retrieved with:
#' \code{config()["fastq"]}
#' @param bam.dir directory where ORFik puts bam file directories,
#'  default: file.path(base.dir, "processed_data"), which is retrieved with:
#'  \code{config()["bam"]}
#' @param reference.dir directory where ORFik puts reference file directories,
#'  default: file.path(base.dir, "references"), which is retrieved with:
#'  \code{config()["ref"]}
#' @param exp.dir directory where ORFik puts experiment csv files,
#'  default: file.path(base.dir, "ORFik_experiments/"), which is retrieved with:
#'  \code{config()["exp"]}
#' @param base.dir base directory for all output directories, default:
#'  "~/Bio_data"
#' @param conf data.frame of complete conf object, default:
#' data.frame(type = c("fastq", "bam", "ref", "exp"),
#'  directory = c(fastq.dir, bam.dir, reference.dir, exp.dir))
#' @return invisible(NULL), file saved to disc
#' @export
#' @examples
#' # Overwrite default config, with new base directory for files
#' #config.save(base.dir = "/media/Bio_data/") # Output files go here instead
#' # of ~/Bio_data
#' ## Dont do this, but for understanding here is how to make a second config
#' #new_config_path <- config_file(query = "ORFik_config_2")
#' #config.save(new_config_path, "/media/Bio_data/raw_data/",
#' # "/media/Bio_data/processed_data", /media/Bio_data/references/)
config.save <- function(file = config_file(),
                        fastq.dir = file.path(base.dir, "raw_data"),
                        bam.dir = file.path(base.dir, "processed_data"),
                        reference.dir = file.path(base.dir, "references"),
                        exp.dir = file.path(base.dir, "ORFik_experiments/"),
                        base.dir = "~/Bio_data",
                        conf = data.frame(type = c("fastq", "bam", "ref", "exp"),
                                          directory = c(fastq.dir, bam.dir, reference.dir, exp.dir))) {
  stopifnot(nrow(conf) == 4)
  stopifnot(colnames(conf) == c("type", "directory"))
  data.table::fwrite(conf, file, col.names = TRUE)
  return(invisible(NULL))
}

#' List genomes created with ORFik
#'
#' Given the reference.folder, list all valid references.
#' An ORFik genome is defined as a folder with a file called output.rds
#' that is a named R vector with names gtf and genome, where the values
#' are character paths to those files inside that folder. This makes sure
#' that this reference was made by ORFik and not some other program.
#' @param reference.folder character path, default:
#' \code{ORFik::config()["ref"]}.
#' @return a data.table with 5 columns:\cr
#'  - character (name of folder)\cr
#'  - logical (does it have a gtf)\cr
#'  - logical (does it have a fasta genome)\cr
#'  - logical (does it have a STAR index)\cr
#'  - logical (only displayed if some are TRUE, does it have protein structure
#'   predictions of ORFs from alphafold etc, in folder called
#'   'protein_structure_predictions')\cr
#'  - logical (only displayed if some are TRUE, does it have gene symbol fst file
#'    from bioMart etc, in file called 'gene_symbol_tx_table.fst')
#' @export
#' @examples
#' ## Run with default config path
#' #list.genomes()
#' ## Run with custom config path
#' list.genomes(tempdir())
#' ## Get the path to fasta genome of first organism in list
#' #readRDS(file.path(config()["ref"], list.genomes()$name, "outputs.rds")[1])["genome"]
list.genomes <- function(reference.folder = ORFik::config()["ref"]) {
  candidates <- list.dirs(reference.folder, recursive = FALSE)
  outputs <- file.path(list.dirs(reference.folder, recursive = FALSE), "outputs.rds")
  valid <- file.exists(outputs)
  if (sum(valid) == 0) {
    message("No valid ORFik-made references in folder!")
    return(data.table())
  }
  candidates <- candidates[valid]
  outputs <- outputs[valid]

  availableGenomes <- data.table()
  for(i in seq_along(outputs)) {
    out <- outputs[i]
    cand <- candidates[i]
    annotation <- readRDS(outputs[i])
    availableGenomes <- rbindlist(list(availableGenomes,
                                       data.table(name = basename(cand),
                                                  gtf = "gtf" %in% names(annotation),
                                                  genome = "genome" %in% names(annotation))))
  }
  indices <- file.path(candidates, "STAR_index", "outputs.rds")
  availableGenomes[, STAR_index := file.exists(indices)]
  protein_structures <- file.path(candidates, "protein_structure_predictions")
  if (any(file.exists(protein_structures))) availableGenomes[, protein_structures := file.exists(protein_structures)]
  gene_symbols <- file.path(candidates, "gene_symbol_tx_table.fst")
  if (any(file.exists(gene_symbols))) availableGenomes[, gene_symbols := file.exists(gene_symbols)]
  availableGenomes[]
  return(availableGenomes)
}
