
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
#' @export
#' @return named character vector of paths for experiment
#' @examples
#' ## Save to default config location
#' #config.exper("Alexaki_Human", "Homo_sapiens_GRCh38_101", c("Ribo-seq", "RNA-seq"))
config.exper <- function(experiment, assembly, type,
                         config = ORFik::config()) {
  # Create fastq dir
  dirs <- file.path(config["fastq"], type, experiment, "")
  # Create bam dir
  dirs <- c(dirs, file.path(config["bam"], type, experiment, ""))
  # Create reference dir
  dirs <- c(dirs, file.path(config["ref"], assembly, ""))
  # Add experiment name
  dirs <- c(dirs, paste(experiment, type, sep = "_"))
  names(dirs) <- c(paste("fastq", type), paste("bam", type),
                   "ref", paste("exp", type))
  return(dirs)
}

#' Read directory config for ORFik experiments
#'
#' Defines a folder for:
#' 1. fastq files (raw data)\cr
#' 2. bam files (processed data)\cr
#' 3. references (organism annotation and STAR index)\cr
#' 4. exp (Location to store and load all \code{\link{experiment}} .csv files)
#' Update or use another config using \code{config.save()} function.
#' @param file file of config for ORFik, default: "~/Bio_data/ORFik_config.csv"
#' @return a named character vector of length 3
#' @export
#' @examples
#' ## Make with default config path
#' #config()
#' ## Load another config (not adviced!)
#' config_location <- "/media/Bio_data/ORFik_config.csv"
#' #config(config_location)
config <- function(file = "~/Bio_data/ORFik_config.csv") {
  default.exp.path <- "~/Bio_data/ORFik_experiments/"
  if(!file.exists(file)) {
    message("--------------------------------")
    message("Setting up config file for ORFik")
    message(paste("Saving to: ", file))
    if (!dir.exists(dirname(file))) {
      dir.create(dirname(file))
    }
    config.save(file, "~/Bio_data/raw_data", "~/Bio_data/processed_data",
                "~/Bio_data/references", default.exp.path)
  }
  dt <- data.table::fread(file, header = TRUE)
  stopifnot(colnames(dt) == c("type", "directory"))
  old_config_format <- nrow(dt) == 3
  if (old_config_format) {
    dt <- rbind(dt, data.table(type = "exp", directory = default.exp.path))
    config.save(file, dt$directory[1], dt$directory[2],
                dt$directory[3], dt$directory[4])
  }

  res <- dt$directory
  names(res) <- dt$type
  return(res)
}

#' Save/update directory config for ORFik experiments
#'
#' Defines a folder for fastq files (raw_data), bam files (processed data) and
#' references (organism annotation and STAR index)
#' @inheritParams config
#' @param fastq.dir directory where ORFik puts fastq file directories,
#'  default: config()["fastq"]
#' @param bam.dir directory where ORFik puts bam file directories,
#'  default: config()["bam"]
#' @param reference.dir directory where ORFik puts reference file directories,
#'  default: config()["ref"]
#' @param exp.dir directory where ORFik puts experiment csv files,
#' default: "~/Bio_data/ORFik_experiments/", which is retrieved with
#' \code{config()["exp"]}
#' @return invisible(NULL), file saved to disc
#' @export
#' @examples
#' ## Save at another config location
#' config_location <- "/media/Bio_data/ORFik_config.csv"
#' #config.save(config_location, "/media/Bio_data/raw_data/",
#' # "/media/Bio_data/processed_data", /media/Bio_data/references/)
config.save <- function(file = "~/Bio_data/ORFik_config.csv",
                        fastq.dir, bam.dir, reference.dir,
                        exp.dir = "~/Bio_data/ORFik_experiments/") {
  conf <- data.frame(type = c("fastq", "bam", "ref", "exp"),
                     directory = c(fastq.dir, bam.dir, reference.dir, exp.dir))
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
