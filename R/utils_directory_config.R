
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
#' 1. fastq files (raw_data)\cr
#' 2. bam files (processed data)\cr
#' 3. references (organism annotation and STAR index)
#'
#' Update or use another config using \code{config.save()} function.
#' @param file file of config for ORFik, default: "~/Bio_data/ORFik_config.csv"
#' @return a named character vector of length 3
#' @export
config <- function(file = "~/Bio_data/ORFik_config.csv") {
  if(!file.exists(file)) {
    message("--------------------------------")
    message("Setting up config file for ORFik")
    message(paste("Saving to: ", file))

    config.save(file, "~/Bio_data/raw_data", "~/Bio_data/processed_data",
                "~/Bio_data/references")
  }
  dt <- data.table::fread(file, header = TRUE)
  if (nrow(dt) != 3) stop("config files must have exactly 3 rows!")

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
#' @return invisible(NULL), file saved to disc
#' @export
config.save <- function(file = "~/Bio_data/ORFik_config.csv",
                        fastq.dir, bam.dir, reference.dir) {
  conf <- data.frame(type = c("fastq", "bam", "ref"),
                     directory = c(fastq.dir, bam.dir, reference.dir))
  data.table::fwrite(conf, file, col.names = TRUE)
  return(invisible(NULL))
}


