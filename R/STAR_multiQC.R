#' Create STAR multiQC plot and table
#'
#' Takes a folder with multiple Log.final.out files
#' from STAR, and create a multiQC report. This is automatically run with
#' STAR.align.folder function.
#' @param folder path to main output folder of STAR run. The folder that contains
#' /aligned/, "/trim/, "contaminants_depletion" etc. To find the LOGS folders in, to
#' use for summarized statistics.
#' @param steps a character, default "auto". Find which steps you did.
#' If manual, a combination of "tr-co-ge". See STAR alignment functions for description.
#' @param plot.ext character, default ".pdf". Which format to save QC plot.
#' Alternative: ".png".
#' @param output.file character, file path, default:
#' file.path(folder, "full_process.csv")
#' @return data.table of main statistics, plots and data saved to disc. Named:
#' "/00_STAR_LOG_plot.pdf" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table merge.data.table
#' @family STAR
#' @export
#' @examples
#' process_dir <- system.file("extdata/test_processing/", package = "ORFik")
#' STAR.allsteps.multiQC(process_dir)
#' STAR.allsteps.multiQC(process_dir, steps = "tr-ge")
STAR.allsteps.multiQC <- function(folder, steps = "auto", plot.ext = ".pdf",
                                  output.file = file.path(folder, "full_process.csv")) {
  if (is.null(steps)) return(data.table())
  if (steps == "auto") {
    steps <- paste(ifelse(dir.exists(file.path(folder, "aligned/")), "ge", ""),
                   ifelse(dir.exists(file.path(folder, "contaminants_depletion/")), "co", ""),
                   ifelse(dir.exists(file.path(folder, "trim/")), "tr", ""), sep = "-")

    if (steps == "--") stop("'folder' is not correct ORFik alignment folder")
  }

  res <- NULL
  if (grepl("ge", steps)) { # If genome alignment done
    aligned <- STAR.multiQC(folder, plot.ext = plot.ext)
    res <- aligned
  }
  if (grepl("co", steps)) { # If contamination depletion was done
    co <- STAR.multiQC(folder, "contaminants_depletion", plot.ext = plot.ext)

    co$sample <- gsub("contaminants_", "", co$sample)
    if (!is.null(res)) {
      if (nrow(co) != nrow(res)) {
        warning("Not equal number of aligned and contaminant depleted log files!
                Skipping contamint statistics for now.")

      } else {
        res_temp <- data.table::merge.data.table(aligned, co, by = "sample",
                                                 suffixes = c("-genome", "-contamination"))
        if (nrow(res_temp) < nrow(res)) {
          warning("Contamination statistics could not be appended safely,
                non default naming of files. Fix them before reruning multiQC for valid output")
        } else res <- res_temp
      }
    } else res <- co

  }

  if (grepl("tr", steps)) {
    tr <- trimming.table(file.path(folder, "trim/"))
    if (!is.null(res)) {
      res_temp <- data.table::merge.data.table(res, tr, by.x = "sample", by.y = "raw_library")
      if (nrow(res_temp) < nrow(res)) {
        warning("Trimming statistics could not be appended safely,
                non default naming of files. Fix them before reruning multiQC for valid output")
      } else res <- res_temp
    } else res <- tr
  }
  colnames(res) <- gsub("\\.x", "", colnames(res))
  message("Final statistics:")


  if (grepl("tr", steps)) {
    if (any(res$`% trimmed` > 40)) {
      warning("A sample lost > 40% of reads during trimming")
    }
    if (grepl("ge", steps)) {
      genome_reads <- if (is.null(res$`total mapped reads #-genome`)) {
        res$`total mapped reads #`
      } else res$`total mapped reads #-genome`

      res$`total mapped reads %-genome vs raw` <-
        round((genome_reads/ res$raw_reads) * 100, 4)
      res$`total mapped reads %-genome vs trim` <-
        round((genome_reads / res$trim_reads) * 100, 4)
      very_low_alignment <- !is.na(res$`total mapped reads %-genome vs trim`) &&
        any(res$`total mapped reads %-genome vs trim` < 3)
      if (very_low_alignment) {
        warning("A sample aligned with < 3%, are you using the correct genome?")
      }
    }
  }

  fwrite(res, output.file)
  res[]
  print(res)
  return(res)
}

#' Create STAR multiQC plot and table
#'
#' Takes a folder with multiple Log.final.out files
#' from STAR, and create a multiQC report
#' @param folder path to LOGS folder of ORFik STAR runs. Can also be the path
#' to the aligned/ (parent directory of LOGS), then it will move into LOG
#' from there. Only if no files with pattern Log.final.out are found in
#' parent directory. If no LOGS folder is found it can check for a folder
#' /aligned/LOGS/ so to go 2 folders down.
#' @param type a character path, default "aligned".
#' Which subfolder to check for. If you want log files for contamination
#' do type = "contaminants_depletion"
#' @param plot.ext character, default ".pdf". Which format to save QC plot.
#' Alternative: ".png".
#' @param log_files character, path to "Log.final.out" STAR files,
#'  default: dir(folder, "Log.final.out", full.names = TRUE)
#' @param simplified_table logical, default TRUE. Subset columns, to
#' the most important ones.
#' @return a data.table with all information from STAR runs,
#'  plot and data saved to disc. Named:
#' "/00_STAR_LOG_plot.pdf" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table melt
#' @family STAR
#' @export
#' @examples
#' #' @examples
#' process_dir <- system.file("extdata/test_processing/", package = "ORFik")
#' STAR.multiQC(process_dir)
#'
STAR.multiQC <- function(folder, type = "aligned", plot.ext = ".pdf",
                         log_files = dir(folder, "Log.final.out", full.names = TRUE),
                         simplified_table = TRUE) {
  if (!(type %in% c("aligned", "contaminants_depletion")))
      stop("type must be either aligned or contaminants_depletion")
  if (!dir.exists(folder)) stop("folder does not specify existing directory!")
  stopifnot(is(log_files, "character"))


  if (length(log_files) == 0) {
    # Go to subfolder directory called /LOGS/ or /aligned/LOGS/
    new_path <- ifelse(dir.exists(file.path(folder, type, "LOGS/")),
                       file.path(folder, type, "LOGS/"),
                       file.path(folder, "LOGS/"))
    return(STAR.multiQC(new_path, type = type, plot.ext = plot.ext))
  }
  message("Runing alignment MultiQC (", type, ")")
  dt_logs <- STAR_read_log_files(log_files)

  # Save table to disc
  fwrite(dt_logs, file.path(folder, "00_STAR_LOG_table.csv"))
  # create plot
  STAR.multiQC_plot(dt_logs, folder, plot.ext)

  if (simplified_table) {
    dt_logs <- dt_logs[, c("sample",
                 "total mapped reads %", "total mapped reads #",
                 "Uniquely mapped reads %","Uniquely mapped reads #",
                 "% of reads multimapped",
                 "# of reads multimapped")]
  }

  return(dt_logs)
}

#' Create trimming table
#'
#' From fastp runs in ORFik alignment process
#' @param trim_folder folder of trimmed files, only reads
#' fastp .json files. Can be NULL if raw_libraries is defined
#' @param raw_libraries character, default:
#'  \code{dir(trim_folder, "\\.json", full.names = TRUE)},
#'  file paths of all json file paths.
#' @param include_adapter logical, default FALSE. If TRUE, will add
#' an extra column: adapter, with adapter found. If not found, it will specify:
#'  "passed".
#' @return a data.table with 6 columns, raw_library (names of library),
#'  raw_reads (numeric, number of raw reads),
#'  trim_reads (numeric, number of trimmed reads),
#'  % trimmed (numeric, percentage of trimmed reads),
#'  raw_mean_length (numeric, raw mean read length),
#'  trim_mean_length (numeric, trim mean read length).
#'  Optional columns:
#'  adapter (character, adapter, if not found "passed")
#' @export
#' @importFrom jsonlite fromJSON
#' @examples
#' # Location of fastp trimmed .json files
#' trimmed_file <- system.file("extdata/test_processing/trim",
#'  "output_template.json", package = "ORFik")
#' trimmed_folder <- dirname(trimmed_file)
#' trimming.table(trimmed_folder)
#' trimming.table(NULL, trimmed_file)
#' trimming.table(NULL, trimmed_file, include_adapter = TRUE)
trimming.table <- function(trim_folder,
                           raw_libraries = dir(trim_folder, "\\.json$", full.names = TRUE),
                           include_adapter = FALSE) {

  raw_library <- raw_libraries
  if (length(raw_library) == 0) stop("fastp .json files not found!,",
                                     " did you change wd or deleted them?")
  raw_reads <- data.table()
  trim_reads <- data.table()
  raw_mean_length <- data.table()
  trim_mean_length <- data.table()
  adapters <- data.table()
  lib <- raw_library[1]
  for (lib in raw_library) { # Read json files
    json <- jsonlite::fromJSON(lib)
    summary <- json$summary
    raw_reads  <- rbind(raw_reads,  summary$before_filtering$total_reads)
    trim_reads <- rbind(trim_reads, summary$after_filtering$total_reads)
    raw_mean_length  <- rbind(raw_mean_length,
                              summary$before_filtering$read1_mean_length)
    trim_mean_length <- rbind(trim_mean_length,
                              summary$after_filtering$read1_mean_length)
    if (include_adapter) {
      adapter <- json$adapter_cutting$read1_adapter_sequence
      if (is.null(adapter)) adapter <- "passed"
      adapters <- rbind(adapters, adapter)
    }
  }
  if (!include_adapter) adapters <- NULL
  raw_data <- cbind(raw_library = basename(raw_library), raw_reads = raw_reads,
                    trim_reads = trim_reads,
                    `% trimmed` = round((1 - (trim_reads / raw_reads)) * 100, 3),
                    raw_mean_length = raw_mean_length, trim_mean_length = trim_mean_length,
                    adapter = adapters)

  raw_data$raw_library <- gsub("report_|\\.json$",
                               x = raw_data$raw_library, replacement = "")
  colnames(raw_data) <- gsub("\\.x", "", colnames(raw_data))

  return(raw_data)
}

STAR.multiQC_plot <- function(dt_f, folder, plot.ext = ".pdf") {
  dt_plot <- melt(dt_f, id.vars = c("sample", "sample_id"))
  sample.col <- if (nrow(dt_f) > 12) {
    dt_plot$sample } else NULL
  dt_plot[value == 0, value := 1]
  gg_STAR <- ggplot(dt_plot,
                    aes(x=sample_id, y = value, group = sample, fill = sample)) +
    geom_bar(aes(color = sample.col), stat="identity", position=position_dodge()) +
    scale_fill_brewer(palette="Paired")+
    ylab("Value (log10)") +
    xlab("Samples") +
    facet_wrap(  ~ variable, scales = "free") +
    scale_y_log10() +
    theme_minimal()

  ggsave(file.path(folder, paste0("00_STAR_LOG_plot", plot.ext)), gg_STAR,
         width = 18, height = 9)
  return(invisible(NULL))
}

STAR_read_log_files <- function(log_files, clean_output = TRUE) {
  # Read log files 1 by 1 (only data column)
  dt <- lapply(log_files, function(file)
    fread(file, sep = c("\t"),  blank.lines.skip = TRUE, fill = TRUE)[,2])
  # Read log files, only 1 (only info column)
  dt_all <- fread(log_files[1], sep = c("\t"),
                  blank.lines.skip = TRUE, fill = TRUE)[,1]
  for (i in dt) {
    dt_all <- cbind(dt_all, as.data.table(i))
  }

  if (clean_output) {
    sample_names <- gsub("_Log.final.out", "", basename(log_files))
    colnames(dt_all) <- c("Info", sample_names)
    dt_all$Info <- gsub(" \\|", "", dt_all$Info)
    dt_all$Info <- gsub("Number", "#", dt_all$Info, ignore.case = TRUE)
    dt_all$Info <- gsub("mapped to multiple loci", "multimapped", dt_all$Info)
    dt_all$Info <- gsub("too many mismatches", "mismatches", dt_all$Info)
    dt_dates <- dt_all[1:3, ]

    dt_data <- dt_all[-c(1:3, grep("READS", dt_all$Info)),]

    dt_temp <- dt_data[,-1]
    for (i in seq(ncol(dt_temp))) {
      dt_temp[,i] <- as.numeric(gsub("%", "", unlist(dt_temp[, i, with = FALSE])))
    }

    dt_f <- dt_temp
    dt_f <- data.table(t(dt_f))
    colnames(dt_f) <- unlist(dt_data[,1])
    dt_f <- cbind(`total mapped reads %` = dt_f$`Uniquely mapped reads %` + dt_f$`% of reads multimapped`,
                  `total mapped reads #` = dt_f$`Uniquely mapped reads #` + dt_f$`# of reads multimapped`,
                  dt_f)
    dt_all <- cbind(sample = factor(sample_names, labels = sample_names,
                                  levels = sample_names, ordered = TRUE),
                  sample_id = factor(sample_names, labels = as.character(seq(length(sample_names))),
                                     levels = sample_names, ordered = TRUE),
                  dt_f)
  }
  return(dt_all)
}
