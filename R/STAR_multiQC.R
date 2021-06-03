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
#' @return data.table of main statistics, plots and data saved to disc. Named:
#' "/00_STAR_LOG_plot.pdf" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table merge.data.table
#' @family STAR
#' @export
STAR.allsteps.multiQC <- function(folder, steps = "auto", plot.ext = ".pdf") {
  if (steps == "auto") {
    steps <- paste(ifelse(dir.exists(file.path(folder, "aligned/")), "ge", ""),
                   ifelse(dir.exists(file.path(folder, "contaminants_depletion/")), "co", ""),
                   ifelse(dir.exists(file.path(folder, "trim/")), "tr", ""), sep = "-")

    if (steps == "--") stop("'folder' is not correct ORFik alignment folder")
  }


  output.file <- file.path(folder, "full_process.csv")
  res <- NULL
  if (1 %in% grep("ge", steps)){
    # If genome alignment done
    aligned <- STAR.multiQC(folder, plot.ext = plot.ext)
    aligned <- aligned[, c("sample", "sample_id",
                           "total mapped reads %", "total mapped reads #",
                           "Uniquely mapped reads %","Uniquely mapped reads #",
                           "% of reads multimapped",
                           "# of reads multimapped")]
    res <- aligned
  }
  if (1 %in% grep("co", steps)) {
    # If contamination depletion was done
    co <- STAR.multiQC(folder, "contaminants_depletion", plot.ext = plot.ext)
    co <- co[, c("sample",
                 "total mapped reads %", "total mapped reads #",
                 "Uniquely mapped reads %","Uniquely mapped reads #",
                 "% of reads multimapped",
                 "# of reads multimapped")]
    co$sample <- gsub("contaminants_", "", co$sample)
    if (!is.null(res)) {
      res <- data.table::merge.data.table(aligned, co, by = "sample",
                                          suffixes = c("-genome", "-contamination"))
    } else res <- co

  }

  if (1 %in% grep("tr", steps)) {
    tr <- ORFik:::trimming.table(file.path(folder, "trim/"))
    if (!is.null(res)) {
      res <- data.table::merge.data.table(res, tr, by.x = "sample", by.y = "raw_library")
    } else res <- tr

  }
  colnames(res) <- gsub("\\.x", "", colnames(res))
  message("Final statistics:")


  if (1 %in% grep("tr", steps)) {
    if (any(res$`% trimmed` > 40)) {
      warning("A sample lost > 40% of reads during trimming")
    }
    if (1 %in% grep("ge", steps)) {
      genome_reads <- if (is.null(res$`total mapped reads #-genome`)) {
        res$`total mapped reads #`
      } else res$`total mapped reads #-genome`

      res$`total mapped reads %-genome vs raw` <-
        round((genome_reads/ res$raw_reads) * 100, 4)
      res$`total mapped reads %-genome vs trim` <-
        round((genome_reads / res$trim_reads) * 100, 4)

      if (any(res$`total mapped reads %-genome vs trim` < 3)) {
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
#' @return a data.table with all information from STAR runs,
#'  plot and data saved to disc. Named:
#' "/00_STAR_LOG_plot.pdf" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table melt
#' @family STAR
#' @export
STAR.multiQC <- function(folder, type = "aligned", plot.ext = ".pdf") {
  if (!(type %in% c("aligned", "contaminants_depletion")))
      stop("type must be either aligned or contaminants_depletion")
  if (!dir.exists(folder)) stop("folder does not specify existing directory!")
  pattern <- "Log.final.out"
  log_files <- dir(folder, pattern, full.names = TRUE)
  if (length(log_files) == 0) {
    # Go to subfolder directory called /LOGS/ or /aligned/LOGS/
    new_path <- ifelse(dir.exists(file.path(folder, type, "LOGS/")),
                       file.path(folder, type, "LOGS/"),
                       file.path(folder, "LOGS/"))
    return(STAR.multiQC(new_path, type = type, plot.ext = plot.ext))
  }
  message("Runing alignment MultiQC")
  # Read log files 1 by 1 (only data column)
  dt <- lapply(log_files, function(file)
    fread(file, sep = c("\t"),  blank.lines.skip = TRUE, fill = TRUE)[,2])
  # Read log files, only 1 (only info column)
  dt_all <- fread(log_files[1], sep = c("\t"),
                  blank.lines.skip = TRUE, fill = TRUE)[,1]
  for (i in dt) {
    dt_all <- cbind(dt_all, as.data.table(i))
  }
  sample_names <- gsub("_Log.final.out", "",dir(folder, pattern))
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
  dt_f <- cbind(sample = factor(sample_names, labels = sample_names,
                                levels = sample_names, ordered = TRUE),
                sample_id = factor(sample_names, labels = as.character(seq(length(sample_names))),
                                   levels = sample_names, ordered = TRUE),
                dt_f)
  # Save table to disc
  fwrite(dt_f, file.path(folder, "00_STAR_LOG_table.csv"))
  # create plot
  dt_plot <- melt(dt_f, id.vars = c("sample", "sample_id"))
  sample.col <- if (nrow(dt_f) > 12) {
    dt_plot$sample } else NULL

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
  return(dt_f)
}

#' Create trimming table
#'
#' From fastp runs in ORFik alignment process
#' @param trim_folder folder of trimmed files
#' @return a data.table with 4 columns, raw_library (names of library),
#'  raw_reads (numeric, number of raw reads),
#'  trim_reads (numeric, number of trimmed reads),
#'  % trimmed (numeric, percentage of trimmed reads)
trimming.table <- function(trim_folder) {

  raw_library <- dir(trim_folder, "\\.json", full.names = TRUE)
  if (length(raw_library) == 0) stop("fastp .json files not found!,",
                                     " did you change wd delete them?")
  raw_reads <- data.table()
  trim_reads <- data.table()
  lib <- raw_library[1]
  for (lib in raw_library) { # Read json files
    a <- data.table::fread(lib, nrows = 7, skip = 2, sep = ":", col.names = c("id", "value"))
    a$value <- as.numeric(gsub(",", "", a$value))
    a$id <- gsub("\\\t", "", a$id)
    # Paired end reads
    b <- data.table::fread(lib, nrows = 7, skip = 13, sep = ":",
                           col.names = c("id", "value"))
    if (length(grep("_reads", b$id[1])) == 0) { # Single end reads
      b <- data.table::fread(lib, nrows = 7, skip = 12, sep = ":",
                             col.names = c("id", "value"))
    }

    b$value <- as.numeric(gsub(",", "", b$value))
    b$id <- gsub("\\\t", "", b$id)
    raw_reads <- rbind(raw_reads, a$value[1])
    trim_reads <- rbind(trim_reads, b$value[1])
  }
  raw_data <- cbind(raw_library = basename(raw_library), raw_reads = raw_reads,
                    trim_reads = trim_reads,
                    `% trimmed` = round((1 - (trim_reads / raw_reads)) * 100, 3))

  raw_data$raw_library <- gsub("report_|\\.json$",
                               x = raw_data$raw_library, replacement = "")

  return(raw_data)
}
