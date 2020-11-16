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
#' @return invisible(NULL), plot and data saved to disc. Named:
#' "/00_STAR_LOG_plot.png" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table melt
#' @family STAR
#' @export
STAR.multiQC <- function(folder, type = "aligned") {
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
    STAR.multiQC(new_path)
    return(invisible(NULL))
  }
  message("Runing alignment MultiQC")
  # Read log files 1 by 1 (only data column)
  dt <- lapply(log_files, function(file)
    fread(file, sep = c("\t"),  blank.lines.skip = TRUE, fill = TRUE)[,2])
  # Read log files 1 by 1 (only info column)
  dt_all <- fread(log_files[1], sep = c("\t"),  blank.lines.skip = TRUE, fill = TRUE)[,1]
  for (i in dt) {
    dt_all <- cbind(dt_all, as.data.table(i))
  }
  sample_names <- gsub("_Log.final.out", "",dir(folder, pattern))
  colnames(dt_all) <- c("Info", sample_names)
  dt_all$Info <- gsub(" \\|", "", dt_all$Info)
  dt_all$Info <- gsub("Number", "#", dt_all$Info)
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
  dt_f <- cbind(sample = factor(sample_names, labels = sample_names,
                                levels = sample_names, ordered = TRUE),
                sample_id = factor(sample_names, labels = as.character(seq(length(sample_names))),
                                   levels = sample_names, ordered = TRUE),
                dt_f)
  # Save table to disc
  fwrite(dt_f, file.path(folder, "00_STAR_LOG_table.csv"))
  # create plot
  dt_plot <- melt(dt_f, id.vars = c("sample", "sample_id"))
  gg_STAR <- ggplot(dt_plot, aes(x=sample_id, y = value, group = sample, fill = sample)) +
    geom_bar(aes(color = sample), stat="identity", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    ylab("STAR alignment statistics value, log10 scaled") +
    xlab("Samples") +
    facet_wrap(  ~ variable, scales = "free") +
    scale_y_log10() +
    theme_minimal()

  ggsave(file.path(folder, "00_STAR_LOG_plot.png"), gg_STAR,
         width = 18, height = 9)
  return(invisible(NULL))
}

#' Create trimming table
#'
#' From fastp runs in ORFik alignment process
#' @param trim_folder folder of trimmed files
#' @return a data.table with 3 columns, raw_library (names of library),
#'  raw_reads (numeric, number of raw reads),
#'  trim_reads (numeric, number of trimmed reads)
trimming.table <- function(trim_folder) {

  raw_library <- dir(trim_folder, "\\.json", full.names = TRUE)

  raw_reads <- data.table()
  trim_reads <- data.table()
  lib <- raw_library[1]
  for (lib in raw_library) { # Read json files
    a <- data.table::fread(lib, nrows = 7, skip = 2, sep = ":", col.names = c("id", "value"))
    a$value <- as.numeric(gsub(",", "", a$value))
    a$id <- gsub("\\\t", "", a$id)
    b <- data.table::fread(lib, nrows = 7, skip = 12, sep = ":", col.names = c("id", "value"))
    b$value <- as.numeric(gsub(",", "", b$value))
    b$id <- gsub("\\\t", "", b$id)
    raw_reads <- rbind(raw_reads, a$value[1])
    trim_reads <- rbind(trim_reads, b$value[1])
  }
  raw_data <- cbind(raw_library = basename(raw_library), raw_reads = raw_reads,
                    trim_reads = trim_reads)

  raw_data$raw_library <- gsub("report_|\\.json$",
                               x = raw_data$raw_library, replacement = "")

  return(raw_data)
}
