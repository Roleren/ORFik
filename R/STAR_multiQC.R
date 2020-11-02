#' Create STAR multiQC plot and table
#'
#' Takes a folder with multiple Log.final.out files
#' from STAR, and create a multiQC report
#' @param folder path to LOGS folder of ORFik STAR runs. Can also be the path
#' to the aligned/ (parent directory of LOGS), then it will move into LOG
#' from there. Only if no files with pattern Log.final.out are found in
#' parent directory. If no LOGS folder is found it can check for a folder
#' /aligned/LOGS/ so to go 2 folders down.
#' @return invisible(NULL), plot and data saved to disc. Named:
#' "/00_STAR_LOG_plot.png" and "/00_STAR_LOG_table.csv"
#' @importFrom data.table melt
#' @family STAR
#' @export
STAR.multiQC <- function(folder) {

  if (!dir.exists(folder)) stop("folder does not specify existing directory!")
  pattern <- "Log.final.out"
  log_files <- dir(folder, pattern, full.names = TRUE)
  if (length(log_files) == 0) {
    # Go to subfolder directory called /LOGS/ or /aligned/LOGS/
    new_path <- ifelse(dir.exists(paste0(folder, "/aligned/LOGS/")),
                       paste0(folder, "/aligned/LOGS/"),
                       paste0(folder, "/LOGS/"))
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
  fwrite(dt_f, paste0(folder, "/00_STAR_LOG_table.csv"))
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

  ggsave(paste0(folder, "/00_STAR_LOG_plot.png"), gg_STAR, width = 18, height = 9)
  return(invisible(NULL))
}
