#' Plot DTEG result
#'
#' For explanation of plot catagories, see \code{\link{DTEG.analysis}}
#' @param dt a data.table with the results from \code{\link{DTEG.analysis}}
#' @param output.dir a character path, default NULL(no save), or a directory
#' to save to a file will be called "DTEG_plot.pdf"
#' @param p.value a numeric, default 0.05 in interval (0,1)
#' or "" to not show.
#' What p-value used for the analysis? Will be shown as a caption.
#' @param plot.title title for plots, usually name of experiment etc
#' @param width numeric, default 6 (in inches)
#' @param height numeric, default 6 (in inches)
#' @param dot.size numeric, default 0.4, size of point dots in plot.
#' @param xlim numeric vector, default c(-5, 5)
#' @param ylim numeric vector, default c(-10, 10)
#' @param relative.name character, Default: "DTEG_plot.pdf".
#' Relative name of file to be saved in folder specified in output.dir.
#' Change to .pdf if you want pdf file instead of png.
#' @return a ggplot object
#' @family TE
#' @export
#' @examples
#' #df.rfp <- read.experiment("Riboseq")
#' #df.rna <- read.experiment("RNAseq")
#' #dt <- DTEG.analysis(df.rfp, df.rna)
#' #DTEG.plot(dt, xlim = c(-2, 2), ylim = c(-2, 2))
DTEG.plot <- function(dt, output.dir = NULL,
                      p.value = 0.05,
                      plot.title = "", plot.ext = ".pdf", width = 6,
                      height = 6, dot.size = 0.4,
                      xlim = c(-5, 5), ylim = c(-10, 10),
                      relative.name = paste0("DTEG_plot", plot.ext)) {
  color.values <- c("black", "orange4", "purple", "darkgreen", "blue", "yellow", "aquamarine")
  regulation.levels <- c("No change", "Translation", "Buffering", "mRNA abundance", "Expression", "Forwarded", "Inverse")
  dt[, Regulation :=
               factor(Regulation,
                      levels = regulation.levels,
                      ordered = TRUE)]
  p.caption <- if (p.value != "") {
    labs(caption = paste("P-value <", p.value))
  } else NULL
  p.title <- if (plot.title != "") {
    ggtitle(label = plot.title)
  } else NULL

  dot.size <- rep(dot.size, nrow(dt))
  dot.size[dt$Regulation != "No change"] <- dot.size[1]*2
  plot.between <- ggplot(data = dt,
                         aes(x = rna, y = rfp, color = Regulation)) +
    geom_point(alpha = 0.5, size = dot.size) +
    scale_color_manual(values = color.values) +
    theme_minimal() +
    geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
    geom_vline(aes(xintercept =  0), alpha = 0.2, color = "red") +
    geom_abline(alpha = 0.2, color = "gray67", linetype = "dashed", intercept = 0, slope = 1) +
    xlab("mRNA (log2 fold change)") +
    ylab("RFP (log2 fold change)") +
    p.title +
    p.caption +
    facet_wrap(~ variable, ncol = 2) +
    xlim(xlim) + ylim(ylim) +
    guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 1.3)))
  plot(plot.between)
  if (!is.null(output.dir)) {
    ggsave(file.path(output.dir, relative.name), plot.between,
           width = width, height = height, dpi = 300)
  }
  return(plot.between)
}

#' Translational efficiency plots
#'
#' Create TE plot of:\cr
#' - Within sample (TE log2 vs mRNA fpkm)\cr
#'
#' @inheritParams DTEG.plot
#' @inheritParams te.table
#' @param dt a data.table with the results from \code{\link{te.table}}
#' @param output.dir a character path, default NULL(no save), or a directory
#' to save to a file will be called "TE_within.pdf"
#' @param height a numeric, width of plot in inches. Default "auto".
#' @return a ggplot object
#' @family TE
#' @export
#' @examples
#' #df.rfp <- read.experiment("Riboseq")
#' #df.rna <- read.experiment("RNAseq")
#' #dt <- te.table(df.rfp, df.rna)
#' #te_rna.plot(dt)
te_rna.plot <- function(dt, output.dir = NULL,
                        filter.rfp = 1, filter.rna = 1,
                        plot.title = "", plot.ext = ".pdf",
                        width = 6, height = "auto",
                        dot.size = 0.4) {

  if (height == "auto") height <- 3+length(unique(dt$variable))
  caption <- paste("Filter: RFP >", filter.rfp, " & mRNA >", filter.rna, "(FPKM)")
  if (nrow(df.rfp) > 1 & nrow(df.rna) == 1)
    caption <- paste(subtitle, "(Single mRNA sample)")
  p.caption <- if (filter.rfp != "") {
    labs(caption = caption)
  } else NULL
  p.title <- if (plot.title != "") {
    ggtitle(label = plot.title)
  } else NULL

  plot <- ggplot(data = dt, aes(x = rna_log10, y = TE_log2)) +
    geom_point(alpha = 0.3, size = dot.size) +
    theme_minimal() +
    geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
    xlab("mRNA FPKM (log10)") +
    ylab("TE (log2)") +
    p.caption +
    p.title +
    ggtitle(label = plot.title) +
    xlim(c(filter.rna, filter.rna + 2.5)) +
    facet_wrap(~ variable, nrow = 1)

  plot(plot)
  if (!is.null(output.dir)) {
    ggsave(file.path(output.dir, paste0("TE_within", plot.ext)), plot,
           width = width, height = height, dpi = 300)
  }
  return(plot)
}

#' Translational efficiency plots
#'
#' Create 2 TE plots of:\cr
#' - Within sample (TE log2 vs mRNA fpkm) ("default")\cr
#' - Between all combinations of samples
#' (x-axis: rna1fpkm - rna2fpkm, y-axis rfp1fpkm - rfp2fpkm)
#'
#' Ribo-seq and RNA-seq must have equal nrows, with matching samples. Only
#' exception is if RNA-seq is 1 single sample. Then it will use that for
#' each of the Ribo-seq samples.
#' Same stages, conditions etc, with a unique pairing 1 to 1. If not you can
#' run collapse = "all". It will then merge all and do combined of all
#' RNA-seq vs all Ribo-seq
#' @param df.rfp a \code{\link{experiment}} of Ribo-seq or 80S from TCP-seq.
#' @param df.rna a \code{\link{experiment}} of RNA-seq
#' @param output.dir directory to save plots, plots will be named
#' "TE_between.pdf" and "TE_within.pdf"
#' @param type which plots to make, default: c("default", "between"). Both plots.
#' @param filter.rfp numeric, default 1. minimum fpkm value to be included in plots
#' @param filter.rna numeric, default 1. minimum fpkm value to be included in plots
#' @param plot.title title for plots, usually name of experiment etc
#' @inheritParams countTable
#' @param width numeric, default 6 (in inches)
#' @param height numeric or character, default "auto", which is:
#' 3 + (ncol(RFP_CDS_FPKM)-2).
#' Else a numeric value of height (in inches)
#' @return a data.table with TE values, fpkm and log fpkm values, library
#' samples melted into rows with split variable called "variable".
#' @export
#' @examples
#' ##
#' # df.rfp <- read.experiment("zf_baz14_RFP")
#' # df.rna <- read.experiment("zf_baz14_RNA")
#' # te.plot(df.rfp, df.rna)
#' ## Collapse replicates:
#' # te.plot(df.rfp, df.rna, collapse = TRUE)
te.plot <- function(df.rfp, df.rna,
                    output.dir = paste0(dirname(df.rfp$filepath[1]),
                                        "/QC_STATS/"),
                    type = c("default", "between"),
                    filter.rfp = 1, filter.rna = 1,
                    collapse = FALSE,
                    plot.title = "", plot.ext = ".pdf",
                    width = 6,
                    height = "auto") {
  RNA_MRNA_FPKM <- countTable(df.rna, "mrna", type = "fpkm", collapse = collapse)
  RNA_MRNA_FPKM <- data.table(id = rownames(RNA_MRNA_FPKM), RNA_MRNA_FPKM)

  RFP_CDS_FPKM <- countTable(df.rfp, "cds", type = "fpkm", collapse = collapse)
  RFP_CDS_FPKM <- data.table(id = rownames(RFP_CDS_FPKM), RFP_CDS_FPKM)

  if (!identical(nrow(RNA_MRNA_FPKM), nrow(RFP_CDS_FPKM)))
    stop("Not equal rows in count tables, did you match rfp and rna from different genome builds?")
  single.rna <- length(bamVarName(df.rna, skip.libtype = TRUE)) == 1
  if (!single.rna & !identical(length(bamVarName(df.rfp, skip.libtype = TRUE)),
                               length(bamVarName(df.rna, skip.libtype = TRUE))))
    stop("Not equal samples in rfp and rna, did you subset or reorder one of the experiments?")
  if ((filter.rfp < 0) | (filter.rna < 0))
    stop("filter value is < 0, not allowed!")

  dt <- data.table::merge.data.table(RNA_MRNA_FPKM, RFP_CDS_FPKM, by = "id")
  txNames <- dt$id; dt$id <- NULL

  #dt <- dt[RNA_MRNA_FPKM > filter.rfp & RFP_CDS_FPKM > filter.rna, ]
  dt <- dt[rowMin(as.matrix(dt)) > max(filter.rna, filter.rfp), ]
  dt.log <- pseudo.transform(dt)
  dt.log10 <- pseudo.transform(dt, log10)
  dt.melt.rna <- melt(dt.log[, colnames(dt.log) %in% colnames(RNA_MRNA_FPKM)[-1], with = FALSE])
  dt.melt.rna.10 <- melt(dt.log10[, colnames(dt.log10) %in% colnames(RNA_MRNA_FPKM)[-1], with = FALSE])
  dt.melt.rfp <- melt(dt.log[, colnames(dt.log) %in% colnames(RFP_CDS_FPKM)[-1], with = FALSE])
  dt.final <- cbind(dt.melt.rfp, dt.melt.rna$value, dt.melt.rna.10$value)
  filter.names <- paste("RFP", "LSU", paste0(df.rfp@experiment, "_"), sep = "|")
  dt.final[, variable := gsub(filter.names, "", variable)]
  dt.final[, variable := gsub("^_|_$", "", variable)]
  colnames(dt.final) <- c("variable", "rfp_log2","rna_log2", "rna_log10")
  dt.final[, LFC_TE := rfp_log2 - rna_log2]

  message(paste("Filter kept", round((nrow(dt) / length(txNames)) *100, 1), "% of transcripts"))

  if (height == "auto") height <- 3+(ncol(RFP_CDS_FPKM)-2)
  subtitle <- paste("Filter: RFP >", filter.rfp, " & mRNA >", filter.rna, "(FPKM)")
  if (nrow(df.rfp) > 1 & nrow(df.rna) == 1)
    subtitle <- paste(subtitle, "(Single mRNA sample)")
  if ("default" %in% type) {
    plot <- ggplot(data = dt.final) +
      geom_point(aes(x = rna_log10, y = LFC_TE), alpha = 0.2) +
      theme_minimal() +
      geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
      xlab("mRNA FPKM (log10)") +
      ylab("TE (log2 fold change)") +
      ggtitle(label = plot.title, subtitle = subtitle) +
      xlim(c(filter.rna, filter.rna + 2.5)) +
      facet_wrap(~ variable, ncol = 1)

    plot(plot)
    ggsave(file.path(output.dir, paste0("TE_within", plot.ext)), plot,
           width = width, height = height, dpi = 300)
  }
  if ("between" %in% type) {
    pairs <- combn.pairs(dt.final$variable)

    dt.between <- data.table()
    for (i in pairs) {
      name <- paste("Comparison:", i[1], "vs", i[2])
      rna <- dt.final[variable == i[1], rna_log2] - dt.final[variable == i[2], rna_log2]
      rfp <- dt.final[variable == i[1], rfp_log2] - dt.final[variable == i[2], rfp_log2]
      dt.between <- rbindlist(list(dt.between, data.table(variable = name, rna, rfp)))
    }
    plot <- ggplot(data = dt.between) +
      geom_point(aes(x = rna, y = rfp), alpha = 0.2) +
      theme_minimal() +
      geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
      geom_vline(aes(xintercept =  0), alpha = 0.2, color = "red") +
      xlab("mRNA (log2 fold change)") +
      ylab("RFP (log2 fold change)") +
      ggtitle(label = plot.title, subtitle) +
      facet_wrap(~ variable, ncol = 2) +
      xlim(c(-5, 5))
    plot(plot)
    ggsave(file.path(output.dir, paste0("TE_between", plot.ext)), plot,
           width = width, height = height, dpi = 300)
  }
  return(dt)
}
