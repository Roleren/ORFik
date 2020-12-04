
#' Translational efficiency plots
#'
#' Create 2 TE plots of:\cr
#' - Within sample (TE log2 vs mRNA fpkm) ("default")\cr
#' - Between all combinations of samples
#' (x-axis: rna1fpkm - rna2fpkm, y-axis rfp1fpkm - rfp2fpkm)
#'
#' Ribo-seq and RNA-seq must have equal nrows, with matching samples.
#' Same stages, conditions etc, with a unique pairing 1 to 1. If not you can
#' run collapse = "all". It will then merge all and do combined of all
#' RNA-seq vs all Ribo-seq
#' @param df.rfp a \code{\link{experiment}} of Ribo-seq
#' @param df.rna a \code{\link{experiment}} of RNA-seq
#' @param output.dir directory to save plots, plots will be named
#' "TE_between.png" and "TE_within.png"
#' @param type which plots to make, default: c("default", "between"). Both plots.
#' @param filter.rfp numeric, default 1. minimum fpkm value to be included in plots
#' @param filter.rna numeric, default 1. minimum fpkm value to be included in plots
#' @inheritParams countTable
#' @return a data.table with TE values, fpkm and log fpkm values, library
#' samples melted into rows with split variable called "variable".
#' @importFrom utils combn
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
                    collapse = FALSE) {
  RNA_MRNA_FPKM <- countTable(df.rna, "mrna", type = "fpkm", collapse = collapse)
  RNA_MRNA_FPKM <- data.table(id = rownames(RNA_MRNA_FPKM), RNA_MRNA_FPKM)

  RFP_CDS_FPKM <- countTable(df.rfp, "cds", type = "fpkm", collapse = collapse)
  RFP_CDS_FPKM <- data.table(id = rownames(RFP_CDS_FPKM), RFP_CDS_FPKM)

  if (!identical(nrow(RNA_MRNA_FPKM), nrow(RFP_CDS_FPKM)))
    stop("Not equal rows in count tables, did you match rfp and rna from different genome builds?")
  if (!identical(nrow(bamVarName(df.rfp, skip.libtype = TRUE)),
                 nrow(bamVarName(df.rna, skip.libtype = TRUE))))
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
  dt.final <- cbind(dt.melt.rna, dt.melt.rfp$value, dt.melt.rna.10$value)
  dt.final[, variable := gsub("RNA_", "", variable)]
  colnames(dt.final) <- c("variable", "rna_log2", "rfp_log2", "rna_log10")
  dt.final[, LFC_TE := rfp_log2 - rna_log2]

  message(paste("Filter kept", round((nrow(dt) / length(txNames)) *100, 1), "% of data"))


  subtitle <- paste("Filter: RFP >", filter.rfp, " & mRNA >", filter.rna, "(FPKM)")
  if ("default" %in% type) {
    plot <- ggplot(data = dt.final) +
      geom_point(aes(x = rna_log10, y = LFC_TE), alpha = 0.2) +
      theme_minimal() +
      geom_hline(aes(yintercept =  0), alpha = 0.2, color = "red") +
      xlab("mRNA FPKM (log10)") +
      ylab("TE (log2 fold change)") +
      ggtitle(label = "", subtitle = subtitle) +
      xlim(c(filter.rna, filter.rna + 2.5)) +
      facet_wrap(~ variable, ncol = 1)

    plot(plot)
    ggsave(file.path(output.dir, "TE_within.png"), plot,
           width = 7, height = 3+(ncol(RNA_MRNA_FPKM)-2), dpi = 300)
  }
  if ("between" %in% type) {
    pairs <- list() # creating compairisons :list of pairs
    my_comparison <- combn(unique(dt.final$variable),1)
    pairs <- list()
    for (i in 1:ncol(my_comparison)) {
      pairs[[i]] <- c(my_comparison[1,i], my_comparison[2,i])
    }
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
      ggtitle(label = "Comparison: 48h vs 24h", subtitle) +
      facet_wrap(~ variable, ncol = 2) +
      xlim(c(-5, 5))
    plot(plot)
    ggsave(file.path(output.dir, "TE_between.png"), plot,
           width = 4, height = 3 + (ncol(RNA_MRNA_FPKM)-2), dpi = 300)
  }
  return(dt)
}
