#' Run differential TE analysis
#'
#' Creates a total of 3 DESeq models (given x is design argument input
#'  (usually stage or condition) and libraryType is RNA-seq and Ribo-seq):\cr
#' 1. Ribo-seq model: design = ~ x (differences between the x groups in Ribo-seq)\cr
#' 2. RNA-seq model: design = ~ x (differences between the x groups in RNA-seq)\cr
#' 3. TE model: design = ~ x + libraryType + libraryType:x
#' (differences between the x and libraryType groups and the interaction between them)\cr
#' Using an equal reimplementation of the deltaTE algorithm (see reference).
#' You need at least 2 groups and 2 replicates per group. The Ribo-seq counts will
#' be over CDS and RNA-seq over mRNAs, per transcript. \cr
#'
#' The respective groups are defined as this
#' (given a user defined p value, shown here as 0.05):\cr
#' 1. Translation -  te.p.adj < 0.05 & rfp.p.adj < 0.05 & rna.p.adj > 0.05\cr
#' 2. mRNA abundance - te.p.adj > 0.05 & rfp.p.adj < 0.05 & rna.p.adj > 0.05\cr
#' 3. Buffering - te.p.adj < 0.05 & rfp.p.adj > 0.05 & rna.p.adj > 0.05\cr\cr
#'
#' Buffering also adds some close group which are split up if you set
#' complex.categories = TRUE (You will then get in addition)
#' See Figure 1 in the reference article for a clear definition of the groups!\cr
#' If you do not need isoform variants, subset to longest isoform per gene
#' either before or in the returned object (See examples).
#' If you do not have RNA-seq controls, you can still use DESeq on Ribo-seq alone.
#' \cr The LFC values are shrunken by lfcShrink(type = "normal").\cr \cr
#' Remember that DESeq by default can not
#' do global change analysis, it can only find subsets with changes in LFC!
#' @inheritParams DTEG.plot
#' @param df.rfp a \code{\link{experiment}} of Ribo-seq or 80S from TCP-seq.
#' @param df.rna a \code{\link{experiment}} of RNA-seq
#' @param design a character vector, default "stage". The columns in the
#' ORFik experiment that represent the comparison contrasts. Usually found
#' in "stage", "condition" or "fraction" column.
#' @param output.dir output.dir directory to save plots,
#' plot will be named "TE_between.png". If NULL, will not save.
#' @param RFP_counts a SummarizedExperiment, default:
#' countTable(df.rfp, "cds", type = "summarized"), all transcripts.
#' Assign a subset if you don't want to analyze all genes.
#' It is recommended to not subset, to give DESeq2 data for variance analysis.
#' @param RNA_counts a SummarizedExperiment, default:
#' countTable(df.rna, "mrna", type = "summarized"), all transcripts.
#' Assign a subset if you don't want to analyze all genes.
#' It is recommended to not subset, to give DESeq2 data for variance analysis.
#' @param batch.effect, logical, default FALSE. If you believe you might have batch effects,
#' set to TRUE, will use replicate column to represent batch effects.
#' Batch effect usually means that you have a strong variance between
#' biological replicates. Check PCA plot on count tables to verify if
#' you need to set it to TRUE.
#' @param complex.categories logical, default FALSE. Seperate into more groups,
#' will add Inverse (opposite diagonal of mRNA abundance) and Expression (only significant mRNA-seq)
#' @references doi: 10.1002/cpmb.108
#' @return a data.table with 9 columns.
#' (log fold changes, p.ajust values, group, regulation status and gene id)
#' @family TE
#' @export
#' @import DESeq2
#' @importFrom data.table rbindlist
#' @examples
#' ## Simple example
#' #df.rfp <- read.experiment("Riboseq")
#' #df.rna <- read.experiment("RNAseq")
#' #dt <- DTEG.analysis(df.rfp, df.rna)
#' ## Restrict DTEGs by log fold change (LFC):
#' ## subset to abs(LFC) < 1.5 for both rfp and rna
#' #dt[abs(rfp) < 1.5 & abs(rna) < 1.5, Regulation := "No change"]
#'
#' ## Only longest isoform per gene:
#' #tx_longest <- filterTranscripts(df.rfp, 0, 1, 0)
#' #dt <- dt[id %in% tx_longest,]
#' ## Convert to gene id
#' #dt[, id := txNamesToGeneNames(id, df.rfp)]
#' ## To get by gene symbol, use biomaRt conversion
DTEG.analysis <- function(df.rfp, df.rna,
                          output.dir = paste0(dirname(df.rfp$filepath[1]),
                                              "/QC_STATS/"),
                          design = "stage", p.value = 0.05,
                          RFP_counts = countTable(df.rfp, "cds", type = "summarized"),
                          RNA_counts = countTable(df.rna, "mrna", type = "summarized"),
                          batch.effect = FALSE,
                          plot.title = "", width = 6,
                          height = 6, dot.size = 0.4,
                          relative.name = "DTEG_plot.png", complex.categories = FALSE) {
  if (!is(df.rfp, "experiment") | !is(df.rna, "experiment"))
    stop("df.rfp and df.rna must be ORFik experiments!")
  if (length(unique(unlist(df.rfp[, design]))) == 1)
    stop("Design column needs at least 2 unique values!")
  if (nrow(df.rfp) < 4)
    stop("Experiment needs at least 4 rows, with minimum 2 per design group!")
  if (p.value > 1 | p.value <= 0)
    stop("p.value must be in interval (0,1]")
  if (!is(RFP_counts, "SummarizedExperiment") |
      !is(RNA_counts, "SummarizedExperiment")) stop("counts must be of type SummarizedExperiment")
  if (nrow(RFP_counts) != nrow(RNA_counts)) stop("counts must have equall number of rows!")

  # Designs
  be <- ifelse(batch.effect, "replicate + ", "")
  te.design <- as.formula(paste0("~ libtype + ", be, design, "+ libtype:", design))
  main.design <- as.formula(paste0("~ ", be, design))
  # TE
  se <- cbind(assay(RFP_counts), assay(RNA_counts))
  colData <- rbind(colData(RFP_counts), colData(RNA_counts))
  combined_se <- SummarizedExperiment(list(counts = se),
                                      rowRanges = rowRanges(RFP_counts),
                                      colData = colData)
  #combined_se <- combined_se[rowMeans(assay(combined_se)) > 1,]
  message("----------------------"); message("Model 1/3: TE")
  message("----------------------")
  combined_DESEQ <- DESeqDataSet(combined_se, design = te.design)
  dds.te <- DESeq2::DESeq(combined_DESEQ)
  resultsNames(dds.te)

  # Ribo
  message("----------------------"); message("Model 2/3: Ribo-seq")
  message("----------------------")
  ddsMat_ribo <- DESeqDataSet(se = RFP_counts, design = main.design)
  ddsMat_ribo <- DESeq(ddsMat_ribo)

  # RNA
  message("----------------------"); message("Model 3/3: RNA-seq")
  message("----------------------")
  ddsMat_rna <- DESeqDataSet(se = RNA_counts, design = main.design)
  ddsMat_rna <- DESeq(ddsMat_rna)
  message("----------------------")
  # Create contrasts
  pairs <- combn.pairs(unlist(df.rfp[, design]))

  # Per contrast
  dt.between <- data.table()
  #dt.between <- bplapply(pairs[1:2], FUN = function(i, dds.te, ddsMat_ribo, ddsMat_rna, design) {
  for(i in pairs) {
    name <- paste("Comparison:", i[1], "vs", i[2])
    message(name)
    # Results
    current.contrast <- c(design, i[1], i[2])
    res_te <- results(dds.te, contrast = current.contrast)

    res_ribo <- results(ddsMat_ribo, contrast=current.contrast)
    suppressMessages(res_ribo <- lfcShrink(ddsMat_ribo, contrast=current.contrast,
                                           res=res_ribo, type = "normal"))

    res_rna <- results(ddsMat_rna, contrast = current.contrast)
    suppressMessages(res_rna <- lfcShrink(ddsMat_rna, contrast = current.contrast,
                                          res = res_rna, type = "normal"))

    # The differential regulation groupings (padj is padjusted)
    both <- which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj < p.value)
    ## The 4 classes of genes
    # Forwarded are non significant in TE, diagonal line
    forwarded <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj < p.value & res_rna$padj < p.value)]
    # These two are the X and Y axis
    exclusive.translation <- rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj > p.value)]
    exclusive.expression <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj > p.value & res_rna$padj < p.value)]
    # These are the remaining groups
    intensified <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] > 0)]] # Called mRNA abundance
    inverse <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] < 0)]]
    buffered <- c(rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj > p.value & res_rna$padj < p.value)])
    # buffered <- c(rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj > p.value & res_rna$padj < p.value)],
    #               buffered)

    n <- rownames(res_te)
    Regulation <- rep("No change", nrow(res_te))
    Regulation[n %in% forwarded] <- "Forwarded" # Old Buffering
    Regulation[n %in% buffered] <- "Buffering"
    Regulation[n %in% inverse] <- "Inverse"
    Regulation[n %in% exclusive.translation] <- "Translation"
    Regulation[n %in% exclusive.expression] <- "Expression"
    Regulation[n %in% intensified] <- "mRNA abundance"

    if (!complex.categories) {
      Regulation[n %in% exclusive.expression] <- "Buffering"
      Regulation[n %in% inverse] <- "Buffering"
      Regulation[n %in% forwarded] <- "Buffering"
    }
    print(table(Regulation))


    dt.between <-
      rbindlist(list(dt.between,
                     data.table(variable = name,
                                Regulation = Regulation,
                                id = rownames(ddsMat_rna),
                                rna = res_rna$log2FoldChange,
                                rfp = res_ribo$log2FoldChange,
                                te = res_te$log2FoldChange,
                                rna.padj = res_rna$padj,
                                rfp.padj = res_ribo$padj,
                                te.padj = res_te$padj
                     )))
  }
  dt.between[, Regulation :=
               factor(Regulation,
                      levels = c("No change", "Translation", "Buffering", "mRNA abundance", "Expression", "Forwarded", "Inverse"),
                      ordered = TRUE)]
  # Plot the result
  message("----------------------"); message("Plotting all results")
  plot <- DTEG.plot(dt.between, output.dir, p.value, plot.title, width, height,
                    dot.size, relative.name = relative.name)
  return(dt.between)
}

#' Create a TE table
#'
#' Creates a data.table with 6 columns, column names are:\cr
#' variable, rfp_log2, rna_log2, rna_log10, TE_log2, id
#' @inheritParams DTEG.analysis
#' @inheritParams countTable
#' @param filter.rfp numeric, default 1. What is the minimum fpkm value?
#' @param filter.rna numeric, default 1. What is the minimum fpkm value?
#' @return a data.table with 6 columns
#' @family TE
#' @export
#' @examples
#' #df.rfp <- read.experiment("Riboseq")
#' #df.rna <- read.experiment("RNAseq")
#' #te.table(df.rfp, df.rna)
te.table <- function(df.rfp, df.rna,
                     filter.rfp = 1, filter.rna = 1,
                     collapse = FALSE) {
  if (!is(df.rfp, "experiment") | !is(df.rna, "experiment"))
    stop("df.rfp and df.rna must be ORFik experiments!")

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
  filtered <- rowMin(as.matrix(dt[,-1])) > max(filter.rna, filter.rfp)
  txNames <- dt$id; dt$id <- NULL

  #dt <- dt[RNA_MRNA_FPKM > filter.rfp & RFP_CDS_FPKM > filter.rna, ]
  dt <- dt[filtered, ]
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
  dt.final[, TE_log2 := rfp_log2 - rna_log2]
  dt.final$id <- rep(txNames[filtered], length(unique(dt.final$variable)))

  message(paste("Filter kept", round((nrow(dt) / length(txNames)) *100, 1), "% of transcripts"))
  return(dt.final)
}
