#' Run differential TE analysis
#'
#' Expression analysis of 1 dimension, usually between conditions of RNA-seq.\cr
#' Using the standardized DESeq2 pipeline flow.\cr
#' Creates a DESeq model (given x is the target.contrast argument)
#'  (usually 'condition' column)\cr
#' 1. RNA-seq model: design = ~ x (differences between the x groups in RNA-seq)\cr
#'
#' #' Analysis is done between each possible
#' combination of levels in the target contrast If target contrast is the condition column,
#' with factor levels: WT, mut1 and mut2 with 3 replicates each. You get comparison
#' of WT vs mut1, WT vs mut2 and mut1 vs mut2. \cr
#' The respective result categories are defined as:
#' (given a user defined p value, shown here as 0.05):\cr
#' Significant -  p-value adjusted < 0.05 (p-value cutoff decided by 'p.value argument)\cr
#' \cr The LFC values are shrunken by lfcShrink(type = "normal").\cr \cr
#' Remember that DESeq by default can not
#' do global change analysis, it can only find subsets with changes in LFC!
#' @inherit DTEG.analysis
#' @param df an \code{\link{experiment}} of usually RNA-seq.
#' @param counts a SummarizedExperiment, default:
#' countTable(df, "mrna", type = "summarized"), all transcripts.
#' Assign a subset if you don't want to analyze all genes.
#' It is recommended to not subset, to give DESeq2 data for variance analysis.
#' @return a data.table with columns:
#' (contrast variable, gene id, regulation status, log fold changes, p.adjust values, mean counts)
#' @export
#' @examples
#' ## Simple example (use ORFik template, then use only RNA-seq)
#' df <- ORFik.template.experiment()
#' df.rna <- df[df$libtype == "RNA",]
#' design(df.rna) # The full experimental design
#' design(df.rna)[1] # Default target contrast
#' #dt <- DEG.analysis(df.rna)
DEG.analysis <- function(df, target.contrast = design[1],
                         design = ORFik::design(df), p.value = 0.05,
                         counts = countTable(df, "mrna", type = "summarized"),
                         batch.effect = TRUE,
                         pairs = combn.pairs(unlist(df[, target.contrast]))) {
  # Get DESeq model
  ddsMat_rna <- DEG_model(df, target.contrast, design, p.value,
                          counts, batch.effect)
  # Do result analysis: per contrast selected
  return(DEG_model_results(ddsMat_rna, target.contrast, pairs, p.value))
}

#' Run differential TE analysis
#'
#' Expression analysis of 2 dimensions, usually Ribo-seq vs RNA-seq.\cr
#' Using an equal reimplementation of the deltaTE algorithm (see reference).\cr
#' Creates a total of 3 DESeq models (given x is the target.contrast argument)
#'  (usually 'condition' column) and libraryType is RNA-seq and Ribo-seq):\cr
#' 1. Ribo-seq model: design = ~ x (differences between the x groups in Ribo-seq)\cr
#' 2. RNA-seq model: design = ~ x (differences between the x groups in RNA-seq)\cr
#' 3. TE model: design = ~ x + libraryType + libraryType:x
#' (differences between the x and libraryType groups and the interaction between them)\cr
#' You need at least 2 groups and 2 replicates per group. By default, the Ribo-seq counts will
#' be over CDS and RNA-seq counts over whole mRNAs, per transcript. \cr
#'
#' Log fold changes and p-values are created from a Walds test on the comparison contrast described bellow.
#' The RNA-seq and Ribo-seq LFC values are shrunken using DESeq2::lfcShrink(type = "normal"). Note
#' that the TE LFC values are not shrunken (as following specifications from deltaTE paper)\cr\cr
#'
#' Analysis is done between each possible
#' combination of levels in the target contrast If target contrast is condition column,
#' with factor levels: WT, mut1 and mut2 with 3 replicates each. You get comparison
#' of WT vs mut1, WT vs mut2 and mut1 vs mut2. \cr
#' The respective result categories are defined as:
#' (given a user defined p value, shown here as 0.05):\cr
#' 1. Translation -  te.p.adj < 0.05 & rfp.p.adj < 0.05 & rna.p.adj > 0.05\cr
#' 2. mRNA abundance - te.p.adj > 0.05 & rfp.p.adj < 0.05 & rna.p.adj > 0.05\cr
#' 3. Buffering - te.p.adj < 0.05 & rfp.p.adj > 0.05 & rna.p.adj > 0.05\cr\cr
#'
#' Buffering will be broken down into sub-categories if you set
#' complex.categories = TRUE
#' See Figure 1 in the reference article for a clear definition of the groups!\cr
#' If you do not need isoform variants, subset to longest isoform per gene
#' either before or in the returned object (See examples).
#' If you do not have RNA-seq controls, you can still use DESeq on Ribo-seq alone.
#' \cr The LFC values are shrunken by lfcShrink(type = "normal").\cr \cr
#' Remember that DESeq by default can not
#' do global change analysis, it can only find subsets with changes in LFC!
#' @inheritParams DTEG.plot
#' @param df.rfp a \code{\link{experiment}} of usually Ribo-seq or 80S from TCP-seq.
#' (the numerator of the experiment, usually having a primary role)
#' @param df.rna a \code{\link{experiment}} of usually RNA-seq.
#' (the denominator of the experiment, usually having a normalizing function)
#' @param target.contrast a character vector, default \code{design[1]}.
#' The column in the ORFik experiment that represent the comparison contrasts.
#' By default: the first design factor of the full experimental design.
#' This is the factor you will do the comparison on. DESeq will normalize
#' the counts based on the full design, but the log fold change values will
#' be based on this contrast only. It is usually the 'condition' column.
#' @param design a character vector, default \code{design(df.rfp)}.
#' The full experiment design. Which factors have more than 1 level.
#' Example: stage column are all HEK293, so it can not be a design factor.
#' The condition column has 2 possible values, WT and mutant, so it is
#' a factor of the experiment. Replicates column is not part of design,
#' that is inserted later with setting \code{batch.effect = TRUE}.
#' Library type 'libtype' column, can also no be part of initial design,
#' it is always added inside the function, after initial setup.
#' @param p.value a numeric, default 0.05 in interval (0,1). Defines adjusted
#' p-value to be used as significance threshold for the result groups. I.e.
#' for exclusive translation group significant subset for p.value = 0.05 means:
#' TE$padj < 0.05 & Ribo$padj < 0.05 & RNA$padj > 0.05.
#' @param output.dir character, default \code{QCfolder(df.rfp)}.
#' output.dir directory to save plots,
#' plot will be named "TE_between". If NULL, will not save.
#' @param RFP_counts a \code{\link{SummarizedExperiment}}, default:
#' \code{countTable(df.rfp, "cds", type = "summarized")},
#' unshifted libraries, all transcript CDSs.
#' If you have pshifted reads and countTables, do:
#' \code{countTable(df.rfp, "cds", type = "summarized", count.folder = "pshifted")}
#' Assign a subset if you don't want to analyze all genes.
#' It is recommended to not subset, to give DESeq2 data for variance analysis.
#' @param RNA_counts a SummarizedExperiment, default:
#' countTable(df.rna, "mrna", type = "summarized"), all transcripts.
#' Assign a subset if you don't want to analyze all genes.
#' It is recommended to not subset, to give DESeq2 data for variance analysis.
#' @param batch.effect, logical, default TRUE. Makes replicate column of the experiment
#' part of the design.\cr
#' If you believe you might have batch effects, keep as TRUE.
#' Batch effect usually means that you have a strong variance between
#' biological replicates. Check out \code{\link{pcaExperiment}} and see if replicates
#' cluster together more than the design factor, to verify if you need to set it to TRUE.
#' @param pairs list of character pairs, the experiment contrasts. Default:
#'  \code{combn.pairs(unlist(df.rfp[, target.contrast])}
#' @param complex.categories logical, default FALSE. Seperate into more groups,
#' will add Inverse (opposite diagonal of mRNA abundance) and Expression (only significant mRNA-seq)
#' @references doi: 10.1002/cpmb.108
#' @return a data.table with columns:
#' (contrast variable, gene id, regulation status, log fold changes, p.adjust values, mean counts)
#' @family DifferentialExpression
#' @export
#' @import DESeq2
#' @importFrom data.table rbindlist
#' @examples
#' ## Simple example (use ORFik template, then split on Ribo and RNA)
#' df <- ORFik.template.experiment()
#' df.rfp <- df[df$libtype == "RFP",]
#' df.rna <- df[df$libtype == "RNA",]
#' design(df.rfp) # The experimental design, per libtype
#' design(df.rfp)[1] # Default target contrast
#' #dt <- DTEG.analysis(df.rfp, df.rna)
#' ## If you want to use the pshifted libs for analysis:
#' #dt <- DTEG.analysis(df.rfp, df.rna,
#' #                    RFP_counts = countTable(df.rfp, region = "cds",
#' #                       type = "summarized", count.folder = "pshifted"))
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
#' ## To flip directionality of contrast pair nr 2:
#' #design <- "condition"
#' #pairs <- combn.pairs(unlist(df.rfp[, design])
#' #pairs[[2]] <- rev(pars[[2]])
#' #dt <- DTEG.analysis(df.rfp, df.rna,
#' #                    RFP_counts = countTable(df.rfp, region = "cds",
#' #                       type = "summarized", count.folder = "pshifted"),
#' #                       pairs = pairs)
DTEG.analysis <- function(df.rfp, df.rna,
                          output.dir = QCfolder(df.rfp),
                          target.contrast = design[1],
                          design = ORFik::design(df.rfp), p.value = 0.05,
                          RFP_counts = countTable(df.rfp, "cds", type = "summarized"),
                          RNA_counts = countTable(df.rna, "mrna", type = "summarized"),
                          batch.effect = FALSE, pairs = combn.pairs(unlist(df.rfp[, design])),
                          plot.title = "", plot.ext = ".pdf", width = 6,
                          height = 6, dot.size = 0.4,
                          relative.name = paste0("DTEG_plot", plot.ext),
                          complex.categories = FALSE) {
  # Input validation
  DTEG_input_validation(df.rfp, df.rna, RFP_counts, RNA_counts, design,
                        target.contrast, p.value)

  # Designs
  be <- ifelse(batch.effect, "replicate + ", "")
  te.design <- as.formula(paste0("~ libtype + ", be,
                                 paste(design, collapse = " + "),
                                 "+ libtype:", target.contrast))
  main.design <- as.formula(paste0("~ ", be, paste(design, collapse = " + ")))
  message("----------------------")
  message("Full exper. design: ", main.design)
  message("Interaction design: ", te.design)
  message("Target -- contrast: ", target.contrast)
  message("----------------------")
  # TE count table (cbind of both)
  se <- cbind(assay(RFP_counts), assay(RNA_counts))
  colData <- rbind(colData(RFP_counts), colData(RNA_counts))
  combined_se <- SummarizedExperiment(list(counts = se),
                                      rowRanges = rowRanges(RFP_counts),
                                      colData = colData)
  ## DESeq models (total: 3)
  # TE
  message("----------------------")
  ddsMat_te <- DEG_DESeq(combined_se, te.design, "Model 1/3: TE")

  # Ribo
  ddsMat_ribo <- DEG_DESeq(RFP_counts, main.design, "Model 2/3: Ribo-seq")

  # RNA
  ddsMat_rna <- DEG_DESeq(RNA_counts, main.design, "Model 3/3: RNA-seq")

  # Do result analysis: per contrast selected
  dt <- DTEG_model_results(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                           target.contrast, pairs, p.value = p.value,
                           complex.categories)

  # Plot the result
  message("----------------------"); message("Plotting all results")
  plot <- DTEG.plot(dt, output.dir, p.value, plot.title, width, height,
                    dot.size, relative.name = relative.name)
  return(dt)
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
#' @family DifferentialExpression
#' @importFrom data.table setkey
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df.rfp <- df[df$libtype == "RFP",]
#' df.rna <- df[df$libtype == "RNA",]
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

  uid <- unique(c(RNA_MRNA_FPKM$id, RFP_CDS_FPKM$id))
  RNA_MRNA_FPKM$uid <- match(RNA_MRNA_FPKM$id, uid)
  RNA_MRNA_FPKM$id <- NULL
  RFP_CDS_FPKM$uid <- match(RFP_CDS_FPKM$id, uid)
  RFP_CDS_FPKM$id <- NULL
  setkey(RNA_MRNA_FPKM, uid)
  setkey(RFP_CDS_FPKM, uid)

  dt <- data.table::merge.data.table(RNA_MRNA_FPKM, RFP_CDS_FPKM, by = "uid")
  txNames <- uid[dt$uid]
  dt$uid <- NULL
  filtered <- rowMin(as.matrix(dt)) > max(filter.rna, filter.rfp)

  #dt <- dt[RNA_MRNA_FPKM > filter.rfp & RFP_CDS_FPKM > filter.rna, ]
  dt <- dt[filtered, ]
  dt.log <- pseudo.transform(dt)
  dt.log10 <- pseudo.transform(dt, log10)
  RNA_MRNA_FPKM$uid <- NULL
  RFP_CDS_FPKM$uid <- NULL
  dt.melt.rna <- melt(dt.log[, colnames(dt.log) %in% colnames(RNA_MRNA_FPKM), with = FALSE])
  dt.melt.rna.10 <- melt(dt.log10[, colnames(dt.log10) %in% colnames(RNA_MRNA_FPKM), with = FALSE])
  dt.melt.rfp <- melt(dt.log[, colnames(dt.log) %in% colnames(RFP_CDS_FPKM), with = FALSE])
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
