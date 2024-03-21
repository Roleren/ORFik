DEG_input_validation <- function(df, counts, design, target.contrast, p.value,
                                 input_check_name = "") {
  message("----------------------")
  message("Input check: ", input_check_name)
  if (!is(df, "experiment"))
    stop("df must be ORFik experiments!")
  if (!is(design, "character") | length(design) == 0) {
    stop("Design must be character of length > 0, don't use formula as input here")
  }
  if (!is(target.contrast, "character") | (length(target.contrast) != 1)){
    stop("target.contrast must be character of length 1, don't use formula as input here")
  }
  if (!(target.contrast %in% design))
    stop("target.contrast must be a column in the 'design' argument")
  if ("rep" %in% design) stop("Do not use rep in design, use batch.effect argument")
  if (length(unique(unlist(df[, target.contrast]))) == 1)
    stop("target.contrast column needs at least 2 unique values, to be able to compare anything!")
  if (nrow(df) < 4)
    stop("Experiment needs at least 4 rows, with minimum 2 per design group!")
  if (p.value > 1 | p.value <= 0)
    stop("p.value must be in interval (0,1]")
  if (!is(counts, "SummarizedExperiment")) stop("counts must be of type SummarizedExperiment")
  return(invisible(NULL))
}

DEG_design <- function(design, target.contrast, batch.effect) {
  be <- ifelse(batch.effect, "replicate + ", "")
  main.design <- as.formula(paste("~", be, paste(design, collapse = " + ")))
  message("----------------------")
  message("Full exp design: ", main.design)
  message("Target contrast: ", target.contrast)
  message("----------------------")
  return(main.design)
}

DEG_DESeq <- function(counts, main.design, message = "Creating DESeq model:") {
  message(message)
  message("----------------------")
  ddsMat <- DESeqDataSet(se = counts, design = main.design)
  ddsMat <- DESeq(ddsMat)
  message("----------------------")
  return(ddsMat)
}

#' Get DESeq2 model without running results
#'
#' This is the preparation step of DESeq2 analysis using ORFik::DEG.analysis.
#' It is exported so that you can do this step in standalone, usually you
#' want to use DEG.analysis directly.
#' @inheritParams DEG.analysis
#' @return a DESeqDataSet object with results stored as metadata columns.
#' @family DifferentialExpression
#' @export
#' @examples
#' ## Simple example (use ORFik template, then use only RNA-seq)
#' df <- ORFik.template.experiment()
#' df.rna <- df[df$libtype == "RNA",]
#' design(df.rna) # The full experimental design
#' target.contrast <- design(df.rna)[1] # Default target contrast
#' #ddsMat_rna <- DEG_model(df.rna, target.contrast)
DEG_model <- function(df,
                      target.contrast = design[1],
                      design = ORFik::design(df),
                      p.value = 0.05,
                      counts = countTable(df, "mrna", type = "summarized"),
                      batch.effect = TRUE) {
  # Input validation
  DEG_input_validation(df, counts, design, target.contrast, p.value)

  # Design
  main.design <- DEG_design(design, target.contrast, batch.effect)

  # DESeq2 model
  ddsMat_rna <- DEG_DESeq(counts, main.design)
}

#' Get DESeq2 model results from DESeqDataSet
#'
#'
#' @param ddsMat_rna a DESeqDataSet object with results stored as metadata columns.
#' @inheritParams DEG.analysis
#' @return a data.table
#' @export
#' @examples
#' ## Simple example (use ORFik template, then use only RNA-seq)
#' df <- ORFik.template.experiment()
#' df.rna <- df[df$libtype == "RNA",]
#' design(df.rna) # The full experimental design
#' target.contrast <- design(df.rna)[1] # Default target contrast
#' #ddsMat_rna <- DEG_model(df.rna, target.contrast)
#' #pairs <- combn.pairs(unlist(df[, target.contrast]))
#' #dt <- DEG_model_results(ddsMat_rna, target.contrast, pairs)
DEG_model_results <- function(ddsMat_rna, target.contrast, pairs,
                              p.value = 0.05) {
  # Do result analysis: per contrast selected
  dt.between <- data.table()
  for(i in pairs) {
    name <- paste("Comparison:", i[1], "vs", i[2])
    message(name)
    # Results
    current.contrast <- c(target.contrast, i[1], i[2])
    res_rna <- results(ddsMat_rna, contrast = current.contrast)
    suppressMessages(res_rna <- lfcShrink(ddsMat_rna, contrast = current.contrast,
                                          res = res_rna, type = "normal"))

    # The differential regulation groupings (padj is padjusted)
    expressed <- which(res_rna$padj < p.value)
    n <- rownames(res_rna)
    Regulation <- rep("No change", nrow(res_rna))
    Regulation[expressed] <- "Significant" # Old Buffering
    print(table(Regulation))

    dt.between <-
      rbindlist(list(dt.between,
                     data.table(contrast = name,
                                Regulation = Regulation,
                                id = rownames(ddsMat_rna),
                                meanCounts = res_rna$baseMean,
                                LFC = res_rna$log2FoldChange,
                                padj = res_rna$padj
                     )))
  }
  dt.between[, Regulation :=
               factor(Regulation,
                      levels = c("No change", "Significant"),
                      ordered = TRUE)]
  message("----------------------")
  return(dt.between)
}

#' Simple Fpkm ratio test DEG
#'
#' If you do not have a valid DESEQ2 experimental setup (contrast), you
#' can use this simplified test
#' @inheritParams DEG.analysis
#' @return a data.table of fpkm ratios
#' @export
#' @examples
#' ## Simple example (use ORFik template, then use only RNA-seq)
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RNA",]
#' #dt <- DEG_model_simple(df)
DEG_model_simple <- function(df,
                             target.contrast = design[1],
                             design = ORFik::design(df),
                             p.value = 0.05,
                             counts = countTable(df, "mrna", type = "summarized"),
                             batch.effect = FALSE) {
  # Design
  main.design <- DEG_design(design[1], target.contrast, batch.effect)
  design_ids <- as.character(df[, target.contrast])
  ids <- rownames(counts)
  counts <- as.data.table(DESeq2::fpkm(DESeq2::DESeqDataSet(counts, ~ 1)))
  # Get LFC matrix
  counts_group <- lapply(unique(design_ids),
                         function(x) rowMeans(counts[,design_ids == x, with = FALSE]))
  setDT(counts_group)
  colnames(counts_group) <- unique(design_ids)
  dt <- data.table(LFC = log2((counts_group[,1][[1]] + 0.0001) / (counts_group[,2][[1]] + 0.0001)),
                   meanCounts = c(rowMeans(counts_group[, c(1,2), with = FALSE])))
  dt[, id := ids]
  dt[, contrast := paste(unique(design_ids)[1:2], collapse = "_vs_")]
  return(dt)
}

DTEG_input_validation <- function(df.rfp, df.rna, RFP_counts, RNA_counts,
                                  design, target.contrast, p.value) {
  message("----------------------")
  message("Input check:")
  if (!is(df.rfp, "experiment") | !is(df.rna, "experiment"))
    stop("df.rfp and df.rna must be ORFik experiments!")

  DEG_input_validation(df.rfp, RFP_counts, design, target.contrast, p.value,
                       "df.rfp")
  DEG_input_validation(df.rna, RNA_counts, design, target.contrast, p.value,
                       "df.rna")

  if (nrow(RFP_counts) != nrow(RNA_counts))
    stop("Count tables must have equall number of rows!")
  message("----------------------")
  return(invisible(NULL))
}

DTEG_model_results <- function(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                               target.contrast, pairs, p.value = 0.05,
                               complex.categories) {
  # Do result analysis: per contrast selected
  dt.between <- data.table()
  for(i in pairs) {
    name <- paste("Comparison:", i[1], "vs", i[2])
    message(name)
    # Results
    current.contrast <- c(target.contrast, i[1], i[2])
    res_te <- results(ddsMat_te, contrast = current.contrast)

    res_ribo <- results(ddsMat_ribo, contrast = current.contrast)
    suppressMessages(res_ribo <- lfcShrink(ddsMat_ribo, contrast=current.contrast,
                                           res=res_ribo, type = "normal"))

    res_rna <- results(ddsMat_rna, contrast = current.contrast)
    suppressMessages(res_rna <- lfcShrink(ddsMat_rna, contrast = current.contrast,
                                          res = res_rna, type = "normal"))

    ## The differential regulation groupings (padj is padjusted)
    both <- which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj < p.value)

    ## The 4 classes of genes
    # Forwarded are non significant in TE, diagonal line
    forwarded <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj < p.value & res_rna$padj < p.value)]
    # These two are the X and Y axis
    exclusive.translation <- rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj > p.value)]
    exclusive.expression <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj > p.value & res_rna$padj < p.value)]
    ## These are the remaining groups
    # Also called mRNA abundance
    intensified <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] > 0)]]
    # Also called inverse mRNA abundance
    inverse <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] < 0)]]
    # Stable protein output
    buffered <- c(rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj > p.value & res_rna$padj < p.value)])

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
                     data.table(contrast = name,
                                Regulation = Regulation,
                                id = rownames(ddsMat_rna),
                                rna = res_rna$log2FoldChange,
                                rfp = res_ribo$log2FoldChange,
                                te = res_te$log2FoldChange,
                                rna.padj = res_rna$padj,
                                rfp.padj = res_ribo$padj,
                                te.padj = res_te$padj,
                                rna.meanCounts = res_rna$baseMean,
                                rfp.meanCounts = res_ribo$baseMean
                     )))
  }
  regulation_levels <-  c("No change", "Translation", "Buffering",
                          "mRNA abundance", "Expression", "Forwarded",
                          "Inverse")
  dt.between[, Regulation :=
               factor(Regulation, levels = regulation_levels, ordered = TRUE)]
  return(dt.between)
}
