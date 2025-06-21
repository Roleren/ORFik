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
#' @param design_ids character vector of contrast group ids, e.g.
#'  c("WT", "WT", "Mutant", "Mutant")
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
                             design_ids = as.character(df[, target.contrast]),
                             batch.effect = FALSE) {
  stopifnot(ncol(counts) == length(design_ids))
  stopifnot(length(unique(design_ids)) >= 2)
  # Design
  main.design <- DEG_design(design[1], target.contrast, batch.effect)

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
  dt[]
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
                               complex.categories, lfcShrinkType = "normal") {
  # Do result analysis: per contrast selected
  dt_all_pairs <- data.table()
  for(pair in pairs) {
    contrast_pair <- c(target.contrast, pair[1], pair[2])
    dt <- DTEG_pair_results(ddsMat_te, ddsMat_ribo, ddsMat_rna, contrast_pair, lfcShrinkType, p.value)
    dt <- DTEG_add_regulation_categories(dt, complex.categories)

    print(table(dt$Regulation))
    dt_all_pairs <- rbindlist(list(dt_all_pairs, dt))
  }

  attr(dt_all_pairs, "p.value") <- p.value
  return(dt_all_pairs)
}

DTEG_pair_results <- function(ddsMat_te, ddsMat_ribo, ddsMat_rna, contrast_vec,
                              lfcShrinkType, p.value, return_type = "data.table") {
  stopifnot(length(contrast_vec) == 3)
  name <- paste("Comparison:", contrast_vec[2], "vs", contrast_vec[3])
  message(name)

  res_te <- results(ddsMat_te, contrast = contrast_vec)

  res_ribo <- results(ddsMat_ribo, contrast = contrast_vec)
  suppressMessages(res_ribo <- lfcShrink(ddsMat_ribo, contrast=contrast_vec,
                                         res=res_ribo, type = lfcShrinkType))

  res_rna <- results(ddsMat_rna, contrast = contrast_vec)
  suppressMessages(res_rna <- lfcShrink(ddsMat_rna, contrast = contrast_vec,
                                        res = res_rna, type = lfcShrinkType))

  if (return_type != "data.table") {
    return(list(res_te, res_ribo, res_rna))
  }
  dt <- data.table(contrast = name,
                   Regulation = rep("No change", nrow(res_te)),
                   id = rownames(ddsMat_rna),
                   rna.lfc = res_rna$log2FoldChange,
                   rfp.lfc = res_ribo$log2FoldChange,
                   te.lfc = res_te$log2FoldChange,
                   rna.padj = res_rna$padj,
                   rfp.padj = res_ribo$padj,
                   te.padj = res_te$padj,
                   rna.meanCounts = res_rna$baseMean,
                   rfp.meanCounts = res_ribo$baseMean
  )

  dt[, rna.sign := ifelse(na_safe(rna.padj, "<=", p.value), TRUE, FALSE)]
  dt[, rfp.sign := ifelse(na_safe(rfp.padj, "<=", p.value), TRUE, FALSE)]
  dt[, te.sign := ifelse(na_safe(te.padj, "<=", p.value), TRUE, FALSE)]
  dt[, all_models_sign := FALSE]
  dt[te.sign & rfp.sign & rna.sign, all_models_sign := TRUE]
  return(dt)
}

na_safe <- function(x, op, threshold) {
  op_func <- match.fun(op)  # Convert string to function

  if (op %in% c("<", "<=")) {
    # For "<" and "<=", NA should be FALSE (assume p.adj = 1)
    !is.na(x) & op_func(x, threshold)
  } else {
    # For ">", ">=", "==", NA should be TRUE (assume non-significance)
    is.na(x) | op_func(x, threshold)
  }
}

#' Append gene symbols to a data.table with tx ids
#'
#' Main use case is to add gene symbols to data.table outputs from ORFik
#' with tx ids only, like the DTEG.analysis etc.
#' @param dt a data.table, must have a id_col with transcript ids
#' @param symbols_dt the data.table with symbols, must have a column
#'  with tx, transcript or value in the name. And only 1 of those!
#' @param extend_id logical, if TRUE, paste together old id from dt,
#' with the symbol id like: tx_id(symbol_id)
#' @param id_col character, default "id". The name of the id column in dt.
#' @return a data.table
#' @export
#' @examples
#' library(data.table)
#' df <- ORFik.template.experiment()
#'
#' cds_names <- names(loadRegion(df, "cds"))
#' dt <- data.table(id = cds_names[-1], LFC = seq(5), p.value = 0.05)
#'
#' symbols_dt <- data.table(ensembl_tx_name = cds_names,
#'  ensembl_gene_id = txNamesToGeneNames(cds_names, df),
#'  external_gene_name = c("ATF4", "AAT1", "ML4", "AST2", "RPL4", "RPL12"))
#' append_gene_symbols(dt, symbols_dt)
#' append_gene_symbols(dt, symbols_dt, extend_id = FALSE)
append_gene_symbols <- function(dt, symbols_dt, extend_id = TRUE,
                                id_col = "id") {
  stopifnot(is(dt, "data.table"))
  stopifnot(!is.null(dt$id))

  dt_with_symbols <- copy(dt)
  if (length(symbols_dt) > 0 & nrow(symbols_dt) > 0) {
    tx_column <- grep("tx|transcript|value", colnames(symbols_dt), ignore.case = TRUE, value = TRUE)
    if (length(tx_column) != 1) {
      if (length(tx_column) == 0) {
        warning("Could not find any column of symbols table with tx, transcript or value in column name, please rename/add")
      } else{
        warning("Found multiple columns of symbols table with tx, transcript or value in column name, please rename/remove")
      }
    }

    dt_with_symbols <- data.table::merge.data.table(dt, symbols_dt, by.x = id_col, by.y = tx_column, sort = FALSE)
    if (extend_id) {
      dt_with_symbols[, id_original := id]
      if (!is.null(dt_with_symbols$external_gene_name)) {
        dt_with_symbols[, id := paste0(id, "(", external_gene_name, ")")]
      } else if (!is.null(dt_with_symbols$label)) {
        dt_with_symbols[, id := paste0(id, "(", sub("-.*", "", label), ")")]
      } else warning("Could not find any column of symbols table with external_gene_name or label in column name, please rename/add")
    }
  }
  dt_with_symbols[]
  return(dt_with_symbols)
}

#' Add regulation categories
#' @noRd
DTEG_add_regulation_categories <- function(dt, complex.categories) {
  # Forwarded are non significant in TE, diagonal line
  dt[!te.sign & rfp.sign & rna.sign, Regulation := "Forwarded"]
  # These two are the X and Y axis
  dt[te.sign & rfp.sign & !rna.sign, Regulation := "Translation"]
  dt[!te.sign & !rfp.sign & rna.sign, Regulation := "Expression"]
  ## These are the remaining groups

  dt[all_models_sign & na_safe(te.lfc * rna.lfc, ">", 0), Regulation := "Intensified"]
  dt[Regulation == "Intensified", Regulation := "mRNA abundance"] #  Rename
  # Also called inverse mRNA abundance
  dt[all_models_sign & na_safe(te.lfc * rna.lfc, "<", 0), Regulation := "Inverse"]
  # Stable protein output
  dt[te.sign & !rfp.sign & rna.sign, Regulation := "Buffering"]


  if (!complex.categories) {
    dt[Regulation == "Expression", Regulation := "Buffering"]
    dt[Regulation == "Inverse", Regulation := "Buffering"]
    dt[Regulation == "Forwarded", Regulation := "Buffering"]
  }
  dt[, all_models_sign := NULL]
  regulation_levels <-  c("No change", "Translation", "Buffering", "mRNA abundance",
                          if(complex.categories) c("Expression", "Forwarded","Inverse") else NULL)

  stopifnot(all(dt$Regulation %in% regulation_levels))
  dt[, Regulation :=
               factor(Regulation, levels = regulation_levels, ordered = TRUE)]
  return(dt)
}

