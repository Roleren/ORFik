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
#'  (usually 'condition' column) and libraryType is RNA-seq and Ribo-seq):\cr\cr
#' ** \strong{The 3 differential sub models} **
#' \itemize{
#'  \item{1. Ribo-seq model : design = ~ x (differences between the x groups in Ribo-seq)}
#'  \item{2. RNA-seq model: design = ~ x (differences between the x groups in RNA-seq)}
#'  \item{3. TE model: design = ~ x + libraryType + libraryType:x
#'  (differences between the x and libraryType groups and the interaction between them)}
#' }
#' You need at least 2 groups and 2 replicates per group. By default, the Ribo-seq counts will
#' be over CDS and RNA-seq counts over whole mRNAs, per transcript. See notes section below
#' for more details.\cr
#'
#' Log fold changes and p-values are created from a Walds test on the comparison contrast described bellow.
#' The RNA-seq and Ribo-seq LFC values are shrunken using DESeq2::lfcShrink(type = "normal"). Note
#' that the TE LFC values are not shrunken (as following specifications from deltaTE paper).
#' The adjusted p-values are created using DESEQ pAdjustMethod = "BH" (Benjamini-Hochberg correction).
#' All other DESEQ2 arguments are default.\cr\cr
#'
#' Analysis is done between each possible
#' combination of levels in the target contrast If target contrast is condition column,
#' with factor levels: WT, mut1 and mut2 with 3 replicates each. You get comparison
#' of WT vs mut1, WT vs mut2 and mut1 vs mut2. \cr
#'
#' The respective result categories are defined through 4 main categories,
#' first some intuition.
#' The number of ribosomes (Ribo-seq) is significantly different between 2
#' contrast elements in the model if the relative counts is
#' statistically higher/lower, for mRNA levels (RNA-seq) it is the same.
#' So TE is then RFP / RNA which is basically how many ribosomes translated
#' per mRNA in the sample, if contrast group 1 has TE of 2, it means 2 ribosomes
#' per mrna fragment, while TE of 4 would be a doubling of 4 ribosomes per mRNA.
#'
#' Mathematically the groups are defined by the p adjusted values as the
#' following (te.sign means na_safe(te.padj < p.value),
#' na_safe is a function where NA values are FALSE for '<=' test
#' and TRUE for '>' test),
#' we also use a helper function:
#' te.sign & rfp.sign & rna.sign, all_models_sign := TRUE.\cr\cr
#' ** \strong{Signicant DTEG Classifications} **
#' \itemize{
#'  \item{No change : None of the below categories}
#'  \item{Translation (only RFP) : te.sign & rfp.sign & !rna.sign}
#'  \item{Expression (only RNA) : !te.sign & !rfp.sign & rna.sign}
#'  \item{mRNA abundance : all_models_sign & na_safe(te.lfc * rna.lfc, ">", 0)}
#'  \item{Inverse (inverse mRNA abundance) : all_models_sign & te.lfc * rna.lfc, "<", 0)}
#'  \item{Buffering (Stable protein output) : te.sign & !rfp.sign & rna.sign}
#'  \item{Forwarded (diagonal bottom left to top right) : !te.sign & rfp.sign & rna.sign}
#' }
#'
#' If complex.categories is FALSE, then Expression, Inverse and forwarded are defined 'Buffering'.
#' mRNA abundance is called"Intensified" in original article
#' For code, of classification, run: View(ORFik:::DTEG_add_regulation_categories).
#' Feel free to redefine the categories as you want them.
#'
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
#' @param complex.categories logical, default FALSE. Separate into more groups,
#' see above for details.
#' @references \url{https://doi.org/10.1002/cpmb.108}
#' @return a data.table with columns:
#' (contrast variable, gene id, regulation status, log fold changes, p.adjust values,
#' mean counts, significant (as logical))
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
#' #dt_with_gene_ids <- append_gene_symbols(dt, symbols(df))
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
                          complex.categories = FALSE,
                          plot_to_console = TRUE) {
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
                    dot.size, relative.name = relative.name,
                    plot_to_console = plot_to_console)
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

#' GO analysis with GOrilla
#'
#' Supports Gene symbols as default, and will produce the best
#' results.
#' You can also use ensembl gene ids, refseq gene ids and Entrez gene
#' ids, but this will give weaker results.
#' @param target_genes a path to a txt file with the target Gene symbols,
#' or presumed to be a character vector of genes (if length > 1).
#' Minimum 10 genes, maximum 1 million.
#' @param background_genes a path to a txt file with the background Gene symbols,
#' or presumed to be a character vector of genes (if length > 1).
#' Minimum 10 genes, maximum 2 million.
#' @param organism organism(df), example "Homo sapiens"
#' @param analysis_name character name, default "test".
#' Used for saved file names and analysis name in GOrilla.
#' @param open_browser = TRUE, open the URL
#' @param pvalue_thresh fixed set numeric, default 0.001, Alternatives: 10e-3 to 10e-11
#' @param db character, default "all". Which GO onthology categories to use,
#' all means process, function and component.
#' Alternatives: "proc", "func" and "comp", if you only want that single category subset.
#' @references https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-48
#' @return a url path to results, will also open your default web browser if open_browser
#' is TRUE.
#' @importFrom httr POST
#' @export
#' @examples
#' target_genes <- system.file("extdata/Homo_sapiens_sample/QC_STATS",
#'   "/DTEG_Comparison_Translation.txt", package = "ORFik")
#' background_genes <- system.file("extdata/Homo_sapiens_sample/QC_STATS",
#'   "/DTEG_Comparison_Background.txt", package = "ORFik")
#' #go_analaysis_gorilla(target_genes, background_genes, "Homo sapiens",
#' #                     analysis_name = "Translation vs background")
go_analaysis_gorilla <- function(target_genes, background_genes,
                                 organism,
                                 analysis_name = paste0("Go_analysis_", organism),
                                 open_browser = TRUE,
                                 pvalue_thresh = 10e-3,
                                 db = "all") {

  gorilla_species <- c(
    "Arabidopsis thaliana"       = "ARABIDOPSIS_THALIANA",
    "Saccharomyces cerevisiae"   = "SACCHAROMYCES_CEREVISIAE",
    "Caenorhabditis elegans"     = "CAENORHABDITIS_ELEGANS",
    "Drosophila melanogaster"    = "DROSOPHILA_MELANOGASTER",
    "Danio rerio (Zebrafish)"    = "DANIO_RERIO",
    "Homo sapiens"               = "HOMO_SAPIENS",
    "Mus musculus"               = "MUS_MUSCULUS",
    "Rattus norvegicus"          = "RATTUS_NORVEGICUS"
  )
  organism <- toupper(gsub("_x_.*", "", gsub(" |-", "_", organism)))
  if (!(organism %in% gorilla_species)) {
    print(gorilla_species)
    stop("Invalid Gorilla organism, see above for valid options")
  }

  file_input <- length(target_genes) == 1 & length(background_genes) == 1
  if (file_input) {
    stopifnot(file.exists(target_genes))
    stopifnot(file.exists(background_genes))
    target_genes <- read.table(target_genes)[[1]]
    background_genes <- read.table(background_genes)[[1]]
  }

  stopifnot(is.character(target_genes) & is.character(background_genes))
  stopifnot(length(target_genes) > 10 & length(background_genes) > 10)
  stopifnot(length(target_genes) < 1e6 & length(background_genes) < 2e6)
  stopifnot(is.logical(open_browser))
  message("- Target genes: ", length(target_genes))
  message("- Background genes: ", length(background_genes))
  message("-- Sending GO analysis request..")
  res <- httr::POST(
    url = "https://cbl-gorilla.cs.technion.ac.il/servlet/GOrilla",
    body = list(
      application = "gorilla",
      species = organism,
      run_mode = "hg",  # 'hg' = two unranked lists
      target_set = paste(target_genes, collapse = "\n"),
      background_set = paste(background_genes, collapse = "\n"),
      db = db,  # or proc/func/comp
      pvalue_thresh = as.character(pvalue_thresh),
      analysis_name = analysis_name,
      output_excel = "on",
      output_revigo = "on"
    ),
    encode = "form"
  )
  if (res$status_code != 200) {
    print(res)
    stop("Gorilla did not resond with OK (200), see above for more info")
  }

  gorilla_url <- res$url
  message("- Complete")
  if (open_browser) browseURL(gorilla_url)
  stable_url <- paste0("https://cbl-gorilla.cs.technion.ac.il/GOrilla/", sub(".*\\?id=", "", gorilla_url), "/GOResults.html")
  attr(gorilla_url, "stable_url") <- stable_url
  return(gorilla_url)
}
