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

  if (!(db %in% c("all", "proc", "func", "comp"))) {
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

#' GO analysis with GOrilla
#'
#' Per contrast, split list into regulation groups and run GOrilla go analysis.
#' If you want to change LFC threshold or p-value cutoff, you reassign
#' Regulation column by LFC and/or p-value.
#' @param dt a data.table of DEG or DTEG results, must also have appended
#' gene symbols column called "external_gene_name"
#' @param output_dir path to save results
#' @inheritParams go_analaysis_gorilla
#' @return a data.table with urls per contrast, this is also saved in
#' output_dir
#' @export
DEG_gorilla <- function(dt, output_dir, organism) {
  stopifnot(is(dt, "data.table"))
  stopifnot(!is.null(dt$external_gene_name))
  stopifnot(nrow(dt) > 20)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  contrast_sel <- unique(dt$contrast)[1]
  all_gorilla_ids <- data.table()
  for (contrast_sel in unique(dt$contrast)) {
    message("- ", contrast_sel)
    res_contrast <- dt[contrast == contrast_sel]

    all_regulations <- unique(res_contrast$Regulation)
    for (reg in all_regulations) {
      message("-- ", reg)
      contrast_basename <- paste0("DTEG_", contrast_sel, "_",
                                  ifelse(reg == "No change", "Background", reg),".txt")

      filename <- file.path(output_dir, contrast_basename)
      fwrite(data.table(unique(res_contrast[Regulation == reg, external_gene_name])), file = filename,
             col.names = FALSE, sep = "\n")
      # if (reg != "No change") {
      #   filename_up <- sub("\\.txt", "_up.txt", filename)
      #   fwrite(data.table(unique(res_contrast[Regulation == reg, external_gene_name])), file = filename_up,
      #          col.names = FALSE, sep = "\n")
      # }
    }
    background_file <- file.path(output_dir,
                                 paste0("DTEG_", contrast_sel, "_", "Background", ".txt"))
    if (file.exists(background_file)) {
      message("- GOrilla analysis..")
      all_significant_regulations <- all_regulations[all_regulations != "No change"]
      reg <- all_significant_regulations[1]
      urls <- lapply(all_significant_regulations, function(reg) {
        message("-- ", reg)
        target_file <- file.path(output_dir,
                                 paste0("DTEG_", contrast_sel, "_",
                                        ifelse(reg == "No change", "Background", as.character(reg)),".txt"))
        name <- gsub("_", " ", paste(sub("\\.txt$", "", basename(target_file)), "vs Background"))
        url <- ORFik:::go_analaysis_gorilla(target_file, background_file, organism,
                                            analysis_name = name, open_browser = FALSE)
        Sys.sleep(3)
        names(url) <- name
        return(url)
      })
      urls_stable <- sapply(urls, function(url) attr(url, "stable_url"))
      urls_unstable <- unlist(urls, use.names = TRUE)
      gorilla_ids <- data.table(id = names(urls_unstable), url = urls_stable)
      fwrite(gorilla_ids, file.path(output_dir, paste0(contrast_sel, "_GOrilla_urls.csv")))
      all_gorilla_ids <- rbindlist(list(all_gorilla_ids, gorilla_ids))
    }
  }
  if (nrow(all_gorilla_ids) > 0)
    fwrite(all_gorilla_ids, file.path(output_dir, paste0("All_contrasts", "_GOrilla_urls.csv")))
  return(all_gorilla_ids)
}
