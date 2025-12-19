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
#' @return a data.table with 2 columns (id: name of contrast), urls: url to html online,
#' this table is also saved in output_dir.
#' @export
#' @family GOrilla
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
        url <- go_analaysis_gorilla(target_file, background_file, organism,
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

#' Copy GOrilla result htmls to local
#'
#' Will retrieve full html, png and xls structure so analysis can be used even
#' when results are deleted online (1 month after creation).\cr
#' Files are saved to disc in directory "./GOrilla_local_html_outputs/", relative
#' to input directory 'gorilla_output_dir'. There is 1 subfolder per
#' analysis url. Open the GOResults.html to view.
#' @param gorilla_output_dir character, directory path to existing results of a
#' DEG_gorilla call. Must contain a file with relative path "./All_contrasts_GOrilla_urls.csv"
#' @param local_html_dir character, output directory,
#' default: \code{file.path(gorilla_output_dir, "GOrilla_local_html_outputs")}
#' @return invisible(NULL), files are saved to disc.
#' @export
#' @family GOrilla
DEG_gorilla_copy_to_local <- function(gorilla_output_dir,
                                      local_html_dir = file.path(gorilla_output_dir, "GOrilla_local_html_outputs")) {
  stopifnot(dir.exists(gorilla_output_dir) && length(gorilla_output_dir) == 1)
  urls <- fread(file.path(gorilla_output_dir, "All_contrasts_GOrilla_urls.csv"))
  message("Creating local copy of GOrilla html structures..")
  for (i in seq_along(urls$url)) {
    gorilla_copy_to_local(urls$url[i], urls$id[i], local_html_dir, i, nrow(urls))
  }
  message("GOrilla htmls also saved localy at location:\n", local_html_dir)
  return(invisible(NULL))
}

gorilla_copy_to_local <- function(url, id, local_html_dir, this_url_index = 1,
                                  total_urls = 1) {
  stopifnot(is.character(url) & length(url) == 1)
  stopifnot(is.character(id) & length(id) == 1)

  message("- ", id, " (", this_url_index, "/", total_urls, ")")
  this_dir <- file.path(local_html_dir, gsub(" ", "_", id))
  dir.create(this_dir, showWarnings = FALSE, recursive = TRUE)

  headers <- c("", "top")
  go_types <- c("PROCESS", "FUNCTION", "COMPONENT")
  html_components <- c(headers, go_types)
  for (go_type in html_components) {
    cat(paste0(ifelse(go_type == "", "RESULTS", go_type)))
    go_type_url <- sub("\\.html$", paste0(go_type, ".html"), url)
    if (go_type == "top") go_type_url <- sub("GOResults", "", go_type_url)
    out_file <- file.path(this_dir, basename(go_type_url))
    cat("(HTML)")
    download.file(go_type_url, out_file)
    if (go_type %in% go_types) {
      for (additional_format in c("png", "xls")) {
        cat(paste0("(", toupper(additional_format), ")"))
        go_type_url <- file.path(dirname(go_type_url), paste0("GO", go_type, ".", additional_format))
        out_file <- file.path(this_dir, basename(go_type_url))
        download.file(go_type_url, out_file)
      }
    }
  }
  cat("\n")
  return(invisible(NULL))
}

#' Load all GOrilla xls files from study
#'
#' Output as data.table with column 'analysis' describing DEG
#' contrast and column 'go_category' describing COMPONENT, FUNCTION or PROCESS.
#' @param gorilla_local_dir path to directory of local GORilla html and xls data
#' @return a data.table
#' @export
DEG_gorilla_local_load_data <- function(gorilla_local_dir) {
  stopifnot(dir.exists(gorilla_local_dir))
  files <- list.files(gorilla_local_dir, pattern = "\\.xls", recursive = TRUE,
                      full.names = TRUE)
  dt <- rbindlist(lapply(files, fread), fill = TRUE, idcol = TRUE)
  analysis_categories <- gsub("_", " ", sub("^DTEG_Comparison:_", "", basename(dirname(files))))
  dt[, analysis := analysis_categories[`.id`]]

  go_categories <- gsub("^GO|\\.xls$", "", basename(files))
  dt[, `.id` := go_categories[`.id`]]
  dt[, p.adjust := `FDR q-value`]
  colnames(dt)[1] <- "go_category"
  dt[, Enrichment := as.numeric(sub("^1,", "", Enrichment))]
  # Parse GeneRatio "b/n" -> numeric; add -log10(adjusted p)
  dt[, GeneRatio := b / n]
  attr(dt, "species") <- DEG_gorilla_local_species(gorilla_local_dir)
  dt[] # Make print work first time
  return(dt)
}

#' Plot GOrilla xls results
#'
#' Inspired by the enrichment plot from the package clusterProfiler.
#' @param dt a data.table with results loaded using
#' \code{DEG_gorilla_local_load_data}
#' @param top_n maximum number of GO terms per analysis / component split. The
#' sorting is done as -GeneRatio, -p.adjust, so it extracts the top 20
#' highest GeneRatio terms by p.adjust values (decreasing sort).
#' @param enrich_cutoff numeric, default 8. Minimum enrichment, set it lower
#' to get more generic groups, higher to get more specific groups.
#' @return a ggplot, it uses facet_wrap to split analysis / components.
#' @export
DEG_gorilla_plot <- function(dt, top_n = 20L, enrich_cutoff = 8) {
  stopifnot(is.numeric(top_n) & is.numeric(enrich_cutoff))
  stopifnot(is(dt, "data.table"))
  cols <- c(
    "go_category", "GO Term", "Description", "P-value", "FDR q-value", "Enrichment",
    "N", "B", "n", "b", "Genes", "analysis", "p.adjust", "GeneRatio"
  )
  stopifnot(all(cols %in% colnames(dt)))
  species <- attr(dt, "species")
  if (!is.null(species)) species <- paste0("(", species, ")")

  plot_dt <- dt[is.finite(p.adjust) & is.finite(GeneRatio) & Enrichment > enrich_cutoff]
  plot_dt <- plot_dt[order(-GeneRatio, -p.adjust), head(.SD, min(top_n, .N)), by = .(analysis, go_category)]
  # order y-axis (GO terms) by GeneRatio (or by negLog10Padj if you prefer)
  # Highest on top:
  setorder(plot_dt, -GeneRatio)  # ascending
  plot_dt[, Description := factor(Description, levels = unique(Description))]
  plot_dt[, Description := factor(Description, levels = rev(levels(Description)))]

  # ggplot dot plot
  g <- ggplot(plot_dt,
              aes(x = GeneRatio, y = Description, size = b, color = p.adjust)) +
    geom_point(alpha = 0.85) +
    scale_size_area(max_size = 12, guide = guide_legend(title = "Gene set size")) +
    scale_color_gradientn(
      colors = rev(c("#4575B4", "#7E3ABF", "#D73027")),  # blue–purple–red
      guide = guide_colorbar(title = "p.adjust"),
      limits = c(min(plot_dt$p.adjust, na.rm = TRUE),
                 max(plot_dt$p.adjust, na.rm = TRUE))
    ) +
    labs(
      x = "Gene ratio", y = NULL,
      title = paste0("GO enrichment"),
      subtitle = paste0(species,
                        "\nTop ", top_n, " terms by adjusted p-value (Enrichment ratio >", enrich_cutoff, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title  = element_text(face = "bold"),
      legend.title = element_text(size = 10)
    ) + facet_wrap(analysis ~ go_category, ncol = 3, scales = "free_y")
  return(g)
}

DEG_gorilla_local_species <- function(gorilla_local_dir) {
  species <- NA_character_
  pattern <- "(PROCESS|FUNCTION|COMPONENT)\\.html"
  file <- list.files(gorilla_local_dir, pattern = pattern, recursive = TRUE,
                     full.names = TRUE)[1]
  if (length(file) == 1) {
    txt <- readLines(file)
    species <- sub(".*Species used:\\s*([^<]+)<.*", "\\1",
                   grep("Species", txt, value = TRUE))
  }
  return(species)
}
