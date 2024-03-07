
categorize_ORFs <- function(orfs_unl, groupings = strtoi(names(orfs_unl)), cds, mrna,
                            verbose = TRUE) {
  stopifnot(is(orfs_unl, "IRanges"))
  stopifnot(is.integer(groupings) & (max(groupings) <= length(mrna)))
  # Prepare CDS
  tx_ORF <- names(mrna)[groupings]
  cds_txcoord <- ranges(unlist(pmapToTranscriptF(cds, mrna)))
  cds_txcoord <- cds_txcoord[chmatch(tx_ORF, names(cds))]

  # Create categories
  is_uORF <- start(orfs_unl) < start(cds_txcoord)
  is_uoORF <- is_uORF & end(orfs_unl) > start(cds_txcoord)
  is_annotated <- orfs_unl == cds_txcoord
  is_NTE <- (start(orfs_unl) < start(cds_txcoord)) & (end(orfs_unl) == end(cds_txcoord))
  is_NTT <- (start(orfs_unl) > start(cds_txcoord)) & (end(orfs_unl) == end(cds_txcoord))
  is_internal <- (start(orfs_unl) > start(cds_txcoord)) & (end(orfs_unl) < end(cds_txcoord))
  is_doORF <- (start(orfs_unl) > start(cds_txcoord))  & (end(orfs_unl) > end(cds_txcoord))
  is_dORF <- is_doORF & (start(orfs_unl) > end(cds_txcoord))
  ORF_type <- rep("a_error", length(orfs_unl))
  ORF_type[is_uORF] <- "uORF"
  ORF_type[is_uoORF] <- "uoORF"
  ORF_type[is_annotated] <- "annotated"
  ORF_type[is_NTE] <- "NTE"
  ORF_type[is_NTT] <- "NTT"
  ORF_type[is_internal] <- "internal"
  ORF_type[is_doORF] <- "doORF"
  ORF_type[is_dORF] <- "dORF"

  if (verbose) print(table(ORF_type))
  if (any(ORF_type == "a_error")) {
    message("Warning some ORFs could not be categorized")
  }
  return(ORF_type)
}

#' Categorize and filter ORFs
#' @inheritParams detect_ribo_orfs
#' @param orfs IRanges of ORF coordinates
#' @param map_to_gr logical, default TRUE (genomic coordinates mapping)
#' @return Categorized and filtered ORFs as IRanges or GRanges depending on map_to_gr,
#' with metacolumn called
#' @noRd
categorize_and_filter_ORFs <- function(orfs, ORF_categories_to_keep,
                                       cds, mrna,
                                       map_to_gr = TRUE) {
  stopifnot(length(ORF_categories_to_keep) > 0)
  stopifnot(all(ORF_categories_to_keep %in% c("uORF", "uoORF", "annotated", "NTE",
                                              "NTT", "internal", "doORF", "dORF", "a_error")))

  orfs_unl <- unlist(orfs, use.names = TRUE)
  groupings <- strtoi(names(orfs_unl))
  names(orfs_unl) <- NULL
  tx_ORF_unique <- names(mrna)[as.integer(names(orfs))]

  ORF_type <- categorize_ORFs(orfs_unl, groupings, cds, mrna)
  mcols(orfs@unlistData)$ORF_type <- ORF_type
  filter_keep <- ORF_type %in% ORF_categories_to_keep
  orfs<- split(orfs@unlistData[filter_keep], groupings[filter_keep])
  ORF_type_keep <- mcols(orfs@unlistData)$ORF_type

  if (map_to_gr) {
    orfs <- pmapFromTranscriptF(orfs, mrna, removeEmpty = TRUE)
  }
  mcols(orfs)$category <- ORF_type_keep
  return(orfs)
}

coveragePerORFStatistics <- function(grl, RFP) {
  cov <- coveragePerTiling(grl, RFP, is.sorted = TRUE, as.data.table = TRUE, withFrames = TRUE)
  cov_stats <- cov[, .(sum = sum(count), median= median(count), mean = mean(count)), by = genes]
  if(any(widthPerGroup(grl, FALSE) == 0)) {
    na.dt <- merge.data.table(cov_stats,
                              data.table(genes = which(!(seq.int(length(grl)) %in% cov_stats$genes))),
                              by = "genes", all = TRUE)
    na.dt[is.na(na.dt)] <- 0
    cov_stats <- na.dt
  }
  orfScores <- orfScore(grl, RFP, TRUE, coverage = cov, stop3 = FALSE)
  cov_stats <- cbind(cov_stats,  orfScores)
  colnames(cov_stats)[5:7] <- c("F1", "F2", "F3")
  return(cov_stats)
}

#' Detect ORFs by Ribosome profiling data
#'
#' Finding all ORFs:
#' 1. Find all ORFs in mRNA using ORFik findORFs, with defined parameters.\cr
#' To create the candidate ORFs (all ORFs returned):\cr
#' Steps (candidate set):\cr
#' Define a candidate search set by these 3 rules:\cr
#'   1.a Allowed ORF type: uORF, NTE, etc (only keep these in candidate list)\cr
#'   1.b Must have at least x reads over whole orf (default 10 reads)\cr
#'   1.c Must have at least x reads over start site (default 3 reads)\cr
#' The total list is defined by these names, and saved according to allowed ORF type/types.\cr
#' To create the prediction status (TRUE/FALSE) per candidate\cr
#'  Steps (prediction status)\cr
#' (UP_NT is a 20nt window upstream of ORF, that stops 2NT before ORF starts) :\cr
#'   1. ORF mean reads per NT > (UP_NT mean reads per NT * 1.3)\cr
#'   2. ORFScore > 2.5\cr
#'   3. TIS total reads + 3 >  ORF median reads per NT\cr
#'   4. Given expression above, a TRUE prediction is defined with the AND operatior: 1. & 2. & 3.
#' \cr\cr
#' In code that is:\cr
#' \code{predicted <- (orfs_cov_stats$mean > upstream_cov_stats$mean*1.3) & orfs_cov_stats$ORFScores > 2.5 &
#'  ((reads_start[candidates] + 3) >  orfs_cov_stats$median)}
#' @inheritParams findORFs
#' @inheritParams outputLibs
#' @param out_folder Directory to save files
#' @param ORF_categories_to_keep options, any subset of: \code{c("uORF", "uoORF", "annotated", "NTE",
#' "NTT", "internal", "doORF", "dORF", "a_error")}.
#' \itemize{
#'  \item{uORF: }{Upstream ORFs (Starting in 5' UTR), not overlapping CDS}
#'  \item{uoORF: }{Upstream ORFs (Starting in 5' UTR), overlapping CDS}
#'  \item{annotated: }{The defined CDS for that transcript}
#'  \item{NTE: }{5' Start codon extension of annotated CDS}
#'  \item{NTT: }{5' Start codon truncation of annotated CDS}
#'  \item{internal: }{Starting inside CDS, ending before CDS ends}
#'  \item{doORF: }{Downstream ORFs (Ending in 3' UTR), overlapping CDS}
#'  \item{dORF: }{Downstream ORFs (Ending in 3' UTR), not overlapping CDS}
#'  \item{a_error: }{Any ORF detect not in the above categories}
#' }
#' @param prefix_result the prefix name of output files to out_folder. Default:
#' \code{paste(c(ORF_categories_to_keep, gsub(" ", "_", organism(df))), collapse = "_")}
#' @param mrna = \code{loadRegion(df, "mrna")}
#' @param cds = \code{loadRegion(df, "cds")}
#' @param libraries the ribo-seq libraries loaded into R as list, default:
#' \code{outputLibs(df, type = "pshifted", output = "envirlist")}
#' @param orf_candidate_ranges IRangesList, =
#'  \code{findORFs(seqs = txSeqsFromFa(mrna, df, TRUE),
#'  longestORF = longestORF, startCodon = startCodon, stopCodon = stopCodon,
#'  minimumLength = minimumLength)}
#' @param export_metrics_table logical, default TRUE. Export table of statistics to file
#' with suffix: "_prediction_table.rds"
#' @param minimum_reads_ORF numeric, default 10, orf removed if less reads overlap whole orf
#' @param minimum_reads_start numeric, default 3, orf removed if less reads overlap start
#' @return invisible(NULL), all ORF results saved to disc
#' @export
#' @examples
#' # Pre requisites
#' # 1. Create ORFik experiment
#' #  ORFik::create.experiment(...)
#' # 2. Create ORFik optimized annotation:
#' # makeTxdbFromGenome(gtf = ORFik:::getGtfPathFromTxdb(df), genome = df@fafile, organism = organism(df), optimize = TRUE)
#' # 3. There must exist pshifted reads, either as default files, or in a relative folder called
#' # "./pshifted/". See ?shiftFootprintsByExperiment
#' # EXAMPLE:
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RFP",][c(1,2),]
#' result_folder <- riboORFsFolder(df, tempdir())
#' results <- detect_ribo_orfs(df, result_folder, c("uORF", "uoORF", "annotated", "NTE"))
#'
#' # Load results of annotated ORFs
#' table <- riboORFs(df[1,], type = "table", result_folder)
#' table # See all statistics
#' sum(table$predicted) # How many were predicted as Ribo-seq ORFs
#' # Load 2 results
#' table <- riboORFs(df[1:2,], type = "table", result_folder)
#' table # See all statistics
#' sum(table$predicted) # How many were predicted as Ribo-seq ORFs
#'
#' # Load GRangesList
#' candidates_gr <- riboORFs(df[1,], type = "ranges_candidates", result_folder)
#' prediction <- riboORFs(df[1,], type = "predictions", result_folder)
#'
#' predicted_gr <- riboORFs(df[1:2,], type = "ranges_predictions", result_folder)
#' identical(predicted_gr[[1]], candidates_gr[[1]][prediction[[1]]])
#' ## Inspect predictions in RiboCrypt
#' # library(RiboCrypt)
#' # Inspect Predicted
#' view <- predicted_gr[[1]][1]
#' #multiOmicsPlot_ORFikExp(view, df, view, leader_extension = 100, trailer_extension = 100)
#' # Inspect not predicted
#' view <- candidates_gr[[1]][!prediction[[1]]][1]
#' #multiOmicsPlot_ORFikExp(view, df, view, leader_extension = 100, trailer_extension = 100)
detect_ribo_orfs <- function(df, out_folder, ORF_categories_to_keep,
                             prefix_result = paste(c(ORF_categories_to_keep, gsub(" ", "_", organism(df))), collapse = "_"),
                             mrna = loadRegion(df, "mrna"), cds = loadRegion(df, "cds"),
                             libraries = outputLibs(df, type = "pshifted", output = "envirlist"),
                             orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df, TRUE), longestORF = longestORF,
                                                            startCodon = startCodon, stopCodon = stopCodon,
                                                            minimumLength = minimumLength),
                             export_metrics_table = TRUE,
                             longestORF = FALSE, startCodon =  startDefinition(1),
                             stopCodon = stopDefinition(1), minimumLength = 0,
                             minimum_reads_ORF = 10, minimum_reads_start = 3) {
  start_timer <- Sys.time()
  message("Finding all candidate ORFs")
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)


  # Filter categories and map to GRanges
  orfs_gr <- categorize_and_filter_ORFs(orf_candidate_ranges,
                                        ORF_categories_to_keep, cds, mrna)
  orf_start_gr <- startSites(orfs_gr, TRUE,TRUE,TRUE)
  orf_stop_gr <- stopSites(orfs_gr, TRUE,TRUE,TRUE)
  # orf_start_overlaps_other <- countOverlaps(orf_start_gr, orfs_gr) - 1

  message("Start Ribo-seq coverage analysis")
  libnames <- name_decider(df, naming = "full")
  symbols <- suppressMessages(symbols(df))
  txdb <- NULL
  if (is(symbols, "try-error") || nrow(symbols) == 0) txdb <- loadTxdb(df)
  ORF_type_keep <- mcols(orfs_gr)$category
  mcols(orfs_gr) <- NULL
  out_file_prefixes <- file.path(out_folder, paste0(prefix_result,
                                                    "_", name(df),
                                                    "_", libnames))
  for (i in seq_along(libraries)) {
    message("- ", libnames[i])
    detect_ribo_orfs_single_cov(orfs_gr, libraries[[i]], out_file_prefixes[i],
                                mrna, txdb = txdb, faFile = df,
                                ORF_type_keep, orf_start_gr, orf_stop_gr,
                                export_metrics_table, symbols,
                                minimum_reads_ORF, minimum_reads_start)

  }
  print(Sys.time() - start_timer)
  message("Done")
  return(NULL)
}

detect_ribo_orfs_single_cov <- function(orfs_gr, RFP, out_file_prefix, mrna,
                                        txdb = NULL, faFile = NULL,
                                        ORF_type_keep = mcols(orfs_gr)$category,
                                        orf_start_gr = startSites(orfs_gr, TRUE,TRUE,TRUE),
                                        orf_stop_gr = stopSites(orfs_gr, TRUE,TRUE,TRUE),
                                        export_metrics_table = TRUE, symbols,
                                        minimum_reads_ORF = 10, minimum_reads_start = 3) {
  stopifnot(!is.null(ORF_type_keep))
  if (ncol(mcols(orfs_gr)) > 0) mcols(orfs_gr) <- NULL
  # Filter candidates by counts
  message("-- Filtering")
  reads_10 <- countOverlapsW(orfs_gr, RFP, weight = "score") > minimum_reads_ORF
  reads_start <- countOverlapsW(orf_start_gr, RFP, weight = "score")
  reads_start_3 <- reads_start > minimum_reads_start

  # Define candidates
  candidates <- reads_10 & reads_start_3
  orfs_cand <- orfs_gr[candidates]
  orfs_kept_percentage <- round((length(orfs_cand) / length(orfs_gr)) * 100, 2)
  message("--- ORFs passed count filters: ", orfs_kept_percentage)
  orf_start_cand <- orf_start_gr[candidates]
  orf_stop_cand <- orf_stop_gr[candidates]
  # Make upstream and downstream window
  upstream_gr <- windowPerGroup(orf_start_cand, mrna, 20, -2)
  downstream_gr <- windowPerGroup(orf_stop_cand, mrna, -2, 20)
  # Calculate coverage for ORF and up/down
  message("-- coverage calculations")
  message("--- ORFs")
  orfs_cov_stats <- coveragePerORFStatistics(orfs_cand, RFP)
  message("--- upstream")
  upstream_cov_stats <- coveragePerORFStatistics(upstream_gr, RFP)
  message("--- downstream")
  downstream_cov_stats <- coveragePerORFStatistics(downstream_gr, RFP)

  #iou <- (orfs_cov_stats$mean + 1) / (upstream_cov_stats$mean + 1)
  predicted <- (orfs_cov_stats$mean > upstream_cov_stats$mean*1.3) & orfs_cov_stats$ORFScores > 2.5 &
    ((reads_start[candidates] + 3) >  orfs_cov_stats$median)

  message("-- Saving ORF and prediction result")
  saveRDS(orfs_cand, paste0(out_file_prefix, "_candidates.rds"))
  saveRDS(predicted, paste0(out_file_prefix, "_prediction.rds"))

  if(export_metrics_table) {
    if (!is.null(faFile)) {
      start_codons <- txSeqsFromFa(startCodons(orfs_cand, TRUE), faFile, TRUE, keep.names = FALSE)
    }
    if (is(symbols, "data.table") && nrow(symbols) > 0 && "ensembl_tx_name" %in% colnames(symbols)) {
      naming <- data.table::merge.data.table(data.table(ensembl_tx_name = names(orfs_cand)),
                                             symbols, by = "ensembl_tx_name", all.x = TRUE, sort = FALSE)
    } else if (!is.null(txdb)) {
      naming <- data.table(gene = txNamesToGeneNames(names(orfs_cand), txdb), tx = names(orfs_cand))
    }
    res <- cbind(naming, type = ORF_type_keep[candidates], predicted, length = widthPerGroup(orfs_cand, FALSE),
                 start_codons, orfs_cov_stats, up = upstream_cov_stats, down = downstream_cov_stats)
    res[, `:=` (down.genes = NULL, up.genes = NULL)]
    setnames(res, "genes", "id")
    saveRDS(res, paste0(out_file_prefix, "_prediction_table.rds"))
  }
  return(invisible(NULL))
}

#' Define folder for prediction output
#'
#' @param df ORFik experiment
#' @param parrent_dir Parrent directory of computed study results, default:
#' resFolder(df)
#' @return a file path (full path)
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RFP",][c(1,2),]
#' riboORFsFolder(df)
#' riboORFsFolder(df, tempdir())
riboORFsFolder <- function(df, parrent_dir = resFolder(df)) {
  dir <- paste0("Ribo_orfs_", gsub(" ", "_", organism(df)))
  return(file.path(parrent_dir, dir))
}

#' Load Predicted translons
#'
#' @param df ORFik experiment
#' @param type default "table", alternatives: c("table", "ranges_candidates",
#'  "ranges_predictions", "predictions")
#' @param folder base folder to check for computed results, default:
#'  riboORFsFolder(df)
#' @return a data.table, GRangesList or list of logical vector depending on input
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[df$libtype == "RFP",][c(1,2),]
#' # riboORFs(df) # Works when you have run prediction
riboORFs <- function(df, type = "table", folder = riboORFsFolder(df)) {
  stopifnot(type %in% c("table", "ranges_candidates", "ranges_predictions", "predictions"))
  stopifnot(is(df, "experiment"))
  stopifnot(dir.exists(folder))
  search_terms <- paste0(c("_prediction_table", "_candidates", "_candidates", "_prediction"), "\\.rds$")
  names(search_terms) <- c("table", "ranges_candidates", "ranges_predictions", "predictions")
  search_term <- search_terms[type]

  file_paths <- list.files(folder, paste(search_term, collapse = "|"),
                           full.names = TRUE)
  libnames <- name_decider(df, naming = "full")
  list <- lapply(seq_along(libnames), function(id) {
    x <- libnames[id]
    file <- grep(pattern = x, x = file_paths, value = TRUE)
    if (length(file) > 1) {
      print(file)
      stop("multiple files found for this library, see above for which")
    } else if (length(file) == 0) {
      print(file_paths)
      stop("Could not find file!")
    }

    rds <- readRDS(file)
    if (type == "table") {
      rds <- cbind(rds, library = x)
      rds[]
    } else if (type == "ranges_predictions") {
      predictions <- riboORFs(df[id, ], "predictions", folder)
      rds <- rds[predictions[[1]]]
    }
    return(rds)
  })
  if (type == "table") list <- rbindlist(list)
  return(list)
}
