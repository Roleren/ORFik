
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
#'   1.b Must have at least x reads (default 10 reads)\cr
#'   1.c Must have at least x reads over start site (default 3 reads)\cr
#' The total list is defined by these names, and saved according to allowed ORF type/types.\cr
#' To create the prediction status (TRUE/FALSE) per candidate\cr
#'  Steps (prediction status)\cr
#' (UP_NT is a 20nt window upstream of ORF, that stops 2NT before ORF starts) :
#'   1. ORF mean reads per NT > (UP_NT mean reads per NT * 1.3)
#'   2. ORFScore > 2.5
#'   3. TIS total reads + 3 >  ORF median reads per NT
#'   4. Given expression above, a TRUE prediction is defined with the AND operatior: 1. & 2. & 3.
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
#' @param mrna = \code{loadRegion(df, "mrna")},
#' @param cds = \code{loadRegion(df, "cds")},
#' @param orf_sequences = \code{findORFs(seqs = txSeqsFromFa(mrna, df, TRUE), longestORF = longestORF,
#'  startCodon = startCodon, stopCodon = stopCodon,
#'  minimumLength = minimumLength)}
#' @param export_metrics_table logical, default TRUE. Export table of statistics to file
#' with suffix: "_prediction_table.rds"
#' @param minimum_reads_ORF numeric, default 10, orf removed if less reads overlap whole orf
#' @param minimum_reads_start numeric, default 3, orf removed if less reads overlap start
#' @return invisible(NULL), all ORF results saved to disc
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
#' result_folder <- file.path(tempdir(), paste0("Ribo_orfs_", organism(df)))
#' results <- detect_ribo_orfs(df, result_folder, c("uORF", "uoORF", "annotated", "NTE"))
#'
#' # Load results of annotated ORFs
#' res_annotated <- file.path(result_folder, "uORF_uoORF_annotated_NTE_Homo sapiens_ORFik_RFP_r1_prediction_table.rds")
#' annotated <- readRDS(res_annotated)
#' annotated # See all statistics
#' sum(annotated$predicted) # How many were predicted as Ribo-seq ORFs
detect_ribo_orfs <- function(df, out_folder, ORF_categories_to_keep,
                             name_of_result = paste(c(ORF_categories_to_keep, organism(df)), collapse = "_"),
                             mrna = loadRegion(df, "mrna"), cds = loadRegion(df, "cds"),
                             orf_sequences = findORFs(seqs = txSeqsFromFa(mrna, df, TRUE), longestORF = longestORF,
                                                      startCodon = startCodon, stopCodon = stopCodon,
                                                      minimumLength = minimumLength),
                             export_metrics_table = TRUE,
                             longestORF = FALSE, startCodon =  startDefinition(1),
                             stopCodon = stopDefinition(1), minimumLength = 0,
                             minimum_reads_ORF = 10, minimum_reads_start = 3) {
  message("Finding all candidate ORFs")
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  stopifnot(length(ORF_categories_to_keep) > 0)
  stopifnot(all(ORF_categories_to_keep %in% c("uORF", "uoORF", "annotated", "NTE",
                                              "NTT", "internal", "doORF", "dORF", "a_error")))

  # New group you can start here
  orfs <- orf_sequences
  orfs_unl <- unlist(orfs, use.names = TRUE)
  groupings <- strtoi(names(orfs_unl))
  names(orfs_unl) <- NULL
  tx_ORF_unique <- names(mrna)[as.integer(names(orfs))]

  ORF_type <- categorize_ORFs(orfs_unl, groupings, cds, mrna)
  mcols(orfs@unlistData)$ORF_type <- ORF_type
  filter_keep <- ORF_type %in% ORF_categories_to_keep
  orfs<- split(orfs@unlistData[filter_keep], groupings[filter_keep])
  ORF_type_keep <- mcols(orfs@unlistData)$ORF_type

  orfs_gr <- pmapFromTranscriptF(orfs, mrna, removeEmpty = TRUE)

  orfs_chr <- seqnamesPerGroup(orfs_gr, FALSE)
  orf_start_gr <- startSites(orfs_gr, T,T,T)
  orf_stop_gr <- stopSites(orfs_gr, T,T,T)
  orf_start_overlaps_other <- countOverlaps(orf_start_gr, orfs_gr) - 1

  message("Start Ribo-seq coverage analysis")
  list <- outputLibs(df, type = "pshifted", output = "envirlist")
  libnames <- bamVarName(df)
  for (i in seq_along(list)) {
    message("- ", libnames[i])
    RFP <- list[[i]]
    out_file_prefix <- file.path(out_folder, paste0(name_of_result,
                                                    "_", name(df),
                                                    "_", libnames[i]))
    reads_10 <- countOverlapsW(orfs_gr, RFP, weight = "score") > minimum_reads_ORF
    reads_start <- countOverlapsW(orf_start_gr, RFP,weight = "score")
    reads_start_3 <- reads_start > minimum_reads_start
    # Define candidates
    candidates <- reads_10 & reads_start_3
    orfs_cand <- orfs_gr[candidates]
    orf_start_cand <- orf_start_gr[candidates]
    orf_stop_cand <- orf_stop_gr[candidates]
    # Make upstream and downstream window
    upstream_gr <- windowPerGroup(orf_start_cand, mrna, 20, -2)
    downstream_gr <- windowPerGroup(orf_stop_cand, mrna, -2, 20)
    # Calculate coverage for ORF and up/down
    orfs_cov_stats <- coveragePerORFStatistics(orfs_cand, RFP)
    upstream_cov_stats <- coveragePerORFStatistics(upstream_gr, RFP)
    downstream_cov_stats <- coveragePerORFStatistics(downstream_gr, RFP)

    #iou <- (orfs_cov_stats$mean + 1) / (upstream_cov_stats$mean + 1)
    predicted <- (orfs_cov_stats$mean > upstream_cov_stats$mean*1.3) & orfs_cov_stats$ORFScores > 2.5 &
      ((reads_start[candidates] + 3) >  orfs_cov_stats$median)

    message("-- Saving ORF and prediction result")
    saveRDS(orfs_cand, paste0(out_file_prefix, "_candidates.rds"))
    saveRDS(predicted, paste0(out_file_prefix, "_prediction.rds"))

    if(export_metrics_table) {
      start_codons <- txSeqsFromFa(startCodons(orfs_cand, TRUE), df, TRUE, keep.names = FALSE)
      res <- cbind(gene = txNamesToGeneNames(names(orfs_cand), df), tx = names(orfs_cand),
                   type = ORF_type_keep[candidates], predicted, length = widthPerGroup(orfs_cand, F),
                   start_codons, orfs_cov_stats, up = upstream_cov_stats, down = downstream_cov_stats)
      res[, `:=` (down.genes = NULL, up.genes = NULL)]
      setnames(res, "genes", "id")
      saveRDS(res, paste0(out_file_prefix, "_prediction_table.rds"))
    }
  }
  message("Done")
  return(NULL)
}
