context("Differential analysis")
library(ORFik)

df <- ORFik.template.experiment()
df.rna <- df[df$libtype == "RNA",]
df.rfp <- df[df$libtype == "RFP",]

make_DTEG_interaction_fixture <- function(libtype_levels = c("RNA", "RFP"),
                                          condition_levels = c("WT", "Mutant"),
                                          add_zero_cases = FALSE) {
  set.seed(1)
  genes <- paste0("gene", seq_len(80))
  base_counts <- sample(80:600, length(genes), replace = TRUE)

  rna_condition_effect <- rep(1, length(genes))
  rfp_condition_effect <- rep(1, length(genes))
  rna_condition_effect[1:20] <- 4
  rfp_condition_effect[1:20] <- 4
  rna_condition_effect[21:40] <- 1
  rfp_condition_effect[21:40] <- 4
  rna_condition_effect[41:60] <- 4
  rfp_condition_effect[41:60] <- 1

  make_count_matrix <- function(condition_effect, libtype) {
    means <- cbind(base_counts, base_counts,
                   base_counts * condition_effect,
                   base_counts * condition_effect)
    matrix(rnbinom(length(means), mu = as.vector(means), size = 30),
           nrow = length(genes),
           dimnames = list(genes,
                           paste0(libtype, "_",
                                  c("WT_1", "WT_2", "Mutant_1", "Mutant_2"))))
  }

  rna_counts <- make_count_matrix(rna_condition_effect, "RNA")
  rfp_counts <- make_count_matrix(rfp_condition_effect * 2, "RFP")
  if (add_zero_cases) {
    rna_counts[61, c("RNA_WT_1", "RNA_WT_2")] <- 0
    rfp_counts[62, c("RFP_WT_1", "RFP_WT_2")] <- 0
    rna_counts[63, c("RNA_Mutant_1", "RNA_Mutant_2")] <- 0
    rfp_counts[64, c("RFP_Mutant_1", "RFP_Mutant_2")] <- 0
  }
  main_coldata <- S4Vectors::DataFrame(
    condition = factor(rep(c("WT", "Mutant"), each = 2),
                       levels = condition_levels),
    row.names = colnames(rna_counts)
  )

  dds_rna <- DESeq2::DESeq(DESeq2::DESeqDataSet(
    SummarizedExperiment::SummarizedExperiment(list(counts = rna_counts),
                                               colData = main_coldata),
    design = ~ condition), fitType = "mean", quiet = TRUE)

  row.names(main_coldata) <- colnames(rfp_counts)
  dds_rfp <- DESeq2::DESeq(DESeq2::DESeqDataSet(
    SummarizedExperiment::SummarizedExperiment(list(counts = rfp_counts),
                                               colData = main_coldata),
    design = ~ condition), fitType = "mean", quiet = TRUE)

  pooled_counts <- cbind(rna_counts, rfp_counts)
  pooled_coldata <- S4Vectors::DataFrame(
    libtype = factor(rep(c("RNA", "RFP"), each = 4),
                     levels = libtype_levels),
    condition = factor(rep(rep(c("WT", "Mutant"), each = 2), times = 2),
                       levels = condition_levels),
    row.names = colnames(pooled_counts)
  )
  dds_te <- DESeq2::DESeq(DESeq2::DESeqDataSet(
    SummarizedExperiment::SummarizedExperiment(list(counts = pooled_counts),
                                               colData = pooled_coldata),
    design = ~ libtype + condition + libtype:condition),
    fitType = "mean", quiet = TRUE)

  list(rna = dds_rna, rfp = dds_rfp, te = dds_te)
}

fixture <- make_DTEG_interaction_fixture()
interaction <- DESeq2::results(fixture$te,
                               name = "libtypeRFP.conditionMutant")
inverse_fixture <- make_DTEG_interaction_fixture(libtype_levels = c("RFP", "RNA"),
                                                 condition_levels = c("Mutant", "WT"),
                                                 add_zero_cases = TRUE)
inverse_interaction <- DESeq2::results(inverse_fixture$te,
                                       name = "libtypeRNA.conditionWT")

sink <- capture.output(dteg_forward <- ORFik:::DTEG_pair_results(
  fixture$te, fixture$rfp, fixture$rna,
  c("condition", "Mutant", "WT"), "normal", 0.05
))
sink <- capture.output(dteg_reverse <- ORFik:::DTEG_pair_results(
  fixture$te, fixture$rfp, fixture$rna,
  c("condition", "WT", "Mutant"), "normal", 0.05
))
sink <- capture.output(dteg_inverse_orientation <- ORFik:::DTEG_pair_results(
  inverse_fixture$te, inverse_fixture$rfp, inverse_fixture$rna,
  c("condition", "WT", "Mutant"), "normal", 0.05
))
sink <- capture.output(dteg_inverse_orientation_forward <- ORFik:::DTEG_pair_results(
  inverse_fixture$te, inverse_fixture$rfp, inverse_fixture$rna,
  c("condition", "Mutant", "WT"), "normal", 0.05
))


test_that("DEG analysis works", {
  sink <- capture.output(dt <- suppressWarnings(DEG.analysis(df.rna)))
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 6)
  expect_equal(unique(as.character(dt$Regulation)), "No change")
})

test_that("DTEG analysis works", {
  sink <- capture.output(dt <- suppressWarnings(DTEG.analysis(df.rfp, df.rna,
                                                              output.dir = NULL,
                                                              plot_to_console = FALSE)))
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 6)
  expect_equal(unique(as.character(dt$Regulation)), "No change")
})

test_that("DTEG TE result uses the libtype-condition interaction", {

  condition_main_effect <- DESeq2::results(fixture$te,
                                           contrast = c("condition", "Mutant", "WT"))

  expect_equal(dteg_forward$te.lfc, interaction$log2FoldChange, tolerance = 1e-8)
  expect_equal(dteg_forward$te.padj, interaction$padj, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(dteg_forward$te.lfc,
                                condition_main_effect$log2FoldChange,
                                tolerance = 1e-8)))
})

test_that("DTEG TE interaction changes sign for reversed contrast", {
  expect_equal(dteg_reverse$te.lfc, -interaction$log2FoldChange, tolerance = 1e-8)
  expect_equal(dteg_reverse$te.padj, interaction$padj, tolerance = 1e-8)
})

test_that("DTEG TE result has expected direction for known TE changes", {
  expect_gt(median(dteg_forward$te.lfc[21:40]), 1)
  expect_lt(median(dteg_forward$te.lfc[41:60]), -1)
})

test_that("DTEG TE result keeps RFP/RNA sign when model encodes RNA/RFP", {
  expect_equal(dteg_inverse_orientation$te.lfc,
               -inverse_interaction$log2FoldChange,
               tolerance = 1e-8)
  expect_equal(dteg_inverse_orientation$te.padj,
               inverse_interaction$padj,
               tolerance = 1e-8)

  expect_equal(dteg_inverse_orientation_forward$te.lfc,
               inverse_interaction$log2FoldChange,
               tolerance = 1e-8)
  expect_equal(dteg_inverse_orientation_forward$te.padj,
               inverse_interaction$padj,
               tolerance = 1e-8)
})

test_that("DTEG TE result avoids ratio infinities for zero-count edge cases", {
  zero_case_rows <- 61:64

  expect_false(any(is.infinite(dteg_inverse_orientation$te.lfc[zero_case_rows])))
  expect_false(any(is.infinite(dteg_inverse_orientation$te.padj[zero_case_rows])))
})

test_that("te.table works", {
  dt <- suppressWarnings(te.table(df.rfp, df.rna))
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 24)
})
