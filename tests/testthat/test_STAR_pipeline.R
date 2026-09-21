# These tests execute the bundled scripts through the public R entry points.
# Stub executables record the arguments they actually receive and create only
# the files needed by subsequent stages; real-tool checks are separate.
star_fixture <- function(paired = FALSE, nested = FALSE) {
  skip_on_os("windows")
  skip_if(Sys.which("python3") == "", "python3 is needed for the tool stub")
  root <- withr::local_tempdir(pattern = "STAR paths ", .local_envir = parent.frame())
  input <- file.path(root, "input")
  dir.create(input)
  library.dir <- if (nested) file.path(input, "nested") else input
  dir.create(library.dir, showWarnings = FALSE)
  reads <- file.path(library.dir, if (paired) c("sample_1.fastq", "sample_2.fastq") else "sample.fastq")
  for (read in reads) writeLines(c("@read", "ACGTACGTACGTACGTACGTACGT", "+", strrep("I", 24)), read)
  index <- file.path(root, "index")
  for (sub in c("genomeDir", "contaminants_genomeDir", "tRNA_genomeDir"))
    dir.create(file.path(index, sub), recursive = TRUE)
  bins <- file.path(root, c("STAR", "fastp"))
  for (bin in bins) {
    file.copy(test_path("fixtures", "star_fastp_stub.py"), bin)
    Sys.chmod(bin, "0755")
  }
  record <- file.path(root, "calls.jsonl")
  withr::local_envvar(c(ORFIK_TEST_CALLS = record, ORFIK_TEST_FAIL = NA),
                     .local_envir = parent.frame())
  list(root = root, input = input, reads = reads, index = index,
       bins = bins, record = record, output = file.path(root, "output $(literal)"))
}
star_calls <- function(f) lapply(readLines(f$record), jsonlite::fromJSON)
star_arg <- function(call, option) {
  at <- match(option, call$args)
  if (is.na(at)) return(NULL)
  call$args[at + 1L]
}
star_run <- function(f, folder = FALSE, ...) {
  common <- list(output.dir = f$output, index.dir = f$index,
                 star.path = f$bins[1], fastp = f$bins[2], max.cpus = 1,
                 multiQC = FALSE, verbose = FALSE, ...)
  if (folder) do.call(STAR.align.folder, c(list(input.dir = f$input), common))
  else do.call(STAR.align.single, c(list(file1 = f$reads[1],
                                        file2 = if (length(f$reads) == 2L) f$reads[2]), common))
}

test_that("both entry points forward both intron modes and default permits discovery", {
  for (folder in c(FALSE, TRUE)) {
    f <- star_fixture()
    star_run(f, folder, steps = "ge")
    expect_identical(star_arg(tail(star_calls(f), 1)[[1]], "--alignIntronMax"), "0")
    star_run(f, folder, steps = "ge", allow.introns = FALSE)
    expect_identical(star_arg(tail(star_calls(f), 1)[[1]], "--alignIntronMax"), "1")
  }
})

test_that("paired folder stages use actual previous outputs and preserve options", {
  for (nested in c(FALSE, TRUE)) {
    f <- star_fixture(paired = TRUE, nested = nested)
    star_run(f, TRUE, paired.end = TRUE, include.subfolders = if (nested) "y" else "n",
             steps = "tr-co-ge", base.correction = TRUE, max.multimap = 7,
             mismatches = 2, keep.contaminants = TRUE, keep.index.in.memory = "noShared")
    calls <- star_calls(f)
    expect_identical(vapply(calls, `[[`, "", "tool"), c("fastp", "STAR", "STAR"))
    expect_true("--correction" %in% calls[[1]]$args)
    expect_identical(star_arg(calls[[2]], "--readFilesIn"),
                     file.path(f$output, "trim/trimmed_sample_1.fastq"))
    expect_identical(star_arg(calls[[3]], "--readFilesIn"),
                     file.path(f$output, "contaminants_depletion/contaminants_sample_1_Unmapped.out.mate1"))
    expect_true(file.path(f$output, "contaminants_depletion/contaminants_sample_1_Unmapped.out.mate2") %in% calls[[3]]$args)
    expect_identical(star_arg(calls[[3]], "--outFilterMultimapNmax"), "7")
    expect_identical(star_arg(calls[[3]], "--outFilterMismatchNmax"), "2")
    expect_identical(star_arg(calls[[2]], "--outSAMtype"), "BAM")
    expect_true(all(vapply(calls[-1], star_arg, "", "--genomeLoad") == "NoSharedMemory"))
    expect_true(file.exists(file.path(f$output, "aligned/LOGS/sample_1_Log.out")))
  }
})

test_that("base correction is optional and validated for both entry points", {
  for (folder in c(FALSE, TRUE)) {
    f <- star_fixture(paired = TRUE)
    extra <- if (folder) list(paired.end = TRUE) else list()
    do.call(star_run, c(list(f, folder, steps = "tr"), extra))
    expect_false("--correction" %in% star_calls(f)[[1]]$args)
    do.call(star_run, c(list(f, folder, steps = "tr", base.correction = TRUE), extra))
    expect_true("--correction" %in% star_calls(f)[[2]]$args)
    expect_error(do.call(star_run, c(list(f, folder, steps = "ge", base.correction = TRUE), extra)), "trimming")
    f <- star_fixture()
    expect_error(star_run(f, folder, steps = "tr", base.correction = TRUE), "paired end")
    expect_error(star_run(f, folder, steps = "tr", base.correction = NA), "TRUE or FALSE")
    expect_error(star_run(f, folder, steps = "ge", allow.introns = NA))
  }
})

test_that("tool failures stop the pipeline and reach the R caller", {
  for (folder in c(FALSE, TRUE)) for (tool in c("fastp", "STAR")) {
    f <- star_fixture()
    withr::local_envvar(c(ORFIK_TEST_FAIL = tool))
    expect_error(star_run(f, folder, steps = "tr-co-ge"), "alignment step failed")
    expect_length(star_calls(f), if (tool == "fastp") 1L else 2L)
    expect_false(file.exists(file.path(f$output, "runCommand.log")))
  }
})

test_that("resuming uses existing intermediate reads and does not rerun trimming", {
  f <- star_fixture()
  star_run(f, TRUE, steps = "tr-co-ge")
  unlink(f$record)
  star_run(f, TRUE, steps = "tr-co-ge", resume = "co")
  expect_length(star_calls(f), 2)
  expect_true(all(vapply(star_calls(f), `[[`, "", "tool") == "STAR"))
  unlink(f$record)
  star_run(f, FALSE, steps = "tr-co-ge", resume = "ge")
  expect_length(star_calls(f), 1)
})

test_that("NoSharedMemory applies to every library and empty folders fail", {
  f <- star_fixture()
  file.copy(f$reads, file.path(f$input, "second.fastq"))
  star_run(f, TRUE, steps = "ge", keep.index.in.memory = "noShared")
  expect_true(all(vapply(star_calls(f), star_arg, "", "--genomeLoad") == "NoSharedMemory"))
  unlink(list.files(f$input, full.names = TRUE))
  expect_error(star_run(f, TRUE, steps = "ge"), "alignment step failed")
})

test_that("tRNA input compression is resolved from the correct preceding stage", {
  f <- star_fixture()
  star_run(f, steps = "tr-tR-ge")
  expect_identical(star_arg(star_calls(f)[[2]], "--readFilesCommand"), "-")
  expect_identical(star_arg(star_calls(f)[[3]], "--readFilesIn"),
                   file.path(f$output, "tRNA_depletion/tRNA_sample_Unmapped.out.mate1"))
})

test_that("step logs identify resolved input, output and executed command", {
  f <- star_fixture(paired = TRUE)
  script <- system.file("STAR_Aligner", "RNA_Align_pipeline_folder.sh", package = "ORFik")
  single <- system.file("STAR_Aligner", "RNA_Align_pipeline.sh", package = "ORFik")
  args <- c(script, "-f", f$input, "-o", f$output, "-g", f$index,
            "-s", "tr-co-ge", "-p", "yes", "-I", single,
            "-S", f$bins[1], "-P", f$bins[2])
  output <- system2("bash", shQuote(args), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"))
  expect_true(any(grepl(paste0("  Input: ", f$output, "/trim/trimmed_sample_1.fastq"), output, fixed = TRUE)))
  expect_true(any(grepl(paste0("  Input: ", f$output, "/contaminants_depletion/contaminants_sample_1_Unmapped.out.mate1"), output, fixed = TRUE)))
  expect_equal(sum(grepl("  Command:", output, fixed = TRUE)), 3)
  expect_true(any(grepl("Splicing: indexed and novel", output, fixed = TRUE)))
})

test_that("homogeneous layout vectors from study metadata work and mixed layouts fail", {
  f <- star_fixture()
  star_run(f, TRUE, steps = "ge", paired.end = c(FALSE, FALSE))
  expect_length(star_calls(f), 1L)
  expect_error(star_run(f, TRUE, steps = "ge", paired.end = c(FALSE, TRUE)), "one layout")
  expect_error(star_run(f, TRUE, steps = "ge", paired.end = NA), "one layout")
})

test_that("bad step order and duplicate output names fail before tools run", {
  f <- star_fixture()
  expect_error(star_run(f, steps = "ge-tr"), "alignment step failed")
  expect_false(file.exists(f$record))
  expect_error(star_run(f, steps = "tr-tr"), "alignment step failed")
  expect_false(file.exists(f$record))
  dir.create(file.path(f$input, "nested"))
  file.copy(f$reads, file.path(f$input, "nested", basename(f$reads)))
  expect_error(star_run(f, TRUE, steps = "ge", include.subfolders = "y"), "alignment step failed")
  expect_false(file.exists(f$record))
})
