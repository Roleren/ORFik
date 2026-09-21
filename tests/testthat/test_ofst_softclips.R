softclip_fixture <- function() data.table::data.table(seqnames="1",strand="-",start=3051114L,
  cigar=c("1S25M","20M","21M","22M","22M3S","23M","23M3S","24M","24M1S","24M2S"),
  score=c(3L,94L,265L,12L,1L,27L,1L,117L,5L,2L))

softclip_files <- function(fun) {
  root <- tempfile("softclip-test-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  paths <- file.path(root,c("a.ofst","b.ofst"))
  for(p in paths) fst::write_fst(softclip_fixture(),p)
  fun(paths,root)
}

test_that("soft clip removal collapses before ordinary and rescue merges without modifying inputs", {
  softclip_files(function(paths,root) {
    hash <- tools::md5sum(paths)
    for(keep in c(FALSE,TRUE)) for(rescue in c(FALSE,TRUE)) {
      got <- suppressMessages(ofst_merge(paths,lib_names=c("a","b"),keep_all_scores=keep,
        remove_softclips=TRUE, filter_target_rows=if(rescue) 6 else 100,
        allow_filtering=FALSE,filter_chunk_rows=14,filter_tmpdir=root))
      expect_equal(nrow(got),6L)
      expect_equal(got$score[match(c("20M","21M","22M","23M","24M","25M"),got$cigar)],
                   2*c(94L,265L,13L,28L,124L,3L))
      if(keep) {expect_equal(got$a,got$score/2);expect_equal(got$b,got$score/2)}
      expect_null(attr(got,"removal_summary"))
    }
    unchanged <- suppressMessages(ofst_merge(paths,keep_all_scores=FALSE))
    expect_equal(nrow(unchanged),10L)
    expect_true(any(grepl("S",unchanged$cigar)))
    expect_identical(tools::md5sum(paths),hash)
    split <- suppressMessages(ofst_merge(paths,keep_all_scores=FALSE,remove_softclips=TRUE,
                                          dt_max_index_size=13,filter_target_rows=100))
    expect_equal(nrow(split),6L)
  })
})

test_that("reference coordinates and coverage stay fixed on both strands", {
  cg <- c("1S25M","2H3S10M2I5M3D4N2=1X4S1H","5H3S10M2S4H","24M2S")
  for(st in c("+","-")) {
    dt <- data.table::data.table(seqnames="1",strand=st,start=100L,cigar=factor(cg),score=1L,
                                qwidth=ORFik:::cigarWidthAlongQuerySpace_compat(cg))
    old <- ORFik:::getGAlignments(data.table::copy(dt))
    clean <- ORFik:::.ofst_remove_softclips(data.table::copy(dt))
    new <- ORFik:::getGAlignments(clean)
    expect_identical(start(old),start(new));expect_identical(end(old),end(new))
    expect_identical(width(old),width(new));expect_identical(njunc(old),njunc(new))
    expect_equal(qwidth(new),ORFik:::cigarWidthAlongQuerySpace_compat(cg,after.soft.clipping=TRUE))
    expect_equal(clean$qwidth,qwidth(new))
    expect_equal(clean$cigar,c("25M","2H10M2I5M3D4N2=1X1H","5H10M4H","24M"))
    seqlengths(old) <- seqlengths(new) <- 500L
    expect_equal(coverage(old),coverage(new))
    expect_identical(as.character(dt$cigar),cg)
  }
})

test_that("paired CIGARs, no CIGAR, empty inputs and malformed flags are handled", {
  pair <- data.table::data.table(seqnames="1",strand="+",start1=10L,start2=30L,
                                 cigar1="2S10M",cigar2="12M3S",score=1L)
  got <- ORFik:::.ofst_remove_softclips(data.table::copy(pair))
  expect_equal(got$cigar1,"10M");expect_equal(got$cigar2,"12M")
  expect_equal(got$start1,pair$start1);expect_equal(got$start2,pair$start2)
  no_cigar <- softclip_fixture()[,cigar:=NULL]
  expect_identical(ORFik:::.ofst_remove_softclips(data.table::copy(no_cigar)),no_cigar)
  expect_equal(nrow(ORFik:::.ofst_remove_softclips(softclip_fixture()[0])),0L)
  bad <- softclip_fixture();bad$cigar[1] <- "5M2S5M"
  expect_error(ORFik:::.ofst_remove_softclips(bad),"internal soft clip")
  bad$cigar[1] <- "5S"
  expect_error(ORFik:::.ofst_remove_softclips(bad),"only soft clips")
  softclip_files(function(paths,root) {
    for(flag in list(NA,1,"TRUE",c(TRUE,FALSE),matrix(TRUE)))
      expect_error(ofst_merge(paths,remove_softclips=flag),"remove_softclips")
    fst::write_fst(pair,paths[1]);fst::write_fst(pair,paths[2])
    got <- suppressMessages(ofst_merge(paths,remove_softclips=TRUE,keep_all_scores=FALSE))
    expect_equal(got$cigar1,"10M");expect_equal(got$cigar2,"12M");expect_equal(got$score,2L)
  })
})

test_that("mergeLibs forwards softclip removal to the master output", {
  softclip_files(function(paths,root) {
    df <- ORFik.template.experiment();df <- df[df$libtype=="RNA",][1:2,]
    out <- file.path(root,"merged")
    suppressMessages(mergeLibs(df,out_dir=out,paths=paths,lib_names_full=c("a","b"),
                               keep_all_scores=FALSE,remove_softclips=TRUE))
    got <- fst::read_fst(file.path(out,"all.ofst"))
    expect_equal(nrow(got),6L);expect_false(any(grepl("S",got$cigar)))
  })
})

test_that("filtered per-input statistics use the same stripped keys", {
  softclip_files(function(paths,root) {
    a <- softclip_fixture()
    b <- data.table::copy(a);b$start <- b$start + 100L
    a <- data.table::rbindlist(list(a,b))
    for(p in paths) fst::write_fst(a,p)
    got <- suppressMessages(ofst_merge(paths,keep_all_scores=FALSE,remove_softclips=TRUE,
      filter_target_rows=6,filter_chunk_rows=28,filter_tmpdir=root,filter_input_summary=TRUE))
    s <- attr(got,"removal_summary")
    expect_true(s$remove_softclips)
    expect_equal(s$rows_before,12)
    expect_equal(s$by_input$rows_before,c(12,12))
    expect_equal(s$by_input$rows_after,c(6,6))
    expect_false(any(grepl("S",got$cigar)))
  })
})
