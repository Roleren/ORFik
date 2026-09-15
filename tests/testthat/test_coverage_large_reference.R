# Large coordinate spans, tiny run-length-encoded objects: never expand these
# Rles with as.numeric()/as.vector() in a test.
large_reference_reads <- function() {
  x <- GAlignments(seqnames=Rle(factor(c("a","b","a","b"),levels=c("a","b","empty"))),
    pos=c(1L,101L,500L,1073741800L),
    cigar=c("22M","11M100N11M","24M","22M"),
    strand=Rle(factor(c("+","*","-","-"),levels=c("+","-","*"))))
  seqlengths(x) <- c(1073741824L,1073741824L,0L)
  mcols(x)$score <- c(2L,3L,5L,7L)
  x
}

test_that("coverage of a >2^31-base reference keeps chromosomes separately encoded", {
  x <- large_reference_reads()
  grl <- grglist(x)
  gr <- unlist(grl,use.names=FALSE)
  w <- rep(as.numeric(mcols(x)$score),lengths(grl))
  for(ignore in c(FALSE,TRUE)) for(size in list(NULL,1L,3L)) {
    got <- suppressMessages(covRleFromGR(x,ignore.strand=ignore,chunk.size=size))
    expect_s4_class(f(got),"SimpleRleList")
    expect_equal(sum(as.double(lengths(f(got)))),2^31)
    expect_identical(seqinfo(got),seqinfo(x))
    for(st in if(ignore) "both" else c("+","-")) {
      keep <- if(ignore) rep(TRUE,length(gr)) else as.character(strand(gr)) %in% c(st,"*")
      expected <- coverage(gr[keep],weight=w[keep])
      actual <- if(st=="-") r(got) else f(got)
      expect_equal(as.list(actual),as.list(expected))
      expect_lt(sum(vapply(as.list(actual),function(z) length(runValue(z)),0L)),30L)
    }
  }
  # Unused chromosomes still contribute their full coordinate lengths.
  empty <- covRleFromGR(x[FALSE])
  expect_s4_class(f(empty),"SimpleRleList")
  expect_equal(sum(as.double(lengths(f(empty)))),2^31)
  expect_equal(sum(vapply(as.list(f(empty)),sum,0)),0)
})

test_that("large-reference covRleList conversion saves and reloads each read length", {
  x <- large_reference_reads()
  root <- tempfile("coverage-large-reference-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  path <- file.path(root,"all.ofst")
  fst::write_fst(data.frame(seqnames=as.character(seqnames(x)),start=start(x),
                            strand=as.character(strand(x)),cigar=cigar(x),score=mcols(x)$score),path)
  for(format in c("rds","qs")) for(force_chunks in c(FALSE,TRUE)) {
    withr::local_options(ORFik.coverage.chunk.size=if(force_chunks) 1L else NULL)
    out <- file.path(root,paste0(format,if(force_chunks) "-chunks" else "-direct"))
    merged <- paste0(out,"-merged")
    suppressMessages(convert_to_covRleList(NULL,in_files=path,out_dir=out,
      out_dir_merged=merged,seq_info=seqinfo(x),format=format,verbose=FALSE))
    reader <- if(format=="rds") readRDS else qs2::qs_read
    filename <- paste0("all.cov",format)
    got <- reader(file.path(out,filename))
    expect_s4_class(got,"covRleList")
    expect_equal(as.integer(names(got@list)),c(22L,24L))
    for(i in seq_along(got@list)) {
      expected <- covRleFromGR(x[readWidths(x)==as.integer(names(got@list)[i])])
      expect_equal(as.list(f(got@list[[i]])),as.list(f(expected)))
      expect_equal(as.list(r(got@list[[i]])),as.list(r(expected)))
      expect_identical(seqinfo(got@list[[i]]),seqinfo(x))
    }
    all <- reader(file.path(merged,filename))
    expect_s4_class(f(all),"SimpleRleList")
    expect_equal(as.list(f(all)),as.list(f(covRleFromGR(x))))
  }
})

test_that("compressed coverage is retained up to the exact cumulative endpoint limit", {
  x <- large_reference_reads()[1:3]
  for(sizes in list(c(1000L,1000L,0L),c(1073741824L,1073741823L,0L))) {
    seqlengths(x) <- sizes
    got <- covRleFromGR(x)
    expect_s4_class(f(got),"CompressedRleList")
    expect_s4_class(r(got),"CompressedRleList")
    expect_equal(as.double(lengths(f(got))),as.double(sizes))
  }
  seqlengths(x) <- c(1073741824L,1073741824L,0L)
  expect_s4_class(f(covRleFromGR(x)),"SimpleRleList")
})
