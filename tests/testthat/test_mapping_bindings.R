test_that("native factories construct valid outputs and never mutate inputs", {
  tx <- GRangesList(a=GRanges("chr1",IRanges(c(100,200),width=20),"+"),
                    b=GRanges("chr2",IRanges(c(400,300),width=20),"-"))
  query <- GRangesList(first=tx[[1]],second=tx[[2]])
  mcols(tx)$annotation <- c("A","B")
  metadata(query) <- list(source="fixture")
  snapshot <- serialize(list(tx,query,ORFik:::.pmap_templates()),NULL)
  result <- pmapToTranscriptF(query,tx)
  expect_true(validObject(result,complete=TRUE))
  expect_identical(as.integer(start(unlist(result))),c(1L,1L))
  expect_identical(as.integer(end(unlist(result))),c(40L,40L))
  expect_identical(names(result),names(query))
  for(drop in c(FALSE,TRUE)) {
    x <- IRanges(c(16,16),c(25,25));names(x)<-c("2","1")
    mapped <- pmapFromTranscriptF(x,tx,drop)
    expect_true(validObject(mapped,complete=TRUE))
    expect_identical(names(mapped),c("b","a"))
    expect_identical(as.character(seqnames(unlist(mapped))),c("chr2","chr2","chr1","chr1"))
  }
  expect_identical(serialize(list(tx,query,ORFik:::.pmap_templates()),NULL),snapshot)
})

test_that("binding validation rejects malformed slots and flags before allocation", {
  tx <- GRangesList(tx=GRanges("chr1",IRanges(100,120),"+"))
  for(flag in list(NA,integer(),TRUE[FALSE],c(TRUE,FALSE),1,"yes")) {
    expect_error(pmapToTranscriptF(tx,tx,ignore.strand=flag),"logical")
    expect_error(pmapFromTranscriptF(IRanges(1,5),tx,removeEmpty=flag),"logical")
  }
  for(name in c("0","-1","2","abc",NA_character_)) {
    x<-IRanges(1,5);names(x)<-name
    expect_error(pmapFromTranscriptF(x,tx),"index")
  }
  bad<-tx;bad@partitioning@end<-2L
  expect_error(pmapToTranscriptF(bad,tx),"partition")
  bad<-tx;bad@unlistData@seqnames@values<-structure(2L,levels="chr1",class="factor")
  expect_error(pmapToTranscriptF(bad,tx),"sequence")
})

test_that("pure interval merging preserves empty intervals and strand separation", {
  tx<-GRangesList(tx=GRanges("chr1",IRanges(c(100,200),width=20),"+"))
  x<-GRangesList(q=GRanges("chr1",IRanges(c(90,100,205),width=3),"+"))
  mapped<-pmapToTranscriptF(x,tx)
  expect_true(validObject(mapped,complete=TRUE))
  expect_identical(as.integer(start(mapped[[1]])),c(1L,26L,0L))
  expect_identical(as.integer(end(mapped[[1]])),c(3L,28L,-1L))
  expect_identical(as.character(strand(mapped[[1]])),c("+","+","*"))
  mapped<-pmapToTranscriptF(ranges(x),tx)
  expect_identical(as.integer(start(mapped[[1]])),c(0L,1L,26L))
  expect_identical(as.integer(end(mapped[[1]])),c(-1L,3L,28L))
})
