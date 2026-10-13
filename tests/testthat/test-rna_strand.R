make_strand_reads <- function() {
 x <- GenomicAlignments::GAlignments(seqnames=Rle(factor(c('chr1','chr1'))),pos=c(11L,101L),cigar=c('5M20N5M','10M'),strand=Rle(factor(c('+','-'),levels=c('+','-','*'))))
 mcols(x)$score <- c(7L,3L);x
}
test_that('RNA reversal preserves alignment spans, splice CIGAR and weights', {
 x<-make_strand_reads();y<-normalize_rna_strand(x,'opposite')
 expect_identical(start(y),start(x));expect_identical(end(y),end(x));expect_identical(cigar(y),cigar(x));expect_identical(mcols(y),mcols(x));expect_identical(seqinfo(y),seqinfo(x))
 expect_equal(as.character(strand(y)),c('-','+'))
 expect_equal(as.character(strand(normalize_rna_strand(x,'alignment'))),c('+','-'))
})
test_that('OFST normalization survives serialization without another flip', {
 x<-normalize_rna_strand(make_strand_reads(),'opposite');file<-tempfile(fileext='.ofst');on.exit(unlink(c(file,ORFik:::rna_strand_sidecar(file))));export.ofst(x,file)
 st<-file.info(file);dir.create(dirname(ORFik:::rna_strand_sidecar(file)),recursive=TRUE,showWarnings=FALSE);saveRDS(list(normalized=TRUE,read_strand='opposite',product_signature=c(size=st$size,mtime=as.numeric(st$mtime))),ORFik:::rna_strand_sidecar(file))
 y<-fimport(file);expect_equal(as.character(strand(y)),c('-','+'));expect_identical(mcols(y)$score,mcols(x)$score)
 expect_equal(ORFik:::rna_strand_resolve(file),'alignment')
})
test_that('Unstranded weighted read counting accepts either gene orientation', {
 x<-normalize_rna_strand(make_strand_reads(),'unstranded');q<-GRanges('chr1',IRanges(c(11,101),width=10),strand=c('-','+'))
 expect_equal(as.numeric(countOverlapsW(q,x,'score')),c(7,3));expect_equal(as.numeric(countOverlapsW(q,x,'score',ignore.strand=FALSE)),c(0,0))
})
test_that('Stale and unknown provenance are rejected', {
 file<-tempfile();writeLines('x',file);sidecar<-ORFik:::rna_strand_sidecar(file);on.exit(unlink(c(file,sidecar)))
 dir.create(dirname(sidecar),recursive=TRUE,showWarnings=FALSE);saveRDS(list(status='unknown'),sidecar);expect_error(ORFik:::rna_strand_resolve(file),'unresolved')
 saveRDS(list(product_signature=c(size=0,mtime=0)),sidecar);expect_error(ORFik:::rna_strand_resolve(file),'stale')
})
test_that('Repeated collapse of one or two alignments never inflates weights', {
 x<-make_strand_reads()
 for(n in 1:2) {
  y<-collapseDuplicatedReads(x[seq_len(n)])
  z<-collapseDuplicatedReads(y)
  expect_equal(sum(mcols(y)$score),sum(mcols(x[seq_len(n)])$score))
  expect_identical(mcols(z)$score,mcols(y)$score)
 }
})

test_that('Opposite RNA orientation preserves strand runs, including unknown strands', {
 x<-rep(make_strand_reads(),c(10000L,20000L));strand(x)<-Rle(factor(c('+','-','*'),levels=c('+','-','*')),c(10000L,19000L,1000L))
 y<-normalize_rna_strand(x,'opposite')
 expect_identical(runLength(strand(y)),runLength(strand(x)))
 expect_equal(as.character(runValue(strand(y))),c('-','+','*'))
 expect_identical(mcols(y),mcols(x));expect_identical(start(y),start(x))
})
