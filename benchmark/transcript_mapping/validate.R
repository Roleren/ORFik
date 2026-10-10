# Run after run.R, in the same session: reuse independent pre-refactor bindings.
old <- list(pmapToTranscriptF=baseline$pmapToTranscriptF,pmapFromTranscriptF=baseline$pmapFromTranscriptF)
new_to <- ORFik::pmapToTranscriptF
new_from <- ORFik::pmapFromTranscriptF
templates <- ORFik:::.pmap_templates()
check <- function(f,g,args) {
 a<-do.call(f,args);b<-tryCatch(do.call(g,args),error=function(e){print(args);stop(e)})
 eq<-all.equal(a,b,check.attributes=TRUE,tolerance=0)
 if(!isTRUE(eq)){print(args);print(eq);print(a);print(b);stop('Binding mismatch')}
 stopifnot(validObject(b,complete=TRUE))
}
set.seed(611L)
checks<-0L
for(iter in seq_len(40L)) {
 n<-sample.int(6L,1L)
 tx<-GRangesList(lapply(seq_len(n),function(i) {
  k<-sample.int(12L,1L);w<-sample.int(25L,k,replace=TRUE)
  s<-cumsum(c(100L,head(w,-1L)+30L));direction<-sample(c('+','-'),1L)
  gr<-GRanges(paste0('chr',1L+i%%2L),IRanges(s,width=w),direction)
  if(direction=='-')gr<-rev(gr)
  gr
 }))
 names(tx)<-paste0('tx',seq_len(n))
 x<-GRangesList(lapply(tx,function(gr) {
  z<-gr[sample.int(length(gr),sample.int(length(gr),1L))];width(z)<-pmin(width(z),5L);z
 }))
 names(x)<-paste0('query',seq_len(n))
 before<-serialize(list(x,tx,templates),NULL)
 for(ignore in c(FALSE,TRUE))for(drop in c(FALSE,TRUE))for(lengths in c(FALSE,TRUE)) {
  check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx,ignore.strand=ignore,reduce.ranges=drop,set.seqlengths=lengths))
  checks<-checks+1L
 }
 w<-widthPerGroup(tx,FALSE);q<-IRanges(pmax(1L,w-3L),w)
 for(drop in c(FALSE,TRUE)) {
  check(old$pmapFromTranscriptF,new_from,list(x=q,transcripts=tx,removeEmpty=drop));checks<-checks+1L
  ix<-sample.int(n,10L,replace=TRUE);z<-q[ix];names(z)<-as.character(ix)
  check(old$pmapFromTranscriptF,new_from,list(x=z,transcripts=tx,removeEmpty=drop));checks<-checks+1L
 }
 stopifnot(identical(before,serialize(list(x,tx,templates),NULL)))
}
# Public schemas, recycling and sort flags.
for(tx in real) {
 for(x in list(tx,ranges(tx),unlist(tx,use.names=FALSE),ranges(unlist(tx,use.names=FALSE)))) {
  # Flat queries paired against a single transcript use reference recycling.
  check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx));checks<-checks+1L
 }
}
scalar<-GRanges('chr1',IRanges(c(100,200,300),width=20),c('+','-','+'))
check(old$pmapToTranscriptF,new_to,list(x=scalar,transcripts=scalar))
# Metadata and empty query groups must follow the established output schema.
x<-GRangesList(a=scalar[1],b=scalar[integer()],c=scalar[3]);tx<-GRangesList(scalar[1],scalar[2],scalar[3]);names(tx)<-c('a','b','c')
check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx))
mcols(x)$tag<-1:3;metadata(x)<-list(source='test');mcols(unlist(x)) # no mutation
check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx))
cat(checks,'full binding comparisons passed, plus scalar/empty/metadata checks.\n')

for(genomic in c(TRUE,FALSE)) {
 tx<-GRangesList(tx=GRanges('chr1',IRanges(c(100,200),width=20),'+'))
 for(s in list(c(90L,90L),c(90L,100L),c(100L,205L,90L))) {
  x<-GRangesList(q=GRanges('chr1',IRanges(s,width=3L),'+'))
  if(!genomic)x<-ranges(x)
  check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx))
 }
}
cat('Unmapped interval reduction comparisons passed.\n')
# Repeated exon coordinates on different chromosomes, optional sorting,
# unnamed references, unknown strands and named IRangesList queries.
tx<-GRangesList(a=GRanges('chr1',IRanges(c(200,100),width=20),'+'),
                b=GRanges('chr2',IRanges(c(100,200),width=20),'-'))
x<-tx
for(ignore in c(FALSE,TRUE))check(old$pmapToTranscriptF,new_to,
  list(x=x,transcripts=tx,ignore.strand=ignore,x.is.sorted=FALSE,tx.is.sorted=FALSE))
names(tx)<-NULL
check(old$pmapToTranscriptF,new_to,list(x=x,transcripts=tx))
tx<-GRangesList(a=GRanges('chr1',IRanges(c(200,100),width=20),'*'))
check(old$pmapToTranscriptF,new_to,list(x=tx,transcripts=tx))
tx<-GRangesList(a=GRanges('chr1',IRanges(c(100,200),width=20),'+'))
x<-IRangesList('1'=IRanges(c(2,25),c(5,28)))
check(old$pmapFromTranscriptF,new_from,list(x=x,transcripts=tx,removeEmpty=TRUE))
cat('Sorting, unnamed references, unknown strands and named IRangesList comparisons passed.\n')
