# Standalone Rcpp/pointer/owned-array comparison; no ORFik data or R package needed.
root<-'benchmark/transcript_mapping'
out<-Sys.getenv('ORFIK_MAPPING_OUTPUT',file.path(root,'results/storage'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(root,'storage.cpp'),cacheDir=file.path(out,'rcpp'))
timings<-list()
for(exons in c(3L,20L,200L)) {
 s<-seq.int(100L,by=30L,length.out=exons)
 args<-list(qs=rep.int(tail(s,1L)+1L,100000L),qe=rep.int(tail(s,1L)+5L,100000L),ts=s,te=s+19L,indices=rep.int(1L,exons),direction='+',offsets=rep.int(0L,100000L))
 expected<-do.call(candidate_storage,c(args,list(mode=0L)))
 for(mode in 1:3) {
  stopifnot(identical(expected,do.call(candidate_storage,c(args,list(mode=mode)))))
  for(i in 1:5)for(method in if(i%%2)c(0L,mode)else c(mode,0L)) {
   gc();elapsed<-system.time(for(j in 1:100)invisible(do.call(candidate_storage,c(args,list(mode=method)))),gcFirst=FALSE)['elapsed']/100
   timings[[length(timings)+1L]]<-data.frame(exons=exons,comparison=mode,method=method,iteration=i,batch=100L,seconds=elapsed)
  }
 }
}
d<-do.call(rbind,timings);write.csv(d,file.path(out,'timings.csv'),row.names=FALSE)
s<-aggregate(seconds~exons+comparison+method,d,median);write.csv(s,file.path(out,'summary.csv'),row.names=FALSE);print(s)
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
