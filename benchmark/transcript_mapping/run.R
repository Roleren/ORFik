# Run from the ORFik repository root after an optimized R CMD INSTALL --preclean.
devtools::load_all('.', compile=FALSE)
root <- 'benchmark/transcript_mapping'
out <- Sys.getenv('ORFIK_MAPPING_OUTPUT', file.path(root,'results'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ns<-asNamespace('ORFik')
baseline<-new.env(parent=ns)
source(file.path(root,'baseline/wrappers.R'),local=baseline)
for(name in c('pmapToTranscripts','pmapFromTranscripts')) {
 text<-paste(readLines(file.path(root,'baseline',paste0(name,'.cpp'))),collapse='\n')
 text<-gsub(paste0(name,'CPP'),paste0('baseline_',name,'CPP'),text,fixed=TRUE)
 source_path<-file.path(out,paste0(name,'.cpp'));writeLines(text,source_path)
 Rcpp::sourceCpp(source_path,env=baseline,cacheDir=file.path(out,'baseline_rcpp'))
 baseline[[paste0(name,'CPP')]]<-baseline[[paste0('baseline_',name,'CPP')]]
}
Rcpp::sourceCpp(file.path(root,'storage.cpp'),cacheDir=file.path(out,'storage_rcpp'))
timings<-list()
measure<-function(label,original,edited,args,batch=30L) {
 a<-do.call(original,args);b<-do.call(edited,args)
 stopifnot(isTRUE(all.equal(a,b,check.attributes=TRUE,tolerance=0)),validObject(b,complete=TRUE))
 for(i in 1:5)for(method in if(i%%2)c('original','edited')else c('edited','original')) {
  f<-if(method=='original')original else edited
  gc();elapsed<-system.time(for(j in seq_len(batch))invisible(do.call(f,args)),gcFirst=FALSE)[['elapsed']]/batch
  timings[[length(timings)+1L]]<<-data.frame(input=label,method=method,iteration=i,batch=batch,seconds=elapsed)
 }
 d<-do.call(rbind,timings);write.csv(d,file.path(out,'timings.csv'),row.names=FALSE)
 s<-aggregate(seconds~input+method,d,median);write.csv(s,file.path(out,'summary.csv'),row.names=FALSE)
 print(s[s$input==label,])
}
cache<-Sys.getenv('DCS_COVERAGE_CACHE',path.expand('~/livemount/forks/dominant_cell_states/.coverage_cache'))
real<-lapply(c('ATF4','FMR1','BOD1L1'),function(gene) {
 file<-list.files(cache,pattern=paste0('^',gene,'__.*\\.qs$'),full.names=TRUE)[1]
 if(is.na(file))stop('Missing external cached annotation for ',gene)
 s<-qs2::qs_read(file)$signature
 tx<-GRangesList(GRanges(s$seqnames,IRanges(s$start,s$end),s$strand));names(tx)<-s$transcript_names;tx
})
if(Sys.getenv("ORFIK_MAPPING_VALIDATE_ONLY") != "1") {
for(n in c(1L,100L,2000L)) {
 tx<-rep(suppressWarnings(do.call(c,real)),n);batch<-if(n==1L)100L else if(n==100L)30L else 5L
 measure(paste0('to_real_',length(tx)),baseline$pmapToTranscriptF,ORFik::pmapToTranscriptF,list(x=tx,transcripts=tx),batch)
 w<-widthPerGroup(tx,FALSE);x<-IRanges(pmax(1L,w-10L),w)
 measure(paste0('from_real_',length(tx)),baseline$pmapFromTranscriptF,ORFik::pmapFromTranscriptF,list(x=x,transcripts=tx,removeEmpty=TRUE),batch)
}
for(exons in c(3L,20L,200L)) {
 s<-seq.int(100L,by=30L,length.out=exons)
 args<-list(qs=rep.int(tail(s,1L)+1L,100000L),qe=rep.int(tail(s,1L)+5L,100000L),ts=s,te=s+19L,indices=rep.int(1L,exons),direction='+',offsets=rep.int(0L,100000L))
 pointer<-function(...)candidate_storage(...,mode=0L)
 for(mode in c(1L,2L,3L)) {
  f<-local({m<-mode;function(...)candidate_storage(...,mode=m)})
  measure(paste0('storage_',exons,'_mode',mode),pointer,f,args,batch=100L)
 }
 native_args<-args;names(native_args)<-c('xStart','xEnd','transcriptStart','transcriptEnd','indices','direction','nExonCumSum')
 measure(paste0('native_to_',exons),baseline$pmapToTranscriptsCPP,ORFik:::pmapToTranscriptsCPP,native_args,100L)
}
}
source(file.path(root,'validate.R'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
cat('All benchmark output comparisons passed.\n')
