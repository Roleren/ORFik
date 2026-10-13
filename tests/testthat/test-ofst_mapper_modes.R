mapper_mode_fixture <- function(collapsed = FALSE, paired = FALSE) {
  root <- tempfile('mapper_modes_'); dir.create(root)
  sam <- file.path(root, 'reads.sam')
  nm <- if (collapsed) c('seq1_x7','seq2_x3','seq3_x5','seq4_x2') else paste0('read',1:4)
  flags <- if (paired) c(99,147,99,147) else c(0,16,0,16)
  lines <- c('@HD\tVN:1.6\tSO:coordinate','@SQ\tSN:chr1\tLN:10000',
    vapply(seq_along(nm),function(i)paste(nm[i],flags[i],'chr1',10L*i,60,
      if (i==2L) '2S8M' else if(i==3L) '5M10N5M' else '10M',
      if(paired)'=' else '*',if(paired)10L*(if(i%%2L)i+1L else i-1L) else 0L,
      0L,'ACGTACGTAC','*',paste0('NH:i:',if(i<=2L)1L else 2L),sep='\t'),character(1)))
  writeLines(lines,sam);Rsamtools::asBam(sam,destination=file.path(root,'reads'),overwrite=TRUE)
  bam <- file.path(root,'reads.bam');stopifnot(file.exists(bam))
  df <- ORFik.template.experiment()[1, ];df@listData$filepath <- bam
  if (paired) {
    sc <- ORFik:::rna_strand_sidecar(bam);dir.create(dirname(sc),recursive=TRUE,showWarnings=FALSE)
    st<-file.info(bam);saveRDS(list(status='stranded',read_strand='opposite',selected_mate=1L,
      product_signature=c(size=st$size,mtime=as.numeric(st$mtime))),sc)
  }
  list(root=root,df=df,bam=bam)
}

test_that('One import produces identical all/unique OFST and weighted diagnostics', {
  for (collapsed in c(FALSE,TRUE)) for (orientation in c('alignment','opposite','unstranded')) {
    f <- mapper_mode_fixture(collapsed);on.exit(unlink(f$root,recursive=TRUE),add=TRUE)
    old <- setNames(file.path(f$root,c('old_all','old_unique')),c('all','unique'))
    new <- setNames(file.path(f$root,c('new_all','new_unique')),c('all','unique'));dist<-list()
    for(mode in names(old)) convert_bam_to_ofst(f$df,out_dir=old[[mode]],only_unique_mappers=mode=='unique',read_strand=orientation,verbose=FALSE)
    convert_bam_to_ofst_modes(f$df,out_dirs=new,read_strand=orientation,after_mode=function(mode,path,distribution)dist[[mode]]<<-distribution)
    for(mode in names(old)) {
      a<-file.path(old[[mode]],'reads.ofst');b<-file.path(new[[mode]],'reads.ofst')
      expect_identical(readBin(a,'raw',file.info(a)$size),readBin(b,'raw',file.info(b)$size))
      x<-fimport(a);expect_equal(sum(dist[[mode]]$count),if(collapsed)sum(mcols(x)$score) else length(x))
      expect_equal(ORFik:::rna_strand_resolve(b),if(orientation=='unstranded')'unstranded' else 'alignment')
    }
  }
})

test_that('Selected mate semantics and successful mode callbacks survive partial failure', {
  f<-mapper_mode_fixture(paired=TRUE);on.exit(unlink(f$root,recursive=TRUE))
  old<-file.path(f$root,'old');convert_bam_to_ofst(f$df,out_dir=old,only_unique_mappers=TRUE,verbose=FALSE)
  out<-c(all=file.path(f$root,'all'),unique=file.path(f$root,'unique'));seen<-character()
  expect_error(convert_bam_to_ofst_modes(f$df,out_dirs=out,after_mode=function(mode,path,distribution){seen<<-c(seen,mode);if(mode=='all')stop('callback failure')}),'callback failure')
  expect_identical(seen,'all');expect_true(file.exists(file.path(out[['all']],'reads.ofst')))
  convert_bam_to_ofst_modes(f$df,modes='unique',out_dirs=out)
  expect_identical(readBin(file.path(old,'reads.ofst'),'raw',1e6),readBin(file.path(out[['unique']],'reads.ofst'),'raw',1e6))
})
