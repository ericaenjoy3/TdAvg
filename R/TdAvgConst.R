#' @include TdAvgClass.R
#' @export chipConst
chipConst<-function(args) {
  dat<-data.frame(fread(args$chip,header=FALSE,col.names=c("chr","start","end","clus"))[,1:4,with=FALSE]);
  dat[,4]=factor(dat[,4])
  return(new("chip",bed=dat))
}

#' @export matlistConst
matlistConst<-function(config) {
  dat<-as.data.frame(fread(config,header=FALSE))
  mat.list<-lapply(dat[,3],function(f)as.matrix(fread(f,header=FALSE)))
  names(mat.list)<-dat[,2]
  ylab<-if (any(grepl("ATAC",dat[,1]))) {
    expression("Read Density (FPKM)")
  } else if (any(grepl("bowtie2",dat[,1]))) {
    expression("Read Density (RPKM)")
  } else {
    expression(paste0(log[2]," (IP/INPUT)"))
  }
  return(new("matlist",ll=mat.list,ylab=ylab))
}