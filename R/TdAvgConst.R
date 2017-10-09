#' @include TdAvgClass.R
#' @title construct an \code{chip} object
#' @name chipConst
#' @rdname chipConst-mehtods
#' @description construct an \code{chip} object from args$chip
#' @param chipF a ChIP-seq bed file
#' @return A \code{chip} object
#' @export chipConst
chipConst <- function(chipF) {
  dat <- data.frame(fread(chipF, header = FALSE,
    col.names = c("chr", "start", "end", "clus"))[, 1:4, with=FALSE]);
  dat[,4] <- as.character(dat[,4])
  return(new("chip", bed = dat))
}

#' @title construct an \code{matlist} object
#' @name matlistConst
#' @rdname matlistConst-methods
#' @description construct an \code{matlist} object according to a BW.config file.
#' @param config A BW.config file
#' @return A \code{matlist} object
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
    expression(paste(log[2],"(IP/INPUT)"))
  }
  return(new("matlist",ll=mat.list,ylab=ylab))
}
