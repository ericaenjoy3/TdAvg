#####rlang .data prevents R CMD check from giving a NOTE about undefined global variables
#' @import ClusterR
#' @import ComplexHeatmap
#' @import EnrichedHeatmap
#' @import RColorBrewer
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom grDevices dev.off png

#' @title An S4 class to represent bed file of chip-seq peaks.
#' @name chip-class
#' @rdname chip-class
#' @description Store cluster sets in a bed object.
#' @slot bed data.frame object
#' @exportClass chip
chip <- setClass(
  "chip",
  slots=c(bed="data.frame"),
  validity=function(object) {
    if(nrow(object@bed)<1) {
      return("Empty data.table was given.")
    }
    if(!is.factor(object@bed[,4])) {
      return("The 4th col of bed cannot be anything else but factor.")
    }
    if(!identical(colnames(object@bed),c("chr","start","end","clus"))) {
      return("Colnames for bed must be \"chr\",\"start\",\"end\", and \"clus\"")
    }
    return(TRUE)
  }
)

#' @title An S4 class to represent information of peaks to keep, to order, and window, step size and coordinates of coverage value sampled.
#' @name info-class
#' @rdname info-class
#' @description Store
#' index of peaks to keep (boolean kpt.idx),
#' order for plotting (integer ord.idx),
#' to use for return index, filter order of peaks (ord.idx) by whether to keep (kpt.idx)
#' width of up/down-stream window (integer w),
#' step size within the window (integer s),
#' coordinates of data collected with reference to window center (numeric points).
#' @slot kpt.idx logical object represeting ChIP-seq peaks to keep.
#' @slot ord.idx integer object represeting order of ChIP-seq peaks to plot in the Tornado plot.
#' @slot w integer object representing number of basepairs flanking the middle of intervals for plotting.
#' @slot s integer object representing every n'th basepairs to sample data values.
#' @slot points numeric object representing coordinates in basepairs with reference to the mid-point of interval for data collected.
#' @exportClass info
info<-setClass(
  "info",
  slots=c(kpt.idx="logical",ord.idx="integer",w="integer",s="integer",points="numeric"),
  validity=function(object){
    if(length(object@kpt.idx)==0) {
      return("kpt.idx is of length 0.")
    }
    if(length(object@ord.idx)==0) {
      return("ord.idx is of length 0.")
    }
    if(length(object@kpt.idx)>0 && length(object@ord.idx)>0 && sum(object@kpt.idx)!=length(object@ord.idx)) {
      return("TRUE kpt.idx does not equal in length to ord.idx.")
    }
    return(TRUE)
  }
)

#' @title An S4 class to represent coverage value and ylab for plotting
#' @name matlist-class
#' @rdname matlist-class
#' @description Store
#' coverage value of all samples in a list (list ll).
#' ylal for plotting.
#' @slot ll list object representing a list of samples' coverage value sampled.
#' @slot ylab expression object for plotting.
#' @exportClass matlist
matlist<-setClass(
  "matlist",
  slots=c(ll="list",ylab="expression"),
  validity=function(object) {
    if(length(object@ll)<1) {
      return("ll was unintialized.")
    }
    if(length(object@ylab)!=1) {
      return("ylab was unintialized.")
    }
    if (any(sapply(object@ll,function(mat)any(dim(mat)==0)))){
      return("Matrices in the list have 0 columns or rows.")
    }
    if(is.null(names(object@ll))) {
      return("Names was not set for list of matrices given.")
    }
    return(TRUE)
  }
)