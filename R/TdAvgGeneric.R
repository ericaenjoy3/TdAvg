#' @include TdAvgClass.R
#' @title order clusters in the 4th col of a bed file according to a configuration file
#' @name orderClus
#' @rdname orderClus-methods
#' @description order clusters in the 4th col of bed file based on a configuration file with 1 column of ordered clusters.
#' @param obj A \code{chip} object.
#' @param orderconfig An order configruation file (character).
#' @return A list of an level-order changed \code{chip} object and an \code{info} object.
#' @export orderClus
setGeneric(name="orderClus",
  def=function(obj,orderconfig){
    standardGeneric("orderClus")
  }
)

#' @title append width and step size to an \code{info} object
#' @name appendInfo
#' @rdname appendInfo-methods
#' @description append width and step size to an \code{info} object
#' @param obj An \code{info} object.
#' @param w The flanking window size in basepairs (integer).
#' @param s The step size for coverage value sampled (integer).
#' @return An updated \code{info} object
#' @export appendInfo
setGeneric(name="appendInfo",
  def=function(obj,w,s){
    standardGeneric("appendInfo")
  }
)

#' @title integer index (peak order) to keep of an \code{info} object
#' @name returnIntIdx
#' @rdname returnIntIdx-methods
#' @description Return integer index (peak order) that is supposed to be kept
#' @param obj An \code{info} object.
#' @return An integer index for ordered peaks to keep
#' @export returnIntIdx
setGeneric(name="returnIntIdx",
  def=function(obj){
    standardGeneric("returnIntIdx")
  }
)

#' @title subset by integer index
#' @name subsetbyIntIdx
#' @rdname subsetbyIntIdx-methods
#' @description Filter peaks to keep or discard.
#' Order of full/subset of original bed rows specified in the integer specified in IntIdx.
#' The function raises an error,
#' if the length of index vector is longer or equal to the number of peaks to be kpt.
#' The function also raises an error,
#' if the max integer position is larger than the lagest indexable postion of peaks to be kpt.
#' @param obj An \code{info} object
#' @param IntIdx An integer/numeric vector specifying index of currently kept peaks to be filtered out from the current point.
#' @param invert when TRUE, specifies selecting the rows, or otherwise FALSE, specifies the complement rows.
#' @return An updated \code{info} object
#' @export subsetbyIntIdx
setGeneric(name="subsetbyIntIdx",
  def=function(obj, IntIdx, invert){
    standardGeneric("subsetbyIntIdx")
  }
)

#' @title subset by boolean index
#' @name subsetbyBoolIdx
#' @rdname subsetbyBoolIdx-methods
#' @description Filter peaks to keep or discard.
#' Set certain TRUE values of kpt.idx of an \code{info} object to FALSE
#' identified by the boolean values specified in BoolIdx.
#' TRUE values in BoolIdx specify of those currently kept peaks to be kept downstream.
#' FLASE values in BoolIdx specify of those currently kept peaks to be filtered out downstream.
#' @param obj An \code{info} object
#' @param BoolIdx A boolean vector of the length as the currently kpt peaks
#' @param invert when TRUE, specifies selecting the rows, or otherwise FALSE, specifies the complement rows.
#' @return An updated \code{info} object
#' @export subsetbyBoolIdx
setGeneric(name="subsetbyBoolIdx",
  def=function(obj, BoolIdx, invert){
    standardGeneric("subsetbyBoolIdx")
  }
)

#' @title orderbyIntIdx
#' @name orderbyIntIdx
#' @rdname orderbyIntIdx-methods
#' @description Order peaks by currently given integer index,
#' which corresponds to the currently sub-indices of currently kept peaks.
#' @param obj An \code{info} object
#' @param IntIdx An integer vector of the same length as the currently kept peaks.
#' @return An updated \code{info} object
#' @export orderbyIntIdx
setGeneric(name="orderbyIntIdx",
  def=function(obj,IntIdx){
    standardGeneric("orderbyIntIdx")
  }
)

#' @title remove small (size<1000) clusters
#' @name rmSmallClus
#' @rdname rmSmallClus-methods
#' @description Set
#' Currently kept peaks (kpt.idx) to FALSE,
#' had the size of cluster is less than 1000.
#' @param info.obj An \code{info} object
#' @param chip.obj A \code{chip} object
#' @return An updated \code{info} object
#' @export rmSmallClus
setGeneric(name="rmSmallClus",
  def=function(info.obj,chip.obj){
    standardGeneric("rmSmallClus")
  }
)

#' @title set currently kept peaks to be discarded if any NA coverage value is reported in any individual samples.
#' @name outliers
#' @rdname outliers-methods
#' @description Set
#' Currently kept peaks (kpt.idx) to FALSE,
#' had coverage values for a peak reported to be NA in any samples.
#' @param info.obj An \code{info} object
#' @param chip.obj A \code{chip} object
#' @param matlist.obj A \code{matlist} object
#' @return An updated \code{info} object
#' @export outliers
setGeneric(name="outliers",
  def=function(info.obj,chip.obj,matlist.obj){
    standardGeneric("outliers")
  }
)

#' @title Update the label of peaks in \code{bed} slot of an \code{chip} object
#' @name updateClus
#' @rdname updateClus-methods
#' @description Update
#' the label of peaks (the 4th column) in \code{bed} slot of an \code{chip} object,
#' according to the integer position (IntIdx),
#' by the new lables (subclus).
#' @param chip.obj A \code{chip} object.
#' @param IntIdx An integer/numeric vector specify a subset or all of labels in the \code{bed} slot of \code{chip} object to be updated by given factor/character vector of subclus.
#' @param subclus A factor/character vector of the same length as IntIdx specifying new labels for peak clusters.
#' @return An updated \code{chip} object
#' @export updateClus
setGeneric(name="updateClus",
  def=function(chip.obj,IntIdx,subclus) {
    standardGeneric("updateClus")
  }
)

#' @title Kmeans subclustering within each clusters of peaks
#' @name clusterR
#' @rdname clusterR-methods
#' @description
#' Find the optimum number of subclusters within each clusters,
#' Kmeans sub-clustering within each clusters,
#' Update the cluster labels (the 4th column) of the \code{bed} slot of a \code{chip} object,
#' by labels of the combinations of cluster and subcluster labels.
#' @param chip.obj A \code{chip} object.
#' @param info.obj An \code{info} object.
#' @param matlist.obj A \code{matlist} object.
#' @param nms A character vector of length one, specifying the directory and the prefix for plotting the optimum number of subclusters within each cluster.
#' @return An updated \code{chip} object.
#' @export clusterR
setGeneric(name="clusterR",
  def=function(chip.obj,info.obj,matlist.obj,nms){
    standardGeneric("clusterR")
  }
)

#' @title Reorder peaks (specified info.obj) by the descending order of coverage signals of individual/grouped samples.
#' @name reord
#' @rdname reord-methods
#' @description Reorder peaks (specified info.obj) by the descending order of coverage signals of individual/grouped samples.
#' Update info.obj by the new peak order.
#' @param info.obj An \code{info} object
#' @param chip.obj A \code{chip} object
#' @param matlist.obj A \code{matlist} object
#' @return An updated \code{info} object
#' @export reord
setGeneric(name="reord",
  def=function(info.obj,chip.obj,matlist.obj){
    standardGeneric("reord")
  }
)

#' @title Update the \code{bed} slot of the \code{chip} object by the order and filteration specified in the \code{info} object
#' @name updateChip
#' @rdname updateChip-methods
#' @description Update the \code{bed} slot of the \code{chip} object,
#' by the order and filteration specified in the \code{info} object.
#' Once updated, the kpt.idx and ord.idx slots of the \code{info} object become obselete for the \code{chip} object.
#' @param chip.obj A \code{chip} object
#' @param info.obj An \code{info} object
#' @return An updated \code{chip} object
#' @export updateChip
setGeneric(name="updateChip",
  def=function(chip.obj,info.obj){
    standardGeneric("updateChip")
  }
)

#' @title Update the \code{ll} slot of the \code{matlist} object by the order and filteration specified in the \code{info} object
#' @name updateMatlist
#' @rdname updateMatlist-methods
#' @description Update the \code{ll} slot of the \code{mat} object,
#' by the order and filteration specified in the \code{info} object.
#' Once updated, the kpt.idx and ord.idx slots of the \code{info} object become obselete for the \code{matlist} object.
#' @param matlist.obj A \code{matlist} object
#' @param info.obj An \code{info} object
#' @return An updated \code{matlist} object
#' @export updateMatlist
setGeneric(name="updateMatlist",
  def=function(matlist.obj,info.obj) {
    standardGeneric("updateMatlist")
  }
)

#' @title Output the \code{bed} slot of the \code{chip} object to file.
#' @name writeChip
#' @rdname writeChip-methods
#' @description
#' Output the \code{bed} slot of the \code{chip} object to file.
#' @param chip.obj A \code{chip} object
#' @param fout A output text file
#' @return invisible of NULL
#' @export writeChip
setGeneric(name="writeChip",
  def=function(chip.obj,fout) {
    standardGeneric("writeChip")
  }
)

#' @title Merge within-group smaples' coverage signals
#' @name mergeRep
#' @rdname mergeRep-methods
#' @description
#' Merge within-group smaples' coverage signals, based on the specification in the repconfig file.
#' @param matlist.obj A \code{matlist} object.
#' @param repconfig A replicates configuration file (1st column: individual sample names as in the 2nd column of BW.config; 2nd column: group names; columns separated by tab).
#' @return A sample-grouped \code{matlist} object
#' @export mergeRep
setGeneric(name="mergeRep",
  def=function(matlist.obj,repconfig){
    standardGeneric("mergeRep")
  }
)

#' @title Plot average coverage signals stratified by clusters and by individual/groupled samples
#' @name avgplot
#' @rdname avgplot-methods
#' @description
#' Plot average coverage signals stratified by clusters and by individual/groupled samples.
#' @param matlist.obj A \code{matlist} object
#' @param chip.obj A \code{chip} object
#' @param info.obj An \code{info} object
#' @param pdffoutFe A pdf graph of average plot stratified by clusters
#' @param pdffoutSm A pdf graph of average plot stratified by individual/groupled samples
#' @return No object
#' @export avgplot
setGeneric(name="avgplot",
  def=function(matlist.obj, chip.obj, info.obj, pdffoutFe, pdffoutSm) {
    standardGeneric("avgplot")
  }
)

#' @title Plot a box plot
#' @name bplot
#' @rdname bplot-method
#' @description Plot a box plot
#' @param matlist.obj A \code{matlist} object
#' @param chip.obj A \code{chip} object
#' @param info.obj An \code{info} object
#' @param pdffoutFe A pdf graph of box/violin plot stratified by clusters
#' @param pdffoutSm A pdf graph of box/violin plot stratified by individual/groupled samples
#' @return No object
setGeneric(name = "bplot",
#' @export bplot
  def = function(matlist.obj, chip.obj, info.obj, pdffoutFe, pdffoutSm) {
    standardGeneric("bplot")
  }
)

#' @title Plot a tornado plot
#' @name tdplot
#' @rdname tdplot-method
#' @description Plot a tornado plot
#' @param matlist.obj A \code{matlist} object
#' @param chip.obj A \code{chip} object
#' @param info.obj An \code{info} object
#' @param pngfout An output png file for the tornado plot
#' @return No object
#' @export tdplot
setGeneric(name="tdplot",
  def=function(matlist.obj, chip.obj, info.obj, pngfout) {
    standardGeneric("tdplot")
  }
)
