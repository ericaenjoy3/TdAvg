#!/usr/bin/env Rscript

###
# Given a query bed file, the query and target config files for matrix series
#
# EL
# Created on 11 July, 2017
# Last modified on 11 July, 2017
###

suppressPackageStartupMessages(library("argparse"))
library(TdAvg)

parser <- ArgumentParser()
parser$add_argument("--chip", type="character", nargs=1, required=TRUE, help="a chipseq bed file")
parser$add_argument("--fs", type="character", nargs="+", required=TRUE, help="a configuration file")
parser$add_argument("--orderconfig", type="character",nargs=1,required=TRUE,help="the configuration file for ordering clusters in the bed, named like '*_order.config'")
parser$add_argument("--repconfig", type="character",nargs="+",required=TRUE,help="the configuration file for grouping replicates")
parser$add_argument("--tf", type="character",nargs=1,required=TRUE,help="the [transcription] factor to be priorised first")
parser$add_argument("--w", type="integer", nargs="+", required=TRUE, help="upstream/downstream window size for each BW.config")
parser$add_argument("--s", type="integer", nargs="+", required=TRUE, help="bin size for each BW.config")
parser$add_argument("--nms", type="character", nargs="+", required=TRUE, help="output directory with file pattern")
parser$add_argument("--kmeans", action="store_true", help="kmeans")
parser$add_argument("--noInd", action="store_true", help="do not plot individual avgplot or TD plots")
parser$add_argument("--noTD", action="store_true", help="do not plot TD plots")
parser$add_argument("--reorder", action="store_true", help="reorder within each clus by signal strength")
parser$add_argument("--rmSmallClus", action="store_true", help="remove cluster size <1000")
args <- parser$parse_args()
stopifnot(length(args$nms)==length(args$fs))
stopifnot(length(args$w)==length(args$s) & length(args$fs)==length(args$w))

# update fs, repconfig, w, s, nms
idx <- c(grep(args$tf, args$nms, ignore.case = TRUE), grep(args$tf, args$nms, ignore.case = TRUE, invert = TRUE))
invisible(sapply(names(args)[sapply(args,length)>1],
  function(m,idx){if(length(args[[m]])>1){args[[m]]<<-args[[m]][idx]}}, idx = idx))

# process chip bed
chip.obj <- chipConst(args$chip)
# create info object
# change cluster order by orderconfig
obj.list <- orderClus(chip.obj,args$orderconfig)
chip.obj <- obj.list$chip.obj
info.obj <- obj.list$info.obj

if (args$rmSmallClus){
  info.obj <- rmSmallClus(info.obj, chip.obj)
}

# process config
for (i in seq_along(args$fs)) {
  matlist.obj<-matlistConst(args$fs[i])
  # add w, s and points info to info.obj
  info.obj <- appendInfo(info.obj, w = as.integer(args$w)[i], s = as.integer(args$s)[i])
  if (i==1) {
    # remove outliers
    info.obj <- outliers(info.obj, chip.obj, matlist.obj)
    # kmeans clustering
    if (isTRUE(args$kmeans)) {
      obj.list <- clusterR(chip.obj, info.obj, matlist.obj, args$nms[i])
      chip.obj <- obj.list$chip.obj
      info.obj <- obj.list$info.obj
    }
    # reorder
    info.obj <- reord(info.obj, chip.obj, matlist.obj)
    # update chip.obj
    chip.obj <- updateChip(chip.obj, info.obj)
    # output clustered bed
    writeChip(chip.obj, paste0(args$nms[i], "_Clus.bed"))
  }
  # update matlist.obj
  matlist.obj <- updateMatlist(matlist.obj, info.obj)
  # combine replicates
  nmatlist.obj <- mergeRep(matlist.obj,args$repconfig[i])
  if (!isTRUE(args$noInd)) {
    # plot average plot
    avgplot(matlist.obj, chip.obj, info.obj,
      pdffoutFe = paste0(args$nms[i], "_Avg.pdf"),
      pdffoutSm = paste0(args$nms[i], "_AvgT.pdf"))
    # plot heatmap
    if (!args$noTD) {
      tdplot(matlist.obj, chip.obj, info.obj,
        pngfout = paste0(args$nms[i], "_TD.png"))
    }
  }
  boxplot(nmatlist.obj, chip.obj, info.obj,
    pdffout = paste0(args$nms[i], "_COMBbox.pdf",
    fout = paste0(args$nms[i], "_COMBbox_pvalue.txt")
  avgplot(nmatlist.obj, chip.obj, info.obj,
    pdffoutFe = paste0(args$nms[i], "_COMBAvg.pdf"),
    pdffoutSm = paste0(args$nms[i], "_COMBAvgT.pdf"))
  tdplot(nmatlist.obj, chip.obj, info.obj,
    pngfout = paste0(args$nms[i], "_COMBTD.png"))
}
