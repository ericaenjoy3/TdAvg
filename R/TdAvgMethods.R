#' @include TdAvgClass.R TdAvgGeneric.R
#' @rdname orderClus-methods
setMethod(f="orderClus",
  signature="chip",
  definition=function(obj, orderconfig){
    ord <- read.table(orderconfig, header = FALSE, as.is = TRUE)[,1]
    ord.idx <- order(match(as.character(obj@bed[,4]), ord))
    info.obj <- new("info", ord.idx = as.integer(ord.idx))
    obj@bed[,4] <- factor(as.character(obj@bed[,4]), levels=ord)
    return(list(chip.obj = obj, info.obj = info.obj))
  }
)

#' @rdname appendInfo-methods
setMethod(f = "appendInfo",
  signature = "info",
  definition = function(obj, w, s){
    obj@w <- w
    obj@s <- s
    obj@points <- seq(-w+s/2, w-s/2, s)
    return(obj)
  }
)

#' @rdname returnIntIdx-methods
setMethod(f = "returnIntIdx",
  signature = "info",
  definition = function(obj) {
    return(obj@ord.idx)
  }
)

#' @rdname subsetbyIdx-methods
setMethod(f = "subsetbyIntIdx",
  signature = c("info", "numeric", "logical"),
  definition = function(obj, IntIdx, invert) {
    stopifnot(length(obj@ord.idx) >= max(IntIdx))
    stopifnot(length(obj@ord.idx) >= length(IntIdx))
    obj@ord.idx <- if (invert) {
      obj@ord.idx[!seq_along(obj@ord.idx) %in% IntIdx]
    } else {
      obj@ord.idx[IntIdx]
    }
    return(obj)
  }
)

#' @rdname subsetbyBoolIdx-methods
setMethod(f = "subsetbyBoolIdx",
  signature = c("info", "logical", "logical"),
  definition = function(obj, BoolIdx, invert) {
    IntIdx <- which(BoolIdx)
    obj <- subsetbyIntIdx(obj, IntIdx, invert)
    return(obj)
  }
)

#' @rdname orderbyIntIdx-methods
setMethod(f="orderbyIntIdx",
  signature="info",
  definition=function(obj, IntIdx) {
    ori.idx <- returnIntIdx(obj)
    stopifnot(length(ori.idx) == length(IntIdx))
    obj@ori.idx <- IntIdex
    return(obj)
  }
)

#' @rdname rmSmallClus-methods
setMethod(f = "rmSmallClus",
  signature = c("info","chip"),
  definition = function(info.obj, chip.obj) {
    ori.idx <- returnIntIdx(info.obj)
    clus.list <- split(chip.obj@bed[ori.idx, 4], chip.obj@bed[ori.idx, 4])
    clus.nms <- names(clus.list)[sapply(clus.list, length) < 1000]
    int.idx <- which(as.character(chip.obj@bed[ori.idx, 4]) %in% clus.nms)
    info.obj <- subsetbyIntIdx(info.obj, int.idx, invert = TRUE)
    validObject(info.obj)
    return(info.obj)
  }
)

#' @rdname outliers-methods
setMethod(f = "outliers",
  signature = c("info", "chip", "matlist"),
  definition = function(info.obj, chip.obj, matlist.obj){
    idx <- returnIntIdx(info.obj)
    outliers <- rowSums(sapply(seq_along(matlist.obj@ll), function(i, matlist.obj, idx)
      rowSums(is.na(matlist.obj@ll[[i]][idx, ])), matlist.obj = matlist.obj, idx = idx))>0
    info.obj <- subsetbyBoolIdx(info.obj, !outliers, invert = FALSE)
    return(info.obj)
  }
)

#' @rdname updateClus-methods
setMethod(f = "updateClus",
  signature = c("chip"),
  definition = function(chip.obj, IntIdx, subclus) {
    stopifnot(length(IntIdx) == length(subclus))
    ss <- as.character(chip.obj@bed[,4])
    ss[IntIdx] <- as.character(subclus)
    ss[!seq_along(ss) %in% IntIdx] <- "NA"
    ss <- factor(ss, levels = unique(subclus))
    chip.obj@bed[,4] <- ss
    return(chip.obj)
  }
)

#' @rdname clusterR-methods
setMethod(f = "clusterR",
  signature = c("chip", "info", "matlist"),
  definition = function(chip.obj, info.obj, matlist.obj, nms) {
    idx <- returnIntIdx(info.obj)
    group <- chip.obj@bed[idx, 4]
    median.col <- floor(median(seq_len(ncol(matlist.obj@ll[[1]]))))
    mat.list <- lapply(matlist.obj@ll, function(mat, median.col, idx){
      mat[idx, median.col]
     }, median.col = median.col, idx = idx)
    mat <- do.call("cbind", mat.list)
    stopifnot(nrow(mat) == length(group))
    if (length(levels(group)) == 1) {
      mat.list <- list(mat);
      names(mat.list) <- levels(group)
    } else {
      mat.list <- split(mat, group)
    }
    mem.list <- lapply(seq_along(mat.list),
      function(i,nms){
        mat <- mat.list[[i]]
        pdf(paste0(nms, "_", names(mat.list)[i], ".pdf"))
        opt<-Optimal_Clusters_KM(mat, max_clusters = min(10,ncol(mat)),
          plot_clusters = TRUE, criterion = 'distortion_fK', fK_threshold = 0.85,
          initializer = 'optimal_init', tol_optimal_init = 0.2)
        dev.off()
        km_mb <- MiniBatchKmeans(mat, clusters = opt, batch_size = 20, num_init = 5,
          max_iters = 100, init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
          verbose = F)
        pr_mb <- predict_MBatchKMeans(mat, km_mb$centroids)
        return(as.numeric(pr_mb))
      }, nms = nms)
    mem <- unsplit(mem.list, group)
    labels <- paste(as.character(group), mem, sep = "_")
    clus <- factor(labels, levels = unique(labels))
    chip.obj <- updateClus(chip.obj, idx, clus)
    info.obj <- orderbyIntIdx(info.obj, order(as.numeric(chip.obj@bed[idx, 4]), mem))
    return(list(chip.obj, info.obj))
  }
)

#' @rdname reord-methods
setMethod(f="reord",
  signature=c("info","chip","matlist"),
  definition=function(info.obj, chip.obj, matlist.obj){
    idx <- returnIntIdx(info.obj)
    rowS <- rowSums(sapply(seq_along(matlist.obj@ll), function(i, matlist.obj, idx){
      rowSums(matlist.obj@ll[[i]][idx,])}, matlist.obj = matlist.obj, idx = idx)) # rowSums
    info.obj <- orderbyIntIdx(info.obj, order(as.numeric(chip.obj@bed[idx, 4]),-rowS))
    return(info.obj)
  }
)

#' @rdname updateChip-methods
setMethod(f="updateChip",
  signature=c("chip", "info"),
  definition=function(chip.obj, info.obj){
    idx <- returnIntIdx(info.obj)
    stopifnot(length(idx)>0)
    chip.obj@bed <- chip.obj@bed[idx,]
    chip.obj@bed[,4] <- droplevels(chip.obj@bed[,4])
    return(chip.obj)
  }
)

#' @rdname updateMatlist-methods
setMethod(f="updateMatlist",
  signature=c("matlist", "info"),
  definition=function(matlist.obj, info.obj){
    idx <- returnIntIdx(info.obj)
    stopifnot(length(idx)>0)
    matlist.obj@ll <- lapply(matlist.obj@ll, function(mat){mat[idx, ]})
    return(matlist.obj)
  }
)

#' @rdname writeChip-methods
setMethod(f = "writeChip",
  signature = c("chip", "character"),
  definition=function(chip.obj, fout) {
    stopifnot(file.info(dirname(fout))$isdir)
    write.table(chip.obj@bed, file = fout, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    return(invisible(NULL))
  }
)

#' @rdname mergeRep-methods
setMethod(f="mergeRep",
  signature=c("matlist"),
  definition=function(matlist.obj,repconfig){
    dat=read.table(repconfig,header=FALSE,as.is=TRUE,sep="\t")
    if (identical(dat[,1],dat[,2])) {return(matlist.obj)}
    cond=dat[,1]
    names(cond)=dat[,2]
    nmat.list=lapply(seq_along(unique(names(cond))), function(i){
      nm=unique(names(cond))[i]
      nmat=Reduce("+",matlist.obj@ll[cond[nm]])/length(cond[nm])
      return(nmat)
    })
    names(nmat.list)=unique(names(cond))
    matlist.obj@ll=nmat.list
    return(matlist.obj)
  }
)

#' @rdname avgplot-methods
setMethod(f="avgplot",
  signature=c("matlist","chip","info"),
  definition=function(matlist.obj, chip.obj, info.obj, pdffoutFe, pdffoutSm){
    size <- info.obj@w
    axis_name <- if (floor(size/1000)>0) {
        c(paste0("-",size/1000," kb"),"0",paste0(size/1000," kb"))
    } else if (size>0) {
        c(paste0("-",size," bp"),"0",paste0(size," bp"))
    }
    dist <- info.obj@points
    mat.list <- matlist.obj@ll
    grp <- chip.obj@bed[,4]
    dd.list <- lapply(1:length(mat.list),function(i){
      dd <- data.frame(mat.list[[i]], grp = grp, sm = names(mat.list)[i]);
      return(dd)
      })
    rm(mat.list)
    dd <- rbindlist(dd.list)
    avg <- dd[, lapply(.SD, mean, na.rm=TRUE), by="grp,sm", .SDcols=colnames(dd)[grep("grp|sm",colnames(dd),invert=TRUE)]]
    setnames(avg, c("grp","sm",dist) )
    avg <- melt(avg,id.vars=c("grp","sm"))
    avg$variable <- as.numeric(as.character(avg$variable))
    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(avg, aes_(x = ~variable, y = ~value, color = ~sm)) + geom_line(size=2)+labs(x = "",y = matlist.obj@ylab)+scale_x_continuous(breaks = c(-size,0,size), labels = axis_name)+theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top")+facet_grid(.~grp)
    ggsave(file = pdffoutFe, plot = p1)
    p2 <- ggplot(avg, aes_(x = ~variable, y = ~value, color = ~grp)) + geom_line(size=2)+labs(x="",y = matlist.obj@ylab)+scale_x_continuous(breaks=c(-size, 0, size), labels = axis_name)+theme(legend.title = element_blank(),panel.spacing = unit(2, "lines"),legend.position = "top")+facet_grid(.~sm)
    ggsave(file = pdffoutSm, plot = p2)
  }
)

#' @rdname tdplot-method
setMethod(f="tdplot",
  signature=c("matlist","chip","info"),
  def=function(matlist.obj, chip.obj, info.obj, pngfout) {
    size<-info.obj@w
    axis_name<- if (floor(size/1000)>0) {
        c(paste0("-",size/1000," kb"),"0",paste0(size/1000," kb"))
    } else if (size>0) {
        c(paste0("-",size," bp"),"0",paste0(size," bp"))
    }
    thresh <- as.numeric(quantile(unlist(matlist.obj@ll),0.95))
    sml <- lapply(matlist.obj@ll,function(mat){mat[mat<0]=0;
      mat[mat>thresh]=thresh;
      attr(mat,"extend") <- 0;
      attr(mat, "upstream_index") <- 1:round(ncol(mat)/2);
      attr(mat, "downstream_index") <- (round(ncol(mat)/2)+1):ncol(mat);
      attr(mat, "target_index") <- integer(0);
      return(mat)})
    names(sml) <- names(matlist.obj@ll)
    grp <- chip.obj@bed[,4]
    options(expressions=500000)
    if (length(unique(as.character(grp)))==1) {
      ht_list <- NULL
    } else {
      col1 <- if(length(levels(grp))<3) {brewer.pal(3, "Dark2")[1:2]} else if (length(levels(grp))>8) {brewer.pal(length(levels(grp)), "Paired")} else {brewer.pal(length(levels(grp)), "Dark2")}
      ht_list <- Heatmap(as.character(grp), col = structure(col1,names = levels(grp)), name = "Groups", show_row_names = FALSE, show_column_names = FALSE, width = unit(3, "mm"),use_raster = FALSE,show_row_dend = FALSE,show_column_dend = FALSE, split = grp, combined_name_fun = NULL)
    }
    rng=range(unlist(sml))
    for (j in 1:length(sml)) {
      if (j==1) {
        ht_list <- ht_list + EnrichedHeatmap(sml[[j]], col=colorRamp2(c(rng[1],rng[2]),
  c("white", "red")), name=names(sml)[j], split=as.numeric(grp), row_title_rot=0,column_title=names(sml)[j], combined_name_fun=NULL, axis_name=axis_name,heatmap_legend_param = list(color_bar="continuous"))
      } else {
        ht_list=ht_list+EnrichedHeatmap(sml[[j]], col=colorRamp2(c(rng[1],rng[2]),
  c("white", "red")), name=names(sml)[j], column_title=names(sml)[j], axis_name=axis_name,heatmap_legend_param = list(color_bar="continuous"))
      }
    }
    png(pngfout,width=max(3000,700*length(sml)),height=3000,res=300)
    draw(ht_list)
    dev.off()
  }
)
