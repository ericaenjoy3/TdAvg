# Modidified from ClusterR::Optimal_Clusters_KMeans
#' @title Find the optimal number of kmeans clusters
#' @name Optimal_Clusters_KM
#' @rdname Optimal_Clusters_KM-methods
#' @description
#' Find the optimal number of kmeans clusters
#' by supplying data (matrix or data.frame), and
#' specifying max_clusters (the maximum number of clusters).
#' The function plots the optimal cluster number plot and returns cluster membership (integers).
#' @param data matrix or data.frame object
#' @param max_clusters numeric/integer vector of length 1
#' @return cluster membership integer vector
#' @export Optimal_Clusters_KM
Optimal_Clusters_KM <- function (data, max_clusters, criterion = "variance_explained",
    fK_threshold = 0.85, num_init = 1, max_iters = 200, initializer = "optimal_init",
    threads = 1, tol = 1e-04, plot_clusters = TRUE, verbose = FALSE,
    tol_optimal_init = 0.3, seed = 1)
{
    if (class(data) == "data.frame")
        data = as.matrix(data)
    if (class(data) != "matrix")
        stop("data should be either a matrix or a data frame")
    if (!is.numeric(max_clusters) || length(max_clusters) !=
        1 || max_clusters < 1)
        stop("max_clusters should be numeric and greater than 0")
    if (!criterion %in% c("variance_explained", "WCSSE", "dissimilarity",
        "silhouette", "distortion_fK", "AIC", "BIC", "Adjusted_Rsquared"))
        stop("available criteria are 'variance_explained', 'WCSSE', 'dissimilarity', 'silhouette', 'distortion_fK', 'AIC', 'BIC' and 'Adjusted_Rsquared'")
    if (num_init < 1)
        stop("the num_init parameter should be greater than 0")
    if (max_iters < 1)
        stop("the max_iters parameter should be greater than 0")
    if (!initializer %in% c("kmeans++", "random", "optimal_init",
        "quantile_init"))
        stop("available initializer methods are 'kmeans++', 'random', 'quantile_init' and 'optimal_init'")
    if (threads < 1)
        stop("threads should be an integer greater than 0")
    if (tol <= 0)
        stop("tol should be a float number greater than 0.0")
    if (!is.logical(plot_clusters))
        stop("the plot_clusters parameter should be either TRUE or FALSE")
    if (!is.logical(verbose))
        stop("the verbose parameter should be either TRUE or FALSE")
    if (tol_optimal_init <= 0)
        stop("tol_optimal_init should be a float number greater than 0.0")
    flag_non_finite = ClusterR:::check_NaN_Inf(data)
    if (!flag_non_finite)
        stop("the data includes NaN's or +/- Inf values")
    vec_out = rep(NA, max_clusters)
    if (verbose) {
        cat("", "\n")
        pb = utils::txtProgressBar(min = 1, max = max_clusters,
            style = 3)
        cat("", "\n")
    }
    for (i in 1:max_clusters) {
        km = ClusterR:::KMEANS_rcpp(data, i, num_init, max_iters, initializer,
            FALSE, threads, FALSE, NULL, tol, 1e-06, tol_optimal_init,
            seed)
        if (criterion == "variance_explained") {
            vec_out[i] = sum(na.omit(as.vector(km$WCSS_per_cluster)))/km$total_SSE
        }
        if (criterion == "WCSSE") {
            vec_out[i] = sum(na.omit(as.vector(km$WCSS_per_cluster)))
        }
        if (criterion == "dissimilarity") {
            eval_km = ClusterR:::evaluation_rcpp(data, as.vector(km$clusters),
                FALSE)
            tmp_dis = mean(na.omit(unlist(lapply(eval_km$INTRA_cluster_dissimilarity,
                mean))))
            vec_out[i] = tmp_dis
        }
        if (criterion == "silhouette") {
            if (i == 1) {
                vec_out[i] = 0
            }
            else {
                eval_km = ClusterR:::evaluation_rcpp(data, as.vector(km$clusters),
                  TRUE)
                tmp_silh = mean(na.omit(unlist(lapply(eval_km$silhouette,
                  mean))))
                vec_out[i] = tmp_silh
            }
        }
        if (criterion == "distortion_fK") {
            vec_out[i] = sum(na.omit(as.vector(km$WCSS_per_cluster)))
        }
        if (criterion == "AIC") {
            m = ncol(km$centers)
            k = nrow(km$centers)
            D = sum(na.omit(km$WCSS_per_cluster))
            vec_out[i] = D + 2 * m * k
        }
        if (criterion == "BIC") {
            m = ncol(km$centers)
            k = nrow(km$centers)
            n = length(km$clusters)
            D = sum(na.omit(km$WCSS_per_cluster))
            vec_out[i] = D + log(n) * m * k
        }
        if (criterion == "Adjusted_Rsquared") {
            vec_out[i] = sum(na.omit(km$WCSS_per_cluster))
        }
        if (verbose) {
            utils::setTxtProgressBar(pb, i)
        }
    }
    if (verbose) {
        close(pb)
        cat("", "\n")
    }
    if (criterion == "Adjusted_Rsquared") {
        vec_out = 1 - (vec_out * (nrow(data) - 1))/(vec_out[1] *
            (nrow(data) - seq(1, max_clusters)))
    }
    if (plot_clusters) {
        tmp_VAL = as.vector(na.omit(vec_out))
        if (length(which(is.na(vec_out))) > 0) {
            x_dis = (1:length(vec_out))[-which(is.na(vec_out))]
            y_dis = vec_out[-which(is.na(vec_out))]
        }
        else {
            x_dis = 1:length(vec_out)
            y_dis = vec_out
        }
        y_MAX = max(tmp_VAL)
        if (criterion %in% c("variance_explained", "WCSSE", "dissimilarity",
            "silhouette", "AIC", "BIC", "Adjusted_Rsquared")) {
            plot(x = x_dis, y = y_dis, type = "l", xlab = "clusters",
                ylab = criterion, col = "blue", lty = 3, axes = FALSE)
            axis(1, at = seq(1, length(vec_out), by = 1))
            if (criterion == "silhouette") {
                axis(2, at = seq(0, y_MAX + 0.05, by = 0.05),
                  las = 1, cex.axis = 0.8)
                abline(h = seq(0, max(as.vector(na.omit(vec_out))),
                  0.05), v = seq(1, length(vec_out), by = 1),
                  col = "gray", lty = 3)
            }
            else {
                tmp_summary = round(summary(y_MAX)[["Max."]])
                out_max_summary = ifelse(tmp_summary == 0, 1,
                  tmp_summary)
                axis(2, at = seq(0, y_MAX + out_max_summary/10,
                  by = out_max_summary/10), las = 1, cex.axis = 0.8)
                abline(h = seq(0, max(as.vector(na.omit(vec_out))),
                  out_max_summary/10), v = seq(1, length(vec_out),
                  by = 1), col = "gray", lty = 3)
            }
            if (criterion %in% c("variance_explained", "Adjusted_Rsquared",
                "dissimilarity", "silhouette")) {
                text(x = 1:length(vec_out), y = vec_out, labels = round(vec_out,
                  2), cex = 0.8, font = 2)
            }
            else {
                text(x = 1:length(vec_out), y = vec_out, labels = round(vec_out,
                  1), cex = 0.8, font = 2)
            }
        }
        if (criterion == "distortion_fK") {
            f_K = ClusterR:::opt_clust_fK(vec_out, ncol(data), fK_threshold)
            fK_vec = as.vector(f_K$fK_evaluation)
            if (length(which(is.na(fK_vec))) > 0) {
                x_fk = (1:length(fK_vec))[-which(is.na(fK_vec))]
                y_fk = fK_vec[-which(is.na(fK_vec))]
            }
            else {
                x_fk = 1:length(fK_vec)
                y_fk = fK_vec
            }
            par(oma = c(0, 2, 0, 0))
            plot(y_fk, type = "l", xlab = "clusters", ylab = "f(K)",
                col = "green", axes = FALSE)
            axis(1, at = x_fk)
            axis(2, at = seq(0, max(y_fk) + 0.1, by = round(summary(y_fk)[["Max."]])/10),
                las = 1, cex.axis = 0.8)
            abline(h = seq(0, max(y_fk), round(summary(y_fk)[["Max."]])/10),
                v = seq(1, length(y_fk), by = 1), col = "gray",
                lty = 3)
            abline(h = fK_threshold, col = "blue", lty = 3)
            mtext("threshold", side = 2, line = 2, at = fK_threshold,
                las = 1, cex = 0.9)
            text(x = x_fk, y = y_fk, labels = round(y_fk, 2),
                cex = 0.8, font = 2)
      	    return(x_fk[which.min(y_fk)]) # added by EL, July 11 2017
        }
    }
    class(vec_out) = "k-means clustering"
    return(vec_out)
}
