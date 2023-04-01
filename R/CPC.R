#' Cluster-Polarization Coefficient
#'
#' Implements clustering algorithms and calculates cluster-polarization coefficient.
#' Contains support for hierarchical clustering, k-means clustering, partitioning
#' around medoids, density-based spatial clustering with noise, and manual assignment
#' of cluster membership.
#'
#' @details
#' \code{type} must take one of six values: \cr
#' \code{"hclust"}: agglomerative hierarchical clustering with \code{\link{hclust}()}, \cr
#' \code{"diana"}: divisive hierarchical clustering with \code{\link{diana}()}, \cr
#' \code{"kmeans"}: k-means clustering with \code{\link{kmeans}()}, \cr
#' \code{"pam"}: k-medoids clustering with \code{\link{pam}()}, \cr
#' \code{"dbscan"}: density-based clustering with \code{\link{dbscan}()}, \cr
#' \code{"manual"}: no clustering is necessary, researcher has specified cluster assignments.
#'
#' For all clustering methods, additional arguments to fine-tune clustering
#' performance, such as the specific algorithm to be used, should be passed to
#' \code{CPC()} and will be inherited by the specified clustering function. In
#' particular, if \code{type = "kmeans"}, using a large number of random starts is
#' recommended. This can be specified with the \code{nstart} argument to
#' \code{\link{kmeans}()}, passed directly to \code{CPC()}.
#'
#' If \code{type = "manual"}, \code{data} must contain a vector identifying cluster
#' membership for each observation, and \code{cols} and \code{clusters} must be
#' defined.
#'
#' @param data a numeric vector or \code{n x k} matrix or data frame. If
#' \code{type = "manual"}, \code{data} must be a matrix containing a vector
#' identifying cluster membership for each observation, to be passed to
#' \code{clusters} argument.
#' @param type a character string giving the type of clustering method to be used.
#' See Details.
#' @param k the desired number of clusters. Required if \code{type} is one of \code{"hclust"},
#' \code{"diana"}, \code{"kmeans"}, or \code{"pam"}.
#' @param epsilon radius of epsilon neighborhood. Required if \code{type = "dbscan"}.
#' @param model a logical indicating whether clustering model output should be
#' returned. Defaults to \code{FALSE}.
#' @param adjust a logical indicating whether the adjusted CPC should be calculated.
#' Defaults to \code{FALSE}. Note that both CPC and adjusted CPC are automatically
#' calculated and returned if \code{model = TRUE}.
#' @param cols columns of \code{data} to be used in CPC calculation. Only used if
#' \code{type = "manual"}.
#' @param clusters column of \code{data} indicating cluster membership for each
#' observation. Only used if \code{type = "manual"}.
#' @param ... arguments passed to other functions.
#'
#' @return If \code{model = TRUE}, \code{CPC()} returns a list with components
#' containing output from the specified clustering function, all sums of squares, the
#' CPC, the adjusted CPC, and associated standard errors. If \code{model = FALSE}, \code{CPC()} returns
#' a numeric vector of length 1 giving the CPC (if \code{adjust = FALSE}) or adjusted CPC (if
#' \code{adjust = TRUE}).
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#' clusters <- matrix(c(rep(1, 25), rep(2, 25)), ncol = 1)
#' data <- cbind(data, clusters)
#'
#' CPC(data[,c(1:2)], "kmeans", k = 2)
#' CPC(data, "manual", cols = 1:2, clusters = 3)
#'
#' @import stats
#' @importFrom cluster pam diana
#' @importFrom dbscan dbscan
#' @export

CPC <- function(data, type, k = NULL, epsilon = NULL, model = FALSE, adjust = FALSE,
                cols = NULL, clusters = NULL, ...) {
  data <- as.matrix(data)
  input <- data[colSums(!is.na(data)) > 0]
  input <- matrix(na.omit(input), ncol = ncol(data))
  cluster <- NULL

  k <- ifelse(!type %in% c("dbscan", "manual"), k, 0)

  if(length(unique(input)) < k){
    warning("More clusters than unique data points; NAs generated")
    return(NA)
  }

  else{
    switch (type,
            dbscan = {
              output_dbscan <- dbscan(x = input, eps = epsilon, ...)
              new_dbscan <- cbind(input, unlist(output_dbscan$cluster))
              new_dbscan <- subset(new_dbscan, new_dbscan[,ncol(new_dbscan)] != 0)

              new_dbscan <- CPCdata.frame(data = new_dbscan,
                                          cols = -ncol(new_dbscan),
                                          clusters = ncol(new_dbscan))
              data_dbscan <- as.matrix(new_dbscan[,-ncol(new_dbscan)])
              data_dbscan <- apply(data_dbscan, 2, as.numeric)
              WSS_dbscan <- c()

              for (i in unique(new_dbscan$cluster)) {
                data_temp <- new_dbscan[new_dbscan$cluster == i,]
                data_temp <- as.matrix(data_temp[,-ncol(new_dbscan)])
                data_temp <- apply(data_temp, 2, as.numeric)
                WSS <- SS(as.matrix(data_temp))
                WSS_dbscan <- c(WSS_dbscan, WSS)
              }

              TSS_dbscan <- SS(as.matrix(data_dbscan))
              TWSS_dbscan <- sum(WSS_dbscan)
              BSS_dbscan <- TSS_dbscan - TWSS_dbscan
              n_i <- nrow(as.matrix(data_dbscan))
              n_j <- ncol(as.matrix(data_dbscan))
              n_k <- length(unique(new_dbscan$cluster))
              CPC <- BSS_dbscan/TSS_dbscan
              CPC_sd <- sqrt((2*(n_j*n_k - n_j)*(n_i - n_j*n_k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_dbscan/TSS_dbscan)*((n_i - n_j)/(n_i - n_j*n_k))
              CPC.adj_sd <- sqrt((2*(n_j*n_k - n_j))/((n_i - n_j*n_k)*(n_i - n_j + 1)))

              if(model){
                list(cluster = output_dbscan$cluster,
                     minPts = output_dbscan$minPts,
                     data = input,
                     WSS = WSS_dbscan,
                     TWSS = TWSS_dbscan,
                     BSS = BSS_dbscan,
                     TSS = TSS_dbscan,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            hclust = {
              input <- apply(input, 2, as.numeric)
              input_dist <- dist(input)
              output_hclust <- hclust(input_dist, ...)
              cut_hclust <- as.data.frame(cutree(output_hclust, k = k))
              colnames(cut_hclust) <- "cluster"
              new_hclust <- cbind(input, cut_hclust)
              WSS_hclust <- c()

              for (i in 1:k) {
                WSS <- SS(new_hclust[new_hclust$cluster == i,])
                WSS_hclust <- c(WSS_hclust, WSS)
              }

              TSS_hclust <- SS(input)
              TWSS_hclust <- sum(WSS_hclust)
              BSS_hclust <- TSS_hclust - TWSS_hclust
              n_i <- nrow(input)
              n_j <- ncol(input)
              CPC <- BSS_hclust/TSS_hclust
              CPC_sd <- sqrt((2*(n_j*k - n_j)*(n_i - n_j*k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_hclust/TSS_hclust)*((n_i - n_j)/(n_i - n_j*k))
              CPC.adj_sd <- sqrt((2*(n_j*k - n_j))/((n_i - n_j*k)*(n_i - n_j + 1)))

              if(model){
                list(merge = output_hclust$merge,
                     height = output_hclust$height,
                     order = output_hclust$order,
                     labels = output_hclust$labels,
                     method = output_hclust$method,
                     call = output_hclust$call,
                     dist.method = output_hclust$dist.method,
                     data = new_hclust,
                     WSS = WSS_hclust,
                     TWSS = TWSS_hclust,
                     BSS = BSS_hclust,
                     TSS = TSS_hclust,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            diana = {
              input <- apply(input, 2, as.numeric)
              input_dist <- dist(input)
              output_diana <- diana(input_dist, ...)
              cut_diana <- as.data.frame(cutree(output_diana, k = k))
              colnames(cut_diana) <- "cluster"
              new_diana <- cbind(input, cut_diana)
              WSS_diana <- c()

              for (i in 1:k) {
                WSS <- SS(new_diana[new_diana$cluster == i,])
                WSS_diana <- c(WSS_diana, WSS)
              }

              TSS_diana <- SS(input)
              TWSS_diana <- sum(WSS_diana)
              BSS_diana <- TSS_diana - TWSS_diana
              n_i <- nrow(input)
              n_j <- ncol(input)
              CPC <- BSS_diana/TSS_diana
              CPC_sd <- sqrt((2*(n_j*k - n_j)*(n_i - n_j*k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_diana/TSS_diana)*((n_i - n_j)/(n_i - n_j*k))
              CPC.adj_sd <- sqrt((2*(n_j*k - n_j))/((n_i - n_j*k)*(n_i - n_j + 1)))

              if(model){
                list(order = output_diana$order,
                     height = output_diana$height,
                     dc = output_diana$merge,
                     diss = output_diana$diss,
                     call = output_diana$call,
                     data = new_diana,
                     WSS = WSS_diana,
                     TWSS = TWSS_diana,
                     BSS = BSS_diana,
                     TSS = TSS_diana,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            kmeans = {
              input <- apply(input, 2, as.numeric)
              output_kmeans <- kmeans(x = input, centers = k, ...)
              cluster_kmeans <- as.data.frame(output_kmeans$cluster)
              colnames(cluster_kmeans) <- "cluster"
              new_kmeans <- cbind(input, cluster_kmeans)

              n_i <- nrow(input)
              n_j <- ncol(input)
              CPC <- output_kmeans$betweenss/output_kmeans$totss
              CPC_sd <- sqrt((2*(n_j*k - n_j)*(n_i - n_j*k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 -
                (output_kmeans$tot.withinss/output_kmeans$totss)*((n_i - n_j)/(n_i - n_j*k))
              CPC.adj_sd <- sqrt((2*(n_j*k - n_j))/((n_i - n_j*k)*(n_i - n_j + 1)))

              if(model){
                list(centers = output_kmeans$centers,
                     size = output_kmeans$size,
                     iter = output_kmeans$iter,
                     ifault = output_kmeans$ifault,
                     data = new_kmeans,
                     WSS = output_kmeans$withinss,
                     TWSS = output_kmeans$tot.withinss,
                     BSS = output_kmeans$betweenss,
                     TSS = output_kmeans$totss,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            pam = {
              input <- apply(input, 2, as.numeric)
              output_pam <- pam(x = input, k = k, ...)
              cluster_pam <- as.data.frame(output_pam$clustering)
              colnames(cluster_pam) <- "cluster"
              new_pam <- cbind(input, cluster_pam)
              WSS_pam <- c()

              for (i in 1:k) {
                WSS <- SS(new_pam[new_pam$cluster == i,])
                WSS_pam <- c(WSS_pam, WSS)
              }

              TSS_pam <- SS(input)
              TWSS_pam <- sum(WSS_pam)
              BSS_pam <- TSS_pam - TWSS_pam
              n_i <- nrow(input)
              n_j <- ncol(input)
              CPC <- BSS_pam/TSS_pam
              CPC_sd <- sqrt((2*(n_j*k - n_j)*(n_i - n_j*k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_pam/TSS_pam)*((n_i - n_j)/(n_i - n_j*k))
              CPC.adj_sd <- sqrt((2*(n_j*k - n_j))/((n_i - n_j*k)*(n_i - n_j + 1)))

              if(model){
                list(medoids = output_pam$medoids,
                     id.med = output_pam$id.med,
                     objective = output_pam$objective,
                     isolation = output_pam$isolation,
                     clusinfo = output_pam$clusinfo,
                     silinfo = output_pam$silinfo,
                     diss = output_pam$diss,
                     call = output_pam$call,
                     data = new_pam,
                     WSS = WSS_hclust,
                     TWSS = TWSS_hclust,
                     BSS = BSS_hclust,
                     TSS = TSS_hclust,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            hclust = {
              input <- apply(input, 2, as.numeric)
              input_dist <- dist(input)
              output_hclust <- hclust(input_dist, ...)
              cut_hclust <- as.data.frame(cutree(output_hclust, k = k))
              colnames(cut_hclust) <- "cluster"
              new_hclust <- cbind(input, cut_hclust)
              WSS_hclust <- c()

              for (i in 1:k) {
                WSS <- SS(new_hclust[new_hclust$cluster == i,])
                WSS_hclust <- c(WSS_hclust, WSS)
              }

              TSS_hclust <- SS(input)
              TWSS_hclust <- sum(WSS_hclust)
              BSS_hclust <- TSS_hclust - TWSS_hclust
              n_i <- nrow(input)
              n_j <- ncol(input)
              CPC <- BSS_hclust/TSS_hclust
              CPC_sd <- sqrt((2*(n_j*k - n_j)*(n_i - n_j*k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_hclust/TSS_hclust)*((n_i - n_j)/(n_i - n_j*k))
              CPC.adj_sd <- sqrt((2*(n_j*k - n_j))/((n_i - n_j*k)*(n_i - n_j + 1)))

              if(model){
                list(merge = output_hclust$merge,
                     height = output_hclust$height,
                     order = output_hclust$order,
                     labels = output_hclust$labels,
                     method = output_hclust$method,
                     call = output_hclust$call,
                     dist.method = output_hclust$dist.method,
                     data = new_hclust,
                     WSS = WSS_hclust,
                     TWSS = TWSS_hclust,
                     BSS = BSS_hclust,
                     TSS = TSS_hclust,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            },
            manual = {
              input <- CPCdata.frame(data = data, cols = cols, clusters = clusters)
              data_manual <- as.matrix(input[, -ncol(input)])
              data_manual <- apply(data_manual, 2, as.numeric)
              WSS_manual <- c()

              for (i in unique(input$cluster)) {
                data_temp <- input[input$cluster == i,]
                data_temp <- as.matrix(data_temp[, -ncol(input)])
                data_temp <- apply(data_temp, 2, as.numeric)
                WSS <- SS(as.matrix(data_temp))
                WSS_manual <- c(WSS_manual, WSS)
              }

              TSS_manual <- SS(as.matrix(data_manual))
              TWSS_manual <- sum(WSS_manual)
              BSS_manual <- TSS_manual - TWSS_manual
              n_i <- nrow(input)
              n_j <- ncol(input)
              n_k <- length(unique(input$cluster))
              CPC <- BSS_manual/TSS_manual
              CPC_sd <- sqrt((2*(n_j*n_k - n_j)*(n_i - n_j*n_k))/(((n_i - n_j)^2)*(n_i - n_j + 1)))
              CPC.adj <- 1 - (TWSS_manual/TSS_manual)*((n_i - n_j)/(n_i - n_j*n_k))
              CPC.adj_sd <- sqrt((2*(n_j*n_k - n_j))/((n_i - n_j*n_k)*(n_i - n_j + 1)))

              if(model){
                list(data = input,
                     WSS = WSS_manual,
                     TWSS = TWSS_manual,
                     BSS = BSS_manual,
                     TSS = TSS_manual,
                     CPC = CPC,
                     CPC_sd = CPC_sd,
                     CPC.adj = CPC.adj,
                     CPC.adj_sd = CPC.adj_sd)
              }

              else{
                if(adjust){
                  CPC.adj
                }

                else{
                  CPC
                }
              }
            }
    )
  }
}
