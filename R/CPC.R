#' Cluster-Polarization Coefficient
#'
#' Implements clustering algorithms and calculates cluster-polarization coefficient.
#' Contains support for hierarchical clustering, k-means clustering, partitioning
#' around medoids, and manual assignment of cluster membership.
#'
#' @details
#' \code{type} must take one of four values: \code{"hclust"} performs agglomerative
#' hierarchical clustering via \code{\link{hclust}()}. \code{"kmeans"}
#' performs k-means clustering via \code{\link{kmeans}()}. \code{"pam"}
#' performs k-medoids clustering via \code{\link{pam}()}. \code{"manual"}
#' indicates that no clustering is necessary and that the researcher has specified
#' cluster assignments.
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
#' @param k the desired number of clusters.
#' @param type a character string giving the type of clustering method to be used.
#' See Details.
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
#' containing output from the specified clustering function, all sums of squares,
#' CPC, and adjusted CPC. If \code{model = FALSE}, \code{CPC()} returns a numeric
#' vector of length 1 giving the CPC (if \code{adjust = FALSE}) or adjusted CPC (if
#' \code{adjust = TRUE}).
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#' clusters <- matrix(c(rep(1, 25), rep(2, 25)), ncol = 1)
#' data <- cbind(data, clusters)
#'
#' CPC(data[,c(1:2)], 2, "kmeans")
#' CPC(data, 2, "manual", cols = 1:2, clusters = 3)
#'
#' @import stats
#' @import cluster
#'
#' @export

CPC <- function(data, k, type, model = FALSE, adjust = FALSE, cols = NULL,
                clusters = NULL, ...) {
  input <- data[colSums(!is.na(as.matrix(data))) > 0]
  input <- as.matrix(na.omit(input))
  cluster <- NULL

  if(length(unique(input)) < k){
    warning("More clusters than unique data points; NAs generated")
    return(NA)
  }

  else{
    switch (type,
            hclust = {
              input <- apply(input, 2, as.numeric)
              input_dist <- dist(input)
              output_hclust <- hclust(input_dist, ...)
              cut_hclust <- as.data.frame(cutree(output_hclust, k = k))
              colnames(cut_hclust) <- "cluster"
              new_hclust <- cbind(input, cut_hclust)
              WSS_hclust <- c()

              for (i in 1:k) {
                WSS <- SS(subset(new_hclust, cluster == i))
                WSS_hclust <- c(WSS_hclust, WSS)
              }

              TSS_hclust <- SS(input)
              TWSS_hclust <- sum(WSS_hclust)
              BSS_hclust <- TSS_hclust - TWSS_hclust
              CPC <- BSS_hclust/TSS_hclust
              CPC.adj <- 1 - (TWSS_hclust/TSS_hclust)*(1/k)

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
                     CPC.adj = CPC.adj)
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

              CPC <- output_kmeans$betweenss/output_kmeans$totss
              CPC.adj <- 1 - (output_kmeans$tot.withinss/output_kmeans$totss)*(1/k)

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
                     CPC.adj = CPC.adj)
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
                WSS <- SS(subset(new_pam, cluster == i))
                WSS_pam <- c(WSS_pam, WSS)
              }

              TSS_pam <- SS(input)
              TWSS_pam <- sum(WSS_pam)
              BSS_pam <- TSS_pam - TWSS_pam
              CPC <- BSS_pam/TSS_pam
              CPC.adj <- 1 - (TWSS_pam/TSS_pam)*(1/k)

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
                     CPC.adj = CPC.adj)
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
              data_manual <- input[,cols]
              WSS_manual <- c()

              for (i in 1:k) {
                data_temp <- subset(input, cluster == i)
                data_temp <- data_temp[,cols]
                WSS <- SS(data_temp)
                WSS_manual <- c(WSS_manual, WSS)
              }

              TSS_manual <- SS(data_manual)
              TWSS_manual <- sum(WSS_manual)
              BSS_manual <- TSS_manual - TWSS_manual
              CPC <- BSS_manual/TSS_manual
              CPC.adj <- 1 - (TWSS_manual/TSS_manual)*(1/k)

              if(model){
                list(data = input,
                     WSS = WSS_manual,
                     TWSS = TWSS_manual,
                     BSS = BSS_manual,
                     TSS = TSS_manual,
                     CPC = CPC,
                     CPC.adj = CPC.adj)
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
