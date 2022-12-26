#' Data Manipulation for CPC Calculation
#'
#' Converts numeric matrix to data frame with necessary format for
#' \code{"manual"} \code{\link{CPC}()} calculation.
#'
#' @param data a numeric \code{n x k} matrix or data frame.
#' @param cols columns in \code{data} to be used for calculating \code{\link{CPC}()}.
#' @param clusters column in \code{data} giving cluster membership.
#'
#' @return Returns a data frame with dimensions identical to those of \code{data}.
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#' clusters <- matrix(c(rep(1, 25), rep(2, 25)), ncol = 1)
#' data <- cbind(data, clusters)
#' CPCdata.frame(data, 1:2, 3)
#'
#' @export

CPCdata.frame <- function(data, cols, clusters) {
  new_data <- as.data.frame(data[,cols])
  new_clusters <- as.data.frame(data[,clusters])
  colnames(new_clusters) <- "cluster"
  na.omit(cbind(new_data, new_clusters))
}
