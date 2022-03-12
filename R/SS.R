#' Sum-of-Squares Calculation
#'
#' Calculates sums of squares for uni- or multi-dimensional numeric data using the
#' distance matrix.
#'
#' @param data a numeric vector or \code{n x k} matrix or data frame.
#' @param ... arguments passed to \code{\link{dist}()}.
#'
#' @return a numeric vector of length 1.
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#' SS(data)
#'
#' @import stats
#'
#' @export

SS <- function(data, ...) {
  sum(as.matrix(Dist(data)^2))/(2*nrow(data))
}
