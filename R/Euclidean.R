#' Euclidean Distance from Dimension Means
#'
#' Calculates two-dimensional Euclidean distance between all points and dimension means.
#'
#' @param data an \code{n x 2} matrix or data frame. 
#'
#' @return Returns a numeric vector of length 1.
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#'
#' Euclidean(data)
#'
#' @import stats
#'
#' @export

Euclidean <- function(data) {
  data <- as.data.frame(na.omit(data))
  colnames(data) <- c("x", "y")
  
  data$x_mean <- mean(data$x)
  data$y_mean <- mean(data$y)
  data$distance <- sqrt((data$x - data$x_mean)^2 + (data$y - data$y_mean)^2)
  
  out <- sum(data$distance)/nrow(data)
  
  return(out)
}
