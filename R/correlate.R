#' Test for Bivariate Correlation
#'
#' Calculates correlation coefficient between two variables and returns a list containing the
#' correlation estimate, its standard error, the p-value of a null-hypothesis significance test, and the
#' number of observations used.
#' 
#' @details
#' Additional arguments to alter the type of null hypothesis significance test, the method used to
#' calculate the correlation coefficient, the confidence level, or other options should be passed to
#' \code{correlate}() and will be inherited by \code{\link{cor.test}()}. Note that unlike
#' \code{\link{cor.test}()}, both arguments \code{x} and \code{y} are required.
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param ... arguments passed to \code{\link{cor.test}()}.
#'
#' @return Returns a list with elements containing the correlation coefficient estimate, its associated
#' standard error, the p-value of a null-hypothesis significance test, and the number of observations
#' used, all as numeric vectors of length 1.
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#'
#' correlate(data[, 1], data[, 2])
#'
#' @import stats
#'
#' @export

correlate <- function(x, y, ...) {
  cor_result <- cor.test(x, y, ...)
  cor_se <- sqrt((1 - cor_result$estimate[["cor"]]^2)/(cor_result$parameter[["df"]]))
  
  out <- list("estimate" = cor_result$estimate[["cor"]],
              "se" = cor_se,
              "p-value" = cor_result$p.value,
              "observations" = cor_result$parameter[["df"]] + 2)
  
  return(out)
}
