#' Multidimensional Difference-in-Means
#'
#' Calculates average Euclidean distance between means in arbitrary dimensions.
#'
#' @param data a numeric vector or \code{n x k} matrix or data frame containing a vector
#' identifying cluster membership for each observation, to be passed to
#' \code{clusters} argument.
#' @param cols columns of \code{data} to be used in difference-in-means calculation.
#' @param clusters column of \code{data} indicating cluster membership for each
#' observation.
#'
#' @return Returns a numeric vector of length 1.
#'
#' @examples
#' data <- matrix(c(rnorm(50, 0, 1), rnorm(50, 5, 1)), ncol = 2, byrow = TRUE)
#' clusters <- matrix(c(rep(1, 25), rep(2, 25)), ncol = 1)
#' data <- cbind(data, clusters)
#'
#' diff_multidim(data, 1:2, 3)
#'
#' @import stats
#'
#' @export

diff_multidim <- function(data, cols, clusters) {
  input <- na.omit(as.matrix(data)[,c(cols, clusters)])
  colnames(input)[ncol(input)] <- "clusters"
  input <- as.data.frame(input)

  means <- c()

  for (i in unique(input$clusters)) {
    name <- paste0("cluster_", i)

    assign(name, apply(subset(input, clusters == i)[, -clusters], 2, as.numeric))

    obs <- nrow(subset(input, clusters == i))

    if (obs > 1) {
      means <- c(means, colMeans(eval(parse(text = name))[, -ncol(t(as.matrix(eval(parse(text = name)))))]))
    }

    else {
      means <- c(means, colMeans(t(as.matrix(t(as.matrix(eval(parse(text = name))[, -ncol(t(as.matrix(eval(parse(text = name)))))]))))))
    }
  }

  diff <- mean(as.numeric(dist(matrix(means, ncol = length(unique(cols)), byrow = TRUE))))

  return(diff)
}
