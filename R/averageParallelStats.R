#' Average parallel statistics
#'
#' Average parallel statistics with consideration of weights and protection from missing values.
#'
#' @param values A list of numeric vectors of the same length, containing statistics from parallel tests.
#' @inheritParams parallelSimes
#'
#' @details
#' If \code{weights} is supplied, a weighted average is computed for each parallel set of tests.
#' If \code{values} contains missing values, they are ignored and will not contribute to the average.
#'
#' @return A numeric vector of (weighted) averages, of length equal to the lengths of the vectors in \code{values}.
#'
#' @author Aaron Lun
#' 
#' @examples
#' averageParallelStats(list(1:10, 2:11))
#' averageParallelStats(list(1:10, 2:11), weights=c(2:1))
#' averageParallelStats(list(1:10, rep(NA, 10)))
#'
#' @seealso
#' \code{\link{averageParallelStats}}, for the parallel counterpart.
#' @export
averageParallelStats <- function(values, weights=NULL) {
    n <- unique(lengths(values))
    if (length(n)!=1){ 
        stop("all vectors in 'values' should have the same length")
    }
    if (!n) {
        return(numeric(0))
    }

    .valid_parallel_weights(weights)
    if (is.null(weights)) {
        weights <- rep(1, length(values))
    } else {
        stopifnot(length(values)==length(weights)) 
        if (is.list(weights)) {
            n2 <- unique(lengths(weights))
            if (!identical(unname(n2), unname(n))) {
                stop("'weights' must have same 'lengths' as 'values'")
            }
        }
    }

    combined <- total.weight <- 0
    for (x in seq_along(values)) {
        cur.weights <- weights[[x]]
        if (!all(cur.weights > 0)) {
            stop("'weights' must be positive")
        }

        cur.values <- values[[x]]
        product <- cur.values * cur.weights

        # Ignoring NA values and their weights.
        failed <- is.na(cur.values)
        product[failed] <- 0
        w <- ifelse(failed, 0, cur.weights)

        combined <- combined + product
        total.weight <- total.weight + w
    }

    combined/total.weight
}
