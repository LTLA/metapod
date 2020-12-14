#' Average grouped statistics 
#'
#' Average grouped statistics with consideration of weights and protection from missing values.
#' 
#' @param values A numeric vector containing statistics for individual tests.
#' @inheritParams groupedSimes
#'
#' @details
#' If \code{weights} is supplied, a weighted average is computed for each group.
#' If \code{values} contains missing values, they are ignored and will not contribute to the average.
#'
#' @return A named numeric vector of (weighted) averages, of length equal to the number of unique levels in \code{grouping}.
#'
#' @author Aaron Lun
#' 
#' @examples
#' grouping <- sample(LETTERS, 20, replace=TRUE)
#' averageGroupedStats(1:20, grouping)
#' averageGroupedStats(1:20, grouping, weights=runif(20))
#'
#' @seealso
#' \code{\link{averageParallelStats}}, for the parallel counterpart.
#'
#' @export
averageGroupedStats <- function(values, grouping, weights=NULL) {
    grouped <- .prepare_grouped_inputs(grouping, list(values=values, weights=weights))
    RLE <- grouped$runs
    values <- grouped$x$values
    weights <- grouped$x$weights

    last <- 0L
    collated <- numeric(length(RLE$value))
    for (i in seq_along(collated)) {
        run.len <- RLE$length[i]

        chosen <- last + seq_len(run.len)
        curvals <- values[chosen]
        keep <- !is.na(curvals)
        curvals <- curvals[keep]

        if (!is.null(weights)) {
            curw <- weights[chosen[keep]]
            if (!all(curw > 0)) {
                stop("'weights' must be all positive")
            }
            meanval <- sum(curvals * curw)/sum(curw)
        } else {
            meanval <- mean(curvals)
        }

        collated[i] <- meanval
        last <- last + run.len
    }

    names(collated) <- RLE$value
    collated
}

