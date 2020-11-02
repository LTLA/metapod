#' Count directions of parallel tests
#'
#' Count the number of parallel tests that are significant and changing in each direction.
#' Each group of tests is defined as corresponding entries across vectors.
#'
#' @param p.values A list of numeric vectors of the same length, containing the p-values to be combined.
#' @param effects A list of numeric vectors of same structure as \code{p.values}, containing the effect sizes for each set of tests.
#' @param p.threshold Numeric scalar defining the adjusted p-value threshold at which test is significant.
#' @param effect.threshold Numeric scalar defining the threshold at which an effect is \dQuote{"up"} or \dQuote{"down"}.
#' @param method String specifying how the multiple testing correction within each group is to be performed.
#' @param log.p Logical scalar indicating whether the \code{p.values} are log-transformed.
#'
#' @return A list with the \code{"up"} and \code{"down"} entries.
#' Both are integer vectors of length equal to the number of groups,
#' containing the number of significant tests changing each each direction in each group.
#'
#' @details
#' We apply a multiple testing correction within each group and define all significant tests as those with adjusted p-values below \code{p.threshold}.
#' We then count the number of tests with effects greater than \code{effect.threshold} (i.e., going up) or less than the threshold (i.e., going down).
#' This allows us to quantify the direction of change across the group into two counts.
#'
#' The output of this function is simpler to interpret than the \code{summarize*Direction} functions.
#' However, \code{p.threshold} has no clear relationship to the significance threshold applied to the combined p-values.
#' Groups with many low p-values will have milder multiplicity corrections, resulting in more non-zero counts for both the \code{up} and \code{down} tests;
#' this often complicates matters by introducing mixed directions of effect, even when the most significant changes in the group are in one direction. 
#'
#' Note that, if \code{log.p=TRUE}, the function will automatically handle the log-transformation of \code{p.threshold}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{summarizeParallelDirection}}, for another way to summarize the overall direction into a single string.
#'
#' \code{\link{countGroupedDirection}}, for the equivalent function based on a grouping factor.
#'
#' @examples
#' p1 <- rbeta(100, 0.5, 1)
#' eff1 <- rnorm(100)
#' p2 <- rbeta(100, 0.5, 1)
#' eff2 <- rnorm(100)
#'
#' (dirs <- countParallelDirection(list(p1, p2), list(eff1, eff2)))
#'
#' @export
countParallelDirection <- function(p.values, effects, p.threshold=0.05, effect.threshold=0, method=c("BH", "holm"), log.p=FALSE) {
    count_parallel_direction(p.values, effects, method=c(BH=0L, holm=1L)[match.arg(method)], 
        pthreshold=p.threshold, ethreshold=effect.threshold, log=log.p)
}
