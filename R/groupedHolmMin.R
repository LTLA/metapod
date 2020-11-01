#' Combine grouped p-values with the minimum Holm approach
#'
#' Combine p-values from grouped tests with the minimum Holm approach.
#' Groups are defined according to unique levels of a grouping factor.
#'
#' @inheritParams groupedSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
#' This contains the minimum Holm p-value for each group, log-transformed if \code{log.p=TRUE}.
#' Each entry is named according to the group.
#' \item \code{representative}, a named integer scalar specifying the representative test for each group.
#' Each entry is named according to the group.
#' \item \code{influential}, a logical vector of length equal to \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
#' }
#'
#' @inherit parallelHolmMin details
#' @inherit parallelHolmMin references
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' g <- sample(10, length(p1), replace=TRUE)
#'
#' # Standard application:
#' out <- groupedHolmMin(p1, g)
#' str(out)
#'
#' # With weights:
#' out <- groupedHolmMin(p1, g, weights=rexp(length(p1)))
#' str(out)
#' 
#' # With log p-values. 
#' out <- groupedHolmMin(log(p1), g, log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{parallelHolmMin}}, for a version that operates on parallel vectors of p-values.
#'
#' \code{\link{groupedWilkinson}}, for a more relaxed version of this test when hypotheses are independent.
#' @export
groupedHolmMin <- function(p.values, grouping, weights=NULL, log.p=FALSE, min.n=1, min.prop=0.5) {
    .grouped_compute(p.values, grouping, weights=weights, log=log.p, min_n=min.n, min_prop=min.prop, FUN=compute_grouped_holm_min)
}
