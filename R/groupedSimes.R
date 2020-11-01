#' Combine grouped p-values with Simes' method
#'
#' Combine p-values from grouped tests with Simes' method.
#' Groups are defined according to unique levels of a grouping factor.
#'
#' @param p.values A numeric vector containing p-values for individual tests.
#' @param grouping A vector or factor of length equal to \code{p.values}, specifying the group to which each test is assigned.
#' @param weights A numeric vector of length equal to \code{p.values}, containing a positive weight for each test.
#' Alternatively \code{NULL}, in which case equal weights are assigned to all tests.
#' @param log.p Logical scalar indicating whether the p-values in \code{p.values} are log-transformed.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
#' This contains the combined Simes p-value for each group, log-transformed if \code{log.p=TRUE}.
#' Each entry is named according to the group.
#' \item \code{representative}, a named integer scalar specifying the representative test for each group.
#' Each entry is named according to the group.
#' \item \code{influential}, a logical vector of length equal to \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
#' }
#'
#' @inherit parallelSimes details
#' @inherit parallelSimes references
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' g <- sample(10, length(p1), replace=TRUE)
#'
#' # Standard application:
#' out <- groupedSimes(p1, g)
#' str(out)
#'
#' # With weights:
#' out <- groupedSimes(p1, g, weights=rexp(length(p1)))
#' str(out)
#' 
#' # With log p-values. 
#' out <- groupedSimes(log(p1), g, log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{parallelSimes}}, for a version that operates on parallel vectors of p-values.
#'
#' @export
groupedSimes <- function(p.values, grouping, weights=NULL, log.p=FALSE) {
    .grouped_compute(p.values, grouping, weights, log.p, FUN=compute_grouped_simes)
}
