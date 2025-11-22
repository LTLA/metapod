#' Combine grouped p-values with Pearson's method
#'
#' Combine p-values from grouped tests with Pearson's method.
#' Groups are defined according to unique levels of a grouping factor.
#'
#' @inheritParams groupedSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
#' This contains the Pearson combined p-value for each group, log-transformed if \code{log.p=TRUE}.
#' Each entry is named according to the group.
#' \item \code{representative}, a named integer vector containing the indices of the representative tests for each group.
#' Each index refers to an entry of \code{p.values} and is named according to its group.
#' \item \code{influential}, a logical vector of length equal to \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
#' }
#'
#' @inherit parallelPearson details
#' @inherit parallelPearson references
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' g <- sample(10, length(p1), replace=TRUE)
#'
#' # Standard application:
#' out <- groupedPearson(p1, g)
#' str(out)
#'
#' # With log p-values. 
#' out <- groupedPearson(log(p1), g, log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{parallelPearson}}, for a version that operates on parallel vectors of p-values.
#'
#' \code{\link{groupedFisher}} and \code{\link{groupedStouffer}}, for different approaches to testing a joint null of independent hypotheses.
#' @export
groupedPearson <- function(p.values, grouping, log.p=FALSE) {
    .grouped_compute(p.values, grouping, NULL, log.p=log.p, FUN=compute_grouped_pearson)
}
