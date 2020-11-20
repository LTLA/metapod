#' Combine grouped p-values with Berger's IUT
#'
#' Combine p-values from grouped tests with Berger's intersection-union test (IUT). 
#' Groups are defined according to unique levels of a grouping factor.
#'
#' @inheritParams groupedSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
#' This contains the IUT p-value for each group, log-transformed if \code{log.p=TRUE}.
#' Each entry is named according to the group.
#' \item \code{representative}, a named integer scalar specifying the representative test for each group.
#' Each index refers to an entry of \code{p.values} and is named according to its group.
#' \item \code{influential}, a logical vector of length equal to \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
#' }
#'
#' @inherit parallelBerger details
#' @inherit parallelBerger references
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' g <- sample(10, length(p1), replace=TRUE)
#'
#' # Standard application:
#' out <- groupedBerger(p1, g)
#' str(out)
#'
#' # With log p-values. 
#' out <- groupedBerger(log(p1), g, log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{parallelBerger}}, for a version that operates on parallel vectors of p-values.
#'
#' @export
groupedBerger <- function(p.values, grouping, log.p=FALSE) {
    .grouped_compute(p.values, grouping, NULL, log.p=log.p, FUN=compute_grouped_berger)
}
