#' Combine parallel p-values with Pearson's method
#'
#' Combine p-values from parallel tests with Pearson's method.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#'
#' @inheritParams parallelSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined Pearson p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' Specifically, this refers to the index of the \emph{vector} of \code{p.values} containing the representative test.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#'
#' @details
#' Here, the joint null hypothesis for each group is that all of the individual null hypotheses are true.
#' Pearson's method combines information from all individual nulls to determine if the joint null should be rejected.
#' Compared to Stouffer's and Pearson's methods, Pearson's method is more sensitive to the largest individual p-value.
#' This method is only applicable to independent tests and no weights are considered.
#'
#' The representative test for each group is defined as the test with the largest p-value, as this has the greatest effect on the combined p-value. 
#' All tests for each group are considered to be influential as increasing any of them (e.g., to unity) would result in a larger combined p-value.
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelPearson(list(p1, p2, p3))
#' str(out)
#'
#' # With log p-values. 
#' out <- parallelPearson(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{groupedPearson}}, for a version that combines p-values based on a grouping factor.
#'
#' \code{\link{parallelFisher}} and \code{\link{parallelStouffer}}, for different approaches to testing a joint null of independent hypotheses.
#'
#' @references
#' Pearson K (1934).
#' On a new method of deternining \dQuote{goodness of fit.}
#' \emph{Biometrika} 26, 425-442.
#'
#' @export
parallelPearson <- function(p.values, log.p=FALSE) {
    .valid_logp(log.p)
    .valid_parallel_pvalues(p.values)

    compute_parallel_pearson(p.values, NULL, log.p)
}
