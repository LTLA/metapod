#' Combine parallel p-values with Fisher's method
#'
#' Combine p-values from parallel tests with Fisher's method.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#'
#' @inheritParams parallelSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined Fisher p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' Specifically, this refers to the index of the \emph{vector} of \code{p.values} containing the representative test.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#'
#' @details
#' The joint null hypothesis for each group is that all of the individual null hypotheses are true.
#' Fisher's method combines information from all individual nulls to determine if the joint null should be rejected.
#' Compared to Stouffer's and Pearson's methods, Fisher's method provides more sensitivity to the smallest individual p-value.
#' This method is only applicable to independent tests and no weights are considered.
#'
#' The representative test for each group is defined as the test with the lowest p-value, as this has the greatest effect on the combined p-value. 
#' All tests for each group are considered to be influential as increasing any of them (e.g., to unity) would result in a larger combined p-value.
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelFisher(list(p1, p2, p3))
#' str(out)
#'
#' # With log p-values. 
#' out <- parallelFisher(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{groupedFisher}}, for a version that combines p-values based on a grouping factor.
#'
#' \code{\link{parallelStouffer}} and \code{\link{parallelPearson}}, for different approaches to testing a joint null of independent hypotheses.
#'
#' @references
#' Fisher RA (1925).
#' \emph{Statistical Methods for Research Workers}.
#' Oliver and Boyd (Edinburgh).
#' @export
parallelFisher <- function(p.values, log.p=FALSE) {
    .valid_logp(log.p)
    .valid_parallel_pvalues(p.values)

    compute_parallel_fisher(p.values, NULL, log.p)
}
