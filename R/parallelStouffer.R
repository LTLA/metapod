#' Combine p-values with Stouffer's Z-score
#'
#' Combine p-values from parallel tests with Stouffer's Z-score method.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#' The function processes all vectors \dQuote{in parallel} - hence the name.
#'
#' @inheritParams parallelSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined Stouffer p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#'
#' @details
#' The joint null hypothesis for each group is that all of the individual null hypotheses are true.
#' Stouffer's method combines information from all individual nulls to determine if the joint null should be rejected.
#' This serves as a compromise between Fisher's method (sensitive to the smallest p-value) and Pearson's method (sensitive to the largest p-value).
#'
#' Stouffer's method is only applicable for independent tests.
#' Weights are supported by scaling the contribution of each individual null to the summed Z-score.
#' In this manner, more highly weighted tests will have a greater effect on the final combined p-value.
#'
#' The representative test for each group is defined as the test with the most negative weighted Z-score, 
#' as this has the greatest effect on the combined p-value when the joint null is rejected. 
#' All tests for each group are considered to be influential as increasing any of them (e.g., to unity) would result in a larger combined p-value.
#'
#' When a group contains both zero and unity p-values, we compare the sum of weights for all zero p-values and all unity p-values.
#' If the former is larger, the combined p-value is set to zero; if the latter, to unity.
#' If they are equal, we pretend that the two sets of p-values \dQuote{cancel out} and contribute nothing to the summed Z-score.
#' This is not entirely rigorous but provides reasonable output in the presence of such boundary values.
#'
#' @inheritSection parallelSimes Handling weights
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelStouffer(list(p1, p2, p3))
#' str(out)
#'
#' # With weights:
#' out <- parallelStouffer(list(p1, p2, p3), weights=5:7)
#' str(out)
#'
#' # With log p-values. 
#' out <- parallelStouffer(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @references
#' Stouffer SA et al. (1949).
#' \emph{The American Soldier, Vol. 1 - Adjustment during Army Life}.
#' Princeton University Press (Princeton).
#'
#' Whitlock MC (2005).
#' Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach.
#' \emph{J. Evol. Biol.} 18, 5:1368-73.
#' 
#' @export
parallelStouffer <- function(p.values, weights=NULL, log.p=FALSE) {
    compute_parallel_stouffer(p.values, weights, log.p)
}
