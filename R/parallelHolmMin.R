#' Combine parallel p-values with the minimum Holm approach
#'
#' Combine p-values from parallel tests with the minimum Holm approach.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#'
#' @inheritParams parallelSimes
#' @param min.n Integer scalar specifying the minimum number of individual nulls to reject when testing the joint null.
#' @param min.prop Numeric scalar in [0, 1], specifying the minimum proportion of individual nulls to reject when testing the joint null.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' Specifically, this refers to the index of the \emph{vector} of \code{p.values} containing the representative test.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final per-group p-value.
#' }
#'
#' @details
#' Here, the joint null hypothesis for each group is that fewer than \eqn{N} individual null hypotheses are false.
#' The joint null is only rejected if \eqn{N} or more individual nulls are rejected; hence the \dQuote{minimum} in the function name.
#'
#' \eqn{N} is defined as the larger of \code{min.n} and the product of \code{min.prop} with the number of tests in the group (rounded up).
#' This allows users to scale rejection of the joint null with the size of the group, while avoiding a too-low \eqn{N} when the group is small.
#' Note that \eqn{N} is always capped at the total size of the group.
#'
#' To compute the combined p-value, we apply the Holm-Bonferroni correction to all p-values in the set and take the \eqn{N}-th smallest value.
#' This effectively recapitulates the step-down procedure where we reject individual nulls until we reach the \eqn{N}-th test.
#' This method works correctly in the presence of dependencies between p-values.
#' 
#' If non-equal weights are provided, they are used to effectively downscale the p-values. 
#' This aims to redistribute the error rate across the individual tests,
#' i.e., tests with lower weights are given fewer opportunities to drive acceptance of the joint null.
#'
#' The representative test for each group is defined as that with the \eqn{N}-th smallest p-value, as this is directly used as the combined p-value.
#' The influential tests for each group are defined as those with p-values no greater than the representative test's p-value.
#' This is based on the fact that increasing them (e.g., by setting them to unity) would result in a larger combined p-value.
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelHolmMin(list(p1, p2, p3))
#' str(out)
#'
#' # With vector-level weights:
#' out <- parallelHolmMin(list(p1, p2, p3), weights=c(10, 20, 30))
#' str(out)
#' 
#' # With log p-values. 
#' out <- parallelHolmMin(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @references
#' Holm S (1979).
#' A simple sequentially rejective multiple test procedure.
#' \emph{Scand. J. Stat.} 6, 65-70.
#'
#' @seealso
#' \code{\link{groupedHolmMin}}, for a version that combines p-values based on a grouping factor.
#'
#' \code{\link{parallelWilkinson}}, for a more relaxed version of this test when hypotheses are independent.
#'
#' @export
parallelHolmMin <- function(p.values, weights=NULL, log.p=FALSE, min.n=1, min.prop=0.5) {
    .valid_logp(log.p)
    .valid_parallel_pvalues(p.values)
    .valid_min_vals(min.n, min.prop)

    compute_parallel_holm_min(p.values, weights, log.p, min.n, min.prop)
}
