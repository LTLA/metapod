#' Combine p-values with Wilkinson's method
#'
#' Combine p-values from parallel tests with Wilkinson's method.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#' The function processes all vectors \dQuote{in parallel} - hence the name.
#'
#' @inheritParams parallelHolmMin
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer vector specifying the representative test in each group.
#' Specifically, it refers to the index of \code{p.values} containing the representative test for each group,
#' i.e., the p-value for the representative test of group \code{i} is \code{p.values[representative[i]]}.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final per-group p-value.
#' }
#'
#' @details
#' Here, the joint null hypothesis for each group is all individual null hypotheses are false.
#' Rejection of the joint null is heavily favored in situations where \eqn{N} or more individual nulls are rejected.
#' This is achieved in Wilkinson's method by considering the \eqn{N}-th order statistic for uniformly distributed p-values.
#' The individual tests are assumed to be independent, and all weights are ignored.
#'
#' \eqn{N} is defined as the larger of \code{min.n} and the product of \code{min.prop} with the number of tests in the group (rounded up).
#' This allows users to scale rejection of the joint null with the size of the group, while avoiding a too-low \eqn{N} when the group is small.
#' Note that \eqn{N} is always capped at the total size of the group.
#'
#' The representative test for each group is defined as that with the \eqn{N}-th smallest p-value, as this is used to compute the combined p-value.
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
#' out <- parallelWilkinson(list(p1, p2, p3))
#' str(out)
#'
#' # With log p-values. 
#' out <- parallelWilkinson(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @seealso
#' \code{\link{groupedWilkinson}}, for a version that combines p-values based on a grouping factor.
#'
#' \code{\link{parallelHolmMin}}, which has a more strict interpretation of \emph{N}, amongst other things.
#'
#' @references
#' Wilkinson B (1951).
#' A statistical consideration in psychological research.
#' \emph{Psychol. Bull.} 48, 156-158.
#'
#' @export
parallelWilkinson <- function(p.values, log.p=FALSE, min.n=1, min.prop=0.5) {
    .valid_logp(log.p)
    .valid_parallel_pvalues(p.values)
    .valid_min_vals(min.n, min.prop)

    compute_parallel_wilkinson(p.values, NULL, log.p, min.n, min.prop)
}
