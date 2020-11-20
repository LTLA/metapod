#' Combine parallel p-values with Berger's IUT
#'
#' Combine p-values from parallel tests with Berger's intersection-union test (IUT).
#' Each group of p-values is defined from the corresponding entries across all vectors.
#'
#' @inheritParams parallelSimes
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the IUT p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' Specifically, this refers to the index of the \emph{vector} of \code{p.values} containing the representative test.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#'
#' @details
#' The joint null hypothesis for each group is that \emph{any} of the individual null hypotheses are true.
#' Berger's IUT will only reject the joint null if all of the individual nulls are rejected.
#' This method is applicable under arbitrary dependency structures.
#' No weights are considered.
#'
#' The representative test for each group is defined as the test with the largest p-value, as this is ultimately used as the IUT p-value.
#' All tests for each group are considered to be influential as increasing any of them (e.g., to unity) would result in a larger combined p-value.
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelBerger(list(p1, p2, p3))
#' str(out)
#'
#' # With log p-values. 
#' out <- parallelBerger(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @references
#' Berger RL and Hsu JC (1996).
#' Bioequivalence trials, intersection-union tests and equivalence confidence sets.
#' \emph{Statist. Sci.} 11, 283-319.
#' 
#' @seealso
#' \code{\link{groupedBerger}}, for the version that combines p-values based on a grouping factor.
#'
#' @export
parallelBerger <- function(p.values, log.p=FALSE) {
    .valid_logp(log.p)
    .valid_parallel_pvalues(p.values)

    compute_parallel_berger(p.values, NULL, log.p)
}
