#' Summarize overall direction of parallel tests
#'
#' Summarize the overall direction of parallel tests in a meta-analysis,
#' based on the influential tests defined by one of the \code{parallel*} functions.
#'
#' @param effects A list of numeric vectors of the same length, containing the effect sizes for each set of tests.
#' Each group of tests is defined as corresponding entries across vectors.
#' @param influential A list of logical vectors with the same structure as \code{effects},
#' specifying the tests that are influential to the final per-group p-value.
#' @param threshold Numeric scalar defining the threshold at which an effect is \dQuote{"up"} or \dQuote{"down"}.
#'
#' @return A character vector of length equal to the number of groups.
#' Each entry can be:
#' \itemize{
#' \item \code{"up"}, if all influential tests have effects above \code{threshold}.
#' \item \code{"down"}, if all influential tests have effects below \code{threshold}.
#' \item \code{"none"}, if all influential tests have effects equal to \code{threshold}.
#' \item \code{"mixed"}, if there are influential tests with effects above and below \code{threshold}.
#' }
#'
#' @details
#' We focus on the direction of effect for the influential tests that actually contribute to a group's final p-value.
#' For example, if we did our meta-analysis using \code{\link{parallelSimes}},we are not particularly concerned about the direction of tests with large p-values.
#' Thus, we just ignore them when summarizing the group's direction in this function.
#' Otherwise, we would unnecessarily obtain a \code{mixed} direction of effect if a test with a large p-value had a weakly opposing effect.
#'
#' Of course, the interpretation of \dQuote{influential} really depends on the choice of meta-analysis strategy.
#' It is also possible that this function reports a single direction when the group really is mixed,
#' e.g., if the tests with the lowest p-values are changing in one direction but tests with weaker but still interesting effects are changing in the other direction.
#' The extent to which this is of interest is left to the discretion of the user. 
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{parallelSimes}} and related \code{parallel*} functions, to obtain \code{influential}.
#'
#' \code{\link{summarizeGroupedDirection}}, for the equivalent function based on a grouping factor.
#'
#' \code{\link{countParallelDirection}}, to count the number of effects in each direction.
#'
#' @examples
#' p1 <- rbeta(100, 0.5, 1)
#' eff1 <- rnorm(100)
#' p2 <- rbeta(100, 0.5, 1)
#' eff2 <- rnorm(100)
#'
#' out <- parallelSimes(list(p1, p2))
#' (dir <- summarizeParallelDirection(list(eff1, eff2), out$influential))
#'
#' @export
summarizeParallelDirection <- function(effects, influential, threshold=0) {
    .valid_parallel_effects(effects)
    .valid_parallel_influential(influential)
    .valid_summary_threshold(threshold)

    i <- summarize_parallel_direction(effects, influential, threshold)
    c("none", "down", "up", "mixed")[i]
}
