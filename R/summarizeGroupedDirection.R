#' Summarize overall direction of grouped tests
#'
#' Summarize the overall direction of grouped tests in a meta-analysis,
#' based on the influential tests defined by one of the \code{grouped*} functions.
#'
#' @param effects A numeric vector containing the effect size for each test.
#' @param grouping A vector of factor of length equal to \code{effects}, specifying the assigned group for each tests.
#'
#' Alternatively, an \link{rle} object where each run corresponds to a group and specifies the entries of \code{effects} belonging to that group.
#' This assumes that \code{effects} is ordered such that all entries in the same group are adjacent to each other.
#' @param influential A logical vector of length equal to \code{effects},
#' indicating whether each test is influential in its assigned group.
#' @inheritParams summarizeParallelDirection
#'
#' @inherit summarizeParallelDirection return
#' @inherit summarizeParallelDirection details
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{groupedSimes}} and related \code{grouped*} functions, to obtain \code{influential}.
#'
#' \code{\link{summarizeParallelDirection}}, for the equivalent function based on parallel tests.
#'
#' @examples
#' p <- rbeta(100, 0.5, 1)
#' eff <- rnorm(100)
#' g <- sample(20, 100, replace=TRUE)
#'
#' out <- groupedSimes(p, g)
#' (dir <- summarizeGroupedDirection(eff, g, out$influential))
#' @export
summarizeGroupedDirection <- function(effects, grouping, influential, threshold=0) {
    gout <- .prepare_grouped_inputs(grouping, list(effects, influential))
    p.values <- gout$x[[1]]
    weights <- gout$x[[2]]
    runs <- gout$runs

    output <- summarize_grouped_direction(effects, runs$lengths, influential, threshold) 
    output <- c("none", "down", "up", "mixed")[output]
    names(output) <- runs$values
    output
}
