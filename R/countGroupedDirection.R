#' Count directions of grouped tests
#'
#' Count the number of grouped tests that are significant and changing in each direction.
#' Each group of tests is defined according to a grouping factor.
#'
#' @param p.values A numeric vector containing p-values for individual tests.
#' @param grouping A vector or factor of length equal to \code{p.values}, specifying the group to which each test is assigned.
#'
#' Alternatively, an \link{rle} object where each run corresponds to a group and specifies the entries of \code{p.values} belonging to that group.
#' This assumes that \code{p.values} is ordered such that all entries in the same group are adjacent to each other.
#' @param effects A numeric vector of length equal to \code{p.values}, containing the effect size for each test.
#' @inheritParams countParallelDirection
#'
#' @inherit countParallelDirection return
#' @inherit countParallelDirection details
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{summarizeGroupedDirection}}, for another way to summarize the overall direction into a single string.
#'
#' \code{\link{countParallelDirection}}, for the equivalent function operating on parallel tests.
#'
#' @examples
#' p <- rbeta(100, 0.5, 1)
#' eff <- rnorm(100)
#' g <- sample(20, 100, replace=TRUE)
#'
#' (dirs <- countGroupedDirection(p, g, eff))
#'
#' @export
countGroupedDirection <- function(p.values, grouping, effects, p.threshold=0.05, effect.threshold=0, method=c("BH", "holm"), log.p=FALSE) {
    stopifnot(is.numeric(p.values))
    stopifnot(is.numeric(effects))
    .valid_count_thresholds(p.threshold, effect.threshold)
    method <- match.arg(method)
    .valid_logp(log.p)

    gout <- .prepare_grouped_inputs(grouping, list(p.values=p.values, effects=effects))
    p.values <- gout$x[[1]]
    effects <- gout$x[[2]]
    runs <- gout$runs

    out <- count_grouped_direction(p.values, runs$lengths, effects, 
        method=c(BH=0L, holm=1L)[method], 
        pthreshold=p.threshold, ethreshold=effect.threshold, log=log.p)
    names(out$up) <- names(out$down) <- runs$values

    out
}
