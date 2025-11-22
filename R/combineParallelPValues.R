#' Combine parallel p-values
#'
#' Combine p-values from parallel hypothesis tests using a variety of meta-analysis methods.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#' The function processes all vectors \dQuote{in parallel} - hence the name.
#' 
#' @inheritParams parallelHolmMin
#' @param method String specifying the method to use to combine p-values.
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the Simes p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer vector specifying the representative test in each group.
#' Specifically, it refers to the index of \code{p.values} containing the representative test for each group,
#' i.e., the p-value for the representative test of group \code{i} is \code{p.values[representative[i]]}.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#'
#' @details
#' \code{min.prop} and \code{min.n} only have an effect for \code{method="wilkinson"} and \code{"holm-min"}.
#'
#' \code{weights} only has an effect for \code{method="simes"}, \code{"holm-min"} and \code{"stouffer"}.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' p1 <- runif(10000)
#' p2 <- runif(10000)
#' p3 <- runif(10000)
#' 
#' fish <- combineParallelPValues(list(p1, p2, p3), method="fisher")
#' hist(fish$p.value)
#' 
#' z <- combineParallelPValues(list(p1, p2, p3), method="stouffer", weights=1:3)
#' hist(z$p.value)
#' 
#' simes <- combineParallelPValues(list(p1, p2, p3), method="simes")
#' hist(simes$p.value)
#' 
#' berger <- combineParallelPValues(list(p1, p2, p3), method="berger")
#' hist(berger$p.value)
#'
#' @export
combineParallelPValues <- function(p.values,
    method=c("simes", "holm-min", "berger", "fisher", "pearson", "wilkinson", "stouffer"),
    weights=NULL, log.p=FALSE, min.n=1, min.prop=0.5) 
{
    method <- match.arg(method)
    all.args <- list(p.values=p.values, weights=weights, log.p=log.p, min.n=min.n, min.prop=min.prop) 

    FUN <- switch(method,
        simes=parallelSimes,
        `holm-min`=parallelHolmMin,
        berger=parallelBerger,
        fisher=parallelFisher,
        pearson=parallelPearson,
        wilkinson=parallelWilkinson,
        stouffer=parallelStouffer
    )

    keep <- names(formals(FUN))
    do.call(FUN, all.args[keep])
}
