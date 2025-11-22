.prepare_grouped_inputs <- function(grouping, x) {
    if (inherits(grouping, "rle")) {
        o <- seq_len(sum(grouping$lengths))
        runs <- grouping
    } else {
        o <- order(grouping)
        if (is.factor(grouping)) {
            tab <- table(grouping) # this should respect the order.
            stopifnot(identical(levels(grouping), names(tab)))
            runs <- list(lengths=as.integer(tab), values=names(tab))
            class(runs) <- "rle"
        } else {
            runs <- rle(grouping[o])
        }

        for (i in seq_along(x)) {
            if (!is.null(x[[i]])) {
                if (length(x[[i]])!=length(grouping)) {
                    stop("lengths of 'grouping' and '", names(x)[i], "' are not the same")
                }
                x[[i]] <- x[[i]][o]
            }
        }
    }
    
    list(order=o, x=x, runs=runs)
}

.grouped_compute <- function(p.values, grouping, weights=NULL, ..., log.p, FUN) {
    .valid_logp(log.p)
    stopifnot(is.numeric(p.values))
    if (!is.null(weights)) {
        stopifnot(is.numeric(weights))
    }

    gout <- .prepare_grouped_inputs(grouping, list(p.values=p.values, weights=weights))
    o <- gout$order
    p.values <- gout$x[[1]]
    weights <- gout$x[[2]]
    runs <- gout$runs

    output <- FUN(p.values, runs$lengths, weights, ..., log=log.p)
    output$representative <- o[output$representative]
    output$influential[o] <- output$influential
    names(output$p.value) <- names(output$representative) <- as.character(runs$values)

    output
}

#' Combine grouped p-values
#'
#' Combine p-values from grouped hypothesis tests using a variety of meta-analysis methods.
#' Each group of p-values is defined as those assigned to the same level of the grouping factor.
#' 
#' @inheritParams groupedHolmMin
#' @param method String specifying the method to use to combine p-values.
#' 
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
#' This contains the combined p-value for each group, log-transformed if \code{log.p=TRUE}.
#' Each entry is named according to the group.
#' \item \code{representative}, a named integer vector containing the indices of the representative tests for each group.
#' Each index refers to an entry of \code{p.values} and is named according to its group.
#' \item \code{influential}, a logical vector of length equal to \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
#' }
#'
#' @inherit combineParallelPValues details
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' p <- runif(10000)
#' g <- sample(100, 10000, replace=TRUE)
#'
#' fish <- combineGroupedPValues(p, g, method="fisher")
#' hist(fish$p.value)
#' 
#' z <- combineGroupedPValues(p, g, method="stouffer", weights=rexp(10000))
#' hist(z$p.value)
#' 
#' simes <- combineGroupedPValues(p, g, method="simes")
#' hist(simes$p.value)
#' 
#' berger <- combineGroupedPValues(p, g, method="berger")
#' hist(berger$p.value)
#'
#' @export
combineGroupedPValues <- function(p.values, grouping, 
    method=c("simes", "holm-min", "berger", "fisher", "pearson", "wilkinson", "stouffer"),
    weights=NULL, log.p=FALSE, min.n=1, min.prop=0.5)
{
    method <- match.arg(method)
    all.args <- list(p.values=p.values, grouping=grouping, weights=weights, log.p=log.p, min.n=min.n, min.prop=min.prop)

    FUN <- switch(method,
        simes=groupedSimes,
        `holm-min`=groupedHolmMin,
        berger=groupedBerger,
        fisher=groupedFisher,
        pearson=groupedPearson,
        wilkinson=groupedWilkinson,
        stouffer=groupedStouffer
    )

    keep <- names(formals(FUN))
    do.call(FUN, all.args[keep])
}
