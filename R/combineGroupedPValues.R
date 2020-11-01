.prepare_grouped_inputs <- function(grouping, x, is.rle=FALSE) {
    if (!is.rle) {
        o <- order(grouping)
        x <- lapply(x, "[", i=o)
        runs <- rle(grouping[o])
    } else {
        o <- seq_len(sum(grouping$lengths))
        runs <- grouping
    }
    
    list(order=o, x=x, runs=runs)
}

.grouped_compute <- function(p.values, grouping, weights=NULL, ..., FUN, is.rle=FALSE) {
    gout <- .prepare_grouped_inputs(grouping, list(p.values, weights), is.rle=is.rle)
    o <- gout$order
    p.values <- gout$x[[1]]
    weights <- gout$x[[2]]
    runs <- gout$runs

    output <- FUN(p.values, runs$lengths, weights, ...)
    names(output$p.value) <- names(output$representative) <- runs$values
    output$influential[o] <- output$influential

    output
}
