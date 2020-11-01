.grouped_compute <- function(p.values, grouping, weights=NULL, log=FALSE, FUN) {
    o <- order(grouping)
    p.values <- p.values[o]
    weights <- weights[o]
    runs <- rle(grouping[o])

    output <- FUN(p.values, runs$length, weights, log)
    names(output$p.value) <- names(output$representative) <- runs$value
    output$influential[o] <- output$influential

    output
}
