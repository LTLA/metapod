# This tests the Stouffer-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-stouffer.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

REF <- function(..., weights=NULL) {
    Q <- qnorm(rbind(...))
    if (is.null(weights)) weights <- rep(1, nrow(Q))
    Q <- colSums(Q * weights)/sqrt(sum(weights))
    pnorm(Q)
}

test_that("parallelStouffer works correctly", {
    pout <- parallelStouffer(list(p1, p2, p3))
    expect_equal(pout$p.value, REF(p1, p2, p3)) 
    expect_equal(pout$representative, max.col(-cbind(p1, p2, p3)))
    expect_true(all(vapply(pout$influential, all, TRUE)))

    parallelTester(p1, p2, p3, FUN=parallelStouffer)

    # Handles ties correctly.
    pout <- parallelStouffer(list(p1, p2, p1))
    expect_equal(pout$p.value, REF(p1, p2, p1))

    pout <- parallelStouffer(list(p1, p1, p1))
    expect_equal(pout$p.value, REF(p1, p1, p1))

    # Works correctly with weights.
    pout <- parallelStouffer(list(p1, p2, p3), c(3,1,2))
    expect_equal(pout$p.value, REF(p1, p2, p3, weights=c(3,1,2)))
    parallelTesterWithWeights(p1, p2, p3, FUN=parallelStouffer)

    # Behaves sensibly at edge cases.
    expect_equal(parallelStouffer(list(0, 0))$p.value, 0)
    expect_equal(parallelStouffer(list(0, 1))$p.value, 1)
    expect_equal(parallelStouffer(list(1, 1))$p.value, 1)
})
