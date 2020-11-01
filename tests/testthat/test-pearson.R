# This tests the Pearson-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-pearson.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

REFFUN <- function(...) {
    x <- cbind(...)
    Y <- rowSums(log(1-x))
    pchisq(-2*Y, df=2*ncol(x), lower.tail=TRUE)
}

test_that("parallelPearson works correctly", {
    pout <- parallelPearson(list(p1, p2, p3))

    expect_equal(pout$p.value, REFFUN(p1, p2, p3))
    expect_equal(pout$representative, max.col(cbind(p1, p2, p3)))
    expect_true(all(vapply(pout$influential, all, TRUE)))

    parallelTester(p1, p2, p3, FUN=parallelPearson)

    # Handles ties correctly.
    pout <- parallelPearson(list(p1, p2, p1))
    expect_equal(pout$p.value, REFFUN(p1, p1, p2))

    pout <- parallelPearson(list(p1, p1, p1))
    expect_equal(pout$p.value, REFFUN(p1, p1, p1))

    # Behaves sensibly at edge cases.
    expect_equal(parallelPearson(list(0, 0))$p.value, 0)
    expect_equal(parallelPearson(list(0, 1))$p.value, 1)
    expect_equal(parallelPearson(list(1, 1))$p.value, 1)
})
