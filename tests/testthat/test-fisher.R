# This tests the Fisher-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-fisher.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

REFFUN <- function(...) {
    x <- cbind(...)
    Y <- rowSums(log(x))
    pchisq(-2*Y, df=2*ncol(x), lower.tail=FALSE)
}

test_that("parallelFisher works correctly", {
    pout <- parallelFisher(list(p1, p2, p3))

    expect_equal(pout$p.value, REFFUN(p1, p2, p3))
    expect_equal(pout$representative, max.col(-cbind(p1, p2, p3)))
    expect_true(all(vapply(pout$influential, all, TRUE)))

    parallelTester(p1, p2, p3, FUN=parallelFisher)

    # Handles ties correctly.
    pout <- parallelFisher(list(p1, p2, p1))
    expect_equal(pout$p.value, REFFUN(p1, p1, p2))

    pout <- parallelFisher(list(p1, p1, p1))
    expect_equal(pout$p.value, REFFUN(p1, p1, p1))

    # Behaves sensibly at edge cases.
    expect_equal(parallelFisher(list(0, 0))$p.value, 0)
    expect_equal(parallelFisher(list(0, 1))$p.value, 0)
    expect_equal(parallelFisher(list(1, 1))$p.value, 1)
})

test_that("groupedFisher works correctly", {
    g <- sample(100, length(p1), replace=TRUE)
    groupedTester(p1, g, pFUN=parallelFisher, gFUN=groupedFisher)
})
