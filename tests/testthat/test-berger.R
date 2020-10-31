# This tests the Berger-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-berger.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

test_that("parallelBerger works correctly", {
    pout <- parallelBerger(list(p1, p2, p3))
    expect_equal(pout$p.value, pmax(p1, p2, p3))
    expect_equal(pout$representative, max.col(cbind(p1, p2, p3)))
    expect_true(all(vapply(pout$influential, all, TRUE)))

    parallelTester(p1, p2, p3, FUN=parallelBerger)

    # Handles ties correctly.
    pout <- parallelBerger(list(p1, p2, p1))
    expect_equal(pout$p.value, pmax(p1, p2))

    pout <- parallelBerger(list(p1, p1, p1))
    expect_equal(pout$p.value, p1)

    # Behaves sensibly at edge cases.
    expect_equal(parallelBerger(list(0, 0))$p.value, 0)
    expect_equal(parallelBerger(list(0, 1))$p.value, 1)
    expect_equal(parallelBerger(list(1, 1))$p.value, 1)
})
