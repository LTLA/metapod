# This tests the Simes-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-simes.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

test_that("parallelSimes works correctly", {
    pout <- parallelSimes(list(p1, p2, p3))
    expect_equal(pout$p.value, apply(cbind(p1, p2, p3), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 
    parallelTester(p1, p2, p3, FUN=parallelSimes)

    # Handles ties correctly.
    pout <- parallelSimes(list(p1, p2, p1))
    expect_equal(pout$p.value, apply(cbind(p1, p2, p1), 1, FUN=function(p) { min(p.adjust(p, method="BH")) })) 

    pout <- parallelSimes(list(p1, p1, p1))
    expect_equal(pout$p.value, p1)

    # Works correctly with weights.
    pout <- parallelSimes(list(p1, p2, p3), c(3,1,2))
    ref <- parallelSimes(rep(list(p1, p2, p3), c(3,1,2)))
    expect_identical(pout$p.value, ref$p.value)
    parallelTesterWithWeights(p1, p2, p3, FUN=parallelSimes)

    # Behaves sensibly at edge cases.
    expect_equal(parallelSimes(list(0, 0))$p.value, 0)
    expect_equal(parallelSimes(list(0, 1))$p.value, 0)
    expect_equal(parallelSimes(list(1, 1))$p.value, 1)
})

test_that("groupedSimes works correctly", {
    g <- sample(100, length(p1), replace=TRUE)
    groupedTester(p1, g, pFUN=parallelSimes, gFUN=groupedSimes)
    groupedTesterWithWeights(p1, g,  pFUN=parallelSimes, gFUN=groupedSimes)
})
