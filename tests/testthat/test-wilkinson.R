# This tests the Holm-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-wilkinson.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

REFFUN <- function(x, r) {
    rp <- apply(x, 1, function(p) sort(p)[r])
    pbeta(rp, r, ncol(x) - r + 1)
}

test_that("parallelWilkinson method works correctly", {
    parallelTester(p1, p2, p3, FUN=parallelWilkinson)

    # Testing alternative min.num, min.prop settings.
    p4 <- runif(1000)    
    together <- cbind(p1, p2, p3, p4)

    for (i in seq_len(4)) {
        pout <- parallelWilkinson(list(p1, p2, p3, p4), min.n=i, min.prop=0)
        expect_equal(pout$p.value, REFFUN(together, r=i))
        expect_equal(pout$representative, apply(together, 1, FUN=function(p) { order(p)[i] })) 
        expect_equal(do.call(rbind, pout$influential), apply(together, 1, FUN=function(p) { 1:4 %in% order(p)[1:i] })) 
    }

    for (p in c(0.2, 0.25, 0.45, 0.65, 0.75, 0.85)) {
        pout <- parallelWilkinson(list(p1, p2, p3, p4), min.prop=p)
        N <- max(ceiling(p * ncol(together)), 1)
        expect_equal(pout$p.value, REFFUN(together, r=N))
        expect_equal(pout$representative, apply(together, 1, FUN=function(p) { order(p)[N] })) 
        expect_equal(do.call(rbind, pout$influential), apply(together, 1, FUN=function(p) { 1:4 %in% order(p)[1:N] })) 
    }

    # Handles ties correctly.
    pout <- parallelWilkinson(list(p1, p2, p1))
    expect_equal(pout$p.value, REFFUN(cbind(p1, p2, p1), 2))

    pout <- parallelWilkinson(list(p1, p1, p1))
    expect_equal(pout$p.value, REFFUN(cbind(p1, p1, p1), 2))

    # Behaves sensibly at edge cases.
    expect_equal(parallelWilkinson(list(0, 0))$p.value, 0)
    expect_equal(parallelWilkinson(list(0, 1))$p.value, 0)
    expect_equal(parallelWilkinson(list(1, 1))$p.value, 1)
})

test_that("groupedWilkinson works correctly", {
    g <- sample(100, length(p1), replace=TRUE)
    groupedTester(p1, g, pFUN=parallelWilkinson, gFUN=groupedWilkinson)
})
