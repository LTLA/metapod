# This tests the Holm-related functions.
# library(metapod); library(testthat); source("setup.R"); source("test-holm-min.R")

set.seed(20000)
p1 <- runif(1000)
p2 <- runif(1000)
p3 <- runif(1000)    

test_that("parallelHolmMin method works correctly", {
    parallelTester(p1, p2, p3, FUN=parallelHolmMin)

    # Testing alternative min.num, min.prop settings.
    p4 <- runif(1000)    
    together <- cbind(p1, p2, p3, p4)

    for (i in seq_len(4)) {
        pout <- parallelHolmMin(list(p1, p2, p3, p4), min.n=i, min.prop=0)
        expect_equal(pout$p.value, apply(together, 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[i] })) 
        expect_equal(pout$representative, apply(together, 1, FUN=function(p) { order(p)[i] })) 
        expect_equal(do.call(rbind, pout$influential), apply(together, 1, FUN=function(p) { 1:4 %in% order(p)[1:i] })) 
    }

    for (p in c(0.2, 0.25, 0.45, 0.65, 0.75, 0.85)) {
        pout <- parallelHolmMin(list(p1, p2, p3, p4), min.prop=p)
        N <- max(ceiling(p * ncol(together)), 1)
        expect_equal(pout$p.value, apply(together, 1, FUN=function(p) { sort(p.adjust(p, method="holm"))[N] })) 
        expect_equal(pout$representative, apply(together, 1, FUN=function(p) { order(p)[N] })) 
        expect_equal(do.call(rbind, pout$influential), apply(together, 1, FUN=function(p) { 1:4 %in% order(p)[1:N] })) 
    }

    # Handles ties correctly.
    pout <- parallelHolmMin(list(p1, p2, p1))
    expect_equal(pout$p.value, apply(cbind(p1, p2, p1), 1, FUN=function(p) { median(p.adjust(p, method="holm")) })) 

    pout <- parallelHolmMin(list(p1, p1, p1))
    expect_equal(pout$p.value, apply(cbind(p1, p1, p1), 1, FUN=function(p) { median(p.adjust(p, method="holm")) })) 

    # Behaves as a special case of lowest-Bonferroni.
    pout <- parallelHolmMin(list(p1, p2, p3), min.prop=0)
    expect_equal(pout$p.value, pmin(p1, p2, p3, 1/3)*3)

    # Handles weights correctly.
    parallelTesterWithWeights(p1, p2, p3, FUN=parallelHolmMin)

    for (i in seq_len(4)) {
        w <- runif(4)
        pout <- parallelHolmMin(list(p1, p2, p3, p4), min.n=i, min.prop=0, weights=w)
        expect_equal(pout$p.value, apply(together, 1, FUN=function(p) { 
            p <- p/w
            o <- order(p)
            wx <- sum(w) - c(0, head(cumsum(w[o]), -1))
            out <- cummax(p[o] * wx)[i]
            min(out, 1)
        }))
        expect_equal(pout$representative, apply(together, 1, FUN=function(p) { order(p/w)[i] })) 
        expect_equal(do.call(rbind, pout$influential), apply(together, 1, FUN=function(p) { 1:4 %in% order(p/w)[1:i] })) 
    }

    # Behaves sensibly at edge cases.
    expect_equal(parallelHolmMin(list(0, 0))$p.value, 0)
    expect_equal(parallelHolmMin(list(0, 1))$p.value, 0)
    expect_equal(parallelHolmMin(list(1, 1))$p.value, 1)
})

test_that("groupedHolmMin works correctly", {
    g <- sample(100, length(p1), replace=TRUE)
    groupedTester(p1, g, pFUN=parallelHolmMin, gFUN=groupedHolmMin)
    groupedTesterWithWeights(p1, g,  pFUN=parallelHolmMin, gFUN=groupedHolmMin)
})
