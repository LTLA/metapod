# This tests the averaging functions.
# library(testthat); library(metapod); source("test-average.R")

test_that("parallel averaging works as expected", {
    p <- list(rnorm(10), rnorm(10))
    out <- averageParallelStats(p)
    expect_equal(out, (p[[1]] + p[[2]])/2)

    # Works for different weight settings.
    p <- list(rnorm(10), rnorm(10))
    out <- averageParallelStats(p, weights=1:2)
    expect_equal(out, (p[[1]] + p[[2]]*2)/3)

    # Works for list weights.
    p <- list(rnorm(10), rnorm(10))
    w <- list(runif(10), runif(10))
    out <- averageParallelStats(p, weights=w)
    expect_equal(out, (p[[1]] * w[[1]] + p[[2]] * w[[2]])/(w[[1]] + w[[2]]))

    # Handles missing values appropriately.
    p <- list(rnorm(10), rnorm(10))
    p[[1]][c(1, 3)] <- NA
    p[[2]][c(2, 3)] <- NA

    out <- averageParallelStats(p)
    expect_equal(out[1], p[[2]][1])
    expect_equal(out[2], p[[1]][2])
    expect_true(is.na(out[3]))
    expect_identical(out[-(1:3)], averageParallelStats(lapply(p, tail, -3)))

    p <- list(rnorm(10), rnorm(10), rnorm(10))
    p[[1]][c(1, 3)] <- NA
    p[[2]][c(2, 3)] <- NA
    w <- list(runif(10), runif(10), runif(10))

    out <- averageParallelStats(p, weights=w)
    expect_equal(out[1], averageParallelStats(p[-1], w[-1])[1])
    expect_equal(out[2], averageParallelStats(p[-2], w[-2])[2])
    expect_equal(out[3], p[[3]][3])
    expect_equal(out[-(1:3)], (Reduce("+", mapply("*", p, w, SIMPLIFY=FALSE))/Reduce("+", w))[-(1:3)])

    # Handles errors properly.
    expect_error(averageParallelStats(list(rnorm(10), 1)), "same length")
    p <- list(rnorm(10), rnorm(10))
    expect_error(averageParallelStats(p, weights=c(-1, 0)), "must be positive")
    expect_error(averageParallelStats(p, weights=1), "not TRUE")
    expect_error(averageParallelStats(p, weights=list(1,1)), "same 'lengths'")
})

test_that("groupe averaging works as expected", {
    p <- rnorm(100)
    g <- sample(LETTERS, length(p), replace=TRUE)

    out <- averageGroupedStats(p, g)
    expect_equal(out, vapply(split(p, g), mean, 0))

    # Handles weighting.
    w <- runif(100)
    out <- averageGroupedStats(p, g, weights=w)
    by.group <- split(seq_along(g), g)
    expect_equal(out, vapply(by.group, function(i) weighted.mean(p[i], w[i]), 0))

    # Handles missing values.
    thrown <- sample(length(p), length(p)/2)
    p2 <- p
    p2[-thrown] <- NA_real_

    out <- averageGroupedStats(p2, g) 
    expect_identical(names(out), sort(unique(g)))

    ref <- averageGroupedStats(p2[thrown], g[thrown])
    common <- names(ref)
    expect_identical(out[common], ref[common])
    expect_true(all(is.na(out[setdiff(names(out), names(ref))])))

    # Missing values and weights interact properly.
    out <- averageGroupedStats(p2, g, weights=w) 
    ref <- averageGroupedStats(p2[thrown], g[thrown], weights=w[thrown])
    common <- names(ref)
    expect_identical(out[common], ref[common])
    expect_true(all(is.na(out[setdiff(names(out), names(ref))])))
})

