# This tests the summarization of the direction.
# library(testthat); library(metapod); source("test-summarize.R")

test_that("summarizeParallelDirection works as expected", {
    influential <- list(
        rbinom(100, 1, 0.5) > 0,
        rbinom(100, 1, 0.5) > 0,
        rbinom(100, 1, 0.5) > 0
    )

    effects <- list(
        rnorm(100),
        rnorm(100),
        rnorm(100)
    )

    out <- summarizeParallelDirection(effects, influential)

    all.I <- do.call(cbind, influential)
    all.E <- do.call(cbind, effects)

    up.E <- rowSums(all.I * (all.E > 0))
    down.E <- rowSums(all.I * (all.E < 0))

    ref <- c("none", "down", "up", "mixed")[1 + (down.E > 0) + (up.E > 0)*2]
    expect_identical(ref, out)

    # Works with edge cases.
    expect_identical(
        summarizeParallelDirection(effects[0], influential[0]),
        character(0)
    )

    # Handles error states.
    expect_error(summarizeParallelDirection(effects, influential[1]), "same structure")
    expect_error(summarizeParallelDirection(effects, lapply(influential, head)), "same structure")
})

test_that("summarizeGroupedDirection works as expected", {
    influential <- rbinom(100, 1, 0.5) > 0
    effects <- rnorm(100)
    g <- sample(50, 100, replace=TRUE)
    out <- summarizeGroupedDirection(effects, grouping=g, influential=influential)

    # Comparing to the parallel method.
    by.g <- split(seq_along(g), g)
    for (x in names(by.g)) {
        current <- by.g[[x]]
        by.g[[x]] <- summarizeParallelDirection(as.list(effects[current]), influential=as.list(influential[current]))
    }

    expect_identical(unlist(by.g), out)

    # Handles error states. 
    expect_error(summarizeGroupedDirection(effects, rle(g[1:10]), influential), "runs.*not the same")
    expect_error(summarizeGroupedDirection(effects[1:10], rle(g), influential[1:10]), "runs.*not the same")
    expect_error(summarizeGroupedDirection(effects[1:10], g[1:10], influential), "not the same")
    expect_error(summarizeGroupedDirection(effects, g, influential[1:10]), "not the same")
})
