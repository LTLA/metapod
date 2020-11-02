# This tests the counting functions.
# library(testthat); library(metapod); source("test-count.R")

test_that("countParallelDirection works as expected", {
    p <- list(
        rbeta(1000, 0.5, 1),
        rbeta(1000, 0.5, 1),
        rbeta(1000, 0.5, 1)
    )

    effects <- list(
        rnorm(1000),
        rnorm(1000),
        rnorm(1000)
    )

    out.bh <- countParallelDirection(p, effects)
    out.holm <- countParallelDirection(p, effects, method="holm")

    outu.bh <- outd.bh <- outu.holm <- outd.holm <- integer(1000)
    for (i in seq_len(1000)) {
        curp <- vapply(p, "[", i=i, 0)
        cure <- vapply(effects, "[", i=i, 0)

        bh.adj <- (p.adjust(curp, method="BH") < 0.05)
        outu.bh[i] <- sum(bh.adj * (cure > 0))
        outd.bh[i] <- sum(bh.adj * (cure < 0))

        holm.adj <- (p.adjust(curp, method="holm") < 0.05)
        outu.holm[i] <- sum(holm.adj * (cure > 0))
        outd.holm[i] <- sum(holm.adj * (cure < 0))
    }

    expect_identical(out.bh$up, outu.bh)
    expect_identical(out.bh$down, outd.bh)
    expect_identical(out.holm$up, outu.holm)
    expect_identical(out.holm$down, outd.holm)

    # Works for log-values.
    lout.bh <- countParallelDirection(lapply(p, log), effects, log=TRUE)
    expect_identical(out.bh, lout.bh)

    lout.holm <- countParallelDirection(lapply(p, log), effects, method="holm", log=TRUE)
    expect_identical(out.holm, lout.holm)

    # Works with NA's.
    p2 <- p
    p2[[1]][] <- NA_real_

    out.na <- countParallelDirection(p2, effects)
    ref <- countParallelDirection(p[-1], effects[-1])
    expect_identical(ref, out.na)

    out <- countParallelDirection(list(NA, NA, NA), list(1, 0, -1))
    expect_identical(out$up, 0L)
    expect_identical(out$down, 0L)

    expect_error(countParallelDirection(p, effects[-1]), "same structure")
})

test_that("countGroupedDirection works as expected", {
    p <- rbeta(100, 0.5, 1)
    effects <- rnorm(100)
    g <- sample(50, 100, replace=TRUE)
    out <- countGroupedDirection(p, grouping=g, effects)

    # Comparing to the parallel method.
    by.g <- split(seq_along(g), g)
    outu <- outd <- integer(length(by.g))
    names(outu) <- names(outd) <- names(by.g)

    for (x in names(by.g)) {
        current <- by.g[[x]]
        ref <- countParallelDirection(as.list(p[current]), effects=as.list(effects[current]))
        outu[x] <- ref$up
        outd[x] <- ref$down
    }

    expect_identical(outu, out$up)
    expect_identical(outd, out$down)

    # Works for log-values.
    lout <- countGroupedDirection(log(p), grouping=g, effects, log=TRUE)
    expect_identical(lout, out)

    # Works for NA values.
    out <- countGroupedDirection(c(NA, NA, NA), grouping=c(1,1,1), c(1, 0, -1))
    expect_identical(out$up[[1]], 0L)
    expect_identical(out$down[[1]], 0L)

    # Handles error states.
    expect_error(countGroupedDirection(p, rle(g[1:10]), effects), "runs.*not the same")
    expect_error(countGroupedDirection(p[1:10], rle(g), effects[1:10]), "runs.*not the same")
    expect_error(countGroupedDirection(p[1:10], g[1:10], effects), "not the same")
    expect_error(countGroupedDirection(p, g, effects[1:10]), "not the same")
})

