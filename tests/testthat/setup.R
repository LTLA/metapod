parallelTester <- function(p1, p2, p3, FUN) {
    all.p <- list(p1, p2, p3)

    ref <- FUN(all.p)
    infmat <- do.call(cbind, ref$influential)
    out <- infmat[cbind(seq_len(nrow(infmat)), ref$representative)]
    expect_true(all(out))

    # Behaves properly upon logging.
    out <- FUN(lapply(all.p, log), log.p=TRUE)
    expect_equal(log(ref$p.value), out$p.value)
    expect_equal(ref$reference, out$reference)
    expect_equal(ref$influential, out$influential)

    # Handles solo inputs.
    for (i in seq_along(all.p)) {
        current <- FUN(all.p[i])
        expect_equal(all.p[[i]], current$p.value)
        expect_true(all(current$representative==1))
        expect_true(all(current$influential[[1]]))
    }

    # Handles empty inputs.
    empty <- FUN(lapply(all.p, "[", i=0))
    expect_equal(empty$p.value, numeric(0))

    # Throws on invalid inputs.
    expect_error(FUN(p1, p2[0]), "length.*should be equal")

    # Handles partial NA values correctly.
    some.na <- sample(length(p1), length(p1)/2)
    p1.na <- p1
    p1.na[some.na] <- NA
    p2.na <- p2
    p2.na[-some.na] <- NA

    ref <- FUN(list(p1, p3))
    ref2 <- FUN(list(p2, p3))
    direct <- FUN(list(p1.na, p2.na, p3))

    refp <- ref$p.value
    refp[some.na] <- ref2$p.value[some.na]
    expect_equal(refp, direct$p.value)

    refr <- c(1L, 3L)[ref$representative]
    refr[some.na] <- ref2$representative[some.na] + 1L
    expect_equal(refr, direct$representative)

    # Handles all-NA values correctly.
    p1.na <- p1
    p2.na <- p2
    p3.na <- p3
    p1.na[1] <- p2.na[1] <- p3.na[1] <- NA_real_
    out <- FUN(list(p1.na, p2.na, p3.na))

    expect_equal(out$p.value[1], NA_real_)
    expect_equal(out$representative[1], 0L)
    expect_false(any(vapply(out$influential, "[", i=1, FALSE)))

    expect_false(any(is.na(out$p.value[-1])))
    expect_false(any(out$representative[-1]==0))
    expect_true(all(Reduce("|", lapply(out$influential, "[", i=-1))))

    TRUE
}

parallelTesterWithWeights <- function(p1, p2, p3, FUN) {
    all.p <- list(p1, p2, p3)
    ref <- FUN(all.p)

    # Checking that weights are actually respected.
    weights <- runif(3)
    out <- FUN(all.p, weights=weights)
    expect_false(isTRUE(all.equal(out,ref)))

    # Checking that list-like weights are handled properly.
    weights2 <- runif(3)
    out <- FUN(lapply(all.p, "[", i=1:10), weights=weights)
    out2 <- FUN(lapply(all.p, "[", i=11:15), weights=weights2)

    lweights <- mapply(c,
        lapply(weights, rep, each=10),
        lapply(weights2, rep, each=5),
        SIMPLIFY=FALSE
    )
    combined <- FUN(lapply(all.p, "[", i=1:15), weights=lweights)

    expect_equal(combined$p.value, c(out$p.value, out2$p.value))
    expect_equal(combined$representative, c(out$representative, out2$representative))

    # Weights and NA's interact correctly.
    weights <- runif(3)
    p1.na <- p1
    p1.na[] <- NA
    out <- FUN(list(p1.na, p2, p3), weights=weights)
    ref <- FUN(list(p2, p3), weights=weights[-1])

    expect_identical(out$p.value, ref$p.value)
    expect_identical(out$representative, ref$representative + 1L)
    expect_false(any(out$influential[[1]]))
    expect_identical(out$influential[-1], ref$influential)
}
