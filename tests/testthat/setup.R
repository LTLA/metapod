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
    expect_error(FUN(list(p1, p2[0])), "same length")

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

    expect_true(is.na(out$p.value[1]))
    expect_true(is.na(out$representative[1]))
    expect_false(any(vapply(out$influential, "[", i=1, FALSE)))

    expect_false(any(is.na(out$p.value[-1])))
    expect_false(any(is.na(out$representative[-1])))
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
    expect_error(FUN(all.p, weights=lweights), "length.*should be equal")

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

    TRUE
}

groupedTester <- function(p, g, gFUN, pFUN) {
    ref <- gFUN(p, g)

    # Logging works.
    lout <- gFUN(log(p), g, log.p=TRUE)
    lout$p.value <- exp(lout$p.value)
    expect_equal(ref, lout)

    # RLE mode works.
    o <- order(g)
    out <- gFUN(p[o], rle(g[o]))
    out$representative[] <- o[out$representative]
    out$influential[o] <- out$influential
    expect_identical(out, ref)

    # Manual looping to compare to parallel function.
    by.g <- split(seq_along(p), g)
    outp <- numeric(length(by.g))
    outrep <- integer(length(by.g))
    outinf <- logical(length(g))

    for (i in seq_along(by.g)) {
        current <- by.g[[i]]
        single <- pFUN(as.list(p[current]))
        outp[i] <- single$p.value
        outrep[i] <- current[single$representative]
        outinf[current] <- unlist(single$influential)
    }

    names(outp) <- names(outrep) <- names(by.g)
    expect_equal(outp, ref$p.value)
    expect_equal(outrep, ref$representative)
    expect_equal(outinf, ref$influential)

    # Robust to NA's.
    thrown <- sample(length(p), length(p)/2)
    p2 <- p
    p2[-thrown] <- NA

    has.na <- gFUN(p2, g)
    nullified <- gFUN(p[thrown], g[thrown])

    commong <- as.character(sort(unique(g[thrown])))
    expect_identical(commong, names(nullified$p.value)) 
    expect_equal(has.na$p.value[commong], nullified$p.value)
    expect_equivalent(has.na$representative[commong], thrown[nullified$representative])
    expect_identical(has.na$influential[thrown], nullified$influential)

    lost <- setdiff(names(has.na$p.value), commong)
    expect_true(all(is.na(has.na$p.value[lost])))
    expect_true(all(is.na(has.na$representative[lost])))
    expect_false(any(has.na$influential[-thrown]))

    # Robust to extreme NA's.
    p2 <- p
    p2[g == g[1]] <- NA
    has.na <- gFUN(p2, g)

    grp <- as.character(g[1])
    expect_true(is.na(has.na$p.value[grp]))
    expect_true(is.na(has.na$representative[grp]))
    expect_false(any(has.na$influential[g == g[1]]))

    TRUE
}

groupedTesterWithWeights <- function(p, g, gFUN, pFUN) {
    w <- rexp(length(p))
    ref <- gFUN(p, g, weights=w)
    out <- gFUN(p, g)
    expect_false(isTRUE(all.equal(out,ref))) # weights actually have an effect.
    
    # Manual looping to compare to parallel function.
    by.g <- split(seq_along(p), g)
    outp <- numeric(length(by.g))
    outrep <- integer(length(by.g))
    outinf <- logical(length(g))

    for (i in seq_along(by.g)) {
        current <- by.g[[i]]
        single <- pFUN(as.list(p[current]), weights=as.list(w[current]))
        outp[i] <- single$p.value
        outrep[i] <- current[single$representative]
        outinf[current] <- unlist(single$influential)
    }

    names(outp) <- names(outrep) <- names(by.g)
    expect_equal(outp, ref$p.value)
    expect_equal(outrep, ref$representative)
    expect_equal(outinf, ref$influential)

    # Weights and NA's interact correctly.
    thrown <- sample(length(p), length(p)/2)
    p2 <- p
    p2[-thrown] <- NA

    has.na <- gFUN(p2, g, weights=w)
    nullified <- gFUN(p[thrown], g[thrown], weights=w[thrown])

    commong <- as.character(sort(unique(g[thrown])))
    expect_identical(commong, names(nullified$p.value)) 
    expect_equal(has.na$p.value[commong], nullified$p.value)
    expect_equivalent(has.na$representative[commong], thrown[nullified$representative])
    expect_identical(has.na$influential[thrown], nullified$influential)

    TRUE
}
