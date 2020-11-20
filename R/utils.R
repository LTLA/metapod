# Note that the repeated definitions are very intentional here; it ensures that
# errors are presented with the correct argument name.

.valid_logp <- function(log.p) {
    stopifnot(isTRUE(log.p) || isFALSE(log.p))
}

.valid_parallel_pvalues <- function(p.values) {
    stopifnot(is.list(p.values), all(vapply(p.values, is.numeric, TRUE)))
}

.valid_count_thresholds <- function(p.threshold, effect.threshold) {
    stopifnot(length(p.threshold)==1L, is.numeric(p.threshold), is.finite(p.threshold))
    stopifnot(length(effect.threshold)==1L, is.numeric(effect.threshold), is.finite(effect.threshold))
}

.valid_summary_threshold <- function(threshold) {
    stopifnot(length(threshold)==1L, is.numeric(threshold), is.finite(threshold))
}

.valid_min_vals <- function(min.n, min.prop) {
    stopifnot(length(min.n)==1L, is.numeric(min.n), is.finite(min.n))
    stopifnot(length(min.prop)==1L, is.numeric(min.prop), is.finite(min.prop))
}

.valid_parallel_effects <- function(effects) {
    stopifnot(is.list(effects), all(vapply(effects, is.numeric, TRUE)))
}

.valid_parallel_influential <- function(influential) {
    stopifnot(is.list(influential), all(vapply(influential, is.logical, TRUE)))
}

.valid_parallel_weights <- function(weights) {
    if (!is.null(weights)) {
        if (is.list(weights)) {
            stopifnot(all(vapply(weights, is.numeric, TRUE)))
        } else {
            stopifnot(is.numeric(weights))
        }
    }
}
