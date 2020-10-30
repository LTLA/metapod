#include "p_berger.h"
#include "p_holm_min.h"
#include "p_lowest_p.h"
#include "p_max_metric.h"
#include "p_simes.h"
#include "p_stouffer.h"
#include "p_fisher.h"

#include "Rcpp.h"
#include <vector>
#include <algorithm>

struct parallel_vectors {
    parallel_vectors() {}

    parallel_vectors(Rcpp::List input) {
        nvectors = input.size();
        vectors.resize(nvectors);
        for (size_t i = 0; i < nvectors; ++i) {
            vectors[i] = Rcpp::NumericVector(input[i]);
        }

        if (nvectors) {
            nelements = vectors.front().size();
            for (size_t p = 1; p < nvectors; ++p) {
                if (nelements != vectors[p].size()) {
                    throw std::runtime_error("p-value vectors should have the same length");
                }
            }
        }
    }

    size_t nvectors;
    size_t nelements = 0;
    std::vector<Rcpp::NumericVector> vectors;
};


class parallel_weight_server {
public:
    parallel_weight_server(size_t nv, size_t ne, Rcpp::RObject weights) : nvectors(nv), nelements(ne), wvecs(nv) {
        if (weights.sexp_type()==VECSXP) {
            wmode = 1;
            wvecs = parallel_vectors(Rcpp::List(weights));
            if (wvecs.nvectors != nvectors || wvecs.nelements != nelements) {
                throw std::runtime_error("lengths of list 'weights' should be equal to lengths of p-value vectors");
            }

        } else if (weights.sexp_type()==REALSXP) {
            wmode = 2;
            wrepeat = Rcpp::NumericVector(weights);
            if (wrepeat.size() != nvectors) {
                throw std::runtime_error("length of vector 'weights' should be equal to number of p-value vectors");
            }
        }
    }

    template <class ITER>
    void prefill(ITER start) {
        // Under these modes, weights are constants and just get filled once.
        if (wmode == 0) {
            std::fill(start, start + nvectors, 1);
        } else if (wmode == 1) {
            std::copy(wrepeat.begin(), wrepeat.end(), start);
        } 
    }

    template <class ITER>
    void fill(size_t i, ITER start) {
        if (wmode == 2) { 
            for (size_t j = 0; j < nvectors; ++j, ++start) {
                *start = wvecs.vectors[j][i]; 
            }
        }
    }

private:
    const size_t nvectors, nelements;
    int wmode = 0;
    Rcpp::NumericVector wrepeat;
    parallel_vectors wvecs;
};


template<class PREP>
Rcpp::List compute_parallel (Rcpp::List pvals, Rcpp::RObject weights, const PREP& pcompute) {
    parallel_vectors pvectors(pvals);
    const size_t np = pvectors.nvectors;
    const size_t nlen = pvectors.nelements;

    parallel_weight_server wserver(np, nlen, weights);
    std::vector<double> tmpweights(np);
    wserver.prefill(tmpweights.begin());
    IndexedPValues pvec(np);

    Rcpp::NumericVector outp(nlen);
    Rcpp::IntegerVector outrep(nlen);

    for (size_t i = 0; i < nlen; ++i) {
        for (size_t j = 0; j < np; ++j) {
            pvec[j].first = pvectors.vectors[j][i];
            pvec[j].second = j;
        }

        wserver.fill(i, tmpweights.begin());

        auto chosen = pcompute(pvec, tmpweights);
        outp[i] = chosen.first;
        outrep[i] = pvec[chosen.second].second;
    }

    return Rcpp::List::create(outp, outrep);
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_simes (Rcpp::List pvals, Rcpp::RObject weights) {
    return compute_parallel(pvals, weights, p_simes());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_berger (Rcpp::List pvals, Rcpp::RObject weights) {
    return compute_parallel(pvals, weights, p_berger());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_holm_min (Rcpp::List pvals, Rcpp::RObject weights, int min_n, double min_prop) {
    return compute_parallel(pvals, weights, p_holm_min(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_lowest_p (Rcpp::List pvals, Rcpp::RObject weights) {
    return compute_parallel(pvals, weights, p_lowest_p());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_max_metric (Rcpp::List pvals, Rcpp::RObject weights, Rcpp::NumericVector metric) {
    return compute_parallel(pvals, weights, p_max_metric(metric));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_stouffer (Rcpp::List pvals, Rcpp::RObject weights) {
    return compute_parallel(pvals, weights, p_stouffer());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_fisher (Rcpp::List pvals, Rcpp::RObject weights) {
    return compute_parallel(pvals, weights, p_fisher());
}
