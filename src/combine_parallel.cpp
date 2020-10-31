#include "p_berger.h"
#include "p_holm_min.h"
#include "p_lowest_p.h"
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

    size_t nvectors = 0;
    size_t nelements = 0;
    std::vector<Rcpp::NumericVector> vectors;
};


class parallel_weight_server {
public:
    parallel_weight_server(size_t nv, size_t ne, Rcpp::RObject weights) : nvectors(nv), nelements(ne) {
        if (!weights.isNULL()) {
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
    }

    template <class ITER>
    void prefill(ITER start) {
        // Under these modes, weights are constants and just get filled once.
        if (wmode == 0) {
            std::fill(start, start + nvectors, 1);
        } else if (wmode == 2) {
            std::copy(wrepeat.begin(), wrepeat.end(), start);
        } 
    }

    template <class ITER>
    void fill(size_t i, ITER start) {
        if (wmode == 1) { 
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
Rcpp::List compute_parallel (Rcpp::List pvals, Rcpp::RObject weights, bool log, const PREP& pcompute) {
    parallel_vectors pvectors(pvals);
    const size_t np = pvectors.nvectors;
    const size_t nlen = pvectors.nelements;

    parallel_weight_server wserver(np, nlen, weights);
    std::vector<double> tmpweights(np);
    wserver.prefill(tmpweights.begin());
    IndexedPValues pvec(np);

    Rcpp::NumericVector outp(nlen);
    Rcpp::IntegerVector outrep(nlen);
    std::vector<Rcpp::LogicalVector> outinf(np);
    for (size_t j = 0; j < np; ++j) {
        outinf[j] = Rcpp::LogicalVector(nlen);
    }

    for (size_t i = 0; i < nlen; ++i) {
        for (size_t j = 0; j < np; ++j) {
            pvec[j].first = pvectors.vectors[j][i];
            pvec[j].second = j;
        }

        wserver.fill(i, tmpweights.begin());

        auto chosen = pcompute(pvec, tmpweights, log);
        outp[i] = chosen.first;
        if (chosen.second==-1) {
            outrep[i] = 0;
        } else {
            size_t rep = pvec[chosen.second].second;
            outrep[i] = rep + 1;
            for (size_t k = 0; k <= rep; ++k) {
                outinf[pvec[k].second][i] = 1;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("p.value")=outp, 
        Rcpp::Named("representative")=outrep,
        Rcpp::Named("influential")=Rcpp::List(outinf.begin(), outinf.end())
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_simes(Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_simes());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_berger(Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_berger());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_holm_min (Rcpp::List pvals, Rcpp::RObject weights, bool log, int min_n, double min_prop) {
    return compute_parallel(pvals, weights, log, p_holm_min(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_lowest_p (Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_lowest_p());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_stouffer (Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_stouffer());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_fisher (Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_fisher());
}
