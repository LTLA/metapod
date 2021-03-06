#include "p_berger.h"
#include "p_holm_min.h"
#include "p_simes.h"
#include "p_stouffer.h"
#include "p_fisher.h"
#include "p_pearson.h"
#include "p_wilkinson.h"

#include "Rcpp.h"
#include <vector>
#include <algorithm>
#include <deque>

class parallel_weight_server {
public:
    parallel_weight_server(size_t nv, size_t ne, Rcpp::RObject weights) : nvectors(nv), nelements(ne) {
        if (!weights.isNULL()) {
            if (weights.sexp_type()==VECSXP) {
                wmode = 1;
                wvecs = parallel_vectors<Rcpp::NumericVector>(Rcpp::List(weights));
                if (wvecs.nvectors != nvectors || wvecs.nelements != nelements) {
                    throw std::runtime_error("lengths of list 'weights' should be equal to lengths of p-value vectors");
                }

            } else { 
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
            for (auto w : wrepeat) {
                if (!R_FINITE(w) || w <= 0) {
                    throw std::runtime_error("all 'weights' must be positive");
                }
            }
        } 
    }

    template <class ITER>
    void fill(size_t i, ITER start) {
        if (wmode == 1) { 
            for (size_t j = 0; j < nvectors; ++j, ++start) {
                double target = wvecs.vectors[j][i]; 
                if (!R_FINITE(target) || target <= 0) {
                    throw std::runtime_error("all 'weights' must be positive");
                }
                *start = target;
            }
        }
    }

private:
    const size_t nvectors, nelements;
    int wmode = 0;
    Rcpp::NumericVector wrepeat;
    parallel_vectors<Rcpp::NumericVector> wvecs;
};

template<class PREP>
Rcpp::List compute_parallel (Rcpp::List pvals, Rcpp::RObject weights, bool log, const PREP& pcompute) {
    parallel_vectors<Rcpp::NumericVector> pvectors(pvals);
    const size_t np = pvectors.nvectors;
    const size_t nlen = pvectors.nelements;

    parallel_weight_server wserver(np, nlen, weights);
    std::vector<double> tmpweights(np);
    wserver.prefill(tmpweights.begin());
    IndexedPValues pvec(np);
    std::deque<size_t> influencers;

    Rcpp::NumericVector outp(nlen);
    Rcpp::IntegerVector outrep(nlen);
    std::vector<Rcpp::LogicalVector> outinf(np);
    for (size_t j = 0; j < np; ++j) {
        outinf[j] = Rcpp::LogicalVector(nlen);
    }

    for (size_t i = 0; i < nlen; ++i) {
        pvec.clear();
        for (size_t j = 0; j < np; ++j) {
            double curp = pvectors.vectors[j][i];
            if (!ISNAN(curp)) {
                pvec.push_back(std::make_pair(curp, j));
            }
        }

        if (pvec.empty()) {
            outp[i] = R_NaReal;
            outrep[i] = NA_INTEGER;
            continue;
        }

        wserver.fill(i, tmpweights.begin());
        influencers.clear();

        auto chosen = pcompute(pvec, tmpweights, log, influencers);
        outp[i] = chosen.first;
        outrep[i] = chosen.second + 1;
        for (auto x : influencers) {
            outinf[x][i] = 1;
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
Rcpp::List compute_parallel_wilkinson(Rcpp::List pvals, Rcpp::RObject weights, bool log, int min_n, double min_prop) {
    return compute_parallel(pvals, weights, log, p_wilkinson(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_holm_min (Rcpp::List pvals, Rcpp::RObject weights, bool log, int min_n, double min_prop) {
    return compute_parallel(pvals, weights, log, p_holm_min(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_stouffer (Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_stouffer());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_fisher (Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_fisher());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_parallel_pearson(Rcpp::List pvals, Rcpp::RObject weights, bool log) {
    return compute_parallel(pvals, weights, log, p_pearson());
}
