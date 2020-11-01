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

template<class PREP>
Rcpp::List compute_grouped (Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log, const PREP& pcompute) {
    Rcpp::NumericVector wvec;
    if (!weights.isNULL()) {
        wvec = Rcpp::NumericVector(weights);
        if (wvec.size()!=pvals.size()) {
            throw std::runtime_error("'weights' and 'pvals' must have the same length");
        }
    } else {
        wvec = Rcpp::NumericVector(pvals.size(), 1);
    }

    IndexedPValues pvec;
    std::deque<size_t> influencers;
    auto pIt = pvals.begin();
    size_t counter = 0;

    Rcpp::NumericVector outp(runs.size());
    Rcpp::IntegerVector outrep(runs.size());
    Rcpp::LogicalVector outinf(pvals.size());
    
    for (size_t g = 0; g < runs.size(); ++g) {
        pvec.clear();
        for (int r = 0; r < runs[g]; ++r, ++pIt, ++counter) {
            if (pIt == pvals.end()) {
                throw std::runtime_error("'sum(runs)' is not the same as 'length(pvals)'");
            }
            if (R_FINITE(*pIt)) {
                pvec.push_back(std::make_pair(*pIt, counter));
            }
        }

        if (pvec.empty()) {
            outp[g] = R_NaReal;
            continue;
        }

        influencers.clear();
        auto out = pcompute(pvec, wvec, log, influencers); 

        outp[g] = out.first;
        outrep[g] = out.second + 1; // 1-based indexing.
        for (auto i : influencers) { 
            outinf[i] = 1;
        }
    }

    if (counter != pvals.size()) {
        throw std::runtime_error("'sum(runs)' is not the same as 'length(pvals)'");
    }

    return Rcpp::List::create(
        Rcpp::Named("p.value")=outp, 
        Rcpp::Named("representative")=outrep,
        Rcpp::Named("influential")=outinf
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_simes(Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log) {
    return compute_grouped(pvals, runs, weights, log, p_simes());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_berger(Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log) {
    return compute_grouped(pvals, runs, weights, log, p_berger());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_wilkinson(Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log, int min_n, double min_prop) {
    return compute_grouped(pvals, runs, weights, log, p_wilkinson(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_holm_min (Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log, int min_n, double min_prop) {
    return compute_grouped(pvals, runs, weights, log, p_holm_min(min_n, min_prop));
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_stouffer (Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log) {
    return compute_grouped(pvals, runs, weights, log, p_stouffer());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_fisher (Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log) {
    return compute_grouped(pvals, runs, weights, log, p_fisher());
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_grouped_pearson(Rcpp::NumericVector pvals, Rcpp::IntegerVector runs, Rcpp::RObject weights, bool log) {
    return compute_grouped(pvals, runs, weights, log, p_pearson());
}
