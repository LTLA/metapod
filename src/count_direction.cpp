#include "Rcpp.h"
#include "utils.h"

void correct_p(IndexedPValues& pvals, int method, bool log) {
    std::sort(pvals.begin(), pvals.end());
    const double total = pvals.size();

    if (method==0) {
        // BH method.
        for (size_t p = 0; p < pvals.size(); ++p) {
            pvals[p].first = multiply(pvals[p].first, total/(p + 1), log);
        }

        double cummin = R_PosInf;
        for (auto pIt = pvals.rbegin(); pIt != pvals.rend(); ++pIt) {
            if (pIt->first > cummin) {
                pIt->first = cummin;
            } else {
                cummin = pIt->first;
            }
        }
    } else {
        // Holm method.
        double cummax = R_NegInf;
        for (size_t p = 0; p < pvals.size(); ++p) {
            pvals[p].first = multiply(pvals[p].first, total - p, log);
            if (pvals[p].first < cummax) {
                pvals[p].first = cummax;
            } else {
                cummax = pvals[p].first;
            }
        }
    }

    return;
}

// [[Rcpp::export(rng=false)]]
Rcpp::List count_parallel_direction(Rcpp::List pvalues, Rcpp::List effects, int method, double pthreshold, double ethreshold, bool log) {
    parallel_vectors<Rcpp::NumericVector> Effects(effects);
    parallel_vectors<Rcpp::NumericVector> PValues(pvalues);
    if (Effects.nvectors != PValues.nvectors || Effects.nelements != PValues.nelements) {
        throw std::runtime_error("'pvalues' should have the same structure as 'effects'");
    }

    IndexedPValues pvals;
    Rcpp::IntegerVector outu(Effects.nelements), outd(Effects.nelements);
    if (log) {
        pthreshold = std::log(pthreshold);
    }

    for (size_t i = 0; i < PValues.nelements; ++i) {
        pvals.clear();
        for (size_t j = 0; j < PValues.nvectors; ++j) {
            if (!ISNAN(PValues.vectors[j][i])) {
                pvals.push_back(std::make_pair(PValues.vectors[j][i], j));
            }
        }
        
        correct_p(pvals, method, log);

        for (auto& p : pvals) {
            if (p.first <= pthreshold) {
                double effector = Effects.vectors[p.second][i];
                if (effector < ethreshold) {
                    ++outd[i];
                } else if (effector > ethreshold) {
                    ++outu[i];
                }
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("down")=outd,
        Rcpp::Named("up")=outu
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List count_grouped_direction(Rcpp::NumericVector pvalues, Rcpp::IntegerVector runs, Rcpp::NumericVector effects, 
    int method, double pthreshold, double ethreshold, bool log) 
{
    if (effects.size() != pvalues.size()) {
        throw std::runtime_error("'effects' and 'pvalues' should have the same length");
    }

    IndexedPValues pvals;
    Rcpp::IntegerVector outu(runs.size()), outd(runs.size());
    if (log) {
        pthreshold = std::log(pthreshold);
    }
    size_t counter = 0;

    for (size_t g = 0; g < runs.size(); ++g) {
        pvals.clear();
        for (int r = 0; r < runs[g]; ++r, ++counter) {
            if (counter >= pvalues.size()) {
                throw std::runtime_error("'sum(runs)' is not the same as 'length(pvalues)'");
            }
            if (!ISNAN(pvalues[counter])) {
                pvals.push_back(std::make_pair(pvalues[counter], counter));
            }
        }

        correct_p(pvals, method, log);

        for (auto& p : pvals) {
            if (p.first <= pthreshold) {
                double effector = effects[p.second];
                if (effector < ethreshold) {
                    ++outd[g];
                } else if (effector > ethreshold) {
                    ++outu[g];
                }
            }
        }
    }

    if (counter != pvalues.size()) {
        throw std::runtime_error("'sum(runs)' is not the same as 'length(pvalues)'");
    }

    return Rcpp::List::create(
        Rcpp::Named("down")=outd,
        Rcpp::Named("up")=outu
    );
}

