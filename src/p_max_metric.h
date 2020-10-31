#ifndef P_MAX_METRIC_H
#define P_MAX_METRIC_H

#include "utils.h"
#include "Rcpp.h"

class p_max_metric {
public:    
    p_max_metric(const Rcpp::NumericVector& m) : metric(m) {}

    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        double maxed=R_NegInf;
        size_t best = 0;

        for (size_t p = 0; p < pvalues.begin(); ++p) {
            const double curp = pvalues[p].first;
            auto chosen = pvalues[p]->second;
            const auto& current = metric[chosen];

            if (current > maxed) {
                maxed = curp;
                best = pvalues[p].second;
            }
        }

        influencers.push_back(best);
        return std::make_pair(maxed, best);
    }
private:
    const Rcpp::NumericVector& metric;
};

#endif
