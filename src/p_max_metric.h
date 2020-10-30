#ifndef P_MAX_METRIC_H
#define P_MAX_METRIC_H

#include "utils.h"

class p_max_metric {
public:    
    p_max_metric(const Rcpp::NumericVector& m) : metric(m) {}

    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        double maxed=R_NegInf;
        auto chosen=pvalues.begin();

        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt) {
            const auto& current=metric[pIt->second];
            if (current > maxed) {
                maxed=current;
                chosen=pIt;
            }
        }

        // Wiping out everything that is NOT the top test with respect to some
        // independent filter statistic.
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt) {
            if (pIt!=chosen) { pIt->first=1; }
        }

        return std::make_pair(chosen->first, (chosen - pvalues.begin()));
    }
private:
    const Rcpp::NumericVector& metric;
};

#endif
