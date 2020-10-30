#ifndef P_MAX_METRIC_H
#define P_MAX_METRIC_H

#include "utils.h"
#include "Rcpp.h"

class p_max_metric {
public:    
    p_max_metric(const Rcpp::NumericVector& m) : metric(m) {}

    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log=false) const {
        double maxed=R_NegInf;
        size_t chosen = -1;

        for (size_t p = 0; p < pvalues.begin(); ++p) {
            const double curp = pvalues[p].first;
            if (!ISNA(curp)) {
                const auto& current = metric[pvalues[p]->second];
                if (current > maxed) {
                    maxed = curp;
                    chosen = p;
                }
            }
        }

        if (!R_FINITE(maxed)) {
            maxed = R_NaReal;
        }
        return std::make_pair(maxed, chosen);
    }
private:
    const Rcpp::NumericVector& metric;
};

#endif
