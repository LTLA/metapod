#ifndef P_HOLM_MIN_H
#define P_HOLM_MIN_H

#include "utils.h"
#include "Rcpp.h"
#include <algorithm>
#include <cmath>

class p_holm_min {
public:    
    p_holm_min (size_t mn, double mp) : min_num(std::max(size_t(1), mn)), min_prop(mp) {}

    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        /* Computing the (weighted) Holm correction. Weights are implemented as
         * scaling factors on the nominal type I error threshold, as described
         * in Holm's original paper.
         */
        double total_weight = 0;

        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt) {
            double& curp = pIt->first;
            const double curw = weights[pIt->second];
            total_weight += curw;
            curp = divide(curp, curw, log);
        }

        size_t index = compute_index(pvalues.size(), min_num, min_prop);
        std::partial_sort(pvalues.begin(), pvalues.begin() + index + 1, pvalues.end());
        double remaining = total_weight;
        double cummax = R_NegInf;

        for (size_t p = 0; p <= index; ++p) { // '<=' is possible, as pvalues is guaranteed to be non-empty.
            double tmp = pvalues[p].first;
            tmp = multiply(tmp, remaining, log);
            tmp = bound_upper(tmp, log);

            if (tmp > cummax) {
                cummax = tmp;
            }

            auto chosen = pvalues[p].second;
            influencers.push_back(chosen);
            remaining -= weights[chosen];
        }

        return std::make_pair(cummax, pvalues[index].second);
    }
private:
    const size_t min_num;
    const double min_prop;
};

#endif
