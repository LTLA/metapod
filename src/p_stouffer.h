#ifndef P_STOUFFER_H
#define P_STOUFFER_H

#include "utils.h"
#include "Rcpp.h"
#include <cmath>

class p_stouffer {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        auto wIt = weights.begin();
        double collated = 0, total_weight = 0;

        // The test with the largest weighted z-score is chosen as the representative.
        size_t lowest_i = 0;
        double lowest_v = 0;

        for (auto pIt = pvalues.begin(); pIt != pvalues.end(); ++pIt, ++wIt) {
            const double to_add = R::qnorm(pIt->first, 0, 1, 1, 0) * (*wIt);
            if (std::abs(to_add) > lowest_v) {
                lowest_v = std::abs(to_add);
                lowest_i = pIt - pvalues.begin();
            }

            collated += to_add;
            total_weight += *wIt;
        }

        collated /= std::sqrt(total_weight);
        const double outp = R::pnorm(collated, 0, 1, 1, 0);
        return std::make_pair(outp, lowest_i);
    }
};

#endif
