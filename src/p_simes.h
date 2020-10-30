#ifndef P_SIMES_H
#define P_SIMES_H

#include "utils.h"
#include <algorithm>

class p_simes {
public:
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        /* Computing the (weighted) FDR and thus the (weighted) Simes value.
         * The weights are implemented as frequency weights, e.g., if you had 2
         * tests with a weight of 10 to 1, you'd consider the one with the
         * higher weight 10 more times to try to reject the global null (i.e.,
         * expanding it in-place in the sorted vector of p-values).
         */
        std::sort(pvalues.begin(), pvalues.end());

        double cumweight = 0;
        auto wIt = weights.begin();
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt, ++wIt) {
            cumweight += *wIt;
            (pIt->first) /= cumweight;
        }

        // Backtracking to create adjusted p-values with a cumulative minimum.
        double curmin=1;
        size_t counter=pvalues.size()-1, minindex=counter;
        for (auto prIt=pvalues.rbegin(); prIt!=pvalues.rend(); ++prIt, --counter) {
            double& current=(prIt->first);
            current *= cumweight; // this is the total weight now.
            if (current < curmin) {
                curmin=current;
                minindex=counter;
            } else {
                current=curmin;
            }
        }

        return std::make_pair(curmin, minindex);
    }
};

#endif
