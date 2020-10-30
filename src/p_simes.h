#ifndef P_SIMES_H
#define P_SIMES_H

#include "utils.h"
#include <algorithm>

class p_simes {
public:
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log=false) const {
        /* Computing the (weighted) FDR and thus the (weighted) Simes value.
         * The weights are implemented as frequency weights, e.g., if you had 2
         * tests with a weight of 10 to 1, you'd consider the one with the
         * higher weight 10 more times to try to reject the global null (i.e.,
         * expanding it in-place in the sorted vector of p-values).
         */
        std::sort(pvalues.begin(), pvalues.end());

        double cumweight = 0;
        auto wIt = weights.begin();
        size_t nonna = 0;
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt, ++wIt) {
            double& curp = pIt->first;
            if (!ISNA(curp)) {
                cumweight += *wIt;
                curp = divide(curp, cumweight, log);
                ++nonna;
            }
        }

        // Backtracking to create adjusted p-values with a cumulative minimum.
        double curmin = 1;
        size_t counter=nonna - 1, minindex=counter;

        for (auto prIt=pvalues.rbegin(); prIt!=pvalues.rend(); ++prIt, --counter) {
            double current = (prIt->first);
            if (!ISNA(current)) {
                current = multiply(current, cumweight, log); // cumweight is the total weight now.
                if (current < curmin) {
                    curmin=current;
                    minindex=counter;
                } else {
                    current=curmin;
                }
            }
        }

        return std::make_pair(curmin, minindex);
    }
};

#endif
