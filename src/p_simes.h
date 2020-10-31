#ifndef P_SIMES_H
#define P_SIMES_H

#include "Rcpp.h"
#include "utils.h"
#include <algorithm>

class p_simes {
public:
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        /* Computing the (weighted) FDR and thus the (weighted) Simes value.
         * The weights are implemented as frequency weights, e.g., if you had 2
         * tests with a weight of 10 to 1, you'd consider the one with the
         * higher weight 10 more times to try to reject the global null (i.e.,
         * expanding it in-place in the sorted vector of p-values).
         */
        std::sort(pvalues.begin(), pvalues.end());

        double cumweight = 0;
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt) {
            double& curp = pIt->first;
            cumweight += weights[pIt->second];
            curp = divide(curp, cumweight, log);
        }

        // Backtracking to create adjusted p-values with a cumulative minimum.
        double curmin = R_PosInf;
        size_t counter = pvalues.size() - 1, minindex = 0;

        for (auto prIt=pvalues.rbegin(); prIt!=pvalues.rend(); ++prIt, --counter) {
            double current = (prIt->first);
            if (current < curmin) {
                curmin=current;
                minindex = counter;
            }
        }

        curmin = multiply(curmin, cumweight, log); // cumweight is the total weight now.
        for (size_t i = 0; i <= minindex; ++i) {
            influencers.push_back(pvalues[i].second);
        }
        return std::make_pair(curmin, pvalues[minindex].second);
    }
};

#endif
