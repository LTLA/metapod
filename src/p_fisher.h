#ifndef P_FISHER_H
#define P_FISHER_H

#include "Rcpp.h"
#include "utils.h"
#include <algorithm>

class p_fisher {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        double collated = 0;
        for (auto pIt = pvalues.begin(); pIt != pvalues.end(); ++pIt) {
            collated += std::log(pIt->first);
        }
        const double outp = R::pchisq(collated, 2 * pvalues.size(), 0, 0);

        // The representative is just chosen as the test with the lowest p-value.
        auto rep = std::min_element(pvalues.begin(), pvalues.end());
        return std::make_pair(outp, rep - pvalues.begin());
    }
};

#endif
