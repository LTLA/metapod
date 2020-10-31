#ifndef P_PEARSON_H
#define P_PEARSON_H

#include "Rcpp.h"
#include "utils.h"
#include <cmath>

class p_pearson {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        double collated = 0;

        // The representative is just chosen as the test with the largest p-value.
        size_t best = 0;
        double largest = R_NegInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            const double curp = pvalues[p].first;
            if (log) {
                collated += std::log(-std::expm1(curp));
            } else {
                collated += std::log1p(-curp);
            }

            auto chosen = pvalues[p].second;
            influencers.push_back(chosen); // everyone contributes to the final p-value.

            if (curp > largest) {
                largest = curp;
                best = chosen;
            }
        }

        const double outp = R::pchisq(-2* collated, 2 * pvalues.size(), 1, log);
        return std::make_pair(outp, best);
    }
};

#endif
