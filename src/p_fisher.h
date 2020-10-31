#ifndef P_FISHER_H
#define P_FISHER_H

#include "Rcpp.h"
#include "utils.h"

class p_fisher {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        double collated = 0;
        double counter = 0;

        // The representative is just chosen as the test with the lowest p-value.
        size_t best = -1;
        double lowest = R_PosInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            const double curp = pvalues[p].first;
            if (!ISNA(curp)) {
                if (log) {
                    collated += curp;
                } else {
                    collated += std::log(curp);
                }
                ++counter;

                auto chosen = pvalues[p].second;
                influencers.push_back(chosen); // everyone contributes to the final p-value.

                if (curp < lowest) {
                    lowest = curp;
                    best = chosen;
                }
            }
        }

        const double outp = (counter ? R::pchisq(collated, 2 * counter, 0, log) : R_NaReal);
        return std::make_pair(outp, best);
    }
};

#endif
