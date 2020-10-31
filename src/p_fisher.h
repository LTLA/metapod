#ifndef P_FISHER_H
#define P_FISHER_H

#include "Rcpp.h"
#include "utils.h"

class p_fisher {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        double collated = 0;

        // The representative is just chosen as the test with the lowest p-value.
        size_t best = 0;
        double lowest = R_PosInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            const double curp = pvalues[p].first;
            if (log) {
                collated += curp;
            } else {
                collated += std::log(curp);
            }

            auto chosen = pvalues[p].second;
            influencers.push_back(chosen); // everyone contributes to the final p-value.

            if (curp < lowest) {
                lowest = curp;
                best = chosen;
            }
        }

        const double outp = R::pchisq(-2* collated, 2 * pvalues.size(), 0, log);
        return std::make_pair(outp, best);
    }
};

#endif
