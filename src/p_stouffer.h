#ifndef P_STOUFFER_H
#define P_STOUFFER_H

#include "utils.h"
#include "Rcpp.h"
#include <cmath>

class p_stouffer {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        double collated = 0, total_weight = 0;

        // The test with the largest weighted z-score is chosen as the representative.
        size_t lowest_i = 0;
        double lowest_v = R_PosInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            double curp = pvalues[p].first;
            auto chosen = pvalues[p].second;
            influencers.push_back(chosen);

            const double curw = weights[chosen];
            const double to_add = R::qnorm(curp, 0, 1, 1, log) * curw;

            if (to_add < lowest_v) {
                lowest_v = to_add;
                lowest_i = chosen;
            }

            collated += to_add;
            total_weight += curw;
        }

        collated /= std::sqrt(total_weight);

        // Protect against NA's due to Inf - Inf, when p-values of 0 and 1 are both present.
        double outp = (ISNAN(collated) ? 1 : R::pnorm(collated, 0, 1, 1, log));
        return std::make_pair(outp, lowest_i);
    }
};

#endif
