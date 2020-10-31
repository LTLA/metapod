#ifndef P_BERGER_H
#define P_BERGER_H

#include "utils.h"
#include "Rcpp.h"

class p_berger {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        size_t best = -1;
        double lowest = R_PosInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            const double curp = pvalues[p].first;
            if (!ISNA(curp)) {
                auto chosen = pvalues[p].second;
                if (curp < lowest) {
                    lowest = curp;
                    best = chosen;
                }
                influencers.push_back(chosen); // every non-NA p-value contributes.
            }
        }

        if (!R_FINITE(best)) {
            lowest = R_NaReal;
        }
        return std::make_pair(lowest, best);
    }
};

#endif
