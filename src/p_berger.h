#ifndef P_BERGER_H
#define P_BERGER_H

#include "utils.h"
#include "Rcpp.h"

class p_berger {
public:    
    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        size_t best = 0;
        double largest = R_NegInf;

        for (size_t p = 0; p < pvalues.size(); ++p) {
            const double curp = pvalues[p].first;
            auto chosen = pvalues[p].second;
            if (curp > largest) {
                largest = curp;
                best = chosen;
            }
            influencers.push_back(chosen); // every non-NA p-value contributes.
        }

        return std::make_pair(largest, best);
    }
};

#endif
