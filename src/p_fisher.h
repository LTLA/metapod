#ifndef P_FISHER_H
#define P_FISHER_H

#include "Rcpp.h"
#include "utils.h"

class p_fisher {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log=false) const {
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

                if (curp < lowest) {
                    lowest = curp;
                    best = p;
                }
            }
        }

        const double outp = (counter ? R::pchisq(collated, 2 * counter, 0, log) : R_NaReal);
        return std::make_pair(outp, best);
    }
};

#endif
