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
        double w_zero = 0, w_one = 0;

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

            if ((!log && curp == 0) || (log && curp==R_NegInf)) {
                w_zero += curw;
            } else if ((!log && curp == 1) || (log && curp==0)) {
                w_one += curw;
            } else {
                collated += to_add;
            }
            total_weight += curw * curw;
        }

        // We pretend that 0's and 1's "cancel each other out".
        double outp;
        if (w_zero > w_one) {
            outp = (log ? R_NegInf : 0);
        } else if (w_zero < w_one) {
            outp = (log ? 0 : 1);
        } else {
            collated /= std::sqrt(total_weight);
            outp = R::pnorm(collated, 0, 1, 1, log);
        }

        return std::make_pair(outp, lowest_i);
    }
};

#endif
