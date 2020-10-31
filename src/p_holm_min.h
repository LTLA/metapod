#ifndef P_HOLM_MIN_H
#define P_HOLM_MIN_H

#include "utils.h"
#include "Rcpp.h"
#include <cmath>

class p_holm_min {
public:    
    p_holm_min (size_t mn, double mp) : min_num(std::max(size_t(1), mn)), min_prop(mp) {}

    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        /* Computing the (weighted) Holm correction. Weights are implemented as
         * scaling factors on the nominal type I error threshold, as described
         * in Holm's original paper.
         */
        double total_weight = 0;
        size_t num_nonna = 0;

        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt) {
            double& curp = pIt->first;
            if (!ISNA(curp)) {
                const double curw = weights[pIt->second];
                total_weight += curw;
                curp = divide(curp, curw, log);
                ++num_nonna;
            }
        }

        std::sort(pvalues.begin(), pvalues.end());

        // Figuring out which index to extract.
        size_t index=std::max(
            min_num, 
            static_cast<size_t>(std::ceil(min_prop * static_cast<double>(num_nonna)))
        );
        index=std::min(index, num_nonna);
        if (index!=0) {
            --index;
        }

        double remaining = total_weight;
        double cummax = 0;
        std::pair<double, size_t> output(R_NaReal, -1);

        for (size_t p = 0; p < pvalues.size(); ++p) {
            double curp = pvalues[p].first;
            if (!ISNA(curp)) {
                curp = multiply(curp, remaining, log);
                curp = bound_upper(curp, log);

                if (curp > cummax) {
                    cummax = curp;
                } else {
                    curp = cummax;
                }

                auto chosen = pvalues[p].second;
                influencers.push_back(chosen);

                if (p == index) {
                    output.first = curp;
                    output.second = chosen;
                    break;
                }

                remaining -= weights[chosen];
            }
        }

        return output;
    }
private:
    const size_t min_num;
    const double min_prop;
};

#endif
