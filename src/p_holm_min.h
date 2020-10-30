#ifndef P_HOLM_MIN_H
#define P_HOLM_MIN_H

#include "utils.h"
#include <algorithm>
#include <cmath>

class p_holm_min {
public:    
    p_holm_min (size_t mn, double mp) : min_num(std::max(size_t(1), mn)), min_prop(mp) {}

    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        /* Computing the (weighted) Holm correction. Weights are implemented as
         * scaling factors on the nominal type I error threshold, as described
         * in Holm's original paper.
         */
        double total_weight = 0;
        auto wIt = weights.begin();
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt, ++wIt) {
            total_weight += *wIt;
            (pIt->first) /= *wIt;
        }

        std::sort(pvalues.begin(), pvalues.end());

        double remaining=total_weight;
        double cummax=0;
        wIt = weights.begin();
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt, ++wIt) {
            double& curp=(pIt->first);
            curp*=remaining;
            remaining -= *wIt;

            if (curp > 1) {
                curp=1;
            }
            if (curp > cummax) {
                cummax=curp;
            } else {
                curp=cummax;
            }
        }

        size_t index=std::max(
            min_num, 
            static_cast<size_t>(std::ceil(min_prop * static_cast<double>(pvalues.size())))
        );
        index=std::min(index, pvalues.size());
        if (index!=0) {
            --index;
        }
        return std::make_pair(pvalues[index].first, index);
    }
private:
    const size_t min_num;
    const double min_prop;
};

#endif
