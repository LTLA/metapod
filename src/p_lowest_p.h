#ifndef P_LOWEST_P_H
#define P_LOWEST_P_H

#include "utils.h"
#include "Rcpp.h"

class p_lowest_p {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        auto lowest = pvalues.begin();
        double curlowest = R_PosInf;

		/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
		 * Weights are defined according to the weighted Bonferroni (see http://arxiv.org/abs/math.ST/0604172,
		 * though some mental arithmetic is needed). These can also be treated as relative frequency weights,
		 * i.e. the total number of tests is rescaled relative to the weight of the current test (so, [10,1] 
		 * weights would consider there to be 1.1 tests for the first one and 11 tests for the second one).
		 */
        double total_weight = 0;
        auto wIt = weights.begin();
        for (auto pIt=pvalues.begin(); pIt!=pvalues.end(); ++pIt, ++wIt) {
            const double curweight = *wIt;
            total_weight += curweight;

            double curp = pIt->first / curweight;
            if (curlowest > curp) {
                curlowest = curp;
                lowest = pIt;
            }
        }

        curlowest *= total_weight;
        if (curlowest > 1) {
            curlowest = 1;
        }
        return std::make_pair(curlowest, lowest - pvalues.begin());
    }
};

#endif
