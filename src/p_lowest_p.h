#ifndef P_LOWEST_P_H
#define P_LOWEST_P_H

#include "utils.h"
#include "Rcpp.h"

class p_lowest_p {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log=false) const {
        size_t lowest = -1;
        double curlowest = R_PosInf;

		/* Computing the Holm p-value for the best window (basically Bonferroni, if we're taking the minimum).
		 * Weights are defined according to the weighted Bonferroni (see http://arxiv.org/abs/math.ST/0604172,
		 * though some mental arithmetic is needed). These can also be treated as relative frequency weights,
		 * i.e. the total number of tests is rescaled relative to the weight of the current test (so, [10,1] 
		 * weights would consider there to be 1.1 tests for the first one and 11 tests for the second one).
		 */
        double total_weight = 0;
        auto wIt = weights.begin();
        for (size_t p = 0; p < pvalues.size(); ++p, ++wIt) { 
            double curp = pvalues[p].first;

            if (!ISNA(curp)) { 
                const double curweight = *wIt;
                total_weight += curweight;
                curp = divide(curp, curweight, log);

                if (curlowest > curp) {
                    curlowest = curp;
                    lowest = p;
                }
            }
        }

        if (R_FINITE(curlowest)) {
            curlowest = multiply(curlowest, total_weight, log);
            curlowest = bound_upper(curlowest, log);
        } else {
            curlowest = R_NaReal;
        }

        return std::make_pair(curlowest, lowest);
    }
};

#endif
