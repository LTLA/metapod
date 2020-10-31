#ifndef P_WILKINSON_H
#define P_WILKINSON_H

#include "Rcpp.h"
#include "utils.h"

class p_wilkinson {
public: 
    p_wilkinson(size_t mn, double mp) : min_num(std::max(size_t(1), mn)), min_prop(mp) {}

    template<class V, class Y>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights, bool log, Y& influencers) const {
        size_t index = compute_index(pvalues.size(), min_num, min_prop);
        std::partial_sort(pvalues.begin(), pvalues.begin() + index + 1, pvalues.end());
        double outp = R::pbeta(pvalues[index].first, index + 1, pvalues.size() - index, 1, log);
        return std::make_pair(outp, pvalues[index].second);
    }
private:
    const size_t min_num;
    const double min_prop;
};

#endif
