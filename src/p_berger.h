#ifndef P_BERGER_H
#define P_BERGER_H

#include "utils.h"
#include <algorithm>

class p_berger {
public:    
    template<class V>
    std::pair<double, size_t> operator()(IndexedPValues& pvalues, const V& weights) const {
        auto out = std::max_element(pvalues.begin(), pvalues.end());
        return std::make_pair(out->first, out - pvalues.begin());
    }
};

#endif
