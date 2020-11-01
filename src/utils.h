#ifndef UTILS_H
#define UTILS_H

#include <deque>
#include <algorithm>
#include <cmath>

typedef std::deque<std::pair<double, int> > IndexedPValues;

inline double divide(double x, double y, bool log) {
    if (log) {
        return x - std::log(y);
    } else {
        return x / y;
    }
}

inline double multiply(double x, double y, bool log) {
    if (log) {
        return x + std::log(y);
    } else {
        return x * y;
    }
}

inline double bound_upper(double x, bool log) {
    if (log) {
        return std::min(x, 0.0);
    } else {
        return std::min(x, 1.0);
    }
}

inline double bound_lower(double x, bool log) {
    if (log) {
        return x;
    } else {
        return std::max(x, 0.0);
    }
}

inline size_t compute_index (size_t ntests, size_t min_num, double min_prop) {
    size_t index=std::max(
        min_num, 
        static_cast<size_t>(std::ceil(min_prop * static_cast<double>(ntests)))
    );

    index = std::min(index, ntests);

    if (index!=0) {
        --index; // zero-indexed.
    }

    return index;
}


#endif
