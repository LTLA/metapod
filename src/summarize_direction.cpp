#include "Rcpp.h"
#include "utils.h"

int choose_direction (int left, int right) {
    return static_cast<int>(left > 0) + static_cast<int>(right > 0) * 2;
}

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector summarize_parallel_direction(Rcpp::List effects, Rcpp::List influential, double threshold) {
    parallel_vectors<Rcpp::NumericVector> Effects(effects);
    parallel_vectors<Rcpp::LogicalVector> Influential(influential);
    if (Effects.nvectors != Influential.nvectors || Effects.nelements != Influential.nelements) {
        throw std::runtime_error("'influential' should have the same structure as 'effects'");
    }

    Rcpp::IntegerVector output(Effects.nelements);
    for (size_t i = 0; i < Effects.nelements; ++i) {
        int left = 0, right = 0;
        for (size_t j = 0; j < Effects.nvectors; ++j) {
            if (Influential.vectors[j][i]) {
                auto effect = Effects.vectors[j][i];
                if (effect < threshold) {
                    ++left;
                } else if (effect > threshold) {
                    ++right;
                }
            }
        }
        output[i] = choose_direction(left, right) + 1;
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector summarize_grouped_direction(Rcpp::NumericVector effects, Rcpp::IntegerVector runs, Rcpp::LogicalVector influential, double threshold) {
    size_t counter = 0;
    Rcpp::IntegerVector output(runs.size());
    if (effects.size() != influential.size()) {
        throw std::runtime_error("'effects' and 'influential' should have the same length");
    }

    for (size_t g = 0; g < runs.size(); ++g) {
        int left = 0, right = 0;
        for (int r = 0; r < runs[g]; ++r, ++counter) {
            if (counter >= effects.size()) {
                throw std::runtime_error("'sum(runs)' is not the same as 'length(effects)'");
            }
            if (influential[counter]) {
                auto effect = effects[counter];
                if (effect < threshold) {
                    ++left;
                } else if (effect > threshold) {
                    ++right;
                }
            }
        }
        output[g] = choose_direction(left, right) + 1;
    }

    if (counter != effects.size()) {
        throw std::runtime_error("'sum(runs)' is not the same as 'length(effects)'");
    }

    return output;
}
