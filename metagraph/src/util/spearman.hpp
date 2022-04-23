#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// Function returns the 1-based rank vector of a set of observations v
template <typename T>
std::vector<double> rankify(const std::vector<T> &v) {
    std::vector<T> sorted = v;
    std::sort(begin(sorted), end(sorted));

    std::vector<double> result(v.size());

    for (size_t i = 0; i < v.size(); i++) {
        const auto lb = std::lower_bound(std::begin(sorted), std::end(sorted), v[i]);
        const auto ub = std::upper_bound(std::begin(sorted), std::end(sorted), v[i]);
        const size_t r = 1 + (lb - std::begin(sorted)), s = ub - lb;

        // Use Fractional Rank formula fractional_rank = r + (s-1)/2
        result[i] = r + (s - 1) * 0.5;
    }

    return result;
}

/* Compute the Pearson correlation coefficient of a and b */
template <typename T>
double pearson(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    T sum_a = 0, sum_b = 0, sum_ab = 0;
    T square_sum_a = 0, square_sum_b = 0;

    for (size_t i = 0; i < a.size(); i++) {
        sum_a = sum_a + a[i];
        sum_b = sum_b + b[i];
        sum_ab = sum_ab + a[i] * b[i];
        square_sum_a = square_sum_a + a[i] * a[i];
        square_sum_b = square_sum_b + b[i] * b[i];
    }

    // compute variances
    T var_a = a.size() * square_sum_a - sum_a * sum_a;
    T var_b = a.size() * square_sum_b - sum_a * sum_b;
    // treat degenerate cases
    if (var_a == 0 && var_b == 0) {
        return 1;
    }
    if (var_a == 0 || var_b == 0) {
        return 0;
    }

    return (a.size() * sum_ab - sum_a * sum_b) / std::sqrt(var_a * var_b);
}

template <typename T, typename U>
double spearman(const std::vector<T> &a, const std::vector<U> &b) {
    std::vector<double> rank1 = rankify(a);
    std::vector<double> rank2 = rankify(b);
    return pearson(rank1, rank2);
}
