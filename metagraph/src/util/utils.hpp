#pragma once

#include "util/multivec.hpp"
#include "util/timer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

namespace ts { // ts = Tensor Sketch

/**
 * Extracts k-mers from a sequence. The k-mer is treated as a number in base alphabet_size and then
 * converted to decimal, i.e. the sequence s1...sk is converted to s1*S^(k-1) + s2*S^(k-2) + ... +
 * sk, where k is the k-mer size.
 * @tparam chr types of elements in the sequence
 * @tparam kmer type that stores a kmer
 * @param seq the sequence to extract kmers from
 * @param kmer_size number of characters in a kmer
 * @param alphabet_size size of the alphabet
 * @return the extracted kmers, as integers converted from base #alphabet_size
 */
template <class chr, class kmer>
std::vector<kmer> seq2kmer(const std::vector<chr> &seq, uint8_t kmer_size, uint8_t alphabet_size) {
    Timer timer("seq2kmer");
    if (seq.size() < (size_t)kmer_size) {
        return std::vector<kmer>();
    }

    std::vector<kmer> result(seq.size() - kmer_size + 1, 0);

    kmer c = 1;
    for (uint8_t i = 0; i < kmer_size; i++) {
        result[0] += c * seq[i];
        c *= alphabet_size;
    }
    c /= alphabet_size;

    for (size_t i = 0; i < result.size() - 1; i++) {
        kmer base = result[i] - seq[i];
        assert(base % alphabet_size == 0);
        result[i + 1] = base / alphabet_size + seq[i + kmer_size] * c;
    }
    return result;
}

template <class T>
T l1_dist(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    T res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        auto el = std::abs(a[i] - b[i]);
        res += el;
    }
    return res;
}


template <class T>
T l2_dist(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    T res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        auto el = std::abs(a[i] - b[i]);
        res += el * el;
    }
    return res;
}


template <class T>
T l1_dist2D_minlen(const Vec2D<T> &a, const Vec2D<T> &b) {
    auto len = std::min(a.size(), b.size());
    T val = 0;
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < a[i].size() and j < b[i].size(); j++) {
            auto el = std::abs(a[i][j] - b[i][j]);
            val += el;
        }
    }
    return val;
}

template <class T>
T l2_dist2D_minlen(const Vec2D<T> &a, const Vec2D<T> &b) {
    auto len = std::min(a.size(), b.size());
    T val = 0;
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < a[i].size() and j < b[i].size(); j++) {
            auto el = (a[i][j] - b[i][j]);
            val += el * el;
        }
    }
    return val;
}


template <class T>
T hamming_dist(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    T diff = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            diff++;
        }
    }
    return diff;
}

template <class seq_type>
int lcs(const std::vector<seq_type> &s1, const std::vector<seq_type> &s2) {
    size_t m = s1.size();
    size_t n = s2.size();
    //        int L[m + 1][n + 1];
    Vec2D<int> L(m + 1, std::vector<int>(n + 1, 0));
    for (size_t i = 0; i <= m; i++) {
        for (size_t j = 0; j <= n; j++) {
            if (i == 0 || j == 0) {
                L[i][j] = 0;
            } else if (s1[i - 1] == s2[j - 1]) {
                L[i][j] = L[i - 1][j - 1] + 1;
            } else {
                L[i][j] = std::max(L[i - 1][j], L[i][j - 1]);
            }
        }
    }
    return L[m][n];
}

template <class seq_type>
size_t lcs_distance(const std::vector<seq_type> &s1, const std::vector<seq_type> &s2) {
    return s1.size() + s2.size() - 2 * lcs(s1, s2);
}

template <class seq_type>
size_t edit_distance(const std::vector<seq_type> &s1, const std::vector<seq_type> &s2) {
    Timer timer("edit_distance");
    const size_t m(s1.size());
    const size_t n(s2.size());

    if (m == 0)
        return n;
    if (n == 0)
        return m;

    auto costs = std::vector<size_t>(n + 1);

    for (size_t k = 0; k <= n; k++)
        costs[k] = k;

    size_t i = 0;
    for (auto it1 = s1.begin(); it1 != s1.end(); ++it1, ++i) {
        costs[0] = i + 1;
        size_t corner = i;

        size_t j = 0;
        for (auto it2 = s2.begin(); it2 != s2.end(); ++it2, ++j) {
            size_t upper = costs[j + 1];
            if (*it1 == *it2) {
                costs[j + 1] = corner;
            } else {
                size_t t(upper < corner ? upper : corner);
                costs[j + 1] = (costs[j] < t ? costs[j] : t) + 1;
            }

            corner = upper;
        }
    }

    size_t result = costs[n];

    return result;
}

template <class T, class = is_u_integral<T>>
T int_pow(T x, T pow) {
    T result = 1;
    for (;;) {
        if (pow & 1)
            result *= x;
        pow >>= 1;
        if (!pow)
            break;
        x *= x;
    }

    return result;
}

std::string
flag_values(char delimiter = ' ', bool skip_empty = false, bool include_flagfile = true);

// If the -o output flag is set, this writes a small shell script <output>.meta containing the
// command line used to generate the output.
void write_output_meta();

// A simple wrapper around std::apply that applies a given lambda on each element of a tuple.
template <typename F, typename T>
void apply_tuple(F &&f, T &tuple_t) {
    std::apply([&](auto &...t) { (f(t), ...); }, tuple_t);
}


// A simple wrapper around std::apply that applies f on pairs of elements of two tuples.
template <typename F, typename T, typename U>
void apply_tuple(F &&f, T &tuple_t, U &tuple_u) {
    std::apply([&](auto &...t) { std::apply([&](auto &...u) { (f(t, u), ...); }, tuple_u); },
               tuple_t);
}


std::pair<double, double> avg_stddev(const std::vector<double> &v);

// v must be sorted.
double median(const std::vector<double> &v);

} // namespace ts
