#pragma once

#include "util/utils.hpp"

#include <cstddef>
#include <random>
#include <vector>


namespace ts {

class Int32Flattener {
  public:
    using sketch_type = std::vector<uint32_t>;

    Int32Flattener(uint32_t flat_dim, uint32_t sketch_dim, uint32_t max_len, uint32_t seed)
        : flat_dim(flat_dim), sketch_dim(sketch_dim), max_len(max_len) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<uint32_t> distribution;
        rand_proj = new2D<uint32_t>(flat_dim, this->max_len * sketch_dim * 2);
        for (auto &v : rand_proj) {
            for (auto &e : v) {
                e = distribution(gen);
            }
        }
    }

    std::vector<uint32_t> flatten(const Vec2D<double> &sketch) {
        Timer timer("Int32Flattener");
        assert(rand_proj.size() == flat_dim);
        std::vector<uint32_t> v(flat_dim, 0);
        for (uint32_t s1 = 0; s1 < flat_dim; s1++) {
            for (uint32_t s2 = 0; s2 < 32; s2++) { // iterate over 32 bits
                size_t j = s1 % sketch_dim;
                double val = 0;
                for (size_t i = 0; i < sketch.size(); i++) {
                    auto bit = rand_proj[s1][i * sketch_dim + j] >> s2; // random bit
                    val += (bit & 1) ? sketch[i][j] : -sketch[i][j];
                }
                v[s1] = (v[s1] << 1) + std::signbit(val); // insert sgn(val) into v[s1]
            }
        }
        return v;
    }

    static double dist(const std::vector<uint32_t> &v1, const std::vector<uint32_t> &v2) {
        Timer timer("Int32Flattener_dist");
        assert(v1.size() == v2.size());
        std::vector<double> d(v1.size());
        double val = 0;
        for (size_t i = 0; i < d.size(); i++) {
            val += __builtin_popcount(v1[i] ^ v2[i]);
        }
        return val;
    }

  private:
    uint32_t flat_dim;
    uint32_t sketch_dim;
    uint32_t max_len;
    Vec2D<uint32_t> rand_proj;
};


class DoubleFlattener {
  public:
    using sketch_type = std::vector<double>;

    DoubleFlattener(uint32_t output_dim,
                    uint32_t input_dim,
                    uint32_t input_max_len,
                    uint32_t seed)
        : flat_dim(output_dim), sketch_dim(input_dim), max_len(input_max_len) {
        std::mt19937 gen(seed);
        std::cauchy_distribution<double> distribution(0, 1.0);
        rand_proj = new2D<double>(this->flat_dim, this->max_len * input_dim * 2);
        for (auto &v : rand_proj) {
            for (double &e : v) {
                e = distribution(gen);
            }
        }
    }

    std::vector<double> flatten(const Vec2D<double> &sketch) {
        Timer timer("DoubleFlattener");
        assert(rand_proj.size() == flat_dim);
        std::vector<double> v(this->flat_dim, 0);
        for (size_t s = 0; s < this->flat_dim; s++) {
            size_t j = s % this->sketch_dim;
            for (size_t i = 0; i < sketch.size(); i++) {
                v[s] += rand_proj[s][i * this->sketch_dim + j] * sketch[i][j];
            }
            v[s] /= (double)(sketch.size()
                             * sketch[0]
                                       .size()); // divide by number of elements to compute the mean
        }

        return v;
    }

    static double dist(const std::vector<double> &v1, const std::vector<double> &v2) {
        Timer timer("DoubleFlattener_dist");
        assert(v1.size() == v2.size());
        std::vector<double> d(v1.size());
        for (size_t i = 0; i < d.size(); i++) {
            d[i] = abs(v1[i] - v2[i]);
        }
        std::sort(d.begin(), d.end());
        return d[d.size() / 2]; // return the median
    }

  private:
    uint32_t flat_dim;
    uint32_t sketch_dim;
    uint32_t max_len;
    Vec2D<double> rand_proj;
};

} // namespace ts
