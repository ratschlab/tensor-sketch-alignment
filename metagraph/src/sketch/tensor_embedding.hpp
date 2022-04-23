#pragma once

#include "immintrin.h" // for AVX
#include "nmmintrin.h" // for SSE4.2
#include "sequence/alphabets.hpp"
#include "sketch//sketch_base.hpp"
#include "util/multivec.hpp"
#include "util/timer.hpp"
#include "util/utils.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

namespace ts { // ts = Tensor Sketch

/**
 * Computes the tensor of subsequence counts for a given sequence as described in
 * https://www.biorxiv.org/content/10.1101/2020.11.13.381814v1.
 * In contrast with TensorSketch, this class does not do any sketching and preserves the full
 * tensor.
 * @tparam seq_type the type of elements in the sequences to be sketched.
 */
template <class seq_type>
class TensorEmbedding : public SketchBase<std::vector<double>, false> {
  public:
    /**
     * @param alphabet_size the number of elements in the alphabet S over which sequences are
     * defined (e.g. 4 for DNA)
     * @param normalize when true the counts will be normalized to relative frequencies with sum 1.
     */
    TensorEmbedding(seq_type alphabet_size,
                    uint32_t t,
                    const std::string &name = "Tensor",
                    bool normalize = true)
        : SketchBase<std::vector<double>, false>(name),
          alphabet_size(alphabet_size),
          t(t),
          normalize(normalize) {}

    void init() {}

    /**
     * Computes the sketch of the given sequence.
     * @param seq the sequence to be sketched
     * @return an array of size alphabet_size^t containing the sequence's sketch
     */
    std::vector<double> compute(const std::vector<seq_type> &seq) {
        // ts[i] contains the counts for subsequences of length i.
        Vec2D<double> ts(t + 1);
        {
            uint32_t num_tmers = 1;
            for (auto &t : ts) {
                t.resize(num_tmers);
                num_tmers *= alphabet_size;
            }
        }

        // The base case is the one empty sequence.
        ts[0][0] = 1;

        // The number of successfully read nucleotides.
        int32_t length = 0;
        for (auto s : seq) {
            // TODO(ragnar): Figure out a nice way to deal with uncertain reads.
            if (s < 0 || s >= alphabet_size)
                continue;
            length += 1;
            for (int32_t i = static_cast<int32_t>(t) - 1; i >= 0; --i) {
                for (size_t j = 0; j < ts[i].size(); ++j) {
                    ts[i + 1][alphabet_size * j + s] += ts[i][j];
                }
            }
        }

        if (normalize) {
            double nchooset = 1;
            for (uint32_t i = 0; i < t; ++i) {
                nchooset = nchooset * (length - i) / (i + 1);
            }

            for (auto &c : ts.back()) {
                c /= nchooset;
            }
        }

        return std::move(ts.back());
    }

    static double dist(const std::vector<double> &a, const std::vector<double> &b) {
        Timer timer("full_tensor_dist");
        return l2_dist(a, b);
    }

  protected:
    /** Size of the alphabet over which sequences to be sketched are defined, e.g. 4 for DNA. */
    const seq_type alphabet_size;

    /** The length of the subsequences considered for sketching. */
    const uint32_t t;

    /** Whether to normalize the counts to relative frequencies with sum 1. */
    const bool normalize;
};

} // namespace ts
