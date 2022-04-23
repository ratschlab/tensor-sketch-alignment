#pragma once

#include "immintrin.h" // for AVX
#include "nmmintrin.h" // for SSE4.2
#include "sketch//sketch_base.hpp"
#include "util/multivec.hpp"
#include "util/timer.hpp"
#include "util/utils.hpp"

#include <cassert>
#include <cmath>
#include <deque>
#include <random>

namespace ts { // ts = Tensor Sketch

/**
 * Computes tensor sketches for a given sequence as described in
 * https://www.biorxiv.org/content/10.1101/2020.11.13.381814v1, with the additional limitation that
 * the subsequences must be made of continuous blocks of a certain size. The adaptation of the
 * recurrence formula for this case is at https://go.grlab.org/tensor_block.
 * For block_size=1, the normal the #TensorBlock sketch is identical with #Tensor sketch.
 * @tparam seq_type the type of elements in the sequences to be sketched.
 */
template <class seq_type>
class TensorBlock : public SketchBase<std::vector<double>, false> {
  public:
    // Tensor sketch output should be transformed if the command line flag is set.
    constexpr static bool transform_sketches = false;

    /**
     * @param block_size only subsequences formed out of block_size continuous elements are sketched
     * @param alphabet_size the number of elements in the alphabet S over which sequences are
     * defined (e.g. 4 for DNA, 20 for protein, etc.)
     * @param sketch_dim the dimension of the embedded (sketched) space, denoted by D in the paper
     * @param subsequence_len the length of the subsequences considered for sketching, denoted by t
     * in the paper
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     */
    TensorBlock(seq_type alphabet_size,
                size_t sketch_dim,
                size_t subsequence_len,
                uint8_t block_size,
                uint32_t seed,
                const std::string &name = "TSB")
        : SketchBase<std::vector<double>, false>(name),
          block_size(block_size),
          alphabet_size(alphabet_size),
          sketch_dim(sketch_dim),
          subsequence_len(subsequence_len),
          rng(seed) {
        assert(block_size > 0 && subsequence_len > 0 && subsequence_len % block_size == 0);
        init();
    }

    void init() {
        hashes = new2D<seq_type>(subsequence_len, alphabet_size);
        signs = new2D<bool>(subsequence_len, alphabet_size);

        std::uniform_int_distribution<seq_type> rand_hash2(0, sketch_dim - 1);
        std::uniform_int_distribution<seq_type> rand_bool(0, 1);

        for (size_t h = 0; h < subsequence_len; h++) {
            for (size_t c = 0; c < alphabet_size; c++) {
                hashes[h][c] = rand_hash2(rng);
                signs[h][c] = rand_bool(rng);
            }
        }
    }

    /**
     * Computes the sketch of the given sequence.
     * @param seq the sequence to be sketched
     * @return an array of size #sketch_dim containing the sequence's sketch
     */
    std::vector<double> compute(const std::vector<seq_type> &seq) {
        Timer timer("tensor_sketch");
        // Tp corresponds to T+, Tm to T- in the paper; Tp[0], Tm[0] are sentinels and contain the
        // initial condition for empty strings; Tp[p], Tm[p] at step i represent the partial sketch
        // when considering hashes h1...hp, over the prefix x1...xi. The final result is then
        // Tp[t]-Tm[t], where t is #sequence_len
        // since the recurrence formula references the T_[1:N-k], i.e. the element situated
        // k=block_size positions behind, we need to always keep the last block_size Tp and Tm
        // matrices. At each iteration we create a new pair of Tp and Tm and then discard the oldest
        // Tp/Tn pair.
        // TODO(ddanciu): use a circular queue on top of vector instead
        std::deque<Vec2D<double>> Tp;
        std::deque<Vec2D<double>> Tm;

        // The number of blocks.
        uint32_t m = subsequence_len / block_size;

        for (uint32_t i = 0; i < block_size; ++i) {
            Tp.push_back(new2D<double>(m + 1, sketch_dim, 0));
            Tm.push_back(new2D<double>(m + 1, sketch_dim, 0));
            // the initial condition states that the sketch for the empty string is (1,0,..)
            Tp.back()[0][0] = 1;
        }

        // the are the "new" Tp and Tm, computed at every iteration and appended to Tp and Tm
        auto nTp = new2D<double>(m + 1, sketch_dim, 0);
        auto nTm = new2D<double>(m + 1, sketch_dim, 0);
        for (uint32_t i = block_size - 1; i < seq.size(); i++) {
            uint32_t block_count = std::min(m, (i + 1) / block_size);
            // must traverse in reverse order, to avoid overwriting the values of Tp and Tm before
            // they are used in the recurrence
            // p must be a multiple of block_size
            for (uint32_t bc = block_count; bc > 0; bc--) {
                uint32_t p = bc * block_size;
                double z = bc / (i + 1.0 - p + bc); // probability that the last index is i
                seq_type r = 0;
                bool s = true;
                for (uint32_t j = 0; j < block_size; ++j) {
                    r += hashes[p - j - 1][seq[i - j]];
                    s = s == signs[p - j - 1][seq[i - j]];
                }
                r %= sketch_dim;
                if (s) {
                    nTp[bc] = this->shift_sum(Tp.back()[bc], Tp[0][bc - 1], r, z);
                    nTm[bc] = this->shift_sum(Tm.back()[bc], Tm[0][bc - 1], r, z);
                } else {
                    nTp[bc] = this->shift_sum(Tp.back()[bc], Tm[0][bc - 1], r, z);
                    nTm[bc] = this->shift_sum(Tm.back()[bc], Tp[0][bc - 1], r, z);
                }
            }
            nTp[0][0] = 1;
            Tp.push_back(std::move(nTp));
            Tm.push_back(std::move(nTm));
            nTp = std::move(Tp.front());
            nTm = std::move(Tm.front());
            for (uint32_t j = 0; j < m + 1; ++j) {
                std::fill(nTp[j].begin(), nTp[j].end(), 0);
                std::fill(nTm[j].begin(), nTm[j].end(), 0);
            }
            Tp.pop_front();
            Tm.pop_front();
        }
        std::vector<double> sketch(sketch_dim, 0);
        for (uint32_t l = 0; l < sketch_dim; l++) {
            sketch[l] = Tp.back()[m][l] - Tm.back()[m][l];
        }

        return sketch;
    }

    /** Sets the hash and sign functions to predetermined values for testing */
    void set_hashes_for_testing(const Vec2D<seq_type> &h, const Vec2D<bool> &s) {
        hashes = h;
        signs = s;
    }

    static double dist(const std::vector<double> &a, const std::vector<double> &b) {
        Timer timer("tensor_sketch_dist");
        return l2_dist(a, b);
    }

  protected:
    /** Computes (1-z)*a + z*b_shift */
    inline std::vector<double> shift_sum(const std::vector<double> &a,
                                         const std::vector<double> &b,
                                         seq_type shift,
                                         double z) {
        assert(a.size() == b.size());
        size_t len = a.size();
        std::vector<double> result(a.size());
        for (uint32_t i = 0; i < a.size(); i++) {
            result[i] = (1 - z) * a[i] + z * b[(len + i - shift) % len];
            assert(result[i] <= 1 + 1e-5 && result[i] >= -1e-5);
        }
        return result;
    }

    /** The size of the block each subsequence is made of. Must be a divisor of #subsequence_len.
     * A block size of 1 means that the class is operating on a character basis. */
    uint8_t block_size;

    /** Size of the alphabet over which sequences to be sketched are defined, e.g. 4 for DNA */
    seq_type alphabet_size;
    /** Number of elements in the sketch, denoted by D in the paper */
    uint32_t sketch_dim;
    /** The length of the subsequences considered for sketching, denoted by t in the paper */
    uint8_t subsequence_len;

    /**
     * Denotes the hash functions h1,....ht:A->{1....D}, where t is #subsequence_len and D is
     * #sketch_dim
     */
    Vec2D<seq_type> hashes;

    /** The sign functions s1...st:A->{-1,1} */
    Vec2D<bool> signs;

    std::mt19937 rng;
};

} // namespace ts
