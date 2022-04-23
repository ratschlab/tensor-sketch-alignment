#pragma once

#include "hash_base.hpp"

#include "util/timer.hpp"
#include "util/utils.hpp"

#include <iostream>
#include <limits>
#include <random>

namespace ts { // ts = Tensor Sketch

/**
 * Naive implementation of weighted min-hash sketching. For more efficient implementations, see
 * https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/36928.pdf and
 * https://www.microsoft.com/en-us/research/wp-content/uploads/2010/06/ConsistentWeightedSampling2.pdf
 *
 * Given a set S, and a sequence s=s1...sn with elements from S, this class computes a vector
 * {hmin_1(s), hmin_2(s), ..., hmin_sketch_size(s)}, where hmin_k(s)=s_i, such that h_k(s_i, #s_i)
 * is the smallest of h_k(s_1, 1..#s_1), h_k(s_2, 1..#s_2), ..., h_k(s_n, 1..#s_n) and
 * h_k:Sx{1..n} -> {1..#set_size} is a random permuation of the elements in S and #s_i denotes the
 * number of occurences of s_i in the sequence s.
 * @tparam T the type of S's elements
 */
template <class T>
class WeightedMinHash : public HashBase<T> {
  public:
    /**
     * Constructs a weighted min-hasher for the given alphabet size which constructs sketches of the
     * given set size, dimension and maximum length.
     * @param set_size the number of elements in S,
     * @param sketch_dim the number of components (elements) in the sketch vector.
     * @param max_len maximum sequence length to be hashed.
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     */
    WeightedMinHash(T set_size,
                    size_t sketch_dim,
                    size_t max_len,
                    HashAlgorithm hash_algorithm,
                    uint32_t seed,
                    const std::string &name = "WMH",
                    size_t kmer_size = 1)
        : HashBase<T>(set_size, sketch_dim, max_len * set_size, hash_algorithm, seed, name, kmer_size),
          max_len(max_len) {}

    std::vector<T> compute(const std::vector<T> &kmers) {
        Timer timer("weighted_minhash");
        std::vector<T> sketch = std::vector<T>(this->sketch_dim);
        if (kmers.empty()) {
            return sketch;
        }

        for (size_t si = 0; si < this->sketch_dim; si++) {
            T min_char = T(0);
            T min_rank = std::numeric_limits<T>::max();
            std::unordered_map<T, uint32_t> cnts;
            for (const auto s : kmers) {
                T r = this->hash(si, s + cnts[s] * this->set_size);
                cnts[s]++;
#ifndef NDEBUG
                assert(cnts[s] != 0); // no overflow
                if (cnts[s] > max_len) {
                    throw std::invalid_argument("Kmer  " + std::to_string(s) + " repeats more than "
                                                + std::to_string(max_len)
                                                + " times. Set --max_len to a higher value.");
                }
#endif

                if (r < min_rank) {
                    min_rank = r;
                    min_char = s;
                }
            }
            sketch[si] = min_char;
        }
        return sketch;
    }

    /**
     * Computes the ordered min-hash sketch for the given sequence.
     * @param sequence the sequence to compute the ordered min-hash for
     * @param k-mer length; the sequence will be transformed into k-mers and the k-mers will be
     * hashed
     * @param number of characters in the alphabet over which sequence is defined
     * @return the ordered min-hash sketch of #sequence
     * @tparam C the type of characters in #sequence
     */
    template <typename C>
    std::vector<T> compute(const std::vector<C> &sequence, uint32_t k, uint32_t alphabet_size) {
        std::vector<T> kmers = seq2kmer<C, T>(sequence, k, alphabet_size);
        std::vector<T> sketch = compute(kmers);
        return sketch;
    }

    static T dist(const std::vector<T> &a, const std::vector<T> &b) {
        Timer timer("weighted_minhash_dist");
        return hamming_dist(a, b);
    }

  private:
    size_t max_len;
};

} // namespace ts
