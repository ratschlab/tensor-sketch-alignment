#pragma once

#include "hash_base.hpp"

#include "util/timer.hpp"
#include "util/utils.hpp"

#include <cstdint>
#include <limits>
#include <random>

namespace ts { // ts = Tensor Sketch

/**
 * Implements min-hash-based sketching, as described in https://en.wikipedia.org/wiki/MinHash.
 * Given a set S, and a sequence s=s1...sn with elements from S, this class computes a vector
 * {hmin_1(s), hmin_2(s), ..., hmin_sketch_size(s)}, where hmin_k(s)=s_i, such that h_k(s_i) is the
 * smallest of h_k(s_1), h_k(s_2), ..., h_k(s_n) and h_k:S->{1..#set_size} is a random permuation of
 * the elements in S.
 * This class assumes that S= {0,1,2....,#set_size}.
 * @tparam T the type of S's elements.
 */
template <class T>
class MinHash : public HashBase<T> {
  public:
    /**
     * Constructs a min-hasher for the given alphabet size which constructs sketches of the set size
     * and sketch dimension.
     * @param set_size the number of elements in S,
     * @param sketch_dim the number of components (elements) in the sketch vector.
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     */
    MinHash(T set_size,
            size_t sketch_dim,
            HashAlgorithm hash_algorithm,
            uint32_t seed,
            const std::string &name = "MH",
            size_t kmer_size = 1)
        : HashBase<T>(set_size, sketch_dim, set_size, hash_algorithm, seed, name, kmer_size) {}

    /**
     * Computes the min-hash sketch for the given kmers.
     * @param kmers kmers extracted from a sequence
     * @return the min-hash sketch of #kmers
     */
    std::vector<T> compute(const std::vector<T> &kmers) {
        Timer timer("minhash");
        std::vector<T> sketch(this->sketch_dim);
        if (kmers.empty()) {
            return sketch;
        }

        for (size_t si = 0; si < this->sketch_dim; si++) {
            T min_char = T(0);
            T min_rank = std::numeric_limits<T>::max();
            for (auto s : kmers) {
                T hash = this->hash(si, s);
                if (hash < min_rank) {
                    min_rank = hash;
                    min_char = s;
                }
            }
            sketch[si] = min_char;
        }
        return sketch;
    }

    /**
     * Computes the min-hash sketch for the given sequence.
     * @param sequence the sequence to compute the min-hash for
     * @param k-mer length; the sequence will be transformed into k-mers and the k-mers will be
     * hashed
     * @param number of characters in the alphabet over which sequence is defined
     * @return the min-hash sketch of sequence
     * @tparam C the type of characters in the sequence
     */
    template <typename C>
    std::vector<T> compute(const std::vector<C> &sequence, uint32_t k, uint32_t alphabet_size) {
        std::vector<T> kmers = seq2kmer<C, T>(sequence, k, alphabet_size);
        return compute(kmers);
    }

    static T dist(const std::vector<T> &a, const std::vector<T> &b) {
        Timer timer("minhash_dist");
        return hamming_dist(a, b);
    }
};
} // namespace ts
