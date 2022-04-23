#pragma once

#include "hash_base.hpp"

#include "util/utils.hpp"

#include <iostream>
#include <random>
#include <string>

namespace ts { // ts = Tensor Sketch

/**
 * Naive implementation of the Ordered MinHash sketching method described in:
 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6612865/
 *
 * @tparam T the type of element in the sequences to be sketched
 */
template <class T>
class OrderedMinHash : public HashBase<T> {
  public:
    /**
     * @param set_size the number of elements in S
     * @param sketch_dim the number of components (elements) in the sketch vector.
     * @param max_len maximum sequence length to be hashed.
     * @param tup_len the sketching will select the tup_len lowest values for each hash function
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     */
    OrderedMinHash(T set_size,
                   size_t sketch_dim,
                   size_t max_len,
                   size_t tup_len,
                   HashAlgorithm hash_algorithm,
                   uint32_t seed,
                   const std::string &name = "OMH",
                   size_t kmer_size = 1)
        : HashBase<T>(set_size, sketch_dim, set_size * max_len, hash_algorithm, seed, name, kmer_size),
          max_len(max_len),
          tup_len(tup_len) {}

    Vec2D<T> compute_2d(const std::vector<T> &kmers) {
        Vec2D<T> sketch(this->sketch_dim);
        if (kmers.size() < tup_len) {
            throw std::invalid_argument("Sequence of kmers must be longer than tuple length");
        }
        for (size_t pi = 0; pi < this->sketch_dim; pi++) {
            std::unordered_map<size_t, uint32_t> counts;
            std::vector<std::pair<T, size_t>> ranks;
            for (size_t i = 0; i < kmers.size(); i++) {
                auto s = kmers[i];
                ranks.push_back({ this->hash(pi, s + this->set_size * counts[s]), i });
                counts[s]++;
#ifndef NDEBUG
                assert(counts[s] != 0); // no overflow
                if (counts[s] > max_len) {
                    throw std::invalid_argument("Kmer  " + std::to_string(s) + " repeats more than "
                                                + std::to_string(max_len)
                                                + " times. Set --max_len to a higher value.");
                }
#endif
            }
            std::sort(ranks.begin(), ranks.end());
            std::vector<size_t> tup;
            for (auto pair = ranks.begin(); pair != ranks.end() && pair != ranks.begin() + tup_len;
                 pair++) {
                tup.push_back(pair->second);
            }
            std::sort(tup.begin(), tup.end()); // sort indices of kmers
            for (auto idx : tup)
                sketch[pi].push_back(kmers[idx]);
        }
        return sketch;
    }

    std::vector<T> compute(const std::vector<T> &kmers) {
        Timer timer("ordered_minhash");
        std::vector<T> sketch;

        Vec2D<T> sketch2D = compute_2d(kmers);
        for (const auto &tuple : sketch2D) {
            T sum = 0;
            for (const auto &item : tuple) {
                sum = sum * this->set_size + item; // TODO: deal with overflows
            }
            sketch.push_back(sum);
        }

        return sketch;
    }

    /**
     * Computes the ordered min-hash sketch for the given sequence.
     * @param sequence the sequence to compute the ordered min-hash for
     * @param k-mer length; the sequence will be transformed into k-mers and the k-mers will be
     * hashed
     * @param number of characters in the alphabet over which sequence is defined
     * @return the ordered min-hash sketch of sequence
     * @tparam C the type of characters in the sequence
     */
    template <typename C>
    std::vector<T> compute(const std::vector<C> &sequence, uint32_t k, uint32_t alphabet_size) {
        return compute(seq2kmer<C, T>(sequence, k, alphabet_size));
    }

    static T dist(const std::vector<T> &a, const std::vector<T> &b) {
        Timer timer("ordered_minhash_dist");
        return hamming_dist(a, b);
    }

  private:
    size_t max_len;
    size_t tup_len;
};

} // namespace ts
