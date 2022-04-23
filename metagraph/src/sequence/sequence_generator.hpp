#pragma once

#include "util/utils.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <utility>

namespace ts { // ts = Tensor Sketch

class SeqGen {
  public:
    SeqGen(uint8_t alphabet_size,
           bool fix_len,
           uint32_t num_seqs,
           uint32_t seq_len,
           uint32_t group_size,
           double max_mutation_rate,
           double min_mutation_rate,
           std::string phylogeny_shape)
        : alphabet_size(alphabet_size),
          fix_len(fix_len),
          num_seqs(num_seqs),
          seq_len(seq_len),
          group_size(group_size),
          max_mutation_rate(max_mutation_rate),
          min_mutation_rate(min_mutation_rate),
          phylogeny_shape(std::move(phylogeny_shape)) {
        assert(group_size >= 2 && "group size<=1 leads to completely independent sequences");
    }

    /**
     * Generate sequences divided into independent groups of size `group_size`, and store
     * ingroup_pairs within each group in `ingroup_pairs`
     * @tparam T character type
     * @tparam C index type
     * @param seqs generated sequences
     * @param pairs sequence ingroup_pairs within each group
     */
    template <class T>
    Vec2D<T> generate_seqs() {
        if (phylogeny_shape == "pair") { // shape=path: implemented as path & group_size=2
            phylogeny_shape = "path";
            group_size = 2;
        }
        Vec2D<T> seqs;
        seqs.reserve(num_seqs);
        while (seqs.size() < num_seqs) {
            Vec2D<T> group;

            // tree-like: g1->g2, add g2 to pool, g1->g3, g2->g4, add g3, g4 to pool
            if (phylogeny_shape == "tree") {
                group = Vec2D<T>(1);
                random_sequence(group[0], seq_len);
                Vec2D<T> children;
                while (group.size() < group_size) {
                    for (auto &seq : group) {
                        std::vector<T> ch;
                        mutate(seq, ch);
                        children.push_back(seq);
                        children.push_back(ch);
                    }
                    std::swap(group, children);
                    children.clear();
                }
            } else if (phylogeny_shape == "path") { // path-like: g0->g1->g2->g3->...
                group = Vec2D<T>(group_size);
                random_sequence(group[0], seq_len);
                for (size_t i = 0; i < group_size - 1; i++) {
                    mutate(group[i], group[i + 1]);
                }
            } else if (phylogeny_shape == "star") { // star-like: g0->g1, g0->g2,g0->g3 ...
                group = Vec2D<T>(1);
                random_sequence(group[0], seq_len);
                for (size_t i = 1; i < group_size; i++) {
                    mutate(group[0], group[i]);
                }
            }

            group.resize(group_size);
            seqs.insert(seqs.end(), group.begin(), group.end());
            if (seqs.size() > num_seqs) {
                seqs.resize(num_seqs);
            }
        }
        return seqs;
    }

    template <class T>
    void ingroup_pairs(std::vector<std::pair<T, T>> &pairs) {
        for (size_t go = 0; go < num_seqs; go += group_size) { // group-offset
            for (size_t i = 0; i < group_size && go + i < num_seqs; i++) { // group-member i
                for (size_t j = i + 1; j < group_size && go + j < num_seqs; j++) { // group-member j
                    pairs.push_back({ go + i, go + j });
                }
            }
        }
    }


  private:
    template <class T>
    void mutate(const std::vector<T> &ref, std::vector<T> &seq) {
        std::uniform_real_distribution<double> unif(min_mutation_rate, max_mutation_rate);
        mutate(ref, seq, unif(gen));
        if (fix_len)
            make_fix_len(seq);
    }

    /**
     * Mutate seq from ref, by mutating each position with the probability = `rate`
     * @tparam T element type in the sequence
     * @param ref
     * @param seq mutated sequence
     * @param rate probability of mutation at each index
     */
    template <class T>
    void mutate(const std::vector<T> &ref, std::vector<T> &seq, double rate) {
        assert((rate >= 0.0) && (rate <= 1.0) && " rate must be strictly in the range [0,1]");
        // probabilities for each index position: no mutation, insert, delete, substitute
        std::discrete_distribution<int> mut { 1 - rate, rate / 3, rate / 3, rate / 3 };
        // the range chosen such that (sub_char+ref % alphabet_size) will different from ref
        std::uniform_int_distribution<T> sub_char(1, alphabet_size - 1);
        // random character from the alphabet
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        for (size_t i = 0; i < ref.size(); i++) {
            switch (mut(gen)) {
                case 0: { // no mutation
                    seq.push_back(ref[i]);
                    break;
                }
                case 1: { // insert
                    seq.push_back(rand_char(gen));
                    i--; // init_tensor_slide_params negate the increment
                    break;
                }
                case 2: { // delete
                    break;
                }
                case 3: { // substitute
                    seq.push_back((sub_char(gen) + ref[i]) % alphabet_size);
                    break;
                }
            }
        }
    }


    template <class T>
    void make_fix_len(std::vector<T> &seq) {
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        if (seq.size() > seq_len) {
            seq = std::vector<T>(seq.begin(), seq.end());
        } else if (seq.size() < seq_len) {
            while (seq.size() < seq_len) {
                seq.push_back(rand_char(gen));
            }
        }
    }

    /**
     * Generate a random sequence of length `len`
     * @tparam T
     * @param seq : the result will be stored in `seq`
     * @param len : length of the random sequence
     */
    template <class T>
    void random_sequence(std::vector<T> &seq, size_t len) {
        seq.resize(len);
        std::uniform_int_distribution<T> rand_char(0, alphabet_size - 1);
        for (uint32_t i = 0; i < len; i++) {
            seq[i] = rand_char(gen);
        }
    }


  private:
    std::mt19937 gen = std::mt19937(341234);

    uint8_t alphabet_size;
    bool fix_len;
    uint32_t num_seqs;
    uint32_t seq_len;
    uint32_t group_size;
    double max_mutation_rate;
    double min_mutation_rate;
    std::string phylogeny_shape;
};

} // namespace ts
