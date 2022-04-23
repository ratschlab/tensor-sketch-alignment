#pragma once

#include <exception>
#include <string>
#include <utility>
namespace ts {

/**
 * A base class for sketch algorithms.
 *
 * @tparam SketchType the type returned by the compute() function.
 * @tparam KmerInput is true when the hash algorithm first should convert the sequence to kmers.
 * */
template <typename SketchType, bool KmerInput>
class SketchBase {
  public:
    // The type that the compute function returns.
    using sketch_type = SketchType;

    // Whether the compute function takes a list of kmers.
    constexpr static bool kmer_input = KmerInput;

    // Whether transformations should be applied to the sketch output of this algorithm.
    constexpr static bool transform_sketches = false;

    // The name of the sketching algorithm.
    const std::string name;

    // If  kmer_input=true, value of kmer_size
    size_t kmer_size = 1;

    explicit SketchBase(std::string name, size_t kmer_size = 1) : name(std::move(name)), kmer_size(kmer_size) {}

    // Must be overridden by implementations.
    // Calling it will initialize the random hashes, overwriting any previous hash functions.
    // May be called multiple times on the same object to reset the state before running it on a new
    // set of sequences.
    void init() { static_assert(!sizeof(SketchType *), "Sketch type should implement init()."); }
};

} // namespace ts
