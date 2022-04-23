#include "sketch/tensor_slide.hpp"

namespace ts {

/**
 * A wrapper class around TensorSlide that flattens the output of TensorSlide::compute() from a 2D
 * vector to a 1D vector. The Flattener type must have a .flatten method that does the conversion.
 * Typically used with a class form sketch/dim_reduce.h.
 *
 * @tparam seq_type the type of elements in the sequences to be sketched.
 * @tparam Flattener: one of the classes from util/dim_reduce.h.
 */
template <class seq_type, class Flattener>
class TensorSlideFlat : public TensorSlide<seq_type> {
    Flattener flattener;

  public:
    using sketch_type = typename Flattener::sketch_type;

    /**
     * @param alphabet_size the number of elements in the alphabet S over which sequences are
     * defined (e.g. 4 for DNA)
     * @param sketch_dim the dimension of the embedded (sketched) space, denoted by D in the paper
     * @param tup_len the length of the subsequences considered for sketching, denoted by t
     * in the paper
     * @param win_len sliding sketches are computed for substrings of size win_len
     * @param stride sliding sketches are computed every stride characters
     * @param flattener the object to use to flatten the compute output.
     * @param seed the seed to initialize the random number generator used for the random hash
     * functions.
     */
    TensorSlideFlat(seq_type alphabet_size,
                    size_t sketch_dim,
                    size_t tup_len,
                    size_t win_len,
                    size_t stride,
                    Flattener flattener,
                    uint32_t seed,
                    const std::string &name = "TSS")
        : TensorSlide<seq_type>(alphabet_size, sketch_dim, tup_len, win_len, stride, seed, name),
          flattener(flattener) {}


    /**
     * Computes sliding sketches for the given sequence.
     * A sketch is computed every #stride characters on substrings of length #window.
     * @return seq.size()/stride sketches of size #sketch_dim
     */
    sketch_type compute(const std::vector<seq_type> &seq) {
        return flattener.flatten(TensorSlide<seq_type>::compute(seq));
    }

    static double dist(const sketch_type &a, const sketch_type &b) { return Flattener::dist(a, b); }
};

} // namespace ts
