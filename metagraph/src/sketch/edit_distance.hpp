#pragma once

#include "sketch/sketch_base.hpp"
#include "util/multivec.hpp"
#include "util/timer.hpp"
#include "util/utils.hpp"

#include <cassert>
#include <cmath>
#include <random>

namespace ts { // ts = Tensor Sketch

template <class seq_type>
class EditDistance : public SketchBase<const std::vector<seq_type> *, false> {
  public:
    explicit EditDistance(const std::string &name = "ED")
        : SketchBase<const std::vector<seq_type> *, false>(name) {
        init();
    }

    void init() {}

    const std::vector<seq_type> *compute(const std::vector<seq_type> &seq) { return &seq; }

    static double dist(const std::vector<seq_type> *a, const std::vector<seq_type> *b) {
        Timer timer("edit_dist");
        return edit_distance(*a, *b);
    }
};

} // namespace ts
