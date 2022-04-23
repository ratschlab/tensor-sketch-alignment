#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace ts {

template <class T>
class transformer {
  public:
    virtual T transform(T val) const = 0;
};

template <class T>
class discretize : public transformer<T> {
  public:
    explicit discretize(size_t num_bins) : num_bins(num_bins) {
        bins = std::vector<double>(num_bins);
        for (size_t b = 0; b < num_bins; b++) {
            bins[b] = std::tan(M_PI * (((double)b + .5) / num_bins - .5));
        }
        bins.push_back(std::numeric_limits<double>::max());
        bins.insert(bins.begin(), -std::numeric_limits<double>::max());
    }


    // bin edges used to discretize the sketch output
    std::vector<double> bins;

    T transform(T val) const override {
        return std::upper_bound(bins.begin(), bins.end(), val) - bins.begin();
    }

  private:
    /** number of bins used to discretize the output*/
    size_t num_bins;
};

template <class T>
class atan_scaler : public transformer<T> {
  public:
    T transform(T val) const override { return atan(val); }
};

} // namespace ts
