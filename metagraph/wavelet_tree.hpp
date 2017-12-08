#ifndef __WAVELET_TREE_HPP__
#define __WAVELET_TREE_HPP__

#include <assert.h>

#define __STDC_FORMAT_MACROS 1
#include <cinttypes>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>


class wavelet_tree {
  public:
    virtual ~wavelet_tree() {};

    virtual uint64_t rank(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t size() const = 0;
    virtual void set(uint64_t id, uint64_t val) = 0;
    virtual uint64_t operator[](uint64_t id) const = 0;
    virtual void insert(uint64_t id, uint64_t val) = 0;
    virtual void remove(uint64_t id) = 0;
    virtual bool deserialise(std::istream &in) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual std::vector<uint64_t> to_vector() const = 0;

    friend inline std::ostream& operator<<(std::ostream &os,
                                           const wavelet_tree &wt);

  protected:
    virtual void print(std::ostream &os) const = 0;
};

inline std::ostream& operator<<(std::ostream &os, const wavelet_tree& wt) {
    wt.print(os);
    return os;
}


class wavelet_tree_stat : public wavelet_tree {
  public:
    explicit wavelet_tree_stat(uint64_t logsigma, uint64_t size = 0, uint64_t value = 0)
            : int_vector_(2 * size + 1, value, 1 << logsigma), n_(size) {}

    template <class Vector>
    wavelet_tree_stat(uint64_t logsigma, const Vector &vector)
            : int_vector_(vector.size(), 0, logsigma), n_(vector.size()) {
        for (uint64_t i = 0; i < n_; ++i) {
            int_vector_[i] = vector[i];
        }
    }

    uint64_t size() const { return n_; }

    bool deserialise(std::istream &in) {
        int_vector_.load(in);
        wwt_.load(in);
        n_ = int_vector_.size();
        // we assume the wwt_ has been build before serialization
        requires_update_ = false;
        return true;
    }

    void serialise(std::ostream &out) const {
        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();
        int_vector_.serialize(out);
        wwt_.serialize(out);
    }

    void insert(uint64_t id, uint64_t val) {
        if (n_ == size()) {
            int_vector_.resize(2 * n_ + 1);
        }
        n_++;
        if (size() > 1)
            std::copy_backward(int_vector_.begin() + id,
                               int_vector_.begin() + n_ - 1,
                               int_vector_.begin() + n_);
        int_vector_[id] = val;
        requires_update_ = true;
    }

    void remove(uint64_t id) {
        if (this->size() > 1)
            std::copy(int_vector_.begin() + id + 1,
                      int_vector_.begin() + n_,
                      int_vector_.begin() + id);
        n_--;
        requires_update_ = true;
    }

    uint64_t rank(uint64_t c, uint64_t i) const {
        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();

        return wwt_.rank(std::min(i + 1, size()), c);
    }

    uint64_t select(uint64_t c, uint64_t i) const {
        assert(i > 0 && size() > 0);

        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();

        assert(i <= rank(c, size() - 1));
        return wwt_.select(i, c);
    }

    uint64_t operator[](uint64_t id) const {
        assert(id < size());
        return int_vector_[id];
    }

    std::vector<uint64_t> to_vector() const {
        std::vector<uint64_t> result(size());
        for (uint64_t i = 0; i < size(); ++i) {
            result[i] = int_vector_[i];
        }
        return result;
    }

    void set(uint64_t id, uint64_t val) {
        requires_update_ = true;
        int_vector_[id] = val;
    }
    

  protected:
    void print(std::ostream& os) const {
        os << wwt_;
    }

  private:
    sdsl::int_vector<> int_vector_;
    sdsl::wt_int<> wwt_;
    bool requires_update_ = true;
    uint64_t n_;

    void init_wt() {
        // std::cout << "Initializing WT index" << std::endl;
        int_vector_.resize(n_);
        sdsl::construct_im(wwt_, int_vector_);
        requires_update_ = false;
    }
};


class wavelet_tree_dyn : public wavelet_tree {
  public:
    explicit wavelet_tree_dyn(uint64_t b) : wavelet_tree_(b) {}

    template <class Vector>
    wavelet_tree_dyn(uint64_t b, const Vector &W_stat, unsigned int parallel = 1)
        : wavelet_tree_(initialize_tree(W_stat, b, parallel), b, W_stat.size()) {}

    wavelet_tree_dyn(uint64_t b, libmaus2::bitbtree::BitBTree<6, 64> *bt, uint64_t n)
        : wavelet_tree_(bt, b, n) {}

    uint64_t size() const {
        return wavelet_tree_.size();
    }

    bool deserialise(std::istream &in) {
        wavelet_tree_.R =
            libmaus2::wavelet::DynamicWaveletTree<6, 64>::loadBitBTree(in);
        const_cast<uint64_t&>(wavelet_tree_.b) =
            libmaus2::util::NumberSerialisation::deserialiseNumber(in);
        wavelet_tree_.n =
            libmaus2::util::NumberSerialisation::deserialiseNumber(in);
        return true;
    }

    void serialise(std::ostream &out) const {
        wavelet_tree_.serialise(out);
    }

    void insert(uint64_t id, uint64_t val) {
        wavelet_tree_.insert(val, id);
    }

    void remove(uint64_t id) {
        wavelet_tree_.remove(id);
    }

    uint64_t rank(uint64_t c, uint64_t i) const {
        return size() > 0
                ? wavelet_tree_.rank(c, std::min(i, size() - 1))
                : 0;
    }

    uint64_t select(uint64_t c, uint64_t i) const {
        assert(i > 0 && size() > 0 && i <= rank(c, size() - 1));
        return wavelet_tree_.select(c, i - 1);
    }

    uint64_t operator[](uint64_t id) const {
        assert(id < size());
        return wavelet_tree_[id];
    }

    void set(uint64_t id, uint64_t val) {
        remove(id);
        insert(id, val);
    }

    std::vector<uint64_t> to_vector() const {
        std::vector<uint64_t> W_stat;

        size_t n = size();

        const size_t b = 4;

        std::vector<uint64_t> offsets(1ull << b, 0);
        std::queue<uint64_t> blocks;
        std::queue<uint64_t> new_blocks;

        blocks.push(n);

        size_t pos = 0;
        uint64_t o = 0;

        for (size_t ib = 0; ib < b; ++ib) {
            while (!blocks.empty()) {
                uint64_t cnt = blocks.front();
                blocks.pop();
                uint64_t epos = pos + cnt;
                for ( ; pos < epos; ++pos) {
                    offsets.at(o) += !get_bit_raw(pos);
                }
                if (ib < b - 1) {
                    new_blocks.push(offsets.at(o));
                    new_blocks.push(cnt - offsets.at(o));
                }
                o++;
            }
            if (ib < b - 1)
                blocks.swap(new_blocks);
        }

        //std::cerr << "R size: " << R->size() << std::endl;

        bool bit;
        std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);
        uint64_t p, co, v, m;
        for (size_t i = 0; i < n; ++i) {
            m = 1ull << (b - 1);
            //v = (uint64_t) W_stat.at(i);
            v = 0;
            o = 0;
            p = i;
            co = 0;
            for (size_t ib = 0; ib < b - 1; ++ib) {
                bit = get_bit_raw(ib * n + p + co);
                if (bit) {
                    v |= m;
                    co += offsets.at(o);
                    p -= upto_offsets.at(o);
                } else {
                    p -= (p - upto_offsets.at(o));
                    upto_offsets.at(o) += 1;
                }
                o = 2 * o + 1 + bit;
                m >>= 1;
            }
            bit = get_bit_raw((b - 1) * n + p + co);
            if (bit) {
                v |= m;
            }
            W_stat.push_back(v);
        }
        return W_stat;
    }

  
  protected:
    void print(std::ostream& os) const {
        os << wavelet_tree_;
    }

  private:
    libmaus2::wavelet::DynamicWaveletTree<6, 64> wavelet_tree_;

    bool get_bit_raw(uint64_t id) const {
        return wavelet_tree_.R->operator[](id);
    }

    template <class Vector>
    libmaus2::bitbtree::BitBTree<6, 64>* initialize_tree(const Vector &W_stat,
                                                         const uint64_t b,
                                                         unsigned int parallel = 1) {
        // return dynamic_cast<libmaus2::bitbtree::BitBTree<6, 64>*>(
        //     makeTree(W_stat, b, parallel)), b, W_stat.size()
        // );

        uint64_t const n = W_stat.size();

        // compute total offsets for the individual bins
        std::vector<uint64_t> offsets((1ull << (b - 1)) - 1, 0);
        #pragma omp parallel num_threads(parallel)
        {
            uint64_t v, m, o;
            #pragma omp for
            for (uint64_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = static_cast<uint64_t>(W_stat[i]);
                o = 0;
                for (uint64_t ib = 1; ib < b; ++ib) {
                    bool const bit = m & v;
                    if (!bit) {
                        #pragma omp critical
                        offsets.at(o) += 1;
                    }
                    o = 2 * o + 1 + bit;
                    m >>= 1;
                }
            }
        }

        auto *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);
        std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);

        #pragma omp parallel num_threads(parallel)
        {
            uint64_t m, v, o, p, co;
            bool bit;
            #pragma omp for
            for (uint64_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat[i];
                o = 0;
                p = i;
                co = 0;
                for (uint64_t ib = 0; ib < b - 1; ++ib) {
                    bit = m & v;
                    if (bit) {
                        #pragma omp critical
                        tmp->setBitQuick(ib * n + p + co, true);
                        co += offsets.at(o);
                        p -= upto_offsets.at(o);
                    } else {
                        p -= (p - upto_offsets.at(o));
                        upto_offsets.at(o) += 1;
                    }
                    //std::cerr << "o: " << o << " offset[o]: " << offsets.at(o) << std::endl;
                    o = 2 * o + 1 + bit;
                    m >>= 1;
                }
                bit = m & v;
                if (bit) {
                    //std::cerr << "b - 1: " << b - 1 << " n: " << n << " p: " << p << " co: " << co << std::endl;
                    #pragma omp critical
                    tmp->setBitQuick((b - 1) * n + p + co, true);
                }
            }
        }

        return tmp;
    }
};


#endif // __WAVELET_TREE_HPP__