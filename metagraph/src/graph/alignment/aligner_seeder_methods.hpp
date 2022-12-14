#ifndef __ALIGNER_SEEDER_METHODS_HPP__
#define __ALIGNER_SEEDER_METHODS_HPP__

#include "alignment.hpp"
#include "common/vectors/bitmap.hpp"
#include <unordered_map>
#include "sequence/alphabets.hpp"
#include "sequence/fasta_io.hpp"
#include "sketch/edit_distance.hpp"
#include "sketch/hash_base.hpp"
#include "sketch/hash_min.hpp"
#include "sketch/hash_ordered.hpp"
#include "sketch/hash_weighted.hpp"
#include "sketch/tensor.hpp"
#include "sketch/tensor_block.hpp"
#include "sketch/tensor_embedding.hpp"
#include "sketch/tensor_slide.hpp"
#include <boost/multiprecision/cpp_int.hpp>

namespace mtg {
namespace graph {
namespace align {

class ISeeder {
  public:
    virtual ~ISeeder() {}

    virtual const DBGAlignerConfig& get_config() const = 0;
    virtual std::vector<Seed> get_seeds() const = 0;
    virtual size_t get_num_matches() const = 0;

    virtual std::vector<Alignment> get_alignments() const {
        std::vector<Alignment> alignments;
        std::vector<Seed> seeds = get_seeds();
        alignments.reserve(seeds.size());
        for (const Seed &seed : seeds) {
            alignments.emplace_back(seed, get_config());
            alignments.back().trim_offset();
        }
        return alignments;
    }
};


class ManualMatchingSeeder : public ISeeder {
  public:
    ManualMatchingSeeder(std::vector<Seed>&& seeds,
                         size_t num_matching,
                         const DBGAlignerConfig &config)
          : config_(config), seeds_(std::move(seeds)), num_matching_(num_matching) {}

    virtual ~ManualMatchingSeeder() {}

    std::vector<Seed> get_seeds() const override { return seeds_; }
    const DBGAlignerConfig& get_config() const override { return config_; }
    size_t get_num_matches() const override final { return num_matching_; }
    std::vector<Seed>& data() { return seeds_; }

  private:
    const DBGAlignerConfig &config_;
    std::vector<Seed> seeds_;
    size_t num_matching_;
};

class ManualSeeder : public ISeeder {
  public:
    ManualSeeder(std::vector<Alignment>&& seeds = {}, size_t num_matching = 0)
        : seeds_(std::move(seeds)), num_matching_(num_matching) {}

    virtual ~ManualSeeder() {}

    std::vector<Seed> get_seeds() const override {
        throw std::runtime_error("Not implemented");
    }

    const DBGAlignerConfig& get_config() const override {
        throw std::runtime_error("Not implemented");
    }

    std::vector<Alignment> get_alignments() const override { return seeds_; }
    size_t get_num_matches() const override final { return num_matching_; }

    std::vector<Alignment>& data() { return seeds_; }

  private:
    std::vector<Alignment> seeds_;
    size_t num_matching_;
};

class ExactSeeder : public ISeeder {
  public:
    typedef DeBruijnGraph::node_index node_index;

    ExactSeeder(const DeBruijnGraph &graph,
                std::string_view query,
                bool orientation,
                std::vector<node_index>&& nodes,
                const DBGAlignerConfig &config);

    virtual ~ExactSeeder() {}

    const DBGAlignerConfig& get_config() const override { return config_; }
    std::vector<Seed> get_seeds() const override;
    size_t get_num_matches() const override final { return num_matching_; }

  protected:
    const DeBruijnGraph &graph_;
    std::string_view query_;
    bool orientation_;
    std::vector<node_index> query_nodes_;
    const DBGAlignerConfig &config_;
    size_t num_matching_;

    size_t num_exact_matching() const;
};

class MEMSeeder : public ExactSeeder {
  public:
    template <typename... Args>
    MEMSeeder(Args&&... args) : ExactSeeder(std::forward<Args>(args)...) {}

    virtual ~MEMSeeder() {}

    std::vector<Seed> get_seeds() const override;

    virtual const bitmap& get_mem_terminator() const = 0;
};

class SketchSeeder : public ISeeder {
    public:
        typedef DeBruijnGraph::node_index node_index;

        SketchSeeder(const DeBruijnGraph &graph,
                     std::string_view query,
                     bool orientation,
                     std::vector<node_index>&& nodes,
                     const DBGAlignerConfig &config);

        virtual ~SketchSeeder() {}

        const DBGAlignerConfig& get_config() const override { return config_; }
        std::vector<Seed> get_seeds() const override;
        std::vector<Alignment> get_alignments() const override;
        size_t get_num_matches() const override final { return 0; } // TODO: Implement this

    protected:
        std::unordered_map<int, std::string> sketches;
        const DeBruijnGraph &graph_;
        std::string_view query_;
        bool orientation_;
        std::vector<node_index> query_nodes_;
        const DBGAlignerConfig &config_;
};

class UniMEMSeeder : public MEMSeeder {
  public:
    template <typename... Args>
    UniMEMSeeder(Args&&... args)
          : MEMSeeder(std::forward<Args>(args)...),
            is_mem_terminus_([&](auto i) {
                                 return graph_.has_multiple_outgoing(i)
                                     || !graph_.has_single_incoming(i);
                             },
                             graph_.max_index() + 1) {
        assert(is_mem_terminus_.size() == graph_.max_index() + 1);
    }

    virtual ~UniMEMSeeder() {}

    const bitmap& get_mem_terminator() const override { return is_mem_terminus_; }

  private:
    bitmap_lazy is_mem_terminus_;
};

template <class BaseSeeder>
class SuffixSeeder : public BaseSeeder {
  public:
    template <typename... Args>
    SuffixSeeder(Args&&... args) : BaseSeeder(std::forward<Args>(args)...) {
        generate_seeds();
    }

    virtual ~SuffixSeeder() {}

    std::vector<Seed> get_seeds() const override { return seeds_; }

  protected:
    void generate_seeds();

  private:
    std::vector<Seed> seeds_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_SEEDER_METHODS_HPP__
