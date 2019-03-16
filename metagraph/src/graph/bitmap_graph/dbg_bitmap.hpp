#ifndef __DBG_BITMAP_HPP__
#define __DBG_BITMAP_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"
#include "bit_vector.hpp"


class DBGSDConstructor;

/**
 * Node-centric de Bruijn graph
 *
 * In canonical mode, for each k-mer in the graph, its
 * reverse complement is stored in the graph as well.
 */
class DBGSD : public DeBruijnGraph {
    friend DBGSDConstructor;

  public:
    // Initialize complete graph
    explicit DBGSD(size_t k, bool canonical_mode = false);

    // Initialize graph from builder
    explicit DBGSD(DBGSDConstructor *builder);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const;
    void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const;

    size_t outdegree(node_index) const {
        // TODO: Complete outdegree for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    size_t indegree(node_index) const {
        // TODO: Complete outdegree for DBGSD.
        throw std::runtime_error("Not implemented");
    }

    template <class... T>
    using Call = typename std::function<void(T...)>;

    // traverse all nodes in graph
    void call_kmers(Call<node_index, const std::string&> callback) const;

    // call paths (or simple paths if |split_to_contigs| is true) that cover
    // exactly all kmers in graph
    void call_sequences(Call<const std::string&> callback,
                        bool split_to_contigs = false) const;

    node_index kmer_to_node(const std::string &kmer) const;

    std::string get_node_sequence(node_index node) const;

    inline size_t get_k() const { return k_; }
    inline bool is_canonical_mode() const { return canonical_mode_; }

    uint64_t num_nodes() const;


    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    bool operator==(const DBGSD &other) const { return equals(other, false); }
    bool operator!=(const DBGSD &other) const { return !(*this == other); }

    bool equals(const DBGSD &other, bool verbose = false) const;

    friend std::ostream& operator<<(std::ostream &out, const DBGSD &graph);

    const std::string &alphabet;

    using Chunk = bit_vector_sd;

  private:
    using Kmer = KmerExtractor2Bit::Kmer64;

    void add_sequence(const std::string &,
                      bit_vector_dyn *) {
        throw std::runtime_error("Not implemented");
    }

    Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                   bool canonical = false) const;

    uint64_t node_to_index(node_index node) const;
    Kmer node_to_kmer(node_index node) const;
    node_index to_node(const Kmer &kmer) const;

    void call_paths(Call<const std::vector<node_index>,
                         const std::vector<uint8_t>&> callback,
                    bool split_to_contigs) const;

    size_t k_;
    bool canonical_mode_;
    KmerExtractor2Bit seq_encoder_;

    bit_vector_sd kmers_;

    static constexpr auto kExtension = ".bitmapdbg";
};

#endif // __DBG_BITMAP_HPP__