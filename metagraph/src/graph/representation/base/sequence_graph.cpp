#include "sequence_graph.hpp"

#include <cassert>
#include <progress_bar.hpp>
#include <sdsl/int_vector.hpp>

#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "graph/representation/canonical_dbg.hpp"
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
#include <queue>

namespace mtg {
namespace graph {
using namespace boost::multiprecision;
typedef DeBruijnGraph::node_index node_index;

static const uint64_t kBlockSize = 1 << 14;
static_assert(!(kBlockSize & 0xFF));


/*************** SequenceGraph ***************/

void SequenceGraph::call_nodes(const std::function<void(node_index)> &callback,
                               const std::function<bool()> &stop_early) const {
    assert(num_nodes() == max_index());

    const auto nnodes = num_nodes();
    for (node_index i = 1; i <= nnodes && !stop_early(); ++i) {
        callback(i);
    }
}

void SequenceGraph::add_extension(std::shared_ptr<GraphExtension> extension) {
    assert(extension.get());
    extensions_.push_back(extension);
}

void SequenceGraph::serialize_extensions(const std::string &filename) const {
    for (const auto &extension : extensions_) {
        extension->serialize(utils::make_suffix(filename, file_extension()));
    }
}

/*************** DeBruijnGraph ***************/

node_index DeBruijnGraph::kmer_to_node(std::string_view kmer) const {
    assert(kmer.size() == get_k());

    node_index node = npos;
    map_to_nodes_sequentially(kmer, [&node](node_index i) { node = i; });
    return node;
}

// Check whether graph contains fraction of nodes from the sequence
bool DeBruijnGraph::find(std::string_view sequence,
                         double discovery_fraction) const {
    if (sequence.length() < get_k())
        return false;

    const size_t num_kmers = sequence.length() - get_k() + 1;
    const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
    const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
    size_t num_kmers_discovered = 0;
    size_t num_kmers_missing = 0;

    map_to_nodes(sequence,
        [&](node_index node) {
            if (node) {
                num_kmers_discovered++;
            } else {
                num_kmers_missing++;
            }
        },
        [&]() { return num_kmers_missing > max_kmers_missing
                        || num_kmers_discovered >= min_kmers_discovered; }
    );

    return num_kmers_missing <= max_kmers_missing;
}

bool DeBruijnGraph::operator==(const DeBruijnGraph &) const {
    throw std::runtime_error("Not implemented");
    return false;
}

void DeBruijnGraph::traverse(node_index start,
                             const char *begin,
                             const char *end,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    assert(start >= 1 && start <= max_index());
    assert(end >= begin);
    if (terminate())
        return;

    for (; begin != end && !terminate(); ++begin) {
        start = traverse(start, *begin);

        if (start == npos)
            return;

        callback(start);
    }
}

void call_sequences_from(const DeBruijnGraph &graph,
                         node_index start,
                         const DeBruijnGraph::CallPath &callback,
                         sdsl::bit_vector *visited,
                         sdsl::bit_vector *discovered,
                         ProgressBar &progress_bar,
                         bool call_unitigs,
                         uint64_t min_tip_size,
                         bool kmers_in_single_form) {
    assert(start >= 1 && start <= graph.max_index());
    assert((min_tip_size <= 1 || call_unitigs)
                && "tip pruning works only for unitig extraction");
    assert(visited);
    assert(discovered);
    assert(!(*visited)[start]);

    std::vector<node_index> queue = { start };
    (*discovered)[start] = true;

    std::vector<node_index> path;
    std::string sequence;

    std::vector<std::pair<node_index, char>> targets;

    // keep traversing until we have worked off all branches from the queue
    while (queue.size()) {
        node_index node = queue.back();
        queue.pop_back();
        assert((*discovered)[node]);
        if ((*visited)[node])
            continue;

        path.resize(0);
        path.push_back(node);
        sequence = graph.get_node_sequence(node);
        start = node;

        // traverse simple path until we reach its tail or
        // the first edge that has been already visited
        while (true) {
            assert(node);
            assert((*discovered)[node]);
            assert(!(*visited)[node]);
            assert(sequence.length() >= graph.get_k());
            (*visited)[node] = true;
            ++progress_bar;

            targets.clear();
            graph.call_outgoing_kmers(node,
                [&](node_index next, char c) { targets.emplace_back(next, c); }
            );

            if (targets.empty())
                break;

            // in call_unitigs mode, all nodes with multiple incoming
            // edges are marked as discovered
            assert(!call_unitigs || graph.has_single_incoming(targets.front().first)
                                 || (*discovered)[targets.front().first]);
            if (targets.size() == 1) {
                const auto &[next, c] = targets[0];

                if ((*visited)[next])
                    break;

                if (!call_unitigs || !(*discovered)[next]) {
                    sequence.push_back(c);
                    path.push_back(next);
                    (*discovered)[next] = true;
                    node = next;
                    continue;
                }
                assert(call_unitigs && graph.indegree(next) >= 2);
                break;
            }

            node_index next_node = DeBruijnGraph::npos;
            //  _____.___
            //      \.___
            for (const auto &[next, c] : targets) {
                if (next_node == DeBruijnGraph::npos
                        && !call_unitigs
                        && !(*visited)[next]) {
                    (*discovered)[next] = true;
                    next_node = next;
                    sequence.push_back(c);
                    path.push_back(next);
                } else if (!(*discovered)[next]) {
                    (*discovered)[next] = true;
                    queue.push_back(next);
                }
            }

            if (next_node == DeBruijnGraph::npos)
                break;

            node = next_node;
        }

        assert(sequence.size() >= graph.get_k());
        assert(path.size());

        if (!call_unitigs
                  // check if long
                  || sequence.size() >= graph.get_k() + min_tip_size - 1
                  // check if not tip
                  || graph.indegree(start) + graph.outdegree(node) >= 2) {

            assert(path == map_to_nodes_sequentially(graph, sequence));

            if (kmers_in_single_form) {
                // get dual path (mapping of the reverse complement sequence)
                std::string rev_comp_seq = sequence;
                reverse_complement(rev_comp_seq.begin(), rev_comp_seq.end());

                auto dual_path = map_to_nodes_sequentially(graph, rev_comp_seq);
                std::reverse(dual_path.begin(), dual_path.end());
                size_t begin = 0;

                assert(std::all_of(path.begin(), path.end(),
                                   [&](auto i) { return (*visited)[i] && (*discovered)[i]; }));

                progress_bar += std::count_if(dual_path.begin(), dual_path.end(),
                                              [&](auto node) { return node && !(*visited)[node]; });

                // Mark all nodes in path as unvisited and re-visit them while
                // traversing the path (iterating through all nodes).
                std::for_each(path.begin(), path.end(),
                              [&](auto i) { (*visited)[i] = false; });

                // traverse the path with its dual and visit the nodes
                for (size_t i = 0; i < path.size(); ++i) {
                    assert(path[i]);
                    (*visited)[path[i]] = true;

                    if (!dual_path[i])
                        continue;

                    // check if reverse-complement k-mer has been traversed
                    if (!(*visited)[dual_path[i]] || dual_path[i] == path[i]) {
                        (*visited)[dual_path[i]] = (*discovered)[dual_path[i]] = true;
                        continue;
                    }

                    // The reverse-complement k-mer has been visited
                    // -> Skip this k-mer and call the traversed path segment.
                    if (begin < i)
                        callback(sequence.substr(begin, i + graph.get_k() - 1 - begin),
                                 { path.begin() + begin, path.begin() + i });

                    begin = i + 1;
                }

                // Call the path traversed
                if (!begin) {
                    callback(std::move(sequence), std::move(path));

                } else if (begin < path.size()) {
                    callback(sequence.substr(begin),
                             { path.begin() + begin, path.end() });
                }

                for (size_t i = 0; i < path.size(); ++i) {
                    if (!dual_path[i])
                        continue;

                    // schedule traversal branched off from each terminal dual k-mer
                    if (i == 0) {
                        graph.adjacent_outgoing_nodes(dual_path[i],
                            [&](auto next) {
                                if (!(*discovered)[next]) {
                                    (*discovered)[next] = true;
                                    queue.push_back(next);
                                }
                            }
                        );
                    } else if (!dual_path[i - 1]) {
                        size_t num_outgoing = 0;
                        graph.adjacent_outgoing_nodes(dual_path[i],
                            [&](node_index next) { num_outgoing++; node = next; }
                        );
                        // schedule only if it has a single outgoing k-mer,
                        // otherwise it's a fork which will be covered in the forward pass.
                        if (num_outgoing == 1) {
                            if (!(*discovered)[node]) {
                                (*discovered)[node] = true;
                                queue.push_back(node);
                            }
                        }
                    }
                }

            } else {
                callback(std::move(sequence), std::move(path));
            }
        }
    }
}

void call_sequences(const DeBruijnGraph &graph,
                    const DeBruijnGraph::CallPath &callback,
                    size_t num_threads,
                    bool call_unitigs,
                    uint64_t min_tip_size,
                    bool kmers_in_single_form) {
    // TODO: port over the implementation from BOSS
    std::ignore = num_threads;

    sdsl::bit_vector discovered(graph.max_index() + 1, true);
    graph.call_nodes([&](auto node) { discovered[node] = false; });
    sdsl::bit_vector visited = discovered;

    ProgressBar progress_bar(visited.size() - sdsl::util::cnt_one_bits(visited),
                             "Traverse graph",
                             std::cerr, !common::get_verbose());

    auto call_paths_from = [&](node_index node) {
        call_sequences_from(graph,
                            node,
                            callback,
                            &visited,
                            &discovered,
                            progress_bar,
                            call_unitigs,
                            min_tip_size,
                            kmers_in_single_form);
    };

    if (call_unitigs) {
        // mark all source and merge nodes (those with indegree 0 or >1)
        //  .____  or  .____  or  ____.___
        //              \___      ___/
        //
        #pragma omp parallel for num_threads(get_num_threads())
        for (uint64_t begin = 0; begin < discovered.size(); begin += kBlockSize) {
            call_zeros(discovered,
                begin,
                std::min(begin + kBlockSize, discovered.size()),
                [&](auto i) { discovered[i] = !graph.has_single_incoming(i); }
            );
        }

        // now traverse graph starting at these nodes
        call_zeros(visited, [&](auto node) {
            assert(discovered[node] == !graph.has_single_incoming(node));
            if (discovered[node])
                call_paths_from(node);
        });

    } else {
        // start at the source nodes (those with indegree == 0)
        //  .____  or  .____
        //              \___
        //
        call_zeros(visited, [&](auto node) {
            if (graph.has_no_incoming(node) && !visited[node])
                call_paths_from(node);
        });
    }

    // then forks
    //  ____.____
    //       \___
    //
    call_zeros(visited, [&](auto node) {
        // TODO: these two calls to outgoing nodes could be combined into one
        if (graph.has_multiple_outgoing(node)) {
            graph.adjacent_outgoing_nodes(node, [&](auto next) {
                if (!visited[next])
                    call_paths_from(next);
            });
        }
    });

    // then the rest (loops)
    call_zeros(visited, call_paths_from);
}

void DeBruijnGraph::call_sequences(const CallPath &callback,
                                   size_t num_threads,
                                   bool kmers_in_single_form) const {
    ::mtg::graph::call_sequences(*this, callback, num_threads, false, 0, kmers_in_single_form);
}

void DeBruijnGraph::call_unitigs(const CallPath &callback,
                                 size_t num_threads,
                                 size_t min_tip_size,
                                 bool kmers_in_single_form) const {
    ::mtg::graph::call_sequences(*this, callback, num_threads, true, min_tip_size, kmers_in_single_form);
}

/**
 * Traverse graph and iterate over all nodes
 */
void DeBruijnGraph
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    sdsl::bit_vector visited(max_index() + 1, false);
    std::stack<std::pair<node_index, std::string>> nodes;

    call_nodes([&](node_index i) {
        if (visited[i])
            return;

        nodes.emplace(i, get_node_sequence(i));
        while (nodes.size()) {
            // FYI: structured binding is a new thing that often
            // leads to issues, so avoid using it.
            // https://bit.ly/2JI6oci
            auto [node, sequence] = std::move(nodes.top());
            nodes.pop();

            while (!visited[node]) {
                visited[node] = true;
                callback(node, sequence);
                sequence.assign(sequence.begin() + 1, sequence.end());

                auto next_node = npos;
                char next_c = '\0';
                call_outgoing_kmers(
                    node,
                    [&, &sequence=sequence](const auto &next, char c) {
                        if (visited[next])
                            return;

                        if (next_node == npos) {
                            std::tie(next_node, next_c) =  std::tie(next, c);
                        } else {
                            nodes.emplace(next, sequence + c);
                        }
                    }
                );

                if (next_node == npos)
                    break;

                node = next_node;
                sequence.push_back(next_c);
            }
        }
    });
}

std::unordered_map<uint64_t, uint64_t> DeBruijnGraph::get_forward_tmers(node_index start, uint32_t t, uint32_t n) const {
    std::unordered_map<uint64_t, uint64_t> output;
    output.reserve(n);
    // TODO: Fix this nasty implementation
    uint32_t kmer_stream = 0;
    uint32_t cleanup = pow(2, 2 * t) - 1;
    node_index node = start;

    // I am k-1 nodes back
    // ASCZXCZXC - node
    // ASCZXCZXC - seed

    // QQQQQQASCZXCZXC - node
    // DZXCZXDDD
    //         DZXCZXDDD
    // Q = k - 1 nodes back
    // QQQQQQQQASCZXCZXC
    // ^ start
    // ASCZXCZXC - seed

    // Go one back
    adjacent_incoming_nodes(node, [&](node_index back_node){
        node = back_node;
    });

    // Go k back
//    for (uint32_t x = 0; x < get_k(); ++x) {
//        adjacent_incoming_nodes(node, [&](node_index back_node){
//            node = back_node;
//        });
//    }
    // QQQQQQQQQASCZXCZXC
    // ASCZXCZXC - seed

    char c;
    uint32_t delta = 0;
    for(uint32_t i = 0; i < n; ++i) { call_outgoing_kmers(node, [&](node_index next, char next_char){ node = next; c = next_char; }); kmer_stream <<= 2; kmer_stream |= ts::char2int(c); kmer_stream &= cleanup; if (i >= t - 1) {
            output[kmer_stream] = delta++;
        }
    }

    return output;
}
void DeBruijnGraph::compute_sketches(uint64_t kmer_word_size,
                                     size_t embed_dim,
                                     size_t tuple_length,
                                     uint32_t n_times_sketch,
                                     uint32_t num_threads,
				     const char* fname,
				     bool load_index) {

    if (load_index) {
    	//read index from file
	index = static_cast<faiss::IndexIDMap2*>(faiss::read_index(fname));
	return;
    }
    uint32_t k = get_k();
    uint32_t total_ram = sizeof(double) * (embed_dim + 1) * static_cast<uint32_t>(std::ceil((double)num_nodes() / (double)get_k()));

    double stride_ratio = 0.1;
    double window_ratio = 0.2;
    uint32_t stride = static_cast<uint32_t>(stride_ratio * (double)k);
    uint32_t window = static_cast<uint32_t>(window_ratio * (double)k);
    uint32_t num_windows = static_cast<uint32_t>(std::ceil((double)(k - window + 1) / (double)stride));
    index_ = (faiss::IndexHNSW*)index_factory(embed_dim * num_windows, "HNSW32", faiss::METRIC_L2);
    index_->hnsw.efSearch = 1000;
    /* index_->hnsw.search_bounded_queue = true; */
    index = new faiss::IndexIDMap2(index_);
    for (uint32_t n_repeat = 0; n_repeat < n_times_sketch; n_repeat++) {
        call_unitigs(
            [&](const std::string &s, const std::vector <uint64_t> &v) {
                std::vector<uint8_t> node_to_int;
                for (unsigned char c: s) {
                    node_to_int.push_back(ts::char2int(c));
                }
                // RAM Used = sizeof(double) * embed_dim * num_nodes / k
                //          + sizeof(uint64_t) * num_nodes / k
                // 64 * (embed_dim * num_nodes + num_nodes) = 64 * (embed_dim + 1)  * num_nodes / k

                uint32_t current_ram_usage = 0;

                uint32_t unitig_length = s.size();
                idx_t n = s.size() - k + 1;
                float* sketch_arr = new float[embed_dim * num_windows  * static_cast<uint32_t>(std::ceil((double)n / (double)k))];
                idx_t* positions = new idx_t[static_cast<uint32_t>(std::ceil((double)n / (double)k))];
                uint32_t pos = 0;
//                std::cout << "Initial kmer: " << get_node_sequence(v[0]) << std::endl;
                for(uint32_t kmer_start = k; kmer_start < unitig_length - k + 1; kmer_start += k) {
//                    std::cout << "Kmer: " << kmer_start << " --  " << get_node_sequence(v[kmer_start]) << std::endl;
//                    std::cout << v[kmer_start] << " " << v[kmer_start - (k - 1)] << std::endl;
                    // Get kmer
                    std::vector<uint8_t> kmer_string(node_to_int.begin() + kmer_start, node_to_int.begin() + kmer_start + k);

                    // Get full sketch
                    std::vector<double> full_sketch;
                    full_sketch.reserve(embed_dim * num_windows);

                    for(uint32_t mmer_start = 0; mmer_start < k - window + 1; mmer_start += stride) {
                        std::vector<uint8_t> mmer_string(kmer_string.begin() + mmer_start, kmer_string.begin() + mmer_start + window);
                        std::vector<double> sketch = ts::Tensor<uint8_t>(kmer_word_size, embed_dim, tuple_length, n_repeat).compute(mmer_string);
                        full_sketch.insert(full_sketch.end(), sketch.begin(), sketch.end());
                    }

                    positions[pos] = v[kmer_start - (k - 1)];
//                    positions[pos] = v[kmer_start];
//                    std::cout << "Adding: " << v[kmer_start] << " " <<get_node_sequence(v[kmer_start]) << " at: " << pos << std::endl;
//                    debugmap[v[kmer_start - (k - 1)]] = v[kmer_start];
                    for(uint32_t j = 0; j < embed_dim * num_windows; ++j) {
                        sketch_arr[pos * embed_dim * num_windows + j] = full_sketch[j];
                    }
                    pos++;
                    current_ram_usage += sizeof(double) * (embed_dim + 1);
                }
                index->train(pos, sketch_arr);
                index->add_with_ids(pos, sketch_arr, positions);
                // debug
//                float* result;
//                index->index->reconstruct(71152, result);
//
//                for(int q = 0; q < embed_dim * num_windows; ++q)
//                    std::cout << result[q] << " ";
//                std::cout << std::endl;
            });
        
    }
    faiss::write_index(index, fname);
    
}

void DeBruijnGraph::print(std::ostream &out) const {
    std::string vertex_header("Vertex");
    vertex_header.resize(get_k(), ' ');

    out << "Index"
        << "\t" << vertex_header
        << std::endl;

    call_nodes([&](node_index i) {
        out << i
            << "\t" << "Node seq:" << get_node_sequence(i) << "\t" << "Node:" << i
            << std::endl;
    });
}

void DeBruijnGraph
::call_source_nodes(const std::function<void(node_index)> &callback) const {
    call_nodes([&](node_index i) {
        if (has_no_incoming(i))
            callback(i);
    });
}


std::ostream& operator<<(std::ostream &out, const DeBruijnGraph &graph) {
    graph.print(out);
    return out;
}

// returns the edge rank, starting from zero
size_t incoming_edge_rank(const SequenceGraph &graph,
                          SequenceGraph::node_index source,
                          SequenceGraph::node_index target) {
    assert(source >= 1 && source <= graph.max_index());
    assert(target >= 1 && target <= graph.max_index());

    size_t edge_rank = 0;
    bool done = false;

    graph.adjacent_incoming_nodes(target, [&](auto node) {
        if (node == source)
            done = true;

        if (!done)
            edge_rank++;
    });

    assert(done && "the edge must exist in graph");

    return edge_rank;
}

std::vector<SequenceGraph::node_index>
map_to_nodes_sequentially(const SequenceGraph &graph, std::string_view sequence) {
    std::vector<SequenceGraph::node_index> nodes;
    nodes.reserve(sequence.size());

    graph.map_to_nodes_sequentially(sequence,
        [&nodes](auto node) { nodes.push_back(node); }
    );

    return nodes;
}

std::vector<SequenceGraph::node_index>
map_to_nodes(const SequenceGraph &graph, std::string_view sequence) {
    std::vector<node_index> nodes;
    nodes.reserve(sequence.size());

    graph.map_to_nodes(sequence, [&](node_index i) { nodes.push_back(i); });

    return nodes;
}

void reverse_complement_seq_path(const SequenceGraph &graph,
                                 std::string &seq,
                                 std::vector<SequenceGraph::node_index> &path) {
    if (const auto *canonical_dbg = dynamic_cast<const CanonicalDBG*>(&graph)) {
        canonical_dbg->reverse_complement(seq, path);
        return;
    }

    reverse_complement(seq.begin(), seq.end());
    path = map_to_nodes_sequentially(graph, seq);
}

} // namespace graph
} // namespace mtg
