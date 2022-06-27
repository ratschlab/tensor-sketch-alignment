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

void DeBruijnGraph::compute_sketches(uint64_t kmer_word_size,
                                     size_t embed_dim,
                                     size_t tuple_length,
                                     size_t m_stride,
                                     uint32_t n_times_sketch,
                                     uint32_t minimizer_window,
                                     uint32_t num_threads) {
    sketch_maps = std::vector<std::unordered_map<key_type, std::vector<node_index>>>(n_times_sketch);
    map_backward = std::unordered_map<node_index, node_index>();
    random_directions = std::vector<std::vector<uint8_t>>(n_times_sketch);
    int ratio = 5;
    int m = (m_stride * get_k()) / ratio;
    std::vector<uint8_t> random_direction(embed_dim * (ratio - 1));
    uint32_t k = get_k();
    key_type discretized_sketch = 0;
    uint32_t delta = minimizer_window;
    call_unitigs(
            [&](const std::string &s, const std::vector <uint64_t> &v) {
                for (uint32_t n_repeat = 0; n_repeat < n_times_sketch; n_repeat++) {
                    // Generate random vector and save it
                    srand(n_repeat);
                    std::generate(random_direction.begin(), random_direction.end(),[]() { return rand() % 2; });
                    random_directions[n_repeat] = random_direction;

                    // Generate sequence of ints
//                    std::vector<uint8_t> node_sequence_to_int(s.size());
//                    for (uint32_t i = 0; i < s.size(); ++i) {
//                        node_sequence_to_int[i] = ts::char2int(s[i]);
//                    }

                    // Sketch
                    ts::TensorSlide <uint8_t> tensor = ts::TensorSlide<uint8_t>(kmer_word_size,
                                                                                embed_dim,
                                                                                tuple_length,
                                                                                m,
                                                                                1,
                                                                                n_repeat);
                    //std::vector <uint64_t> sanity_sketches = tensor.compute_discretized(node_sequence_to_int);
                    // Break a unitig into multiple pieces
                    std::vector<std::vector<uint64_t>> m_sketches_per_thread(num_threads);
                    #pragma omp parallel for num_threads(num_threads)
                    for(uint32_t thread = 0; thread < num_threads; ++thread) {
                        uint32_t start = thread * (s.size() / num_threads) + (thread > 0) * (-k + 1);
                        uint32_t end = (thread + 1) * (s.size() / num_threads) + (thread == num_threads - 1) * (s.size() % num_threads);

                        std::string my_piece_of_s = s.substr(start, end - start + 1);
                        std::vector<uint8_t> my_node_sequence_to_int(my_piece_of_s.size());
                        for (uint32_t i = 0; i < my_piece_of_s.size(); ++i) {
                            my_node_sequence_to_int[i] = ts::char2int(my_piece_of_s[i]);
                        }

                        std::vector <uint64_t> my_m_sketches = tensor.compute_discretized(my_node_sequence_to_int);
                        m_sketches_per_thread[thread] = my_m_sketches;
                    }
                    // Merge it all
                    std::vector<uint64_t> m_sketches(s.size() - m + 1);
                    for(uint32_t thread = 0; thread < num_threads; ++thread) {
                        if (thread == 0)
                            std::move(m_sketches_per_thread[thread].begin(), m_sketches_per_thread[thread].end(), m_sketches.begin());
                        else
                            std::move(m_sketches_per_thread[thread].begin(), m_sketches_per_thread[thread].end(), m_sketches.begin() + m_sketches_per_thread[thread - 1].size());
                    }

                    // Compute minimizers
                    // Take each kmer in [0, delta], [delta, 2delta] etc
                    // Instead of 0, we take [k - 1, delta], [delta, 2delta], [3delta, v.size()] etc
                    uint32_t num_delta_windows = v.size() / delta;
                    for (uint32_t delta_window = 0; delta_window < num_delta_windows; ++delta_window) {
                        uint32_t start_kmer = (delta_window == 0) ? (k - 1) : delta_window * delta;
                        uint32_t end_kmer = (delta_window == num_delta_windows - 1) ? v.size() : (delta_window + 1) * delta;

                        // To save the minimizer
                        uint32_t min_distance = 1e8;

                        // Get argmin (distances in this window)
                        for(uint32_t kmer = start_kmer; kmer < end_kmer; ++kmer) {
                            int64_t distance = 0;
                            uint32_t bit_index = 0;
                            // Build the long thing (kmer from concatted mmers)
                            for (unsigned long mmer = kmer; mmer < kmer + m_stride * m; mmer += (m / m_stride)) {
                                auto m_sketch = m_sketches[mmer];
                                // go through the bits and compare with the random_direction
                                for (uint32_t bit = 0; bit < embed_dim; ++bit, ++bit_index) {
                                    int8_t sketch_bit = ((m_sketch >> bit) & 1);
                                    int8_t random_bit = (random_direction[bit_index]);
                                    //distance += (sketch_bit * random_bit);
                                    distance += (sketch_bit - random_bit) * (sketch_bit - random_bit);
                                }
                            }
                            // Now I have the sketch and distance, check if minimal
                            if (distance < min_distance) {
                                min_distance = distance;
                                if (min_distance == 0)
                                    break; // small optimization
                            }
                        }

                        // Collect minimizers
                        for(uint32_t kmer = start_kmer; kmer < end_kmer; ++kmer) {
                            int64_t distance = 0;
                            discretized_sketch = 0;
                            uint32_t bit_index = 0;

                            // Build the long thing (kmer from concatted mmers)
                            for (unsigned long mmer = kmer; mmer < kmer + m_stride * m; mmer += (m / m_stride)) {
                                auto m_sketch = m_sketches[mmer];
                                discretized_sketch <<= embed_dim;
                                discretized_sketch |= m_sketch;
                                // go through the bits and compare with the random_direction
                                for (uint32_t bit = 0; bit < embed_dim; ++bit, ++bit_index) {
                                    int8_t sketch_bit = ((m_sketch >> bit) & 1);
                                    int8_t random_bit = (random_direction[bit_index]);
//                                    distance += (sketch_bit * random_bit);
                                    distance += (sketch_bit - random_bit) * (sketch_bit - random_bit);
                                }
                            }
                            if (distance == min_distance) {
                                sketch_maps[n_repeat][discretized_sketch].push_back(v[kmer - (k - 1)]);
                                map_backward[v[kmer]] = v[kmer - (k - 1)];
                            }
                        }
                    }
                }
            });
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
