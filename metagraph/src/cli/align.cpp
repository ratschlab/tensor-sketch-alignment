#include "align.hpp"

#include <tsl/ordered_set.h>
#include <random>
#include <unordered_set>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/graph_extensions/node_rc.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include <filesystem>
#include <boost/multiprecision/cpp_int.hpp>
namespace mtg {
namespace cli {

using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace boost::multiprecision;
using mtg::seq_io::kseq_t;
using mtg::common::logger;


DBGAlignerConfig initialize_aligner_config(const Config &config) {
    assert(config.alignment_num_alternative_paths);

    DBGAlignerConfig c = {
        .num_alternative_paths = config.alignment_num_alternative_paths,
        .min_seed_length = config.alignment_min_seed_length,
        .max_seed_length = config.alignment_max_seed_length,
        .max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus,
        .min_path_score = config.alignment_min_path_score,
        .xdrop = config.alignment_xdrop,
        .min_exact_match = config.alignment_min_exact_match,
        .max_nodes_per_seq_char = config.alignment_max_nodes_per_seq_char,
        .max_ram_per_alignment = config.alignment_max_ram,
        .rel_score_cutoff = config.alignment_rel_score_cutoff,
        .gap_opening_penalty = static_cast<int8_t>(-config.alignment_gap_opening_penalty),
        .gap_extension_penalty = static_cast<int8_t>(-config.alignment_gap_extension_penalty),
        .left_end_bonus = config.alignment_end_bonus,
        .right_end_bonus = config.alignment_end_bonus,
        .forward_and_reverse_complement = !config.align_only_forwards,
        .chain_alignments = config.alignment_chain,
        .post_chain_alignments = config.alignment_post_chain,
        .alignment_edit_distance = config.alignment_edit_distance,
        .alignment_match_score = config.alignment_match_score,
        .alignment_mm_transition_score = config.alignment_mm_transition_score,
        .alignment_mm_transversion_score = config.alignment_mm_transversion_score,
        .score_matrix = DBGAlignerConfig::ScoreMatrix{},
        .embed_dim = config.embed_dim,
        .tuple_length = config.tuple_length,
        .stride = config.stride,
        .n_times_sketch = config.n_times_sketch,
        .minimizer_window = config.minimizer_window
    };

    c.set_scoring_matrix();

    c.print_summary();

    return c;
}

void map_sequences_in_file(const std::string &file,
                           const DeBruijnGraph &graph,
                           const Config &config,
                           const Timer &timer,
                           ThreadPool *thread_pool = nullptr,
                           std::mutex *print_mutex = nullptr) {
    // TODO: multithreaded
    std::ignore = std::tie(thread_pool, print_mutex);

    std::unique_ptr<std::ofstream> ofile;
    if (config.outfbase.size())
        ofile = std::make_unique<std::ofstream>(config.outfbase);

    std::ostream *out = ofile ? ofile.get() : &std::cout;

    Timer data_reading_timer;

    seq_io::read_fasta_file_critical(file, [&](kseq_t *read_stream) {
        logger->trace("Sequence: {}", read_stream->seq.s);

        if (config.query_presence
                && config.alignment_length == graph.get_k()) {

            bool found = graph.find(read_stream->seq.s,
                                    config.discovery_fraction);

            if (!config.filter_present) {
                *out << found << "\n";

            } else if (found) {
                *out << ">" << read_stream->name.s << "\n"
                            << read_stream->seq.s << "\n";
            }

            return;
        }

        assert(config.alignment_length <= graph.get_k());

        std::vector<DeBruijnGraph::node_index> nodes;
        if (config.alignment_length == graph.get_k()) {
            nodes = map_to_nodes(graph, read_stream->seq.s);

        } else {
            // TODO: make more efficient
            // TODO: canonicalization
            const DBGSuccinct &dbg = static_cast<const DBGSuccinct&>(graph);
            if (dbg.get_mode() == DeBruijnGraph::PRIMARY)
                logger->warn("Sub-k-mers will be mapped to unwrapped primary graph");

            for (size_t i = 0; i + config.alignment_length <= read_stream->seq.l; ++i) {
                nodes.emplace_back(DeBruijnGraph::npos);
                dbg.call_nodes_with_suffix_matching_longest_prefix(
                    std::string_view(read_stream->seq.s + i, config.alignment_length),
                    [&](auto node, auto) {
                        if (nodes.back() == DeBruijnGraph::npos)
                            nodes.back() = node;
                    },
                    config.alignment_length
                );
            }
        }

        size_t num_discovered = std::count_if(nodes.begin(), nodes.end(),
                                              [](const auto &x) { return x > 0; });

        const size_t num_kmers = nodes.size();

        if (config.query_presence) {
            const size_t min_kmers_discovered =
                num_kmers - num_kmers * (1 - config.discovery_fraction);
            if (config.filter_present) {
                if (num_discovered >= min_kmers_discovered)
                    *out << ">" << read_stream->name.s << "\n"
                                << read_stream->seq.s << "\n";
            } else {
                *out << (num_discovered >= min_kmers_discovered) << "\n";
            }
            return;
        }

        if (config.count_kmers) {
            std::sort(nodes.begin(), nodes.end());
            size_t num_unique_matching_kmers = 0;
            auto prev = DeBruijnGraph::npos;
            for (auto i : nodes) {
                if (i != DeBruijnGraph::npos && i != prev)
                    ++num_unique_matching_kmers;

                prev = i;
            }
            *out << fmt::format("{}\t{}/{}/{}\n", read_stream->name.s,
                                num_discovered, num_kmers, num_unique_matching_kmers);
            return;
        }

        // mapping of each k-mer to a graph node
        for (size_t i = 0; i < nodes.size(); ++i) {
            assert(i + config.alignment_length <= read_stream->seq.l);
            *out << fmt::format("{}: {}\n", std::string_view(read_stream->seq.s + i,
                                                             config.alignment_length),
                                            nodes[i]);
        }

    }, config.forward_and_reverse);

    logger->trace("File {} processed in {} sec, current mem usage: {} MB, total time {} sec",
                  file, data_reading_timer.elapsed(), get_curr_RSS() / 1e6, timer.elapsed());
}

std::string sequence_to_gfa_path(const std::string &seq,
                                 const size_t seq_id,
                                 const DeBruijnGraph &graph,
                                 const tsl::ordered_set<uint64_t> &is_unitig_end_node,
                                 const Config *config) {
    auto path_nodes = map_to_nodes_sequentially(graph, seq);

    std::string nodes_on_path;
    std::string cigars_on_path;
    const size_t overlap = graph.get_k() - 1;

    for (size_t i = 0; i < path_nodes.size() - 1; ++i) {
        if (config->output_compacted && !is_unitig_end_node.count(path_nodes[i])) {
            continue;
        }
        nodes_on_path += fmt::format("{}+,", path_nodes[i]);
        cigars_on_path += fmt::format("{}M,", overlap);
    }
    uint64_t last_node_to_print = path_nodes.back();
    // We need to print the id of the last unitig even in the case that
    // this unitig is not completely covered by the query sequence.
    while (config->output_compacted && !is_unitig_end_node.count(last_node_to_print)) {
        uint64_t unique_next_node;
        graph.adjacent_outgoing_nodes(
            last_node_to_print,
            [&](uint64_t node) { unique_next_node = node; }
        );
        last_node_to_print = unique_next_node;
    }
    nodes_on_path += fmt::format("{}+,", last_node_to_print);

    // Remove right trailing comma.
    nodes_on_path.pop_back();
    if (cigars_on_path.size()) {
        cigars_on_path.pop_back();
    }
    return fmt::format("P\t{}\t{}\t{}\n", seq_id, nodes_on_path, cigars_on_path);
}

void gfa_map_files(const Config *config,
                   const std::vector<std::string> &files,
                   const DeBruijnGraph &graph) {
    logger->trace("Starting GFA mapping:");

    tsl::ordered_set<uint64_t> is_unitig_end_node;

    graph.call_unitigs(
        [&](const auto &, const auto &path) {
            is_unitig_end_node.insert(path.back());
        },
        get_num_threads()
    );

    std::ofstream gfa_file(utils::remove_suffix(config->outfbase, ".gfa", ".path") + ".path.gfa");

    for (const std::string &file : files) {
        logger->trace("Loading sequences from FASTA file {} to append GFA paths.", file);

        std::vector<string> seq_queries;
        seq_io::FastaParser fasta_parser(file, false);
        for (const seq_io::kseq_t &kseq : fasta_parser) {
            seq_queries.push_back(kseq.seq.s);
        }
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic) shared(gfa_file)
        for (size_t i = 0; i < seq_queries.size(); ++i) {
            std::string path_string_gfa
                = sequence_to_gfa_path(seq_queries[i], i + 1, graph, is_unitig_end_node, config);
            #pragma omp critical
            gfa_file << path_string_gfa;
        }
    }
}

std::string format_alignment(const std::string &header,
                             const AlignmentResults &paths,
                             const DeBruijnGraph &graph,
                             const Config &config) {
    std::string sout;
    if (!config.output_json) {
        sout += fmt::format("{}\t{}", header, paths.get_query());
        if (paths.empty()) {
            sout += fmt::format("\t*\t*\t{}\t*\t*\t*\n", config.alignment_min_path_score);
        } else {
            for (const auto &path : paths) {
                sout += fmt::format("\t{}", path);
            }
            sout += "\n";
        }

    } else {
        Json::StreamWriterBuilder builder;
        builder["indentation"] = "";

        bool secondary = false;
        for (size_t i = 0; i < paths.size(); ++i) {
            const auto &path = paths[i];

            Json::Value json_line = path.to_json(graph.get_k(), secondary, header);
            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
            secondary = true;
        }

        if (paths.empty()) {
            Json::Value json_line = Alignment().to_json(graph.get_k(), secondary, header);
            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
        }
    }

    return sout;
}
using kmer_type = uint64_t;
using seq_type = uint8_t;

string mutate(string s, int mutation_rate, std::vector<char> alphabet) {
    std::string mutated_string = s;
    for(uint32_t i = 0; i < s.size(); ++i) {
        if (rand() % 100 < mutation_rate) {
            // mutate
            char new_c = alphabet[rand() % alphabet.size()];
            while(new_c == s[i]) {
                new_c = alphabet[rand() % alphabet.size()];
            }

            mutated_string[i] = new_c;
        }
    }
    return mutated_string;
}

void generate_sequences(const DeBruijnGraph &graph,
                        size_t min_path_size,
                        size_t max_path_size,
                        size_t num_paths,
                        int mutation_rate,
                        std::vector<char> alphabet,
                        std::vector<std::string>& spellings,
                        std::vector<std::vector<uint64_t>>& paths) {
    std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dis(1, graph.num_nodes());

    while (spellings.size() < num_paths) {
        uint64_t root_node = dis(gen);
        std::string root_node_seq = graph.get_node_sequence(root_node);
        std::vector <uint64_t> nodes;
        std::string spelling;
        bool hit_dummy = false;

        // If root node contains a $, then just keep looking
        while (root_node_seq.find("$") != string::npos) {
            root_node = dis(gen);
            root_node_seq = graph.get_node_sequence(root_node);
        }

        nodes.push_back(root_node);
        spelling = root_node_seq;
        while (nodes.size() < max_path_size && !hit_dummy && graph.outdegree(nodes.back()) > 0) {
            std::vector<uint64_t> kmers;
            std::vector<uint8_t> next_char;
            graph.call_outgoing_kmers(
                    nodes.back(),
                    [&](uint64_t target, char c) {
                        if(c == '$' || target == DeBruijnGraph::npos) {
                            hit_dummy = true;
                        } else {
                            kmers.push_back(target);
                            next_char.push_back(c);
                        }
                    });
            if (hit_dummy)
                break;
            auto pos = rand() % kmers.size();
            nodes.push_back(kmers[pos]);
            spelling += next_char[pos];
        }

        if (nodes.size() >= min_path_size && nodes.size() <= max_path_size) {
            auto mutated_string = mutate(spelling, mutation_rate, alphabet);
            spellings.push_back(mutated_string);
            paths.push_back(nodes);

            // debug
//           std::cout << spelling << " --- " << mutated_string << std::endl;
//           for(int x = 0; x < nodes.size(); ++x) {
//               std::cout << spelling.substr(x,graph.get_k()) << "--" << mutated_string.substr(x,graph.get_k()) << " " << nodes[x] << "\n";
//           }
//           std::cout << std::endl;
        }

    }
    std::unordered_set<string> unique(spellings.begin(), spellings.end());
    fprintf(stderr, "Unique sequences: %lu/%lu\n", unique.size(), num_paths);
}


int align_to_graph(Config *config) {
    assert(config);
    assert(config->infbase.size());
    // initialize graph
    auto graph = load_critical_dbg(config->infbase);
//    graph->print(std::cout);
    
    fprintf(stderr, "Number of nodes: %llu\n", graph->max_index());
    // initialize alphabet
    ts::init_alphabet("dna4");



    std::vector<std::string> spellings;
    std::vector<std::vector<uint64_t>> paths;
    std::ofstream out(config->output_path);
    if (config->experiment) {
        if (config->min_path_size > graph->max_index()) {
            logger->error("min_path_size = {} is larger than graph->max_index() = {}",
                          config->min_path_size,
                          graph->max_index());
            exit(1);
        }
        generate_sequences(*graph,
                           config->min_path_size,
                           config->max_path_size,
                           config->num_query_seqs,
                           config->mutation_rate,
                           {'A', 'T', 'G', 'C'},
                           spellings,
                           paths);
        for (uint32_t i = 0; i < config->num_query_seqs; ++i) {
            std::string header = ">Q" + std::to_string(i);
            out << header << std::endl;
            out << spellings[i] << std::endl;
        }
        out.close();
    }

    const auto &files = config->fnames;

    if (utils::ends_with(config->outfbase, ".gfa")) {
        gfa_map_files(config, files, *graph);
        return 0;
    }

    // For graphs which still feature a mask, this speeds up mapping and allows
    // for dummy nodes to be matched by suffix seeding
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    if (dbg_succ) {
        dbg_succ->reset_mask();
        if (dbg_succ->get_mode() == DeBruijnGraph::PRIMARY) {
            auto node_rc = std::make_shared<NodeRC>(*dbg_succ);
            if (node_rc->load(config->infbase)) {
                logger->trace("Loaded the adj-rc index (adjacent to reverse-complement nodes)");
                dbg_succ->add_extension(node_rc);
            } else {
                logger->warn("adj-rc index missing or failed to load. "
                             "Alignment speed will be significantly slower. "
                             "Use metagraph transform to generate an adj-rc index.");
            }
        }
    }

    Timer timer;
    ThreadPool thread_pool(get_num_threads());
    std::mutex print_mutex;

    if (config->map_sequences) {
        if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
            graph = std::make_shared<CanonicalDBG>(graph);
            logger->trace("Primary graph wrapped into canonical");
        }

        if (!config->alignment_length) {
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length > graph->get_k()) {
            logger->warn("Mapping to k-mers longer than k is not supported");
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length != graph->get_k() && !dbg_succ) {
            logger->error("Matching k-mers shorter than k only supported for succinct graphs");
            exit(1);
        }

        logger->trace("Map sequences against the de Bruijn graph with k={}",
                      graph->get_k());
        logger->trace("Length of mapped k-mers: {}", config->alignment_length);

        for (const auto &file : files) {
            logger->trace("Map sequences from file {}", file);

            map_sequences_in_file(file, *graph, *config, timer, &thread_pool, &print_mutex);
        }

        thread_pool.join();

        return 0;
    }

    DBGAlignerConfig aligner_config = initialize_aligner_config(*config);

    // compute sketches

    if (config->seeder == "sketch") {
        graph->compute_sketches(aligner_config.kmer_word_size,
                                aligner_config.embed_dim,
                                aligner_config.tuple_length,
                                aligner_config.stride,
                                aligner_config.n_times_sketch,
                                aligner_config.minimizer_window);
    }
    std::unique_ptr<AnnotatedDBG> anno_dbg;
    if (config->infbase_annotators.size()) {
        assert(config->infbase_annotators.size() == 1);
        anno_dbg = initialize_annotated_dbg(graph, *config);
    }

    for (const auto &file : files) {
        logger->trace("Align sequences from file {}", file);
        seq_io::FastaParser fasta_parser(file, config->forward_and_reverse);

        Timer data_reading_timer;

        std::unique_ptr<std::ofstream> ofile;
        if (config->outfbase.size())
            ofile = std::make_unique<std::ofstream>(config->outfbase);

        std::ostream *out = ofile ? ofile.get() : &std::cout;

        const uint64_t batch_size = config->query_batch_size_in_bytes;

        auto it = fasta_parser.begin();
        auto end = fasta_parser.end();

        size_t num_batches = 0;

        std::unordered_map<std::string, std::vector<Seed>> forward_query_seeds;
        std::unordered_map<std::string, std::vector<Seed>> rc_query_seeds;
        std::unordered_map<std::string, double> explored_nodes_per_kmer_per_query;
        std::mutex stats_mutex;

        while (it != end) {
            uint64_t num_bytes_read = 0;

            // Read a batch to pass on to a thread
            typedef std::vector<IDBGAligner::Query> SeqBatch;
            SeqBatch seq_batch;
            num_bytes_read = 0;
            for ( ; it != end && num_bytes_read <= batch_size; ++it) {
                std::string header
                    = config->fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                        && it->comment.l
                            ? utils::join_strings({ it->name.s, it->comment.s },
                                                  config->fasta_anno_comment_delim, true)
                            : std::string(it->name.s);
                seq_batch.emplace_back(std::move(header), it->seq.s);
                num_bytes_read += it->seq.l;
            }

            ++num_batches;
            thread_pool.enqueue([&,graph,batch=std::move(seq_batch)]() {
                // Make a dummy shared_ptr
                auto aln_graph
                    = std::shared_ptr<DeBruijnGraph>(std::shared_ptr<DeBruijnGraph>{}, graph.get());

                // Wrap it in CanonicalDBG if needed. This way, each thread gets its
                // own wrapper (and more importantly, its own local cache).
                if (aln_graph->get_mode() == DeBruijnGraph::PRIMARY) {
                    aln_graph = std::make_shared<CanonicalDBG>(aln_graph);
                    logger->trace("Primary graph wrapped into canonical");
                    // If backwards traversal on DBGSuccinct will be needed, then
                    // add a cache to speed it up.
                    if (dbg_succ)
                        aln_graph->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));
                }

                std::unique_ptr<IDBGAligner> aligner;
                if (anno_dbg) {
                    aligner = std::make_unique<LabeledAligner<>>(*aln_graph, aligner_config,
                                                                 anno_dbg->get_annotator());
                } else if (config->seeder == "default"){
                    logger->trace("Using default seeder");
                    aligner = std::make_unique<DBGAligner<>>(*aln_graph, aligner_config);
                } else if (config->seeder == "sketch") {
                    logger->trace("Using sketch seeder");
                    logger->trace("Sketch size (aligner): {}", aligner_config.embed_dim);
                    logger->trace("Times to sketch: {}", config->n_times_sketch);
                    aligner = std::make_unique<DBGAligner<SketchSeeder, DefaultColumnExtender, LocalAlignmentLess>>(*aln_graph, aligner_config);
                }
                aligner->align_batch(batch,
                    [&](const std::string &header, AlignmentResults&& paths) {
                        const auto &res = format_alignment(header, paths, *graph, *config);
                        std::lock_guard<std::mutex> lock(print_mutex);
                        *out << res;
                    }
                );

                // Merge into the big map
                std::lock_guard<std::mutex> lock1(stats_mutex);
                {
                    forward_query_seeds.merge(aligner->forward_query_seeds);
                    rc_query_seeds.merge(aligner->rc_query_seeds);
                    explored_nodes_per_kmer_per_query.merge(aligner->explored_nodes_per_kmer_per_query);
                }
            });
        };

        thread_pool.join();

        float avg_time = data_reading_timer.elapsed() / config->num_query_seqs;

        // Compute the recall now
        int recalled_paths = 0;
        double precision = 0.0;
        int n_precision_gt_zero = 0;
        for(int seq_i = 0; seq_i < config->num_query_seqs; ++seq_i) {
            std::string header = "Q" + std::to_string(seq_i);

            // Get query sequence and path
            std::string query_seq = spellings[seq_i];
            std::vector<uint64_t> path = paths[seq_i];

            // Get matched seeds (fwd and bwd)
            auto fwd_seeds = forward_query_seeds[header];
            auto rc_seeds = rc_query_seeds[header]; //unused
            precision += explored_nodes_per_kmer_per_query[header];
            if(explored_nodes_per_kmer_per_query[header] > 0.0)
                n_precision_gt_zero++;
            int recalled = 0;
            // Check if the forward sequence recalled
            for(auto seed : fwd_seeds) {
                auto seed_nodes = seed.get_nodes();

                int match_start = seed.get_clipping();
                int num_matched = seed_nodes.size();

                for(int i = 0; i < num_matched; ++i) {
                    std::string kmer = query_seq.substr(match_start + i, graph->get_k());
                    uint64_t node = path[match_start + i];
                    if (config->seeder == "sketch") {
                        if (std::count(seed_nodes.begin(), seed_nodes.end(), graph->map_backward[node])) {
                            recalled++;
                            break;
                        }
                    }
                    else {
                        // Default seeder
                        if (std::count(seed_nodes.begin(), seed_nodes.end(), node)) {
                            recalled++;
                            break;
                        }
                    }
                }

                if (recalled > 0) {
                    break;
                }
            }
            recalled_paths += (recalled > 0);
        }

        if (n_precision_gt_zero > 0.0)
          std::cout << "{"
                    << "\"recall\":" << (float)recalled_paths / (float)config->num_query_seqs << ","
                    << "\"precision\":" << precision/ (double)n_precision_gt_zero << ","
                    << "\"avg_time\":" << avg_time
                    << "}";
        else
          std::cout << "{"
                    << "\"recall\":" << (float)recalled_paths / (float)config->num_query_seqs << ","
                    << "\"precision\":" << 0 << ","
                    << "\"avg_time\":" << avg_time
                    << "}";
        logger->trace("File {} processed in {} sec, "
                      "num batches: {}, batch size: {} KB, "
                      "current mem usage: {} MB, total time {} sec",
                      file, data_reading_timer.elapsed(), num_batches, batch_size / 1e3,
                      get_curr_RSS() / 1e6, timer.elapsed());
    }

    return 0;
}

} // namespace cli
} // namespace mtg