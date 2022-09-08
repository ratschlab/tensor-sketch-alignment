#include "align.hpp"
#include "seqgen.hpp"

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
using namespace mtg::graph::seqgen;
using namespace boost::multiprecision;
using mtg::seq_io::kseq_t;
using mtg::common::logger;

using kmer_type = uint64_t;
using seq_type = uint8_t;

string mutate(string s, int mutation_rate, std::vector<char> alphabet) {
    std::string mutated_string;
    uint32_t counter = 0;
    uint32_t thresh1 = 33;
    uint32_t thresh2 = 66;
    uint32_t nmutations = 0;
    while (counter < s.size()) {
        if (rand() % 100 < mutation_rate) {
            nmutations += 1;
            uint32_t chance = rand() % 100;
            if (chance <= thresh1) {
                // substitution
                char new_c = alphabet[rand() % alphabet.size()];
                while (new_c == s[counter]) {
                    new_c = alphabet[rand() % alphabet.size()];
                }
                mutated_string += new_c;
                counter++;
            } else if (chance > thresh1 && chance <= thresh2) {
                // deletion
                counter+=1;
            } else if (chance > thresh2) {
                // insertion
                for(int i = 0; i < 1; ++i) {
                    char new_c = alphabet[rand() % alphabet.size()];
                    mutated_string += new_c;
                }
            }
        } else {
            // nothing
            mutated_string += s[counter];
            counter++;
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
                        std::vector<std::string>& reference_spellings,
                        std::vector<std::string>& mutated_spellings,
                        std::vector<std::vector<uint64_t>>& paths) {
    std::mt19937 gen(1);
    std::uniform_int_distribution<uint64_t> dis(1, graph.num_nodes());

    while (mutated_spellings.size() < num_paths) {
        uint64_t root_node = dis(gen);
        std::string root_node_seq = graph.get_node_sequence(root_node);
        std::vector <uint64_t> nodes;
        std::string reference_spelling;
        bool hit_dummy = false;

        // If root node contains a $, then just keep looking
        while (root_node_seq.find("$") != string::npos) {
            root_node = dis(gen);
            root_node_seq = graph.get_node_sequence(root_node);
        }

        nodes.push_back(root_node);
        reference_spelling= root_node_seq;
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
            reference_spelling += next_char[pos];
        }

        if (nodes.size() >= min_path_size && nodes.size() <= max_path_size) {
            reference_spellings.push_back(reference_spelling);
            std::string mutated_spelling = mutate(reference_spelling, mutation_rate, alphabet);
            mutated_spellings.push_back(mutated_spelling);
            paths.push_back(nodes);
        }

    }
    std::unordered_set<string> unique(reference_spellings.begin(), reference_spellings.end());
    fprintf(stderr, "Unique sequences: %lu/%lu\n", unique.size(), num_paths);
}


int generate_sequences(Config *config) {
    assert(config);
    assert(config->infbase.size());
    // initialize graph
    auto graph = load_critical_dbg(config->infbase);
    
    fprintf(stderr, "Number of nodes: %llu\n", graph->max_index());
    // initialize alphabet
    ts::init_alphabet("dna4");
    std::vector<std::string> reference_spellings;
    std::vector<std::string> mutated_spellings;
    std::vector<std::vector<uint64_t>> paths;
    std::ofstream reference_out;
    std::ofstream mutated_out;
    reference_out = std::ofstream(config->output_path + "reference_" + std::to_string(config->mutation_rate) + ".fa");
    mutated_out = std::ofstream(config->output_path + "mutated_" + std::to_string(config->mutation_rate) + ".fa");
    if (config->min_path_size > graph->max_index()) {
        logger->error("min_path_size = {} is larger than graph->max_index() = {}",
                      config->min_path_size,
                      graph->max_index());
        exit(1);
    }
    srand(0);
    generate_sequences(*graph,
                       config->min_path_size,
                       config->max_path_size,
                       config->num_query_seqs,
                       config->mutation_rate,
                       {'A', 'T', 'G', 'C'},
                       reference_spellings,
                       mutated_spellings,
                       paths);
    for (uint32_t i = 0; i < config->num_query_seqs; ++i) {
        std::string header = ">Q" + std::to_string(i);
        reference_out << header << std::endl;
        mutated_out << header << std::endl;

        reference_out << reference_spellings[i] << std::endl;
        mutated_out << mutated_spellings[i] << std::endl;

        // debug
//            auto nodes = paths[i];
//            std::cout << header << std::endl;
//            std::cout << reference_spellings[i] << "\n" << mutated_spellings[i]<< std::endl;
//            for(int x = 0; x < nodes.size(); ++x) {
//                std::cout << reference_spellings[i].substr(x,graph->get_k()) << "--" << mutated_spellings[i].substr(x, graph->get_k()) << " " << nodes[x] << "\n";
//            }
//            std::cout << std::endl;
    }
    mutated_out.close();
    reference_out.close();
    return 0;
}

} // namespace cli
} // namespace mtg
