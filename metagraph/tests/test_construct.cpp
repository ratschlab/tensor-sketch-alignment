#include <stdio.h>
#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "dbg_succinct_merge.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F,
                                 Config::StateType state) {
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << *graph->last;
    EXPECT_EQ(last, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    ostr << *(graph->W);
    EXPECT_EQ(W, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    for (size_t i = 0; i < graph->F.size(); ++i) {
        ostr << graph->F[i] << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state;
}


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F) {
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::STAT);
}


TEST(Construct, GraphDefaultConstructor) {
    DBG_succ *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ;
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ();
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        DBG_succ graph;
    });

    ASSERT_NO_THROW({
        DBG_succ graph();
    });
}

TEST(DBGSuccinct, EmptyGraph) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ");
    delete graph;
}

TEST(DBGSuccinct, SwitchState) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::STAT);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    delete graph;
}

TEST(DBGSuccinct, AddSequenceFast) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    std::vector<std::string> sequences;
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        sequences.push_back(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(3, sequences);

    //test graph construction
    test_graph(graph, "00011101101111111111111",
                      "00131124434010141720433",
                      "0 3 11 13 17 22 ");
}

TEST(Construct, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    std::vector<std::string> sequences;
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        sequences.push_back(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(3, sequences);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    EXPECT_EQ(outgoing_edges[1], graph->outgoing(1, DBG_succ::encode('$')));

    for (size_t i = 1; i < graph->get_W().size(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBG_succ::encode('$')) {
            EXPECT_EQ(outgoing_edges[i], graph->outgoing(i, graph->get_W(i)))
                << "Edge index: " << i;

            EXPECT_EQ(
                graph->succ_last(i),
                graph->incoming(graph->outgoing(i, graph->get_W(i)),
                                graph->get_node_first_value(i))
            );
            for (TAlphabet c = 0; c < DBG_succ::alph_size; ++c) {
                uint64_t node_idx = graph->incoming(i, c);
                if (node_idx) {
                    EXPECT_EQ(
                        graph->succ_last(i),
                        graph->outgoing(node_idx, graph->get_node_last_value(i))
                    );
                }
            }
        }

        //test FM index property
        EXPECT_TRUE(graph->get_last(graph->fwd(i)));
        if (graph->get_W(i)) {
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_value(graph->fwd(i)));
        }
    }

    delete graph;
}

TEST(DBGSuccinct, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    std::vector<std::string> sequences;
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        sequences.push_back(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(3, sequences);

    graph->serialize(test_dump_basename);

    DBG_succ loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";
    EXPECT_FALSE(DBG_succ() == loaded_graph);

    delete graph;
}

TEST(DBGSuccinct, AddSequenceSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A'));
        EXPECT_EQ(k + 1, graph.num_nodes());
        EXPECT_EQ(k + 2, graph.num_edges());
    }
}

TEST(DBGSuccinct, AddSequenceBugRevealingTestcase) {
    DBG_succ graph(1);
    graph.add_sequence("CTGAG", false);
}

TEST(DBGSuccinct, AddSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        graph.add_sequence("GAAC");
        EXPECT_EQ(11u, graph.num_nodes());
        EXPECT_EQ(14u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("AACG");
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC");
        graph.add_sequence("GACT");
        graph.add_sequence("ACTA");
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC", true);
        graph.add_sequence("AACG", true);
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(7u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(7u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAACT", true);
        graph.add_sequence("AAATG", true);
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequenceAnyKmerSize) {
    for (size_t k = 1; k < 10; ++k) {
        {
            DBG_succ graph(k);
            graph.add_sequence("AAAC", true);
            graph.add_sequence("AACG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", true);
            graph.add_sequence("AAATG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
        }
    }
}

void test_pred_kmer(const DBG_succ &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::deque<TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(), DBG_succ::encode);
    EXPECT_EQ(expected_idx, graph.pred_kmer(kmer)) << kmer_s << std::endl << graph;
}

TEST(DBGSuccinct, PredKmer) {
    {
        DBG_succ graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("A");

        test_pred_kmer(graph, "ACGCG", 3);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 3);
        test_pred_kmer(graph, "NNNNN", 3);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AC");

        test_pred_kmer(graph, "ACGCG", 4);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 4);
        test_pred_kmer(graph, "NNNNN", 4);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(2);
        graph.add_sequence("ATAATATCC");
        graph.add_sequence("ATACGC");
        graph.add_sequence("ATACTC");
        graph.add_sequence("ATACTA");
        graph.add_sequence("CATT");
        graph.add_sequence("CCC");
        graph.add_sequence("GGGC");
        graph.add_sequence("GGTGTGAC");
        graph.add_sequence("GGTCT");
        graph.add_sequence("GGTA");

        test_pred_kmer(graph, "$$", 4);
        test_pred_kmer(graph, "$A", 5);
        test_pred_kmer(graph, "$T", 26);
        test_pred_kmer(graph, "AT", 29);
        test_pred_kmer(graph, "TT", 35);
        test_pred_kmer(graph, "NT", 35);
        test_pred_kmer(graph, "TN", 35);
    }
}

TEST(DBGSuccinct, PredKmerRandomTest) {
    for (size_t k = 1; k < 8; ++k) {
        DBG_succ graph(k);

        for (size_t p = 0; p < 10; ++p) {
            size_t length = rand() % 400;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, true);
        }

        auto all_kmer_str = utils::generate_strings("ACGT", k);
        for (size_t i = 1; i < k; ++i) {
            auto kmer_str_suffices = utils::generate_strings("ACGT", i);
            for (size_t j = 0; j < kmer_str_suffices.size(); ++j) {
                all_kmer_str.push_back(std::string(k - i, '$')
                                        + kmer_str_suffices[j]);
            }
        }

        for (const auto &kmer_str : all_kmer_str) {
            std::deque<TAlphabet> kmer(kmer_str.size());
            std::transform(kmer_str.begin(), kmer_str.end(),
                           kmer.begin(), DBG_succ::encode);

            uint64_t lower_bound = graph.pred_kmer(kmer);

            EXPECT_FALSE(
                utils::colexicographically_greater(
                    graph.get_node_seq(lower_bound), kmer
                )
            ) << graph
              << "kmer: " << kmer_str << std::endl
              << "lower bound: " << lower_bound << std::endl
              << "which is: " << graph.get_node_str(lower_bound) << std::endl;

            if (lower_bound < graph.get_W().size() - 1) {
                EXPECT_TRUE(
                    utils::colexicographically_greater(
                        graph.get_node_seq(lower_bound + 1), kmer
                    )
                );
            }
        }
    }
}

TEST(DBGSuccinct, FindSequence) {
    for (size_t k = 1; k < 10; ++k) {
        SequenceGraph *graph = new DBG_succ(k);

        graph->add_sequence(std::string(100, 'A'));
        EXPECT_EQ(k + 2, graph->find(std::string(k, 'A')));
        EXPECT_EQ(k + 2, graph->find(std::string(2 * k, 'A')));
        EXPECT_EQ(DBG_succ::npos, graph->find(std::string(k - 1, 'A')));

        delete graph;
    }
}

TEST(DBGSuccinct, Traversals) {
    for (size_t k = 1; k < 10; ++k) {
        SequenceGraph *graph = new DBG_succ(k);

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        auto it = graph->find(std::string(k, 'A'));
        ASSERT_EQ(k + 3, it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it + 1, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it + 1, 'A'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse_back(it + 1, 'G'));

        delete graph;
    }
}