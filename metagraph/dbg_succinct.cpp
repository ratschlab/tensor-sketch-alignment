#include "dbg_succinct.hpp"

/**
 * This class contains a succinct representation of the de bruijn graph
 * following ideas and suggestions presented here:
 * http://link.springer.com/chapter/10.1007/978-3-642-33122-0_18
 *
 * There is also conceptual code available at
 * https://code.google.com/p/csalib/downloads/list
 * that has been used as a reference for this implementation.
 */

#include <assert.h>
#include <pthread.h>
#include <vector>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <cstdio>
#include <parallel/algorithm>

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include <libmaus2/util/NumberSerialisation.hpp>
#include <zlib.h>

#include "serialization.hpp"


const std::string DBG_succ::alphabet = "$ACGTNX$ACGTNXn";
const size_t DBG_succ::alph_size = 7;

const TAlphabet kCharToNucleotide[128] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5
};


DBG_succ::DBG_succ(size_t k, bool sentinel)
      : F(alph_size, 0), k_(k), p_(0) {

    last->insertBit(0, false);
    W->insert(0, 0);

    if (sentinel) {
        last->insertBit(1, true);
        W->insert(0, 0);
        for (size_t j = 1; j < alph_size; j++) {
            F[j] = 1;
        }
        p_ = 1;
    }
    state = Config::DYN;
}

DBG_succ::~DBG_succ() {
    delete W;
    delete last;
}

/**
* Given a pointer to a graph structures G1 and G2, the function compares their elements to the
* each other. It will perform an element wise comparison of the arrays W, last and
* F and will only check for identity. If any element differs, the function will return
* false and true otherwise.
*/
bool DBG_succ::operator==(const DBG_succ &other) const {
    // compare size
    if (W->size() != other.W->size()) {
        std::cerr << "sizes of graphs differ" << std::endl;
        std::cerr << "1: " << W->size() << std::endl;
        std::cerr << "2: " << other.W->size() << std::endl;
        return false;
    }

    if (F.size() != other.F.size()) {
        std::cerr << "sizes of F arrays differ" << std::endl;
        std::cerr << "1: " << F.size() << std::endl;
        std::cerr << "2: " << other.F.size() << std::endl;
        return false;
    }

    // compare last
    for (size_t i = 0; i < W->size(); ++i) {
        if (get_last(i) != other.get_last(i)) {
            std::cerr << "last differs at position " << i << std::endl;
            std::cerr << "1: last[" << i << "] = " << get_last(i)  << std::endl;
            std::cerr << "2: last[" << i << "] = " << other.get_last(i) << std::endl;
            return false;
        }
    }

    // compare W
    for (size_t i = 0; i < W->size(); ++i) {
        if (get_W(i) != other.get_W(i)) {
            std::cerr << "W differs at position " << i << std::endl;
            std::cerr << "1: W[" << i << "] = " << get_W(i)  << std::endl;
            std::cerr << "2: W[" << i << "] = " << other.get_W(i) << std::endl;
            return false;
        }
    }

    // compare F
    for (size_t i = 0; i < F.size(); ++i) {
        if (get_F(i) != other.get_F(i)) {
            std::cerr << "F differs at position " << i << std::endl;
            std::cerr << "1: F[" << i << "] = " << get_F(i) << std::endl;
            std::cerr << "2: F[" << i << "] = " << other.get_F(i) << std::endl;
            return false;
        }
    }

    return true;
}

void DBG_succ::serialize(const std::string &outbase) const {
    // write Wavelet Tree
    std::ofstream outstream(outbase + ".W.dbg");
    W->serialise(outstream);
    outstream.close();

    // write last array
    outstream.open(outbase + ".l.dbg");
    last->serialise(outstream);
    outstream.close();

    // write F values and k
    outstream.open(outbase + ".F.dbg");
    outstream << ">F" << std::endl;
    for (size_t i = 0; i < F.size(); ++i) {
        outstream << F.at(i) << "\n";
    }
    outstream << ">k\n"
              << k_ << "\n"
              << ">p\n"
              << p_ << "\n"
              << ">s\n"
              << state << "\n";
    outstream.close();
}

bool DBG_succ::load(const std::string &infbase) {
    // if not specified in the file, the default for loading is dynamic
    state = Config::DYN;

    // load F and k and p
    F.resize(0);
    F.reserve(alph_size);

    try {
        std::ifstream instream(infbase + ".F.dbg");
        char mode = 0;
        std::string cur_line;

        while (std::getline(instream, cur_line)) {
            if (cur_line[0] == '>') {
                if (cur_line.length() < 2)
                    return false;
                mode = cur_line[1];
                continue;
            }
            switch (mode) {
                case 'F':
                    F.push_back(std::stoul(cur_line));
                    break;
                case 'k':
                    k_ = std::stoul(cur_line);
                    break;
                case 'p':
                    p_ = std::stoul(cur_line);
                    break;
                case 's':
                    state = static_cast<Config::StateType>(std::stoul(cur_line));
                    break;
                default:
                    return false;
            }
        }
        instream.close();

        if (F.size() != alph_size)
            return false;

        // load W and last arrays
        delete W;
        delete last;
        switch (state) {
            case Config::DYN:
                W = new wavelet_tree_dyn(4);
                last = new bit_vector_dyn();
                break;
            case Config::STAT:
                W = new wavelet_tree_stat(4);
                last = new bit_vector_stat();
                break;
            case Config::CSTR:
                assert(false && "Must not happen");
                return false;
        }
        std::ifstream instream_W(infbase + ".W.dbg");
        std::ifstream instream_l(infbase + ".l.dbg");
        return W->deserialise(instream_W) && last->deserialise(instream_l);
    } catch (...) {
        return false;
    }
}

//
//
// QUERY FUNCTIONS
//
//

/**
 * Uses the object's array W, a given position i in W and a character c
 * from the alphabet and returns the number of occurences of c in W up to
 * position i.
 */
uint64_t DBG_succ::rank_W(uint64_t i, TAlphabet c) const {

    // deal with  border conditions
    if (i <= 0)
        return 0;
    return W->rank(c, std::min(i, W->size() - 1)) - (c == 0);
}

/**
 * Uses the array W and gets a count i and a character c from
 * the alphabet and returns the position of the i-th occurence of c in W.
 */
uint64_t DBG_succ::select_W(uint64_t i, TAlphabet c) const {

    // deal with  border conditions
    if (i <= 0)
        return 0;

    return i + (c == 0) <= W->rank(c, W->size() - 1)
                ? W->select(c, i + (c == 0))
                : W->size();
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the last index of a character c preceding in W[1..i].
 */
uint64_t DBG_succ::pred_W(uint64_t i, TAlphabet c) const {
    return select_W(rank_W(i, c), c);
}

/**
 * This is a convenience function that returns for array W, a position i and
 * a character c the first index of a character c in W[i..N].
 */
uint64_t DBG_succ::succ_W(uint64_t i, TAlphabet c) const {
    return select_W(rank_W(i - 1, c) + 1, c);
}

/**
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
uint64_t DBG_succ::rank_last(uint64_t i) const {
    // deal with  border conditions
    if (i <= 0)
        return 0;
    return last->rank1(i);
}

/**
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
uint64_t DBG_succ::select_last(uint64_t i) const {
    // deal with  border conditions
    if (i <= 0)
        return 0;
    // for some reason the libmaus2 select is 0 based ...
    return std::min(last->select1(i), last->size());
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
uint64_t DBG_succ::pred_last(uint64_t i) const {
    return select_last(rank_last(i));
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
uint64_t DBG_succ::succ_last(uint64_t i) const {
    return select_last(rank_last(i - 1) + 1);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character.
 */
uint64_t DBG_succ::bwd(uint64_t i) const {
    // get value of last position in node i
    TAlphabet c = get_node_end_value(i);
    // get the offset for the last position in node i
    uint64_t o = F[c];
    // compute the offset for this position in W and select it
    //usleep(500000);
    return select_W(rank_last(i) - rank_last(o), c);
}

/**
 * This functions gets a position i reflecting the r-th occurence of the corresponding
 * character c in W and returns the position of the r-th occurence of c in last.
 */
uint64_t DBG_succ::fwd(uint64_t i) const {
    // get value of W at position i
    TAlphabet c = (*W)[i] % alph_size;
    // get the offset for position c
    uint64_t o = F[c];
    // get the rank of c in W at position i
    uint64_t r = rank_W(i, c);
    // select the index of the position in last that is rank many positions after offset
    return select_last(rank_last(o) + r);
}


/**
 * Using the offset structure F this function returns the value of the last
 * position of node i.
 */
TAlphabet DBG_succ::get_node_end_value(uint64_t i) const {
    if (i == 0)
        return 0;
    for (size_t j = 0; j < F.size(); j++) {
        if (F[j] >= i)
            return j - 1;
    }
    return F.size() - 1;
}


/**
 * Given index of node i, the function returns the
 * first character of the node.
 */
TAlphabet DBG_succ::get_node_begin_value(uint64_t i) const {
    if (i == 1)
        return 0;

    for (size_t j = 0; j < k_ - 1; ++j) {
        i = bwd(succ_last(i));
        if (i == 1)
            return 0;
    }
    return get_node_end_value(i);
}

/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the outgoing edge with label c.
 */
uint64_t DBG_succ::outgoing_edge_idx(uint64_t i, TAlphabet c) const {

    if (i > W->size())
        return 0;
    if (i == 0)
        return 0;
    std::pair<uint64_t, uint64_t> R = get_equal_node_range(i);

    uint64_t j1 = pred_W(R.second, c);
    uint64_t j2 = pred_W(R.second, c + alph_size);
    uint64_t j = (j1 < j2) ? j2 : j1;
    if (j < R.first || j >= W->size())
        return 0;

    return j;
}



/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the node the edge is pointing to.
 */
uint64_t DBG_succ::outgoing(uint64_t i, TAlphabet c) const {
    if (i > W->size())
        return 0;
    if (i == 0)
        return 0;

    uint64_t j = outgoing_edge_idx(i, c);
    if (j == 0)
        return 0;
    j = fwd(j);
    if (j == 0 || j == W->size())
        return 0;
    return j;
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the node the incoming edge belongs to.
 */
uint64_t DBG_succ::incoming(uint64_t i, TAlphabet c) const {
    if (i == 1)
        return 0;
    c %= alph_size;
    TAlphabet d = get_node_end_value(i);
    uint64_t x = bwd(i);
    if (get_node_begin_value(x) == c) {
        return succ_last(x);
    }
    uint64_t y = succ_W(x + 1, d);
    while (true) {
        x = succ_W(x+1, d + alph_size);
        if (x >= y) {
            break;
        }
        if (get_node_begin_value(x) == c) {
            return succ_last(x);
        }
    }
    return 0;
}

/**
 * Given a node index i, this function returns the number of outgoing
 * edges from node i.
 */
uint64_t DBG_succ::outdegree(uint64_t i) const {
    return (i < W->size()) ? succ_last(i) - pred_last(i - 1) : 0;
}


/**
 * Given a node index i, this function returns the number of incoming
 * edges to node i.
 */
uint64_t DBG_succ::indegree(uint64_t i) const {
    if (i < 2)
        return 0;
    uint64_t x = bwd(succ_last(i));
    TAlphabet d = get_node_end_value(i);
    uint64_t y = succ_W(x + 1, d);
    return 1 + rank_W(y, d + alph_size) - rank_W(x, d + alph_size);
}



/**
 * Given a string str and a maximal number of edit operations
 * max_distance, this function returns all nodes with labels at most
 * max_distance many edits away from str.
 */
std::vector<HitInfo> DBG_succ::index_fuzzy(const std::string &str,
                                           uint64_t max_distance) const {

    std::vector<HitInfo> result;
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits;
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits2;
    uint64_t rl;
    uint64_t ru;

    // walk through pattern, thereby collecting possible partial matches
    // once the end of the pattern is reached, add match to results

    // init match/mismatch to first pattern position
    TAlphabet s = encode(str[0]);
    for (TAlphabet b = 1; b < 5; ++b) {
        rl = succ_last(F[b] + 1);
        ru = F[b + 1];
        //std::cout << "pushing: rl " << rl << " ru " << ru << " str_pos 1 max_distance " << (uint64_t) (b != s) << std::endl;
        //std::cout << "s " << s << " b " << b << std::endl;
        std::vector<uint64_t> tmp;
        hits.push({ rl, ru, 1, 1, static_cast<uint64_t>(b != s),
                    std::string(1, decode(b)), tmp });

        // opening/extending a gap in the pattern starting with the first position
        if (max_distance > 0) {
            for (size_t p = 1; p < str.length() - 1; ++p) {
                TAlphabet ss = encode(str[p]);
                if ((p + (b != ss)) > max_distance)
                    break;
                hits.push({ rl, ru, p + 1, 1, p + (b != ss),
                            std::string(p, 'd') + std::string(1, decode(b)), tmp });
                //std::cout << "a) adding '-'" << std::endl;
            }
        }
    }

    // walk through pattern thereby extending all partial hits
    while (hits.size() > 0) {
        while (hits.size() > 0) {
            HitInfo curr_hit(hits.top());
            hits.pop();
            //std::cout << "loaded: rl " << curr_hit.rl << " ru " << curr_hit.ru << " dist " << curr_hit.distance << std::endl;

            if (curr_hit.str_pos < str.length()) {

                // opening/extending a gap in the graph, leaving current pattern position unmatched
                if (curr_hit.distance < max_distance) {
                    hits2.push({ curr_hit.rl, curr_hit.ru, curr_hit.str_pos + 1,
                                 curr_hit.graph_pos, curr_hit.distance + 1,
                                 curr_hit.cigar + 'd', curr_hit.path });
                    //std::cout << "b) " << curr_hit.cigar << " adding '-'" << std::endl;
                }

                s = encode(str[curr_hit.str_pos]);

                // has the number of matches exceeded the node length?
                // there are three possible scenarios for extension of the path:
                //  1) pattern is shorter than the node length --> get an interval of matching nodes
                //  2) pattern length exactly mathces the node length --> there is one correponding node
                //  3) pattern is longer than the node length --> we append to a path
                if (curr_hit.graph_pos >= k_) {
                //    std::cout << "push back tp path " << curr_hit.rl << std::endl;
                    curr_hit.path.push_back(curr_hit.rl);
                }

                // iterate through all possible extensions of current position
                for (TAlphabet b = 1; b < 5; ++b) {
                    if (curr_hit.distance <= max_distance) {

                        // we cannot afford any more mismatches
                        if ((curr_hit.distance + (b != s)) > max_distance)
                            continue;

                        // re-define range of nodes to check for outgoing nodes
                        rl = std::min(succ_W(pred_last(curr_hit.rl - 1) + 1, b),
                                      succ_W(pred_last(curr_hit.rl - 1) + 1, b + alph_size));
                        ru = std::max(pred_W(curr_hit.ru, b),
                                      pred_W(curr_hit.ru, b + alph_size));

                        // the current range in W does not contain our next symbol
                        if ((rl >= W->size()) || (ru >= W->size()) || (rl > ru))
                            continue;

                        // update the SA range with the current symbol b
                        rl = outgoing(rl, b);
                        ru = outgoing(ru, b);

                        // range is empty
                        if ((rl == 0) && (ru == 0))
                            continue;

                        // add hit for extension in next step
                        hits2.push({ rl, ru, curr_hit.str_pos + 1,
                                     curr_hit.graph_pos + 1, curr_hit.distance + (b != s),
                                     curr_hit.cigar + decode(b), curr_hit.path });

                        // opening/extending a gap in the pattern, leaving current graph position unmatched
                        // --> choose any available mismatching next edge
                        if (b != s) {
                            hits2.push({ rl, ru, curr_hit.str_pos,
                                         curr_hit.graph_pos + 1, curr_hit.distance + 1,
                                         curr_hit.cigar + 'i', curr_hit.path });
                        }
                    }
                }
            } else {
                // collect results
                //std::make_pair(curr_hit.rl < curr_hit.ru ? curr_hit.ru : curr_hit.rl, curr_hit.cigar));
                result.push_back(curr_hit);
            }
        }
        hits.swap(hits2);
    }

    return result;
}


/**
 * Given a node label s, this function returns the index
 * of the corresponding node or the closest predecessor, if no node
 * with the sequence is not found.
 */
uint64_t DBG_succ::colex_upper_bound(std::deque<TAlphabet> str) const {
    // get first
    std::deque<TAlphabet>::iterator it = str.begin();
    TAlphabet s1 = *it % alph_size;

    // init range

    // lower bound
    uint64_t rl = succ_last(F.at(s1) + 1);
    // upper bound
    uint64_t ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->size() - 1);
    while (rl > ru && s1 > 0) {
        s1--;
        rl = succ_last(F.at(s1) + 1);
        ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->size() - 1);
    }
    /*if (s1 == 0) {
        s1 = *it % alph_size + 1;
        rl = succ_last(F.at(s1) + 1);
        ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);
        while (rl > ru && s1 < alph_size) {
            s1++;
            rl = succ_last(F.at(s1) + 1);
            ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);
        }
    }*/
    //std::cerr << "s: " << s1 << " rl: " << rl << " ru: " << ru << std::endl;

    it++;
    uint64_t pll, puu;
    bool before = false;
    //fprintf(stderr, "char: %i rl: %i ru: %i\n", (int) s1, (int) rl, (int) ru);
    // update range iteratively while scanning through s
    for (; it != str.end(); it++) {
        s1 = *it % alph_size;
        pll = this->pred_last(rl - 1) + 1;
        puu = this->succ_last(ru);
        before = false;

        //std::cerr << "s: " << s1 << " rl: " << rl << " ru: " << ru << " pll: " << pll << std::endl;

        rl = std::min(succ_W(pll, s1), succ_W(pll, s1 + alph_size));
        ru = std::max(pred_W(puu, s1), pred_W(puu, s1 + alph_size));
        if (rl > puu) {
            rl = std::max(pred_W(pll, s1), pred_W(pll, s1 + alph_size));
            if (rl == 0) {
                rl = std::min(succ_W(pll, s1), succ_W(pll, s1 + alph_size));
                if (rl >= W->size()) {
                    s1--;
                    while (s1 > 0) {
                        rl = std::max(pred_W(W->size() - 1, s1), pred_W(W->size() - 1, s1));
                        if (rl < W->size())
                            break;
                        s1--;
                    }
                    if (s1 == 0) {
                        s1 = (*it % alph_size) + 1;
                        before = true;
                        while (s1 < alph_size) {
                            rl = std::min(succ_W(1, s1), succ_W(1 + alph_size, s1));
                            if (rl < W->size())
                                break;
                            s1++;
                        }
                    }
                } else {
                    s1 = (*W)[rl];
                    before = true;
                }
            } else {
                s1 = (*W)[rl];
            }
        }

        if (ru == 0)
            ru = rl;

        //std::cerr << "endloop - rl: " << rl << " s1: " << s1 << " ru: " << ru << std::endl;
        assert(rl <= ru);

        rl = outgoing(rl, s1);
        ru = outgoing(ru, s1);
    }
    return pred_last(rl) - before;
}



/**
 * Given a position i, this function returns the boundaries of the interval
 * of nodes identical to node i (ignoring the values in W).
 */
std::pair<uint64_t, uint64_t> DBG_succ::get_equal_node_range(uint64_t i) const {
    return std::make_pair(pred_last(i - 1) + 1, succ_last(i));
}

/**
 * Given index i of a node and a value k, this function
 * will return the k-th last character of node i.
 */
std::pair<TAlphabet, uint64_t> DBG_succ::get_minus_k_value(uint64_t i, uint64_t k) const {
    for (; k > 0; --k)
        i = bwd(succ_last(i));
    return std::make_pair(get_node_end_value(i), bwd(succ_last(i)));
}



/**
 * This function gets two node indices and returns whether the
 * node labels share a k-1 suffix.
 */
bool DBG_succ::compare_node_suffix(uint64_t i1, uint64_t i2) const {
    for (size_t ii = 0; ii < k_ - 1; ++ii) {
        if (get_node_end_value(i1) != get_node_end_value(i2)) {
            return false;
        }
        i1 = bwd(succ_last(i1));
        i2 = bwd(succ_last(i2));
    }
    return true;
}

bool DBG_succ::compare_node_suffix(TAlphabet *ref, uint64_t i2) const {
    TAlphabet *i1 = &ref[k_ - 1];
    for (size_t ii=0; ii < k_ - 1; ii++) {
        if (*i1 != get_node_end_value(i2)) {
            return false;
        }
        i1 = &ref[k_ - 2 - ii];
        i2 = bwd(succ_last(i2));
    }
    return true;
}

/**
 * This function returns true if node i is a terminal node.
 */
bool DBG_succ::is_terminal_node(uint64_t i) const {
    for (size_t ii = 0; ii < k_ - 1; ii++) {
        if (get_node_end_value(i) % alph_size != 0) {
            return false;
        }
        i = bwd(i);
    }
    return true;
}

/**
* Given a node index k_node, this function returns the k-mer sequence of the
* node in a deque data structure.
*/
std::deque<TAlphabet> DBG_succ::get_node_seq(uint64_t k_node) const {
    std::deque<TAlphabet> ret;

    for (uint64_t curr_k = 0; curr_k < this->k_; ++curr_k) {
        auto k_val = get_minus_k_value(k_node, 0);
        ret.push_front(k_val.first);
        k_node = k_val.second;
    }

    return ret;
}

/**
* Given a node index k_node, this function returns the k-mer sequence of the
* node as a string.
*/
std::string DBG_succ::get_node_str(uint64_t k_node) const {
    std::string node_string(k_, 0);
    auto node_encoding = get_node_seq(k_node);
    std::transform(node_encoding.begin(), node_encoding.end(), node_string.begin(),
                   [](TAlphabet c) { return decode(c); });
    return node_string;
}

std::vector<uint64_t> DBG_succ::align(const std::string &sequence,
                                      uint64_t alignment_length) const {
    std::vector<uint64_t> indices;

    if (alignment_length == 0)
        alignment_length = k_;

    std::vector<HitInfo> curr_result;
    for (uint64_t i = 0; i < sequence.size() - alignment_length + 1; ++i) {
        std::string kmer(sequence.data() + i, sequence.data() + i + alignment_length);
        indices.push_back(this->index(kmer, kmer.size()));
    }

    return indices;
}

std::vector<std::vector<HitInfo>> DBG_succ::align_fuzzy(const std::string &sequence,
                                                        uint64_t alignment_length,
                                                        uint64_t max_distance) const {
    std::vector<std::vector<HitInfo>> hit_list;

    if (alignment_length == 0) {

    } else {
        alignment_length = alignment_length < 2 ? 2 : alignment_length;
        alignment_length = alignment_length < sequence.size()
                                    ? alignment_length
                                    : sequence.size();
        for (uint64_t i = 0; i < sequence.size() - alignment_length + 1; ++i) {
            std::string kmer(sequence.data() + i, sequence.data() + i + alignment_length);
            hit_list.push_back(this->index_fuzzy(kmer, max_distance));
        }
    }
    return hit_list;
}

/*
 * Returns the number of nodes on the current graph.
 */
uint64_t DBG_succ::num_nodes() const {
    return rank_last(last->size() - 1);
}

/*
 * Return the number of edges in the current graph.
 */
 //TODO: Check this!!!
uint64_t DBG_succ::num_edges() const {
    return W->size() - 1;
}

/**
 * This function gets a value of the alphabet c and updates the offset of
 * all following values by +1 is positive is true and by -1 otherwise.
 */
void DBG_succ::update_F(TAlphabet c, bool positive) {
    for (TAlphabet i = c + 1; i < F.size(); i++) {
        F[i] += 2 * static_cast<int>(positive) - 1;
    }
}

/**
 * This function gets a local range in W from lower bound l
 * to upper bound u and swaps the inserted element to the
 * righ location.
 */
//TODO: this function can be improved
//TODO: fix_order_in_W_range
void DBG_succ::sort_W_locally(uint64_t l, uint64_t u) {
    assert(state == Config::DYN);

    for (uint64_t s = u; s > l; --s) {
        auto first = get_W(s - 1);
        auto second = get_W(s);
        if ((second % alph_size) < (first % alph_size)) {
            W_set_value(s - 1, second);
            W_set_value(s, first);
        }
    }
    for (uint64_t s = l; s < u; ++s) {
        auto first = get_W(s);
        auto second = get_W(s + 1);
        if ((first % alph_size) > (second % alph_size)) {
            W_set_value(s + 1, first);
            W_set_value(s, second);
        }
    }
}

/**
 * This is a convenience function to replace the value at
 * position i in W with val.
 */
void DBG_succ::W_set_value(uint64_t i, TAlphabet val) {
    W->set(i, val);
}

TAlphabet DBG_succ::encode(char s) {
    return kCharToNucleotide[static_cast<int>(s)];
}

char DBG_succ::decode(TAlphabet c) {
    assert(c < alphabet.size());
    return alphabet[c];
}

void DBG_succ::switch_state(Config::StateType state) {

    //std::cerr << "switching state from " << this->state << " to " << state << std::endl;
    if (this->state == state)
        return;

    switch (state) {
        case Config::CSTR: {
            this->state = Config::CSTR;
            break;
        }
        case Config::DYN: {
            if (this->state == Config::CSTR) {
                delete W;
                W = new wavelet_tree_dyn(4, W_stat);
                W_stat.clear();

                delete last;
                last = new bit_vector_dyn(last_stat);
                last_stat.clear();

                //delete bridge;
                //if (bridge_stat.size()) {
                //    bridge = new bit_vector_dyn(bridge_stat);
                //    bridge_stat.clear();
                //} else {
                //    bridge = NULL;
                //}
            } else {
                wavelet_tree *W_new = new wavelet_tree_dyn(4, *W);
                delete W;
                W = W_new;

                bit_vector *last_new = new bit_vector_dyn(*last);
                delete last;
                last = last_new;
            }
            this->state = Config::DYN;
            break;
        }
        case Config::STAT: {
            wavelet_tree *W_new = new wavelet_tree_stat(4, *W);
            delete W;
            W = W_new;

            bit_vector *last_new = new bit_vector_stat(*last);
            delete last;
            last = last_new;

            this->state = Config::STAT;
            break;
        }
    }
}


void DBG_succ::print_state() const {
    std::cout << "Index" << "\t" << "L"
                         << "\t" << "Vertex"
                         << "\t" << "W" << std::endl;

    for (uint64_t i = 1; i < W->size(); i++) {
        std::cout << i << "\t" << get_last(i)
                       << "\t" << get_node_str(i)
                       << "\t" << decode(get_W(i))
                               << (get_W(i) > alph_size
                                       ? "-"
                                       : "")
                               << (i == this->p_
                                      ? "<"
                                      : "")
                               << std::endl;
    }
}

/*
 * Returns the sequence stored in W and prints the node
 * information in an overview.
 * Useful for debugging purposes.
 */
void DBG_succ::print_seq() const {

    const uint64_t linelen = 80;

    for (uint64_t start = 1; start < W->size(); start += linelen) {
        uint64_t end = start + linelen < W->size() ? start + linelen : W->size();

        for (uint64_t i = start; i < end; i++) {
            if (i % 10 == 0)
                fprintf(stdout, "%d", static_cast<int>((i / 10) % 10));
            else
                fprintf(stdout, " ");
        }
        fprintf(stdout, "\n");

        for (uint64_t i = start; i < end; i++) {
            if (get_W(i) >= alph_size)
                fprintf(stdout, "|");
            else
                fprintf(stdout, " ");
        }
        fprintf(stdout, "\n");

        for (uint64_t i = start; i < end; i++) {
            fprintf(stdout, "%c", decode(get_W(i)));
        }
        fprintf(stdout, "\n");

        std::vector<std::string> nodes;
        for (uint64_t i = start; i < end; i++) {
            if (p_ == i)
                fprintf(stdout, "*");
            else
                fprintf(stdout, " ");
            nodes.push_back(get_node_str(i));
        }
        fprintf(stdout, "\n");

        for (size_t l = 0; l < k_; l++) {
            for (uint64_t i = start; i < end; i++) {
                fprintf(stdout, "%c", nodes[i - start][k_ - l - 1]);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");

        for (uint64_t i = start; i < end; i++) {
            fprintf(stdout, "%d", get_last(i));
        }
        fprintf(stdout, "\n\n");

        for (uint64_t i = start; i < end; ++i) {
            fprintf(stdout, "%d", static_cast<int>(outdegree(i)));
        }
        fprintf(stdout, "\n");

        for (uint64_t i = start; i < end; ++i) {
            fprintf(stdout, "%d", static_cast<int>(indegree(i)));
        }
        fprintf(stdout, "\n\n");
    }
}

void DBG_succ::print_adj_list(const std::string &filename) const {
    std::ofstream of;
    std::ostream &outstream = filename != ""
                                       ? (of = std::ofstream(filename))
                                       : std::cout;

    for (uint64_t edge = 1; edge < W->size(); ++edge) {
        outstream << rank_last(succ_last(edge))
                  << "\t"
                  << rank_last(outgoing(edge, get_W(edge)))
                  << "\n";
    }
}

///////////////
// Construct //
///////////////

// add a full sequence to the graph
void DBG_succ::add_sequence(const std::string &seq) {
    /*
    std::vector<TAlphabet> ckmer(k_, encode('$'));

    for (size_t i = 0; i < seq.length(); ++i) {
        // print the process
        if (i > 0 && i % 1'000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10'000 == 0)
                verbose_cout(i, " - edges ", num_edges(), " / nodes ", num_nodes(), "\n");
        }

        TAlphabet c = encode(seq[i]);
        append_pos(c, &ckmer[i], 0);
        ckmer.push_back(c);
    }

    verbose_cout("edges ", num_edges(), " / nodes ", num_nodes(), "\n");
    */

    // Padding of the input genome / read
    if (W->size() == 2) {
        for (size_t j = 0; j < k_; j++) {
            append_pos(encode('X'));
        }
    }

    // Iterate over input sequence and enumerae all k-mers.
    std::vector<TAlphabet> ckmer(k_ + 1, encode('X'));
    ckmer[std::min(seq.length(), k_)] = encode('$');

    for (size_t i = 0; i < seq.length(); ++i) {
        if (i > 0 && i % 1'000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10'000 == 0)
                verbose_cout(i, " - edges ", num_edges(), " / nodes ", num_nodes(), "\n");
        }

        TAlphabet c = encode(seq[i]);
        append_pos(c, ckmer.data(), 0);
        memmove(ckmer.data(), &ckmer[1], sizeof(TAlphabet) * (k_ - 1));
        ckmer[k_ - 1] = c;
    }
    assert(get_W(0) == encode('$'));

    // Padding after sequence to get back into default state.
    if (get_W(0) == 0) {
        for (size_t j = 0; j < k_; j++) {
            append_pos(encode('X'), ckmer.data(), 0);
            memmove(ckmer.data(), &ckmer[1], sizeof(TAlphabet) * (k_ - 1));
            ckmer[k_ - 1] = encode('X');
        }
        append_pos(encode('$'), ckmer.data());
    }

    verbose_cout("edges ", num_edges(), " / nodes ", num_nodes(), "\n");
}


bool equal_encodings(const char first, const char second) {
    return DBG_succ::encode(first) == DBG_succ::encode(second);
}

void DBG_succ::add_sequence_fast(const std::string &seq,
                                 bool add_bridge, unsigned int parallel,
                                 std::string suffix) {
    std::vector<uint32_t> label_id;

    // ther is nothing to parse
    if (!seq.size()) {
        return;
    }

    char *bridge = (char*)malloc(k_ + 2);
    memset(bridge, 'X', k_);
    bridge[k_] = seq[0];
    bridge[k_ + 1] = 0;

    size_t i = 0;
    //std::cout << "Loading next sequence with " << parallel << " threads\n";
    if (add_bridge) {
        for (i = 0; i < std::min(k_, seq.length()); ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           bridge + k_ - suffix.length(),
                           equal_encodings)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k_ + 1), DBG_succ::encode));
            }
            memmove(bridge, bridge + 1, k_);
            bridge[k_] = (i + 1 < seq.length()) ? seq[i + 1] : 'X';
        }
    }
    if (k_ < seq.length()) {
        #pragma omp parallel num_threads(parallel)
        {
            std::vector<KMer> kmer_priv;
            #pragma omp for nowait
            for (i = 0; i < seq.length() - k_; ++i) {
                if (std::equal(suffix.begin(), suffix.end(),
                               seq.c_str() + i + k_ - suffix.length(),
                               equal_encodings)) {
                    kmer_priv.push_back(KMer::from_string(
                        std::string(seq.c_str() + i, k_ + 1),
                        DBG_succ::encode
                    ));
                }
            }
            #pragma omp critical
            kmers.insert(kmers.end(),
                std::make_move_iterator(kmer_priv.begin()),
                std::make_move_iterator(kmer_priv.end())
            );
        }
        memcpy(bridge, seq.c_str() + seq.length() - k_, k_);
        bridge[k_] = 'X';
    }
    if (add_bridge) {
        for (i = 0; i < k_; ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           bridge + k_ - suffix.length(),
                           equal_encodings)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k_ + 1), DBG_succ::encode));
            }
            memmove(bridge, bridge + 1, k_);
            bridge[k_] = 'X';
        }
    }
    free(bridge);
}

void DBG_succ::construct_succ(unsigned int parallel) {

    // parallel sort of all kmers
    omp_set_num_threads(std::max(static_cast<int>(parallel), 1));
    __gnu_parallel::sort(kmers.begin(),kmers.end());

    auto last = std::unique(kmers.begin(), kmers.end());
    kmers.erase(last, kmers.end()); 

    //DEBUG: output kmers in current bin
    /*
    std::cerr << "\n";
    for (size_t i=0;i<kmers.size();++i) {
        char* curseq = kmer_to_s(kmers[i].first, alphabet, alph_size);
        std::cerr << kmers[i].first << "\t" << curseq+1 << " " << curseq[0] << " " << kmers[i].second << "\n";
        free(curseq);
    }
    */

    size_t curpos = W_stat.size();
    W_stat.resize(W_stat.size() + kmers.size());

    std::vector<uint8_t> last_stat_safe(last_stat.size() + kmers.size(), true);
    std::copy(last_stat.begin(), last_stat.end(), last_stat_safe.begin());

    #pragma omp parallel num_threads(parallel)
    {
        #pragma omp for nowait
        for (size_t i = 0; i < kmers.size(); ++i) {
            //set last
            if (i + 1 < kmers.size()) {
                bool dup = KMer::compare_kmer_suffix(kmers[i], kmers[i + 1]);
                if (dup) {
                    last_stat_safe[curpos + i] = false;
                }
            }
            //set W
            uint8_t curW = kmers[i][0];
            if (curW == 127) {
                std::string curseq = kmers[i].to_string(alphabet);
                std::cerr << "Failure decoding kmer " << i << "\n" << kmers[i] << "\n" << curseq << "\n";
                exit(1);
            }
            if (!curW && curpos+i)
                p_ = curpos + i;
            if (i) {
                for (size_t j = i - 1; KMer::compare_kmer_suffix(kmers[j], kmers[i], 1); --j) {
                    //TODO: recalculating W is probably faster than doing a pragma for ordered
                    if (kmers[j][0] == curW) {
                        curW += alph_size;
                        break;
                    }
                    if (!j)
                        break;
                }
            }
            W_stat[curpos+i] = curW;
        }
    }
    last_stat.resize(last_stat_safe.size());
    std::copy(last_stat_safe.begin(), last_stat_safe.end(), last_stat.begin());

    size_t lastlet = 0;

    for (size_t i = 0; i < kmers.size(); ++i) {
        char cF = alphabet[kmers[i][k_]];
        if (cF != alphabet[lastlet]) {
            for (lastlet++; lastlet < alph_size; lastlet++) {
                F[lastlet] = curpos + i - 1;
                if (alphabet[lastlet] == cF) {
                    break;
                }
            }
        }
    }
    kmers.clear();
}

/**
 * This function takes a character c and appends it to the end of the graph
 * sequence given that the corresponding note is not part of the graph yet.
 */
uint64_t DBG_succ::append_pos(uint64_t c, TAlphabet *ckmer, uint64_t i) {

    // check that the last position of the graph is indeed a terminal
    assert(get_W(p_) == 0);

    TAlphabet c_p = get_node_end_value(p_);
    // get range of identical nodes (without W) pos current end position
    std::pair<uint64_t, uint64_t> R;
    //std::pair<uint64_t, uint64_t> R = get_equal_node_range(p_);
    //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

    bool append = i ? false : true;
    if (!i) {
        R.second = succ_last(p_);
    } else {
        R = get_equal_node_range(p_);
    }
    TAlphabet &p = i ? R.first : p_;

    // get position of first occurence of c in W after p
    uint64_t next_c = succ_W(p, c);
    // check if c is part of range
    bool exist_c = (next_c <= R.second);
    if (!exist_c) {
        // get position of first occurence of c- in W after p
        next_c = succ_W(p, c + alph_size);
        // check if c- is part of range
        exist_c = (next_c <= R.second);
    }

    /**
     * if the character already exists in the range, we delete the terminal symbol
     * at p, insert c at fwd(next_c) and update p.
     */
    if (exist_c) {
        if (!append) {
            //if not appending, then this has already been observed and we only need to delete the extra pointer
            for (i = 0; get_W(p) == 0; ++i, ++p) {
                if (p == p_)
                    continue;

                assert(i < 2);
                assert(p < R.second);

                W->remove(p);
                last->deleteBit(p);
                update_F(c_p, false);
                if (p <= p_) {
                    p_--;
                    assert(get_W(p) == 0);
                }
                if (p <= next_c)
                    next_c--;
            }
            return fwd(next_c);
        }
        uint64_t p_new = fwd(next_c);
        // remove old terminal symbol
        last->deleteBit(p);
        W->remove(p);
        // adapt position if altered by previous deletion
        p_new -= (p < p_new);
        // insert new terminal symbol
        // we have to insert 0 into last as the node already existed in the range
        // and the terminal symbol is always first
        last->insertBit(p_new, false);
        W->insert(p_new, 0);
        // update new terminal position
        p = p_new;
        // take care of updating the offset array F
        update_F(c_p, false);
        //assert(get_node_end_value(p) == c);
        update_F(c, true);
    } else {
        /**
         * We found that c does not yet exist in the current range and now have to
         * figure out if we need to add c or c- to the range.
         * To do this, we check if there is a previous position j1 with W[j1] == c
         * whose node shares a k-1 suffix with the current node. If yes, we add c-
         * instead of c.
         */
        if (!append) {
            //we need to insert a new pointer
            if (p == p_)
                p++;
            if (get_W(p)) {
                //if no placeholder exists
                W->insert(p, 0);
                update_F(c_p, true);
                last->insertBit(p, false);
                R.second++;
                if (p <= p_) {
                    p_++;
                    assert(get_W(p_) == 0);
                }
            }
        }
        // get position of last occurence of c before p (including p - 1)
        uint64_t last_c = pred_W(p - 1, c);
        // if this position exists
        if (last_c > 0) {
            uint64_t x = fwd(last_c);
            assert((*last)[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

            // check, if there are any c or c- symbols following after position p
            uint64_t next_c = succ_W(p + 1, c);
            uint64_t next_cm = succ_W(p + 1, c + alph_size);
            // there is no c between p and next_cm and next_cm is a c- ==> we should add a c-
            // all nodes between W[i] = c and W[j] = c- share a common suffix of length k-1
            bool minus1 = (next_cm < next_c);
            // check, if we share a k-1 suffix with last_c
            if (!minus1) {
                if (ckmer) {
                    minus1 = compare_node_suffix(ckmer, last_c);
                } else {
                    minus1 = compare_node_suffix(p, last_c);
                }
            }

            // adding a new node can influence following nodes that share a k-1 suffix with the
            // new node -> need to adapt the respektive cc to a cc-
            bool minus2 = false;
            if (next_c < W->size()) {
                if (ckmer) {
                    minus2 = compare_node_suffix(ckmer, next_c);
                } else {
                    minus2 = compare_node_suffix(p, next_c);
                }
                if (minus2) {
                    W_set_value(next_c, get_W(next_c) + alph_size);
                }
            }

            W_set_value(p, minus1 ? c + alph_size : c);
            // after we are done, assert that the order within the range we created
            // is still valid within W
            if (p - R.second > 0) {
                sort_W_locally(p, R.second);
            }

            // if minus1 is true, we share a k-1 suffix with the node at
            // last_c and thus need to adapt our insertion position by -1,
            // as we would like to insert before it. Otherwise we insert directly after
            // it as we are now sorted after it.
            if (minus1) {
                p = x;
                last->insertBit(x, false);
                W->insert(x, 0);
            } else if (minus2) {
                p = x + 1;
                last->insertBit(x + 1, false);
                W->insert(x + 1, 0);
            // no node shares a k-1 suffix with last_c and thus the new node comes after
            // the forward of last_c (as the current node came after last_c as well)
            } else {
                p = x + 1;
                last->insertBit(x + 1, true);
                W->insert(x + 1, 0);
            }
        } else {
            uint64_t x = F[c] + 1;
            uint64_t next_c = succ_W(p + 1, c);
            bool minus = false;
            if (next_c < W->size()) {
                if (ckmer) {
                    minus = compare_node_suffix(ckmer, next_c);
                } else {
                    minus = compare_node_suffix(p, next_c);
                }
            }
            W_set_value(p, c);
            if (p - R.second > 0) {
                sort_W_locally(p, R.second);
            }
            p = x;
            if (minus) {
                W_set_value(next_c, (*W)[next_c] + alph_size);
                last->insertBit(x, false);
            } else {
                last->insertBit(x, true);
            }
            W->insert(x, 0);
        }
        if (p < p_ || (!append && p == p_)) {
            p_++;
            assert(get_W(p_) == 0);
        }
        update_F(c, true);
    }
    // update sorting at new location of p
    // with this we assert that $ is always inserted at the first position
    // of a range of equal nodes --> this will help us to prevent multiple insertions
    // of already existing nodes
    R = get_equal_node_range(p);
    if (R.second - R.first > 0) {
        sort_W_locally(R.first, R.second);
        p = R.first;
        if (!append && p == p_)
            p++;
        /*
        while ((*(W))[p] != 0)
            (p)--;
        */
        assert(get_W(p) == 0);
    }
    return p;
}

/** This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's arrays. In almost all cases this will not produce a valid graph and
 * should only be used as a helper in the parallel merge procedure.
 */
void DBG_succ::append_graph(const DBG_succ &G) {
    verbose_cout("    adding ", G.W->size(), " edges\n");

    // handle last and W
    for (size_t j = 1, curr_pos = W->size(); j < G.W->size(); ++j, ++curr_pos) {
        W->insert(curr_pos, G.get_W(j));
        last->insertBit(curr_pos, G.get_last(j));
    }

    verbose_cout("new total edges: ", W->size(), "\n");

    // handle F
    assert(F.size() == G.F.size());
    for (size_t j = 0; j < F.size(); ++j) {
        F.at(j) += G.F.at(j);
    }
}

/**
 * This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's static containers last_stat and W_stat. In almost all cases
 * this will not produce a valid graph and should only be used as a helper in the
 * parallel merge procedure.
 */
void DBG_succ::append_graph_static(const DBG_succ &G) {
    verbose_cout("    adding ", G.W->size(), " edges\n");

    assert(dynamic_cast<wavelet_tree_dyn*>(G.W));
    auto G_W_stat = dynamic_cast<wavelet_tree_dyn*>(G.W)->to_vector();

    W_stat.insert(W_stat.end(), G_W_stat.begin() + 1, G_W_stat.end());

    for (size_t i = 1; i < G.W->size(); ++i) {
        last_stat.push_back(G.get_last(i));
    }

    verbose_cout("new total edges: ", W->size(), "\n");

    // handle F
    assert(F.size() == G.F.size());
    for (size_t j = 0; j < F.size(); ++j) {
        F.at(j) += G.F.at(j);
    }
}

// Given an edge list, remove them from the graph.
void DBG_succ::remove_edges(const std::set<uint64_t> &edges) {
    uint64_t shift = 0;

    for (const auto &edge : edges) {
        assert(edge >= shift);
        uint64_t edge_id = edge - shift;

        uint64_t d = W->operator[](edge_id);
        if (d < alph_size) {
            //fix W array
            uint64_t next = edge_id + 1;
            uint64_t j = succ_W(next, d);
            for (uint64_t i = next; i < j; ++i) {
                if (W->operator[](i) == d + alph_size) {
                    W_set_value(i, d);
                    break;
                }
            }
        }
        W->remove(edge_id);
        update_F(get_node_end_value(edge_id), false);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        if (get_last(edge_id) && (edge >= shift + 1)
                              && !get_last(edge_id - 1)) {
            last->deleteBit(edge_id - 1);
        } else {
            last->deleteBit(edge_id);
        }
        if (edge_id <= p_) {
            p_--;
            assert(W->operator[](p_) == 0);
        }
        shift++;
    }
}

void DBG_succ::add_sink(unsigned int parallel, std::string suffix) {
    add_sequence_fast(std::string(k_ - 1, alphabet[alph_size - 1]) + "$$",
                      false, parallel, suffix);
    add_sequence_fast("$", true, parallel, suffix);
}