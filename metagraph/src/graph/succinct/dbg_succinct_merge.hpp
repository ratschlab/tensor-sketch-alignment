#ifndef __MERGE_HPP__
#define __MERGE_HPP__

#include <cstdint>
#include <vector>
#include <string>

#include "dbg_succinct.hpp"
#include "dbg_succinct_chunk.hpp"


namespace merge {

    /**
     * Given a list of graph structures, this functions
     * integrates all of them into a new graph G.
     */
    DBG_succ* merge(const std::vector<const DBG_succ*> &graphs,
                    bool verbose = false);

    DBG_succ::Chunk* merge_blocks_to_chunk(const std::vector<const DBG_succ*> &graphs,
                                           size_t chunk_idx,
                                           size_t num_chunks,
                                           size_t num_threads,
                                           size_t num_bins_per_thread,
                                           bool verbose = false);

} // namespace merge

#endif // __MERGE_HPP__