#ifndef __SEQGEN_GRAPH_HPP__
#define __SEQGEN_GRAPH_HPP__

#include <memory>

namespace mtg {

namespace graph {
namespace seqgen {
struct DBGAlignerConfig;
} // namespace align
} // namespace graph

namespace cli {
class Config;

int generate_sequences(Config *config);

} // namespace cli
} // namespace mtg

#endif // __SEQGEN_GRAPH_HPP__
