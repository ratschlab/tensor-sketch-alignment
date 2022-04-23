#include "sketch/hash_base.hpp"

namespace ts {

HashAlgorithm parse_hash_algorithm(const std::string &name) {
    if (name == "uniform") {
        return HashAlgorithm::uniform;
    }
    if (name == "crc32") {
        return HashAlgorithm::crc32;
    }
    if (name == "murmur") {
        return HashAlgorithm::murmur;
    }
    assert(false);
}

} // namespace ts
