#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

// Defines the alphabets that TensorSketch can operate on
namespace ts {

extern std::function<uint32_t(uint8_t c)> char2int;
extern const char *alphabet;
extern uint8_t alphabet_size;
extern uint8_t bits_per_char;

void init_alphabet(const std::string &alphabet_str);

} // namespace ts
