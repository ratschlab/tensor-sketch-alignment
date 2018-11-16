#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <string>
#include <vector>
#include <map>
#include <unordered_map>


void serialize_number(std::ostream &out, uint64_t number);

uint64_t load_number(std::istream &in);

template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number = sizeof(T) * 8);

uint64_t load_number_vector_size(std::istream &in);

template <typename T>
std::vector<T> load_number_vector(std::istream &in);

void serialize_number_number_map(std::ostream &out,
                                 const std::map<uint32_t, uint32_t> &M);

std::map<std::uint32_t, uint32_t> load_number_number_map(std::istream &in);

void serialize_string_number_map(std::ostream &out,
                                 const std::unordered_map<std::string, uint32_t> &M);

std::unordered_map<std::string, uint32_t> load_string_number_map(std::istream &in);


#endif // __SERIALIZATION_HPP__