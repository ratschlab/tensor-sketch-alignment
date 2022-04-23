#pragma once

#include "sequence/alphabets.hpp"
#include "util/utils.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

namespace ts { // ts = Tensor Sketch

/**
 * Represents the contents of a single Fasta file.
 * All the sequences in a file should be treated as a single assembly and should be sketched as a
 * whole.
 */
template <typename seq_type>
struct FastaFile {
    /** The name of the file. */
    std::string filename;
    /** The leading comment before each sequence. Always has the same length as sequences. */
    std::vector<std::string> comments;
    /** The sequences in the file. */
    std::vector<std::vector<seq_type>> sequences;
};

/**
 * Reads a fasta file and returns its contents.
 * @tparam seq_type type used for storing a character of the fasta file, typically uint8_t
 */
template <typename seq_type>
FastaFile<seq_type> read_fasta(const std::string &file_name, const std::string &input_format) {
    FastaFile<seq_type> f;

    if (!std::filesystem::exists(file_name)) {
        std::cerr << "Input file does not exist: " << file_name << std::endl;
        std::exit(1);
    }

    std::ifstream infile(file_name);
    if (!infile.is_open()) {
        std::cout << "Could not open " + file_name << std::endl;
        std::exit(1);
    }

    f.filename = std::filesystem::path(file_name).filename();

    std::string line;
    std::vector<seq_type> seq;
    while (std::getline(infile, line)) {
        if (line[0] == '>') {
            if (!seq.empty()) {
                f.sequences.push_back(std::move(seq));
                seq.clear();
            }
            // Drop the leading '>'.
            f.comments.emplace_back(line.begin() + 1, line.end());
        } else if (!line.empty()) {
            if (input_format == "fasta") {
                for (char c : line) {
                    seq.push_back(char2int(c));
                }
            } else if (input_format == "csv") {
                std::stringstream ss(line);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    seq.push_back(std::stoi(item, 0, 16));
                }
                f.comments.push_back("seq" + std::to_string(f.sequences.size()));
                f.sequences.push_back(std::move(seq));
                seq.clear();
            } else {
                std::cerr << "Invalid input foramt: " << input_format << std::endl;
                exit(1);
            }
        }
    }
    if (!seq.empty()) {
        f.sequences.push_back(std::move(seq));
        seq.clear();
    }
    if(f.sequences.size() != f.comments.size()) {
        std::cerr << "Invalid fasta file: " << file_name << std::endl;
        std::exit(1);
    }
    return f;
}

/**
 * Reads all .fasta and .fna files in the given directory and returns them.
 * @tparam seq_type type used for storing a character of the fasta file, typically uint8_t
 */
template <typename seq_type>
std::vector<FastaFile<seq_type>> read_directory(const std::string &directory_name) {
    if (!std::filesystem::exists(directory_name)) {
        std::cerr << "Input directory does not exist: " << directory_name << std::endl;
        std::exit(1);
    }
    std::vector<FastaFile<seq_type>> files;

    // Handle the case where the argument is a single file as well.
    if (std::filesystem::is_regular_file(directory_name)) {
        files.emplace_back(read_fasta<seq_type>(directory_name, "fasta"));
    } else {
        for (const auto &f : std::filesystem::directory_iterator(directory_name)) {
            const std::filesystem::path ext = f.path().extension();
            if (ext == ".fa" || ext == ".fna" || ext == ".fasta") {
                files.emplace_back(read_fasta<seq_type>(f.path(), "fasta"));
            }
        }
    }

    return files;
}

template <class seq_type>
void write_fasta(const std::string &file_name, const Vec2D<seq_type> &sequences, bool Abc = false) {
    std::ofstream fo(file_name);
    fo << "#" + std::to_string(random()) << std::endl;
    fo << "# " << flag_values(' ') << std::endl;
    for (uint32_t si = 0; si < sequences.size(); si++) {
        fo << ">s" << si << "\n";
        auto &seq = sequences[si];
        for (auto &c : seq) {
            if (Abc) {
                fo << (char)(c + (int)'A');
            } else {
                fo << c << ",";
            }
        }
        fo << "\n\n";
    }
}

} // namespace ts
