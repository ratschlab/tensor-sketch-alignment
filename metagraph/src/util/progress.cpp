//
// Created by amir on 1/9/21.
//
#include "util/progress.hpp"
#include <iomanip>
#include <iostream>

namespace ts {

size_t progress_bar::it;
size_t progress_bar::total;
size_t progress_bar::bar_len;
size_t progress_bar::bar_step;

void progress_bar::init(size_t total_iterations, size_t len) {
    progress_bar::it = 0;
    progress_bar::total = total_iterations;
    progress_bar::bar_len = len;
    progress_bar::bar_step = 0;
}

void progress_bar::iter() {
#pragma omp critical
    {
        ++it;
        auto step = (it * bar_len) / total;
        while (step > bar_step) {
            if (bar_step > 0)
                std::cerr << "\b\b\b\b";
            ++bar_step;
            std::cerr << "#" << std::setw(3) << (int)(100.0 * it / total) << "%" << std::flush;
        }
        if (it == total) {
            std::cerr << "\033[2K\r" << std::flush;
        }
    }
}

}; // namespace ts
